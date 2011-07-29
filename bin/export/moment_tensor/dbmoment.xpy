"""
dbmoment.py

@author   Gert-Jan van den Hazel <hazelvd@knmi.nl>, Rob Newman <robertlnewman@gmail.com>
@notes    Look for initialized comments - they define questions to be answered
          e.g. # RLN (2011-07-28): Why this hard-coded value here?
"""

from optparse import OptionParser
import re
import numpy as np
import math as math
from datetime import datetime
from collections import defaultdict
from time import gmtime, time

import antelope.stock as stock
import antelope.datascope as antdb

# RLN (2011-07-27): Setting global variables (or defaults) is bad programming in Python. Comment out.
# Set defaults
'''
pfname  = 'dbmoment'
chan_to_use = 'LH*'
model_type = 'v'
isoflag = 5
dw = False
use_incidence = 0
trim_value = 0.6
statmax = 20
event_db = ''
green_db = ''
wave_db = ''
resp_db = ''
filters = {}
mag_filters = {}
always = True
verbose = False
debug = False
flag = 0
'''

def configure():
    """Function to configure parameters 
    from command-line and the pf-file
    Check command line options
    """
    usage = "Usage: dbmoment [-vV] [-p pfname] orid"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="verbose output")
    parser.add_option("-V", "--debug", action="store_true", dest="debug", default=False, help="debug application")
    parser.add_option("-p", "--pf", action="store", dest="pf", type="string", help="parameter file path")
    (options, args) = parser.parse_args()
    if options.verbose:
        verbose = True
        debug = False
    if options.debug:
        debug = True
        verbose = True
    # RLN (2011-07-29): There should be a better way to search all the PFPATHs for the parameter file!
    if options.pf:
        if os.path.isfile(options.pf):
            pfname = options.pf
        else:
            print "Parameter file %s does not exist. Exiting" % options.pf
    elif os.path.isfile('pf/dbmoment.pf'): # Use local pf if available
        pfname = 'pf/dbmoment.pf'
    elif os.path.isfile('%s/data/pf/dbmoment.pf' % os.environ['ANTELOPE']): # Last option, use default Antelope pf path
        pfname = 'pf/dbmoment.pf'
    else:
        print "You have not defined a valid parameter file"
    if len(args) != 1:
        print usage;
        sys.exit(-1)
    else:
        orid = args[0]
    return pfname, orid, verbose, debug

class MomentTensor():
    """Class for building moment tensors
    """

    def __init__(self, pfname, orid, verbose=False, debug=False):
        """Initialize"""
        self.pfname = pfname
        self.orid = orid
        self.verbose = verbose
        self.debug = debug

    def logmt(self, flag, message):
        """Function to handle log output, which depends 
        on the verbosity setting. Prints messages with a 
        timestamp. Whenever an error occurs it also exits.
        """
        curtime = stock.strtime(time())
        if not flag:
            return
        elif flag < 3:
            print curtime, message
        else:
            print curtime, 'ERROR:',message,'--> exiting'
            sys.exit(-1)

    def parse_pf(self):
        """Parse the parameter file
        and assign to self vars
        """
        if self.debug:
            self.logmt(1, 'Parse parameter file %s' % self.pfname)
        self.chan_to_use = stock.pfget_string(self.pfname, 'chan_to_use')
        self.trim_value = stock.pfget_string(self.pfname, 'trim_value')
        if stock.pfget_string(self.pfname, 'isoflag') == 1: 
            self.isoflag = 6
        else:
            self.isoflag = 5
        self.model_type = stock.pfget_string(self.pfname, 'model_type')
        self.statmax = stock.pfget_int(self.pfname,'statmax')
        self.use_inc = stock.pfget_string(self.pfname, 'use_inc')
        self.event_db = stock.pfget_string(self.pfname,'event_db')
        try:
            self.wave_db = stock.pfget_string(self.pfname, 'wave_db')
        except:
            self.wave_db = False
        self.green_db = stock.pfget_string(self.pfname, 'green_db')
        self.green_pf = stock.pfget_string(self.pfname, 'green_pf')
        try:
            self.resp_db = stock.pfget_string(self.pfname, 'resp_db')
        except:
            self.resp_db = False
        self.mag_filters = stock.pfget_arr(self.pfname,'mag_filters')
        if stock.pfget_string(self.pfname,'distance_weighting') == 'on':
            self.dw = True
        else:
            self.dw = False
        if not self.wave_db:
            self.wave_db = self.event_db
        if not self.resp_db:
            self.resp_db = self.event_db

    def get_view_from_db(self):
        """Open event database and get 
        event parameters for given origin-id
        Returns event parameters or 
        exit when an error occurs
        """
        if self.debug:
            self.logmt(1, 'Get canned view from database')
        if not os.path.isfile(self.event_db):
            self.logmt(3, 'Database (%s) does not exist' % self.event_db)
        evdb = antdb.dbopen(self.event_db, 'r')
        evdb.lookup(table='origin')
        self.logmt(1, 'Processing origin: %s' % self.orid)
        evdb.subset('orid == %s' % self.orid)
        if evdb.nrecs() == 0:
            self.logmt(3, 'Orid %s does not exist in origin table' % self.orid)
        evdb.join('netmag')
        if evdb.nrecs() == 0:
            self.logmt(3, 'Could not join netmag table for orid %s' % self.orid)
        evparams = {}
        evdb[3] = 0 # Just get the first record
        for f in evdb.query('dbTABLE_FIELDS'):
            try: 
                evparams[f] = evdb.getv(f)[0]
            except:
                self.logmt(3, 'Could not find field %s in join of origin & netmag tables for orid %s' % (f, self.orid))
        evdb.join('assoc')
        evdb.join('arrival')
        evdb.sort('delta')
        evdb.subset('iphase=~/.*P.*|.*p.*/')
        if evdb.nrecs() == 0:
            self.logmt(3, 'No arrivals for selected origin %s' % self.orid)
        lower_mags = []
        upper_mags = []
        filters = []
        for line in self.mag_filters:
            splitline = line.split()
            if not len(splitline) == 3:
                lower_magnitude, filter = splitline
            else:
                lower_magnitude, upper_magnitude, filter = splitline
            lower_mags.append(lower_magnitude)
            upper_mags.append(upper_magnitude)
            filters.append(filter)
        min_mag = lower_mags[0] # RLN (2011-07-28): Assumes order correct in pf
        max_mag = upper_mags[-1] # RLN (2011-07-28): Assumes order correct in pf
        filter = filters[-1] # RLN (2011-07-28): Assumes all filters are the same
        if float(evparams['magnitude']) > float(min_mag) and float(evparams['magnitude']) < float(max_mag):
            self.filter_string = filter.replace('_', ' ')
        return evdb, evparams
 
    def get_chan_data(self, dbptr):
        """Opens waveform database and returns 
        trace objects based on sta_chan. Applies 
        calibration, splice, filter, decimation 
        and rotation of channels
        """
        if self.verbose or self.debug:
            self.logmt(1, 'Get channel data (trace objects)')
        try:
            wvdb = antdb.dbopen(self.wave_db, 'r')
        except:
            self.logmt(3, 'Could not open waveform database %s' % self.wave_db)
        wvdb.lookup(table='wfdisc')
        wvdb.subset('samprate>=0.9') # RLN (2011-07-28): Why this hard-coded value here?
        try:
            wvdb.subset('chan =~/%s/' % self.chan_to_use)
        except:
            self.logmt(3, 'No records found in wfdisc table (%s.wfdisc)' % self.wave_db) 
        numrec = dbptr.nrecs()
        if numrec > int(self.statmax):
            numrec = int(self.statmax)
        self.logmt(1, 'Processing %s stations for orid %s' % (numrec, self.orid))
        stachan_traces = {}
        counter = 0
        for i in range(dbptr.nrecs()):
            dbptr[3] = i
            sta, at, esaz = dbptr.getv('sta','arrival.time','esaz')
            st = at - 12 
            et = at + 12
            wvstadb = antdb.dbsubset(wvdb,'sta=~/^%s$/' % sta)
            self.logmt(1, 'Looking for channels with a sample rate >= 1 Hz for sta = %s' % sta)
            try:
                chans = self.choose_chan(wvstadb)
            except:
                self.logmt(1, 'No channels found with a sample-rate >= 1 Hz for sta = %s' % sta) # RLN (2011-07-28): Why hard-coded sample rate here?
            self.logmt(1, 'Channels found: %s %s %s' % (chans[0], chans[1], chans[2]))
            wvstadb.subset('chan =~ /%s|%s|%s/' % (chans[0], chans[1], chans[2]))
            resample = 0
            for j in range(wvstadb.nrecs()):
                wvstadb[3] = j
                if wvstadb.getv('samprate')[0] != 1:
                    logmt(0, 'Samplerate higher than 1.0 --> resample')
                    resample = 1
                break
            if resample == 1:
                logmt(0, 'Resampling to be implemented, skipping station %s for now' % sta) # RLN (2011-07-28): What does this mean?
                continue
            vang = 0.0 # RLN (2011-07-28): Why hard code the vang? You can get this from the db?
            if self.use_inc == 1:
                vang = dbptr.getv('vang')[0]
            trace = wvstadb.load_css(st, et)
            trace.apply_calib()
            trace.splice()
            trace.filter(self.filter_string)
            rotchan = ('R', 'T', 'Z') # RLN (2011-07-28): Why hard-code this? Should be defined at the config phase
            trace.rotate(esaz, vang, rotchan)
            trace.subset('chan =~ /R|T|Z/') # (2011-07-28): Why hard-code this? Should be defined at the config phase
            foundchans = []
            for j in range(trace.nrecs()):
                trace[3] = j
                sta, chan, ns, sr = trace.getv('sta', 'chan', 'nsamp', 'samprate')
                stachan = '%s_%s' % (sta, chan)
                ns_tra = 0
                ns_req = 0
                samples = trace.data()
                for k in range(len(samples)):
                    ns_tra += ns
                    ns_req += int((et-st)*sr)
                if not ns_tra == ns_req:
                    logmt(0, 'Incorrect number of samples for %s, skipping %s' % (stachan, sta))
                    for delchan in foundchans:
                        if stachan_traces[delchan]:
                            del stachan_traces[delchan]
                            logmt(0, 'Deleting channel %s' % delchan)
                            counter -= 1
                stachan_traces[stachan] = samples  
                self.logmt(1, 'Trace extracted for %s' % stachan)
                foundchans.append(stachan)
                counter += 1
            if counter == self.statmax*3:
                break
        return stachan_traces

    def choose_chan(self, dbptr):
        """Get all channels for a given station. 
        Returns the channels with samplerate 
        closest to and bigger than 1
        """
        if self.debug:
            self.logmt(1, 'Choosing channels')
        view = antdb.dbsort(dbptr, 'samprate', unique=True)
        for i in range(view.nrecs()):
            view[3] = i
            chan = view.getv('chan')
            chanview = antdb.dbsubset(dbptr,'chan =~ /%s.*/' % chan[0][:2])
            chanview.sort('chan', unique=True)
            channels = []
            if chanview.nrecs() == 3:
                for j in range(chanview.nrecs()):
                    chanview[3] = j
                    chan = chanview.getv('chan')[0]
                    channels.append(chan)
            else:
                self.logmt(1, 'Did not find 3 components for %s, trying higher sample rate' % chan)
                channels = []
                continue
            if len(channels) == 3:
                break
        return channels

    def construct_data_matrix(self, stachan_traces, this_dict):
        """Construct data matrix from dictionary 
        containing trace_objects per sta_chan
        Returns matrix s, which is a 3D matrix, 
        containing all data per channel
        """
        if self.debug:
            self.logmt(1, 'Constructing data matrix from trace objects')
        numsta = 0
        numchan = 0
        for stachan, tr in sorted(stachan_traces.items()):
            sta, chan = stachan.split('_')
            if self.dw == True:
                self.logmt(1, 'Apply distance weighting')
            if chan == 'R':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['R'][numsta][j] = tr[j]
            if chan == 'T':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['T'][numsta][j] = tr[j]
            if chan == 'Z':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['Z'][numsta][j] = tr[j]
            if numchan == 3:
                numchan = 0
                numsta += 1
        return this_dict

    def construct_data_matrix_hard_coded(self, stachan_traces, this_dict):
        """Construct data matrix from files
        containing trace_objects per sta_chan
        Returns matrix s, which is a 3D matrix, 
        containing all data per channel
        """
        if self.debug:
            self.logmt(1, 'Constructing data matrix from hard-coded files in db/data/data1-3')
        for j in range(3): # Forced to be 3
            file_suffix = j+1
            data = open('db/data/data%d' % file_suffix, 'r')
            samples = []
            for line in data:
                tmp = line.split()
                for samp in tmp:
                    samples.append(float(samp))
            data.close()
            # RLN (2011-07-28): 
            # len(samples) = 600
            # Three channels of 200 samples as per Dregers format
            for i in range(200):
                this_dict['T'][j][i] = samples.pop(0)
            for i in range(200):
                this_dict['R'][j][i] = samples.pop(0)
            for i in range(200):
                this_dict['Z'][j][i] = samples.pop(0)
        return this_dict

    def cross_cor(self, a, b):
        """Calculate cross correlation between data and accompanying 
        Green's function. Returns the max_cor and the shift
        """
        if self.debug:
            self.logmt(1, 'Calculate cross correlation between data and Greens function')
        xcor = np.correlate(a.values(), b.values(), 'full')
        pr = len(xcor)/2 
        maxval = 0
        maxcoef = 0
        for j in range(pr):
            if abs(xcor[j+pr]) > maxval:
                maxval = abs(xcor[j+pr])
                maxcoef = j+1
        return maxval, maxcoef

    def get_greens_functions(self, dbptr):
        """Opens Green's database and returns a 
        dictionary which describes the path to 
        the Green's function per station. Based on azimuth
        and distance. Some settings required from 
        parameter-file for the Green's function creation.

        RLN (2011-07-28): Need to get this working
        NOT CURRENTLY IN USE
        """
        tmp = {}
        if not os.path.isfile(self.green_pf):
            if not os.path.isfile('%s/data/%s' % (os.environ['ANTELOPE'], self.green_pf)):
                green_pf_error = True
        if green_pf_error:
            self.logmt(1, 'Cannot open up Greens functions parameter file %s' % self.green_pf)
            return False
        ddist = stock.pfget_string(self.green_pf,'ddist')
        dazim = stock.pfget_string(self.green_pf,'dazim')
        ddip = stock.pfget_string(self.green_pf,'ddip')
        ddepth = stock.pfget_string(self.green_pf,'ddepth')
        if not os.path.isfile(self.green_db):
            logmt(3, 'Database (%s) does not exist' % self.green_db)
        gdb = antdb.dbopen(green_db, 'r')
        gdb.lookup(table='moment_tensor_greensfuncs')
        for i in range(dbptr.nrecs()):
            dbptr[3] = i
            sta, delta, seaz = evdb.getv('sta', 'delta', 'seaz')
            ddist = float(ddist)
            dazim = float(dazim)
            expr = 'delta > %s && delta <= %s ' % (delta-ddist/2, delta+ddist/2)
            expr += '&& azimuth > %s && azimuth <= %s ' % (seaz-dazim/2, seaz+dazim/2)
            try:
                subgdb = gdb.subset(expr)
            except:
                self.logmt(1, 'No Greens function found for %s' % sta)
            else:
                subgdb[3] = 0
                file = subgdb.extfile()       
                tmp[sta] = file
        return tmp

    def get_greens_functions_hard_coded(self, this_dict):
        """Hard coded way of getting synthetics (Green's
        function). Remove before production.
        """
        if self.debug:
            self.logmt(1, 'Constructing Greens function matrix from hard-coded file in db/data/green')
        # Container list for values
        greens = []
        # Open up the Green function file
        green = open('db/data/green', 'r')
        for line in green:
            # RLN (2011-07-28): Every 12 values?
            for j in range(len(line)/12):
                # RLN (2011-07-28): Append 12 values to the list. So we one giant list. Why?
                greens.append(line[j*12:j*12+12]) 
        green.close()

        # RLN (2011-07-28): Some hard-coded translations to populate Greens function matrix 'g'
        for i in range(3):
            for k in range(8):
                for j in range(200):
                    this_dict[k][i][j] = float(greens[j + k*200])
                    if k in [5,6,7]:
                        this_dict[k][i][j] *= -1
            if self.isoflag == 5:
                for k in [8,9]:
                    for j in range(200):
                        this_dict[k][i][j] = 0.0
            if self.isoflag == 6:
                for k in [8,9]:
                    for j in range(200):
                        this_dict[k][i][j] = float(greens[j + k*200])
        return this_dict

    def get_time_shift(self, data_dict, greens_dict):
        """Get the time shift and return a list
        """
        if self.debug:
            self.logmt(1, 'Get the time shift and return a list')
        timeshift = []
        for i in range(len(data_dict['T'])):
            shift = 0
            xcor  = 0
            if self.cross_cor(data_dict['T'][i], greens_dict[0][i])[0] > xcor:
                xcor  = self.cross_cor(data_dict['T'][i], greens_dict[0][i])[0]
                shift = self.cross_cor(data_dict['T'][i], greens_dict[0][i])[1]
            if self.cross_cor(data_dict['T'][i], greens_dict[1][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['T'][i], greens_dict[1][i])[0]
                shift = self.cross_cor(data_dict['T'][i], greens_dict[1][i])[1]
            if self.cross_cor(data_dict['R'][i], greens_dict[2][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['R'][i], greens_dict[2][i])[0]
                shift = self.cross_cor(data_dict['R'][i], greens_dict[2][i])[1]
            if self.cross_cor(data_dict['R'][i], greens_dict[3][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['R'][i], greens_dict[3][i])[0]
                shift = self.cross_cor(data_dict['R'][i], greens_dict[3][i])[1]
            if self.cross_cor(data_dict['R'][i], greens_dict[4][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['R'][i], greens_dict[4][i])[0]
                shift = self.cross_cor(data_dict['R'][i], greens_dict[4][i])[1]
            if self.cross_cor(data_dict['Z'][i], greens_dict[5][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['Z'][i], greens_dict[5][i])[0]
                shift = self.cross_cor(data_dict['Z'][i], greens_dict[5][i])[1]
            if self.cross_cor(data_dict['Z'][i], greens_dict[6][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['Z'][i], greens_dict[6][i])[0]
                shift = self.cross_cor(data_dict['Z'][i], greens_dict[6][i])[1]
            if self.cross_cor(data_dict['Z'][i], greens_dict[7][i])[0] > xcor:
                xcor = self.cross_cor(data_dict['Z'][i], greens_dict[7][i])[0]
                shift = self.cross_cor(data_dict['Z'][i], greens_dict[7][i])[1]
            timeshift.append(shift)
        return timeshift

    def matrix_AIV_B(self, dict_s, dict_g, list_az, this_timeshift):
        """Construct matrices AIV and B using 
        the data and Green's function matrices
        Return AIV and B for further processing 
        (inversion)
        """
        if self.debug:
            self.logmt(1, 'Construct matrices AIV and B using the data and Greens function matrices')
        AJ = defaultdict(dict) # RLN (2011-07-29): What is AJ? 
        trim = 0
        cnt1 = cnt2 = cnt3 = 0
        trim  = int(len(dict_s['T'][0]) * float(self.trim_value))
        for i in range(len(dict_s)):
            cnt1 = cnt2 = cnt3
            cnt2 += len(dict_s['T'][0]) - trim
            cnt3 += 2*len(dict_s['T'][0]) - 2*trim
            list_az[i] *= math.pi/180
            for j in range(len(dict_s['T'][0])-trim):
                AJ[0][cnt1] =  math.sin(2*list_az[i])*dict_g[0][i][j]/2
                AJ[1][cnt1] = -math.sin(2*list_az[i])*dict_g[0][i][j]/2
                AJ[2][cnt1] = -math.cos(2*list_az[i])*dict_g[0][i][j]
                AJ[2][cnt2] = -math.sin(2*list_az[i])*dict_g[2][i][j]
                AJ[2][cnt3] = -math.sin(2*list_az[i])*dict_g[5][i][j]
                AJ[3][cnt1] = -math.sin(list_az[i])*dict_g[1][i][j]
                AJ[3][cnt2] =  math.cos(list_az[i])*dict_g[3][i][j]
                AJ[3][cnt3] =  math.cos(list_az[i])*dict_g[6][i][j]
                AJ[4][cnt1] =  math.cos(list_az[i])*dict_g[1][i][j]
                AJ[4][cnt2] =  math.sin(list_az[i])*dict_g[3][i][j]
                AJ[4][cnt3] =  math.sin(list_az[i])*dict_g[6][i][j]
                if self.isoflag == 5:
                    AJ[0][cnt2] = (dict_g[4][i][j])/2 - (math.cos(2*list_az[i])*dict_g[2][i][j])/2
                    AJ[0][cnt3] = (dict_g[7][i][j])/2 - (math.cos(2*list_az[i])*dict_g[5][i][j])/2
                    AJ[1][cnt2] = (dict_g[4][i][j])/2 + (math.cos(2*list_az[i])*dict_g[2][i][j])/2
                    AJ[1][cnt3] = (dict_g[7][i][j])/2 + (math.cos(2*list_az[i])*dict_g[5][i][j])/2
                if self.isoflag == 6:
                    AJ[0][cnt2] = (dict_g[4][i][j])/6 - (math.cos(2*list_az[i])*dict_g[2][i][j])/2 + (dict_g[8][i][j])/3
                    AJ[0][cnt3] = (dict_g[7][i][j])/6 - (math.cos(2*list_az[i])*dict_g[5][i][j])/2 + (dict_g[9][i][j])/3
                    AJ[1][cnt2] = (dict_g[4][i][j])/6 + (math.cos(2*list_az[i])*dict_g[2][i][j])/2 + (dict_g[8][i][j])/3
                    AJ[1][cnt3] = (dict_g[7][i][j])/6 + (math.cos(2*list_az[i])*dict_g[5][i][j])/2 + (dict_g[9][i][j])/3
                    AJ[5][cnt1] = 0.0
                    AJ[5][cnt2] = (dict_g[8][i][j])/3  - (dict_g[4][i][j])/3
                    AJ[5][cnt3] = (dict_g[9][i][j])/3 - (dict_g[7][i][j])/3
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1
        AIV = defaultdict(dict) 
        for i in range(self.isoflag):
            for j in range(self.isoflag):
                AIV[i][j] = 0.0
        for i in range(5):
            for j in range(5):
                for k in range(cnt3):
                    AIV[i][j] += AJ[i][k]*AJ[j][k]
        B = defaultdict(dict) 
        for i in range(self.isoflag):
            B[i][0] = 0.0
        cnt1 = cnt2 = cnt3 = 0
        tmp = defaultdict(dict) 
        for i in range(len(dict_s)):
            cnt1 = cnt2 = cnt3
            cnt2 += len(dict_s['T'][0])-trim
            cnt3 += 2*(len(dict_s['T'][0])-trim)
            for j in range(len(dict_s['T'][0])-trim):
                tmp[cnt1] = dict_s['T'][i][j+this_timeshift[i]]
                tmp[cnt2] = dict_s['R'][i][j+this_timeshift[i]]
                tmp[cnt3] = dict_s['Z'][i][j+this_timeshift[i]]
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1
        for i in range(self.isoflag):
            for j in range(cnt3):
                B[i][0] += AJ[i][j]*tmp[j]
        if len(AIV) != self.isoflag or len(AIV[0]) != self.isoflag:
            self.logmt(3, 'Matrix AIV has dimension [%s,%s] should be [%s,%s]' % (len(AIV), len(AIV), self.isoflag, self.isoflag))
        elif len(B) != self.isoflag:
            self.logmt(3, 'Matrix AIV has dimension [%s,i%s] should be [%s,1]' % (len(B), len(B[0]), self.isoflag))
        elif self.verbose or self.debug:
            self.logmt(1, 'Matrices AIV and B created and have correct dimensions')
        return AIV, B

    def swap(self, a, b):
        """Remap list
        """
        tmp = a
        a = b
        b = tmp
        return(a,b)

    def dyadic(self, v, n1, n2, c):
        """Compute the dyadic matrix 
        of vector v
        """
        if self.debug:
            self.logmt(1, 'Compute the dyadic matrix of vector v')
        tmp = np.matrix([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        for i in range(3):
            for j in range(3):
                tmp[i,j] = v[i, n1]*v[j, n2]*c
        return tmp

    def det_solution_vec(self, dict_AIV, dict_B):
        """Solve the inversion problem, 
        returning the moment tensor
        """
        if self.debug:
            self.logmt(1, 'Solve the inversion problem and return the moment tensor')
        ipiv = defaultdict(dict)
        for i in range(self.isoflag):
            ipiv[i] = 0
        for i in range(self.isoflag):
            big = 0.0
            for j in range(self.isoflag):
                if ipiv[j] != 1:
                    for k in range(self.isoflag):
                        if ipiv[k] == 0:
                            if abs(dict_AIV[j][k]) >= big:
                                    big = abs(dict_AIV[j][k])
                                    irow = j
                                    icol = k
                        elif ipiv[k] > 1:
                            self.logmt(3, 'det_solution_vec(): ERROR... 1')
            ipiv[icol] += 1
            if not irow == icol:
                for l in range(self.isoflag):
                    (dict_AIV[irow][l],dict_AIV[icol][l]) = swap(dict_AIV[irow][l],dict_AIV[icol][l])
                for l in range(1):
                    (dict_B[irow][l],dict_B[icol][l]) = swap(dict_B[irow][l],dict_B[icol][l])
            if dict_AIV[icol][icol] == 0.0:
                self.logmt(3, 'det_solution_vec(): ERROR... 2')
            pivinv = 1.0/dict_AIV[icol][icol]
            dict_AIV[icol][icol] = 1.0
            for l in range(self.isoflag):
                dict_AIV[icol][l] *= pivinv
            for l in range(1):
                dict_B[icol][l] *= pivinv
            for h in range(self.isoflag):
                if h != icol:
                    dum = dict_AIV[h][icol]
                    dict_AIV[h][icol] = 0.0
                    for l in range(self.isoflag):
                        dict_AIV[h][l] -= dict_AIV[icol][l]*dum
                    for l in range(1):
                        dict_B[h][l] -= dict_B[icol][l]*dum
        if self.isoflag == 6:
            M = np.matrix([[dict_B[0][0], dict_B[2][0], dict_B[3][0]], [dict_B[2][0], dict_B[1][0], dict_B[4][0]], [dict_B[3][0], dict_B[4][0], dict_B[5][0]]])
        if self.isoflag == 5:
            M = np.matrix([[dict_B[0][0], dict_B[2][0], dict_B[3][0]], [dict_B[2][0], dict_B[1][0], dict_B[4][0]], [dict_B[3][0], dict_B[4][0], -(dict_B[0][0]+dict_B[1][0])]])
        
        """
        # Coded in the original code by Doug Dreger, however, AIV is not used further in the program, so not needed???
        indc = defaultdict(dict)
        indr = defaultdict(dict)
        for i in range(self.isoflag-1,0,-1):
            if indr[i] != indc[i]:
                for j in range(self.isoflag):
                    (dict_AIV[j][indr[i]],dict_AIV[j][indxc[i]]) = swap(dict_AIV[j][indr[i]],dict_AIV[j][indxc[i]])
        """
        if len(M) != 3:
            self.logmt(3, 'Wrong dimension returned for matrix M after inversion')
        elif self.debug:
            self.logmt(1, 'Moment tensor matrix M[3,3] created')
        return M

    def decompose_moment_tensor(self, matrix_M):
        """Decompose moment tensor into eigenvector/values.
        Calculate moment tensor parameters from eigenvalues and vectors. 
        Return M0, Mw, strike, slip, rake and precentages of present
        source characteristics
        """
        if self.debug:
            self.logmt(1, 'Decompose moment tensor into eigenvector/values')
        matrix_M *= -1.0e+20
        trace = 0
        for i in range(3):
            trace += matrix_M[i,i]
        trace /= 3
        
        for i in range(3):
            matrix_M[i,i] -= trace
        miso = np.matrix([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
        eval, evec = np.linalg.eig(matrix_M)
        
        for i in (0,1):
            k = i
            p = eval[i]
            for j in (1,2):
                if abs(eval[j]) < abs(p):
                    k = j
                    p = eval[j]
            if k != i:
                eval[k] = eval[i]
                eval[i] = p
                for j in range(3):
                    p=evec[j,i]
                    evec[j,i] = evec[j,k]
                    evec[j,k] = p
        
        f = -eval[0]/eval[2]
        c = eval[2]*(1-2*f)
        a2a2 = self.dyadic(evec, 2, 2, c)
        c *= -1
        a1a1 = self.dyadic(evec, 1, 1, c)
        mdc = a2a2+ a1a1
        c = 2*eval[2]*f
        a2a2 = self.dyadic(evec, 2, 2, c)
        c = -eval[2]*f
        a1a1 = self.dyadic(evec, 1, 1, c)
        a0a0 = self.dyadic(evec, 0, 0, c)
        mclvd = a2a2+a1a1+a0a0
        
        dd = []
        for i in range(3):
            dd.append(abs(eval[i]))
        for i in range(3):
            for j in range(i,3):
                if dd[j] < dd[i]:
                    dd[i],dd[j] = swap(dd[i],dd[j])
        eps = dd[0]/dd[2]
        pcdc = 100*(1-2*eps)
        pcclvd = 200*eps
        
        for i in range(3):
            if evec[2,i] < 0:
                for j in range(3):
                    evec[j,i] *= -1
        
        azimuth = []
        plunge = [] 
        for i in range(3):
            if evec[1,i] == 0 and evec[0,i] == 0:
                azimuth.append(0.0)
            else:
                tmp = math.degrees(math.atan2(evec[1,i],evec[0,i]))
                if tmp < 0:
                    tmp += 360
                azimuth.append(tmp)
            r = math.sqrt(evec[0,i]*evec[0,i] + evec[1,i]*evec[1,i])
            if evec[2,i] == 0 and r == 0:
                plunge.append(0.0)
            else:
                plunge.append(math.atan2(evec[2,i],r))
        
        axis = []
        axis.append('N')
        if eval[1] > eval[2]:
            axis.append('T')
            axis.append('P')
        else:
            axis.append('P')
            axis.append('T')
        
        p = []
        t = []
        for i in range(3):
            if axis[i] == 'P':
                for j in range(3):
                    p.append(evec[j,i])
            elif axis[i] == 'T':
                for j in range(3):
                    t.append(evec[j,i])
        con = 1/math.sqrt(2)
        tmp1 = []
        tmp2 = []
        for i in range(3):
            tmp1.append(con*(t[i]+p[i]))
            tmp2.append(con*(t[i]-p[i]))
        u  = np.matrix([[tmp1[0],tmp1[1],tmp1[2]],[tmp2[0],tmp2[1],tmp2[2]]])
        nu = np.matrix([[tmp2[0],tmp2[1],tmp2[2]],[tmp1[0],tmp1[1],tmp1[2]]])
        
        dip = []
        slip = []
        strike = []
        for i in range(2):
            dip.append(math.acos(-nu[i,2]))
            if nu[i,0] == 0 and nu[i,1] == 0:
                strike.append(0)
            else:
                strike.append(math.atan2(-nu[i,0],nu[i,1]))
        for i in range(2):
            sstr = math.sin(strike[i])
            cstr = math.cos(strike[i])
            sdip = math.sin(dip[i])
            cdip = math.cos(dip[i])
            if abs(sdip) > 0:
                lamb = math.asin(-u[i,2])/math.sin(dip[i])
            else:
                if u[i,2] > 0:
                    arg = 1
                else:
                    arg = -1
                if arg < 0:
                    lamb = math.pi
                else:
                    lamb = 0
            slamb = math.sin(lamb)
            cdsl = cdip*slamb
            if abs(sstr) > abs(cstr):
                clamb = (u[i,1]+cdsl*cstr)/sstr
            else:
                clamb = (u[i,0]-cdsl*sstr)/cstr
            if slamb == 0 and clamb == 0:
                slip.append(0)
            else:
                slip.append(math.atan2(slamb,clamb))
            if dip[i] > math.pi/2:
                dip[i] = math.pi - dip[i]
                strike[i] += math.pi
                slip[i] = 2*math.pi - slip[i]
            if strike[i] < 0:
                strike[i] += 2*math.pi
            if slip[i] > math.pi:
                slip[i] -= 2*math.pi
        
        m0 = abs(miso[0,0]) + abs(eval[2])
        Mw = math.log10(m0)/1.5 - 10.7
        pciso = abs(miso[0,0])/m0

        if self.verbose or self.debug:
            self.logmt(1, 'Decomposition of moment tensor succesful')

        return m0, Mw, strike, dip, slip, pcdc, pcclvd, pciso

    def fitcheck(self, dict_g, dict_s, list_W, matrix_M, m0, timeshift, az):
        """Calculate the variance, variance reduction
        and flag bad stations
        """
        if self.debug:
            self.logmt(1, 'Calculate the variance, variance reduction & flag bad stations')
        matrix_M /= -1.0e+20
        cnt = 0
        wsum = etot = var = dtot = dvar = 0
        svar = []
        sdpower = []
        for i in range(len(dict_s['T'])):
            dpower = 0
            e = 0
            trimrange = int(len(dict_s['T'][0])*float(self.trim_value))
            for j in range(trimrange):
                etmp  = dict_s['T'][i][j+timeshift[i]] 
                etmp -= matrix_M[0,0]*0.5*dict_g[0][i][j]*math.sin(2*az[i]) 
                etmp += matrix_M[1,1]*0.5*dict_g[0][i][j]*math.sin(2*az[i]) 
                etmp += matrix_M[0,1]*dict_g[0][i][j]*math.cos(2*az[i]) 
                etmp += matrix_M[0,2]*dict_g[1][i][j]*math.sin(az[i]) 
                etmp -= matrix_M[1,2]*dict_g[1][i][j]*math.cos(az[i])
                
                e += etmp*etmp
                
                etmp  = dict_s['R'][i][j+timeshift[i]] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[4][i][j] - 0.5*dict_g[2][i][j]*math.cos(2*az[i]) + dict_g[8][i][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[4][i][j] + 0.5*dict_g[2][i][j]*math.cos(2*az[i]) + dict_g[8][i][j]/3)
                etmp -= matrix_M[2,2]*dict_g[8][i][j]/3 
                etmp += matrix_M[0,1]*dict_g[2][i][j]*math.sin(2*az[i]) 
                etmp -= matrix_M[0,2]*dict_g[3][i][j]*math.cos(az[i])
                etmp -= matrix_M[1,2]*dict_g[3][i][j]*math.sin(az[i])
                
                e += etmp*etmp
                
                etmp  = dict_s['Z'][i][j+timeshift[i]] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[7][i][j] - 0.5*dict_g[5][i][j]*math.cos(2*az[i]) + dict_g[9][i][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[7][i][j] + 0.5*dict_g[5][i][j]*math.cos(2*az[i]) + dict_g[9][i][j]/3)
                etmp -= matrix_M[2,2]*dict_g[9][i][j]/3 
                etmp += matrix_M[0,1]*dict_g[5][i][j]*math.sin(2*az[i]) 
                etmp -= matrix_M[0,2]*dict_g[6][i][j]*math.cos(az[i])
                etmp -= matrix_M[1,2]*dict_g[6][i][j]*math.sin(az[i])
                
                e += etmp*etmp
                dpower += dict_s['T'][i][j+timeshift[i]]*dict_s['T'][i][j+timeshift[i]]
                dpower += dict_s['R'][i][j+timeshift[i]]*dict_s['R'][i][j+timeshift[i]]
                dpower += dict_s['Z'][i][j+timeshift[i]]*dict_s['Z'][i][j+timeshift[i]]
                cnt += 1
            
            wsum += list_W[i]
            etot += e
            var += list_W[i]*e
            dtot += dpower
            dvar += list_W[i]*dpower
            e /= dpower
            svar.append((1.0 - e)*100.0)
            sdpower.append(dpower)
        pvar = etot/(3*cnt - self.isoflag - 1.0)
        etot /= dtot
        pvred = (1-etot)*100
        var /= wsum
        dvar /= wsum
        var /= dvar
        var = (1-var)*100
        return pvar, pvred, var, svar, sdpower

    def quality_check(self, vr):
        """Check the quality of the result
        """
        if vr < 100:
            qlt = 4
        if vr < 80:
            qlt = 3
        if vr < 60:
            qlt = 2
        if vr < 40:
            qlt = 1
        if vr < 20:
            qlt = 0
        return qlt

    def write_results(self, orid, strike, dip, rake):
        """Write out results to database
        """
        if self.debug:
            self.logmt(1, 'Write out moment tensor for orid %s to database table %s.moment' % (orid, self.event_db))
        moment = antdb.dbopen(self.event_db, 'r+')
        try:
            moment.lookup(table='moment')
        except Exception, e:
            self.logmt(3, 'write_results(): ERROR in lookup: %s' % e)
        else:
            try:
                moment.addv('orid', int(orid),
                    'str1', math.degrees(strike[0]),
                    'dip1', math.degrees(dip[0]),
                    'rake1', math.degrees(rake[0]),
                    'str2', math.degrees(strike[1]),
                    'dip2', math.degrees(dip[1]),
                    'rake2', math.degrees(rake[1]))
            except Exception, e:
                self.logmt(3, 'write_results(): ERROR in adding row: %s' % e)
            else:
                self.logmt(1, 'write_results(): Successfully wrote out moment tensor to database')
        return

def main():
    """Get the moment tensor solution
    for a particular origin
    """
    pf, orid, verbose, debug = configure()
    my_mt = MomentTensor(pf, orid, verbose, debug)
    my_mt.parse_pf()
    evdbptr, evparams = my_mt.get_view_from_db()
    stachan_traces = my_mt.get_chan_data(evdbptr)

    """
    Declare empty matrices (for data and Green's functions)

    Currently both s and g are populated from pre-existing
    data files for the test case (method name + '_hard_coded()'). 
    We need to make these database driven
    """
    s = defaultdict(lambda: defaultdict(defaultdict))
    g = defaultdict(lambda: defaultdict(defaultdict))
    # s = my_mt.construct_data_matrix(stachan_traces, s)
    s = my_mt.construct_data_matrix_hard_coded(stachan_traces, s)
    if len(s) != 0:
        print 'Data matrix S created --> %s stations used' % len(s)

    # Read the Green's function from database based on distance and depth.....
    # greens_funcs = my_mt.get_greens_functions(evdbptr)

    # RLN (2011-07-28): green_funcs dict is not currently 
    # used anywhere, so why is it here? Possibly relates to 
    # the commented out dict above?
    green_funcs  = {}

    g = my_mt.get_greens_functions_hard_coded(g)
    if len(g[1]) != len(s['T']):
        print "Number of Green's functions (%d) does not match the number of stations (%d)" % (len(g[1]), len(s['T']))
    else:
        print 'Matrix S (data) and G (synthetics) created'

    timeshift = my_mt.get_time_shift(s, g)
    if len(timeshift) != len(s['T']):
        print 'Error in correlating synthetic and measured data'
    else:
        print 'Timeshift between data and synthetics computed for %s stations' % len(timeshift)

    # RLN (2011-07-28): Where does this come from without explanation? Dregers code?
    az = [10, 40, 50]

    # RLN (2011-07-28): What is W?
    W = []
    for i in range(len(s['T'])):
        W.append(1.0)

    AIV = defaultdict(dict)
    B = defaultdict(dict) 

    AIV, B = my_mt.matrix_AIV_B(s, g, az, timeshift)
    M = my_mt.det_solution_vec(AIV, B)

    gfscale = -1.0e+20
    gfscale = -1.0
    strike = []
    dip = []
    rake = []

    m0, Mw, strike, dip, rake, pcdc, pcclvd, pciso = my_mt.decompose_moment_tensor(M)

    svar = []
    sdpower = []

    E, VR, VAR, svar, sdpower = my_mt.fitcheck(g, s, W, M, m0, timeshift, az)

    qlt = my_mt.quality_check(VR)

    # Append results to moment database-table. If db.moment does not exists, create it
    evdbptr.close()
    my_mt.write_results(orid, strike, dip, rake)

    print 'M0      = %s' % m0
    print 'Mw      = %s' % Mw 
    print 'Strike  = %s' % math.degrees(strike[0])
    print 'Rake    = %s' % math.degrees(rake[0])
    print 'Dip     = %s' % math.degrees(dip[0])
    print 'Strike2 = %s' % math.degrees(strike[1])
    print 'Rake2   = %s' % math.degrees(rake[1])
    print 'Dip2    = %s' % math.degrees(dip[1])
    print 'Pdc     = %s' % pcdc
    print 'Pclvd   = %s' % pcclvd
    print 'Piso    = %s' % pciso
    for j in range(len(svar)):
        print 'Stat(%d) = %s  %s' % (j,svar[j],sdpower[j])
    print 'Var     = %s' % E
    print 'VR      = %s (UNWEIGHTED)' % VR
    print 'VR      = %s (WEIGHTED)' % VR
    print 'Var/Pdc = ', E/pcdc
    print 'Quality = %s' % qlt

    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)

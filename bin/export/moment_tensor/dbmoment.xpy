"""
dbmoment.py

@authors  Gert-Jan van den Hazel <hazelvd@knmi.nl>
          Juan Reyes <jreyes1108@gmail.com>
          Rob Newman <robertlnewman@gmail.com>
          Matt Koes <mattkoes@uvic.ca>
@notes    Look for initialized comments - they define questions to be answered
          e.g. # !!! BUG: Why this hard-coded value here? RLN (2011-07-28)
               # !!! FIX: This is a hacked solution. RLN (2011-07-28)
               # ??? Why is this here? RLN (2011-07-28)
               # !!! NOTE: An explanation. RLN (2011-07-28)
               (http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html#docstrings-comments)
"""

from optparse import OptionParser
from pprint import pprint
import re
import math as math
from datetime import datetime
from collections import defaultdict
from time import gmtime, time
from PIL import Image, ImageDraw, ImageFont

# ANTELOPE
import antelope.stock as stock
import antelope.datascope as antdb

# NUMPY
try:
    import numpy as np
except ImportError:
    print "Import Error: Do you have NumPy installed correctly?"

# MATPLOTLIB
try:
    import matplotlib as mpl
except ImportError:
    print "Import Error: Do you have MatplotLib installed correctly?"
else:
    mpl.use('tkagg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# PYLAB
"""
try:
    from pylab import *
except ImportError:
    print "Import Error: Do you have Pylab installed correctly?"
"""
# FOR MATTS LEAST SQRS ANALYSIS
import scipy.linalg as linalg

# FKRPROG
import moment_tensor.fkrprog as fkr

# OBSPY
try:
    from obspy.imaging.beachball import Beachball
except ImportError:
    print "Import Error: Do you have ObsPy installed correctly?"

# {{{ Global functions
def configure():
    """Function to configure parameters 
    from command-line and the pf-file
    Check command line options
    Use Antelope's built-in PFPATH to
    determine the paths to search for the
    parameter file
    """
    usage = "Usage: dbmoment [-vV] [-p pfname] orid"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="verbose output")
    parser.add_option("-V", "--debug", action="store_true", dest="debug", default=False, help="debug application")
    parser.add_option("-p", "--pf", action="store", dest="pf", type="string", help="parameter file path")
    (options, args) = parser.parse_args()
    if options.verbose:
        verbose = True
    else:
        verbose = False
    if options.debug:
        debug = True
    else:
        debug = False
    if not options.pf:
        pfs_tuple = list(stock.pffiles('dbmoment'))
        pfs_tuple.reverse() # Reverse order to search local pf dir first
        for p in pfs_tuple:
            if os.path.isfile(p):
                pfname = p
                break
        print "Used PFPATH to determine which parameter file to use and found '%s'" % pfname
    else:
        if not os.path.isfile(options.pf):
            print "Command line defined parameter file '%s' does not exist. Exiting" % options.pf
            sys.exit(-1)
        else:
            pfname = options.pf
    if len(args) != 1:
        print usage;
        sys.exit(-1)
    else:
        orid = args[0]
    return pfname, orid, verbose, debug

def determine_verbosity(verbose=False, debug=False):
    """Determine the verbosity"""
    verbosity = 0
    if verbose:
        verbosity = 1
    if debug:
        verbosity = 2
    return verbosity

def parse_pf(pfname, verbosity=0):
    """Parse the parameter file
    and assign to vars that will be
    used throughout the script
    """
    if verbosity > 1:
        print 'Parse parameter file %s' % pfname
    try:
        chan_to_use = stock.pfget_string(pfname, 'chan_to_use')
    except:
        chan_to_use = 'LH.*'
    trim_value = stock.pfget_string(pfname, 'trim_value')
    if stock.pfget_string(pfname, 'isoflag') == 1: 
        isoflag = 6
    else:
        isoflag = 5
    model_name = stock.pfget_string(pfname, 'model_name')
    model_type = stock.pfget_string(pfname, 'model_type')
    statmax = stock.pfget_int(pfname,'statmax')
    clip_values = stock.pfget_arr(pfname, 'clip_values')
    use_inc = stock.pfget_string(pfname, 'use_inc')
    event_db = stock.pfget_string(pfname,'event_db')
    try:
        wave_db = stock.pfget_string(self.pfname, 'wave_db')
    except:
        wave_db = False
    green_db = stock.pfget_string(pfname, 'green_db')
    try:
        resp_db = stock.pfget_string(pfname, 'resp_db')
    except:
        resp_db = False
    mag_filters = stock.pfget_arr(pfname, 'mag_filters')
    mt_images_dir = stock.pfget_string(pfname, 'mt_images_dir')
    ttfont = stock.pfget_string(pfname, 'ttfont')
    obspy_beachball = stock.pfget_arr(pfname, 'obspy_beachball')
    if stock.pfget_string(pfname, 'distance_weighting') == 'on':
        distance_weighting = True
    else:
        distance_weighting = False
    if not wave_db:
        wave_db = event_db
    if not resp_db:
        resp_db = event_db

    return (chan_to_use, trim_value, isoflag, model_name,
            model_type, statmax, clip_values, use_inc, event_db, 
            wave_db, green_db, resp_db, 
            mag_filters, mt_images_dir, ttfont, obspy_beachball, 
            distance_weighting)

def logmt(flag, message):
    """Function to handle log output, which depends 
    on the verbosity setting. Prints messages with a 
    timestamp. Whenever an error occurs it also exits.
    """
    curtime = stock.strtime( time() )
    if not flag:
        return
    elif flag < 3:
        print curtime, message
    else:
        print curtime, 'ERROR:',message,'--> exiting'
        sys.exit(-1)

# }}}

class MomentTensor():
    """Class for building moment tensors
    and doing the inversion
    """
    # {{{

    def __init__(self, distance_weighting, isoflag, trim_value, mt_images_dir, verbosity=0):
        """Initialize"""
        self.distance_weighting = distance_weighting
        self.isoflag = isoflag
        self.trim_value = trim_value
        self.mt_images_dir = mt_images_dir
        self.verbosity = verbosity

    def construct_data_matrix(self, stachan_traces):
        """Construct data matrix from dictionary 
        containing trace_objects per sta_chan
        Returns rearranged 3D matrix, 
        containing all data per channel
        """
        this_dict = defaultdict(lambda: defaultdict(defaultdict))

        if self.verbosity > 0:
            logmt(1, 'Constructing data matrix from trace objects')
        numsta = 0
        numchan = 0
        for stachan, tr in sorted(stachan_traces.items()):
            sta, chan = stachan.split('_')

            if self.distance_weighting == True:
                logmt(1, 'Apply distance weighting')

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

    def plot_cross_cor(self, a, b, xcor, components):
        """Plot the cross 
        correlation for 
        debugging purposes
        """
        stacode, component, shift, maxval = components
        title_string = "Station %s. Filtered data for %s. Shift = %s" % (stacode, component, shift)
        fig = plt.figure()
        fig.subplots_adjust(left=0.04, right=0.96, hspace=0.3)
        ax = fig.add_subplot(311)
        ax.set_title(title_string)
        ax.plot(self.normalize(a), 'b-')
        ax.plot(self.normalize(b), 'g-')
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.legend(('Station', 'Synthetic'), 'upper right', prop={'size':'10'})
        # Print shifted time series
        bx = fig.add_subplot(312)
        bx.set_title('Shifted time series')
        if shift < 0:
            try:
                [a.insert(0, 0) for x in range(abs(shift))]
            except Exception, e:
                logmt(1, "   - plot_cross_cor(): Could not insert values at start of real data. Error: %s" % e)
                print a
        elif shift > 0:
            try:
                [b.insert(0, 0) for x in range(abs(shift))]
            except Exception, e:
                logmt(1, "   - plot_cross_cor(): Could not insert values at start of synthetic data. Error: %s" % e)
                print b
        bx.plot(self.normalize(a), 'r-')
        bx.plot(self.normalize(b), 'c-')
        plt.setp(bx.get_yticklabels(), visible=False)
        bx.legend(('Station', 'Synthetic'), 'upper right', prop={'size':'10'})
        # Print xcor
        cx = fig.add_subplot(313)
        cx.set_title('Cross correlation results (max cross-corr val = %s)' % maxval)
        cx.plot(xcor)
        plt.setp(cx.get_yticklabels(), visible=False)
        cx.set_ylabel('X-cor.')
        cx.set_xlabel('Point shift (in samples)')
        if self.verbosity > 1:
            plt.show()
        cross_corr_outfile = '%s/%s_%s.png' % (self.mt_images_dir, stacode, component)
        plt.savefig(cross_corr_outfile)
        return True
 
    def cross_cor(self, a, b, component=False, stacode=False):
        """Calculate cross correlation 
        between data and accompanying 
        Green's function. Returns the 
        max_cor and the shift
        """
        xcor = np.correlate(a, b, 'full')
        maxval = np.amax(xcor)
        maxshift = np.argmax(xcor)
        shift = maxshift - (len(xcor)-1)/2
        # If shift < 0, STATION data must move to the right
        # If shift > 0, STATION data must move to the left
        if self.verbosity > 1:
            logmt(1, '  - cross_cor(): Component: %s, max val: %s, timeshift (in samples): %s' % (component, maxval, shift))
        if self.verbosity > 0:
            self.plot_cross_cor(a, b, xcor, [stacode, component, shift, maxval])
        return maxval, shift

    def normalize(self, data_as_list):
        """Determine the 
        normalization factor
        for all the data
        """
        try:
            normalizer = max([abs(x) for x in data_as_list])
        except Exception, e:
            logmt(1, "  - Exception encountered: %s" % e)
        else:
            if self.verbosity > 1:
                logmt(1, "  - Normalization factor: %s" % normalizer)
            return [x / normalizer for x in data_as_list]

    def get_time_shift(self, data_dict, greens_dict, stacode, size):
        """Get time shift for each station
        using cross correlation of station 
        data points and the Green's Function
        data points. Make sure both matrices
        are the same size for the comparison.
        Return a list of tuples (stacode, timeshift).

        !!! NOTE: They don't have to be the same
                  size for a successful cross 
                  correlation
        """
        if self.verbosity > 0:
            logmt(1, '  - Get the time shift for all components & return the max cross-cor and shift in samples')
        '''
        !!! NOTE: Loop over each stations data points.
                  Remember that the Green's Function 
                  matrix is just two dimensional - not 
                  replicating the number of stations
                  (3D matrix). RLN (2011-12-02)
        '''
        size = len(greens_dict['XDS'])
        max_shift = 0
        max_xcor = 0
        for c in greens_dict:
            if c == 'REX' or c == 'ZEX':
                continue
            else:
                if c[:1] == 'T':
                    data_key = 'T'
                elif c[:1] == 'X':
                    data_key = 'R'
                elif c[:1] == 'Z':
                    data_key = 'Z'
                this_xcor, this_shift = self.cross_cor(data_dict[data_key][0:size], greens_dict[c][0:size], c, stacode)
                if this_xcor > max_xcor:
                    max_xcor = this_xcor
                    max_shift = this_shift
        if self.verbosity > 0:
            logmt(1, '  - Cross-correlation value: %s, timeshift value: %s' % (max_xcor, max_shift))
        return max_xcor, max_shift

    def convertdict2list(self, A_dict, b_dict):
        ''' 
        Convert a dictionary input to a list ouput
        IN:
           - A_dict = dictionary to be converted to a MxN list (M columns, N rows)
           - b_dict = dictionary to be converted to a Mx1 list (M columns, 1 row )
        OUT: 
           - list A
           - list b
           
        Comments: 
        For the dbmoment inversion script A_dict = AJ (the greens function matrix)
        and b_dict = the waveform data

        A_dict is a dictionary
        b_dict is the station data and is a list (ordered) and looks like:
        b_dict = [ 
                stacode: {
                    'T': [ list_of_data_points ],
                    'R': [ list_of_data_points ],
                    'Z': [ list_of_data_points ]
                    },
                stacode: {
                    'T': [ list_of_data_points ],
                    'R': [ list_of_data_points ],
                    'Z': [ list_of_data_points ]
            ]
        '''
        # A list
        A = []
        tmp = [v for k,v in A_dict.items()]
        for i in range(len(tmp)):
            tmp2 = tmp[i]
            A.append([subv for subv in tmp2.itervalues()])

        # Rotate the matrix. Do we need to do this?
        # Adash = np.array(A)
        # newA = Adash.T.tolist()
        # newA = Adash.T.tolist()

        # print "\n\nA:"
        # print A
        # print "\n\nAdash:"
        # print Adash
        # print "\n\nnewA:"
        # print newA

        # print 'The converted list (from dicitonary) is A =',A
        # b list
        b = []
        for v in b_dict:
            for subk, subv in v.items():
                b.extend(subv)
        # 
        # for i in range(len(tmp2)):
        #     tmp2 = tmp[i]
        #     b += [v for k,v in tmp2.items()]
        
        # print 'The converted list (from dicitonary) is b =',b
        # return(Adash,b)
        return(A, b)
	
    def matrix_AIV_B_SCIPY(self, dict_s, dict_g, ev2sta, size):
        """MATTS INVERSION USING SCIPY
        """
        if self.verbosity > 0:
            logmt(1, 'INVERSION ROUTINE: Construct matrices AIV and B using the data and Greens function matrices')

        AJ = defaultdict(dict)
        trim = 0
        cnt1 = cnt2 = cnt3 = 0

        if self.verbosity > 0:
            # logmt(1, ' - Timeshift: %s' % timeshift)
            # logmt(1, ' - Azimuthal distances tuple: %s' % list_az)
            logmt(1, ' - Number of stations used in calculation: %s' % len(dict_s))

        # Iterate over number of stations
        for i in range(len(ev2sta)):
            if self.verbosity > 0:
                logmt(1, '  - Working on inversion for station (%s)' % ev2sta[i][0])
            # Allocate the dictionary space for each section
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            # Convert azimuth to radians
            this_azimuth = ev2sta[i][1] * math.pi/180
            '''
            !!! NOTE: len(dict_s['T'][0]) is the number of 
                      the datapoints for the first station. 
                      Use var size instead because of timeshift.
                      RLN (2011-12-06)

            '''
            for j in range(size - trim):
                # Mxx term
                AJ[0][cnt1] =  math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Myy term
                AJ[1][cnt1] = -math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Mxy term
                AJ[2][cnt1] = -math.cos(2*this_azimuth)*dict_g[i]['TSS'][j]
                AJ[2][cnt2] = -math.sin(2*this_azimuth)*dict_g[i]['XSS'][j]
                AJ[2][cnt3] = -math.sin(2*this_azimuth)*dict_g[i]['ZSS'][j]

                # Mxz term
                AJ[3][cnt1] = -math.sin(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[3][cnt2] =  math.cos(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[3][cnt3] =  math.cos(this_azimuth)*dict_g[i]['ZDS'][j]

                # Myz term
                AJ[4][cnt1] =  math.cos(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[4][cnt2] =  math.sin(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[4][cnt3] =  math.sin(this_azimuth)*dict_g[i]['ZDS'][j]

                # Vary the other values depending on isoflag value
                if self.isoflag == 5:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                if self.isoflag == 6:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Explosion-related values
                    AJ[5][cnt1] = 0.0
                    AJ[5][cnt2] = (dict_g[i]['REX'][j])/3  - (dict_g[i]['XDD'][j])/3
                    AJ[5][cnt3] = (dict_g[i]['ZEX'][j])/3 - (dict_g[i]['ZDD'][j])/3
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1
        if self.verbosity > 0:
            logmt(1, ' - Created matrix AJ with length: %s' % len(AJ))
            logmt(1, ' - Final counts: cnt1=%s, cnt2=%s, cnt3=%s' % (cnt1, cnt2, cnt3))
            logmt(1, ' - Now apply the timeshift if needed')


        # APPLY MATTS CONVERSION ROUTINE
        listA, listB = self.convertdict2list(AJ, dict_s)
        # print "\n\nFIRST ELEMENT OF listA"
        # print listA[0]
        print "Length of listA = %s" % len(listA)
        print "Length of listA = %s" % len(listA[0])
        # print "\n\nlistB"
        # print listB
        print "Length of listB = %s" % len(listB)

        try:
            x = linalg.lstsq(listA, listB)
        except ValueError, e:
            logmt(3, '  - ValueError: %s' % e)
        else:
           print x
           x = x[0]
           print 'linalg.lstsq solution x =',x
           print 'the size of x is',x.shape		
        exit()

    def matrix_AIV_B(self, dict_s, dict_g, ev2sta, size):
        """INVERSION ROUTINE

        Construct matrices AIV and B using 
        the data and Green's function matrices
        Return AIV and B for further processing 
        This is to normalize the two signals
        """
        if self.verbosity > 0:
            logmt(1, 'INVERSION ROUTINE: Construct matrices AIV and B using the data and Greens function matrices')

        AJ = defaultdict(dict)
        trim = 0
        cnt1 = cnt2 = cnt3 = 0

        if self.verbosity > 0:
            # logmt(1, ' - Timeshift: %s' % timeshift)
            # logmt(1, ' - Azimuthal distances tuple: %s' % list_az)
            logmt(1, ' - Number of stations used in calculation: %s' % len(dict_s))

        # Iterate over number of stations
        for i in range(len(ev2sta)):
            if self.verbosity > 0:
                logmt(1, '  - Working on inversion for station (%s)' % ev2sta[i][0])
            # Allocate the dictionary space for each section
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            # Convert azimuth to radians
            this_azimuth = ev2sta[i][1] * math.pi/180
            '''
            !!! NOTE: len(dict_s['T'][0]) is the number of 
                      the datapoints for the first station. 
                      Use var size instead because of timeshift.
                      RLN (2011-12-06)

            '''
            for j in range(size - trim):
                # Mxx term
                AJ[0][cnt1] =  math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Myy term
                AJ[1][cnt1] = -math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Mxy term
                AJ[2][cnt1] = -math.cos(2*this_azimuth)*dict_g[i]['TSS'][j]
                AJ[2][cnt2] = -math.sin(2*this_azimuth)*dict_g[i]['XSS'][j]
                AJ[2][cnt3] = -math.sin(2*this_azimuth)*dict_g[i]['ZSS'][j]

                # Mxz term
                AJ[3][cnt1] = -math.sin(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[3][cnt2] =  math.cos(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[3][cnt3] =  math.cos(this_azimuth)*dict_g[i]['ZDS'][j]

                # Myz term
                AJ[4][cnt1] =  math.cos(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[4][cnt2] =  math.sin(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[4][cnt3] =  math.sin(this_azimuth)*dict_g[i]['ZDS'][j]

                # Vary the other values depending on isoflag value
                if self.isoflag == 5:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                if self.isoflag == 6:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Explosion-related values
                    AJ[5][cnt1] = 0.0
                    AJ[5][cnt2] = (dict_g[i]['REX'][j])/3  - (dict_g[i]['XDD'][j])/3
                    AJ[5][cnt3] = (dict_g[i]['ZEX'][j])/3 - (dict_g[i]['ZDD'][j])/3
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1
        if self.verbosity > 0:
            logmt(1, ' - Created matrix AJ with length: %s' % len(AJ))
            logmt(1, ' - Final counts: cnt1=%s, cnt2=%s, cnt3=%s' % (cnt1, cnt2, cnt3))
            logmt(1, ' - Now apply the timeshift if needed')

        '''
        !!! NOTE: Make placeholder AIV 
                  (or reset it like Dreger?)
                  RLN (2011-12-02)
        !!! NOTE: Attempt to replace this with scipy.linalg.lstsq
                  http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html#scipy.linalg.lstsq
                  RLN & MK (2012-02-07)
        '''

        AIV = defaultdict(dict) 
        for i in range(self.isoflag):
            for j in range(self.isoflag):
                AIV[i][j] = 0.0
        # Compute AtA
        for i in range(5):
            for j in range(5):
                for k in range(cnt3):
                    AIV[i][j] += AJ[i][k] * AJ[j][k]

        B = defaultdict(dict) 
        for i in range(self.isoflag):
            B[i][0] = 0.0

        cnt1 = cnt2 = cnt3 = 0
        tmp = defaultdict(dict) 

        # Iterate over the number of stations
        for i in range(len(ev2sta)):
            if self.verbosity > 0:
                logmt(1, '  - Applying timeshift to station (%s): Trim: %s, timeshift: %s' % (ev2sta[i][0], trim, ev2sta[i][3]))
            '''
            !!! FIX: Direct copy of Dregers code - don't need to do this
                     RLN (2011-12-08)
            Pre-populate dictionary of values
            |cnt1           |cnt2        |cnt3       |
            |*****----------|*****-------|*****------|
            '''
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            for j in range(size - trim):
                correction_iter = j + ev2sta[i][3]
                try:
                    tmp[cnt1] = dict_s[i]['T'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s with index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['T'])
                try:
                    tmp[cnt2] = dict_s[i]['R'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s with index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['R'])
                try:
                    tmp[cnt3] = dict_s[i]['Z'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s at index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['Z'])
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1

        # Calculate Righthand Side
        for i in range(self.isoflag):
            for j in range(cnt3):
                B[i][0] += AJ[i][j] * tmp[j]

        if len(AIV) != self.isoflag or len(AIV[0]) != self.isoflag:
            logmt(3, 'Matrix AIV has dimension [%s,%s] should be [%s,%s]' % (len(AIV), len(AIV), self.isoflag, self.isoflag))
        elif len(B) != self.isoflag:
            logmt(3, 'Matrix AIV has dimension [%s,i%s] should be [%s,1]' % (len(B), len(B[0]), self.isoflag))
        elif self.verbosity > 0:
            logmt(1, 'Matrices AIV and B created and have correct dimensions')

        return AIV, B

    def swap(self, a, b):
        """Remap list
        """
        tmp = a
        a = b
        b = tmp
        return(a, b)

    def dyadic(self, v, n1, n2, c):
        """Calculate the dyadic 
        matrix of eigenvectors v
        """
        if self.verbosity > 0:
            logmt(1, ' - Compute the dyadic matrix of vector v')
        tmp = np.matrix([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        for i in range(3):
            for j in range(3):
                tmp[i,j] = v[i, n1]*v[j, n2]*c
        return tmp

    def determine_solution_vector(self, dict_AIV, dict_B):
        """Determine the solution vector
        Solve the inversion problem, 
        returning the moment tensor
        From Dreger:
            Call Bobs MT decomposition routines
            The minus one is needed to map Helmbergers convention into Aki's
            Jost and Hermann (1989) state that AKI's convention is -1*LANGSTONS

        !!! NOTE: Matt thinks we can replace this with a library in SciPy
            RLN (2012-02-07)
        """
        if self.verbosity > 0:
            logmt(1, 'Solve the inversion problem and return the moment tensor')
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
                            logmt(3, 'determine_solution_vector(): ERROR... 1: GAUSSJ: Singular Matrix-1')
            ipiv[icol] += 1

            if not irow == icol:
                for l in range(self.isoflag):
                    (dict_AIV[irow][l], dict_AIV[icol][l]) = swap(dict_AIV[irow][l], dict_AIV[icol][l])
                for l in range(1):
                    (dict_B[irow][l], dict_B[icol][l]) = swap(dict_B[irow][l], dict_B[icol][l])

            if dict_AIV[icol][icol] == 0.0:
                logmt(3, 'determine_solution_vector(): ERROR... 2: GAUSSJ: Singular Matrix-2')
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
            logmt(3, 'Wrong dimension returned for matrix M after inversion')
        elif self.verbosity > 0:
            logmt(1, ' - Moment tensor matrix M[3,3] created')
        logmt(1, '\n\t\tMxx=%s\n\t\tMxy=%s\n\t\tMxz=%s\n\t\tMyy=%s\n\t\tMyz=%s\n\t\tMzz=%s\n' % (M[0,0], M[0,1], M[0,2], M[1,1], M[1,2], M[2,2]))
        return M

    def decompose_moment_tensor(self, matrix_M):
        """Decompose moment tensor into eigenvector/values.
        Calculate moment tensor parameters from eigenvalues and vectors. 
        Return M0, Mw, strike, slip, rake and precentages of present
        source characteristics
        From Dreger source: fmap_subs_linux.f
        Mxx = matrix_M[0, 0]
        Mxy = matrix_M[0, 1]
        Mxz = matrix_M[0, 2]
        Myy = matrix_M[1, 1]
        Myz = matrix_M[1, 2]
        Mzz = matrix_M[2, 2]
        """
        if self.verbosity > 0:
            logmt(1, 'Decompose moment tensor into eigenvector/values')
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

        if self.verbosity > 0:
            logmt(1, ' - Decomposition of moment tensor succesful!')

        return m0, Mw, strike, dip, slip, pcdc, pcclvd, pciso

    def fitcheck(self, dict_g, dict_s, matrix_M, m0, ev2sta, size):
        """Calculate the variance, 
        variance reduction
        and flag bad stations
        """
        if self.verbosity > 0:
            logmt(1, 'Calculate the variance, variance reduction & flag bad stations')
        matrix_M /= -1.0e+20
        cnt = 0
        wsum = etot = var = dtot = dvar = 0
        svar = []
        sdpower = []
        # Loop over the number of stations
        for i in range(len(ev2sta)):
            dpower = 0
            e = 0
            trimrange = int(int(size) * float(self.trim_value))
            for j in range(trimrange):
                correction_iter = j + ev2sta[i][3]
                etmp  = dict_s[i]['T'][correction_iter] 
                etmp -= matrix_M[0,0]*0.5*dict_g[i]['TSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp += matrix_M[1,1]*0.5*dict_g[i]['TSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp += matrix_M[0,1]*dict_g[i]['TSS'][j]*math.cos(2*ev2sta[i][1]) 
                etmp += matrix_M[0,2]*dict_g[i]['TDS'][j]*math.sin(ev2sta[i][1]) 
                etmp -= matrix_M[1,2]*dict_g[i]['TDS'][j]*math.cos(ev2sta[i][1])
                
                e += etmp*etmp
                
                etmp  = dict_s[i]['R'][correction_iter] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[i]['XDD'][j] - 0.5*dict_g[i]['XSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['REX'][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[i]['XDD'][j] + 0.5*dict_g[i]['XSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['REX'][j]/3)
                etmp -= matrix_M[2,2]*dict_g[i]['REX'][j]/3 
                etmp += matrix_M[0,1]*dict_g[i]['XSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp -= matrix_M[0,2]*dict_g[i]['XDS'][j]*math.cos(ev2sta[i][1])
                etmp -= matrix_M[1,2]*dict_g[i]['XDS'][j]*math.sin(ev2sta[i][1])
                
                e += etmp*etmp
                
                etmp  = dict_s[i]['Z'][correction_iter] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[i]['ZDD'][j] - 0.5*dict_g[i]['ZSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['ZEX'][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[i]['ZDD'][j] + 0.5*dict_g[i]['ZSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['ZEX'][j]/3)
                etmp -= matrix_M[2,2]*dict_g[i]['ZEX'][j]/3 
                etmp += matrix_M[0,1]*dict_g[i]['ZSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp -= matrix_M[0,2]*dict_g[i]['ZDS'][j]*math.cos(ev2sta[i][1])
                etmp -= matrix_M[1,2]*dict_g[i]['ZDS'][j]*math.sin(ev2sta[i][1])
                
                e += etmp*etmp
                dpower += dict_s[i]['T'][correction_iter]*dict_s[i]['T'][correction_iter]
                dpower += dict_s[i]['R'][correction_iter]*dict_s[i]['R'][correction_iter]
                dpower += dict_s[i]['Z'][correction_iter]*dict_s[i]['Z'][correction_iter]
                cnt += 1
            
            # wsum += list_W[i]
            wsum += 1
            etot += e
            # var += list_W[i]*e
            var += 1*e
            dtot += dpower
            # dvar += list_W[i]*dpower
            dvar += 1*dpower
            e /= dpower
            svar.append((ev2sta[i][0], (1.0 - e)*100.0))
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
    # }}}

class Event():
    """Class for extracting per event 
    data from the database
    and creating the plots
    """
    # {{{

    def __init__(self, orid, event_db, mt_images_dir, verbosity=0):
        """Initialize"""
        self.orid = orid
        self.event_db = event_db
        self.mt_images_dir = mt_images_dir
        self.verbosity = verbosity

        if not os.path.exists(self.mt_images_dir):
            if self.verbosity > 0:
                logmt(1, 'Images dir (%s) does not exist. Try to create...' % self.mt_images_dir)
            try:
                os.makedirs(self.mt_images_dir, 0775)
            except Exception, e:
                logmt(3, 'Moment tensor images dir (%s) does not exist and cannot be created! Exception: %s' % (self.mt_images_dir, e))

    def extract_data(self, mag_filters):
        """Open event database and get 
        event data for given origin-id
        Returns event parameters or 
        exit when an error occurs
        """
        if self.verbosity > 0:
            logmt(1, 'Get canned view from database')
        if not os.path.isfile(self.event_db):
            logmt(3, 'Database (%s) does not exist' % self.event_db)
        try:
            evdb = antdb.dbopen(self.event_db, 'r')
        except Exception, e:
            logmt(3, 'Could not open database (%s). Exception: %s' % (self.event_db, e))
        else:
            evdb.lookup(table='origin')
            if self.verbosity > 0:
                logmt(1, 'Processing origin #: %s' % self.orid)
            evdb.subset('orid == %s' % self.orid)
            if evdb.nrecs() == 0:
                logmt(3, 'Orid (%s) does not exist in origin table' % self.orid)
            evdb.join('netmag', outer=True)
            if evdb.nrecs() == 0:
                logmt(3, 'Could not join netmag table for orid (%s)' % self.orid)
            elif evdb.nrecs() > 1:
                logmt(1, 'More than one entry (%s) in the netmag table for orid (%s)' % (evdb.nrecs(), self.orid))
                evdb[3] = 0 # Force to first record
                magid = evdb.getv('magid') 
                logmt(1, 'Subset to get unique entry. Using magid (%s).' % magid)
                evdb.subset('magid == %d' % magid)
            evdb[3] = 0 # Should only be one record now
            evparams = {}
            for f in evdb.query('dbTABLE_FIELDS'):
                try: 
                    evparams[f] = evdb.getv(f)[0]
                except Exception, e:
                    logmt(3, 'Could not find field (%s) in join of origin & netmag tables for orid %s. Exception: %s' % (f, self.orid, e))
            evdb.join('assoc')
            evdb.join('arrival')
            evdb.join('site')
            evdb.sort('delta')
            evdb.subset('iphase=~/.*P.*|.*p.*/')
            '''
            !!! NOTE: This is regional moment
                      tensor. Will fail if stations
                      are closer than 1 deg (clipping)
                      or further away than 10 deg
                      (teleseismic) to the event.
                      RLN (2011-12-13)
            '''
            evdb.subset('delta >= 1 && delta < 10') # REGIONAL moment tensor,
            if evdb.nrecs() == 0:
                logmt(3, 'No arrivals for selected origin %s' % self.orid)
            elif self.verbosity > 0:
                logmt(1, 'There are %d records between 1 & 10 deg that match this origin arrival and iphase' % evdb.nrecs())
            mag_filter_list = []
            for line in mag_filters:
                splitline = line.split()
                if len(splitline) == 3:
                    lower_magnitude, timepad, filter = splitline
                    upper_magnitude = False
                elif len(splitline) == 4:
                    lower_magnitude, upper_magnitude, timepad, filter = splitline
                else:
                    logmt(3, 'Filter format not recognized')
                mag_filter_list.append({'lowermag':lower_magnitude, 'uppermag':upper_magnitude, 'timepad':timepad, 'filter':filter})
            for mfl in mag_filter_list:
                if float(evparams['magnitude']) >= float(mfl['lowermag']) and float(evparams['magnitude']) < float(mfl['uppermag']):
                    filter_string = (mfl['filter']).replace('_', ' ')
                    filter_timepad = int(mfl['timepad'])
            if not filter_string and not filter_timepad:
                logmt(1, 'Magnitude %s not within magnitude bounds defined in the parameter file:' % evparams['magnitude'])
                print mag_filter_list
                logmt(3, '')
            return evdb, evparams, filter_string, filter_timepad

    def get_stations_and_orientations(self, dbptr):
        """Return a list of 
        all stations that
        recorded the event.
        Split into eight 
        sections:

                 0
             NNW | NNE
           WNW \ | / ENE
        270 ----------- 90
           WSW / | \ ESE
             SSW | SSE  
                180
        """
        sta_list = {'NNE':[], 'ENE':[], 'ESE':[], 'SSE':[], 'SSW':[], 'WSW':[], 'WNW':[], 'NNW':[]}
        for i in range(dbptr.nrecs()):
            dbptr[3] = i
            sta, esaz, depth, at, ev_lat, ev_lon, site_lat, site_lon = dbptr.getv('sta', 'esaz', 'depth', 'origin.time', 'lat', 'lon', 'site.lat', 'site.lon')
            distance_deg = dbptr.ex_eval('distance(%s, %s, %s, %s)' % (ev_lat, ev_lon, site_lat, site_lon))
            distance_km = dbptr.ex_eval('deg2km(%s)' % (distance_deg))
            if esaz >= 0 and esaz < 45:
                sta_list['NNE'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 45 and esaz < 90:
                sta_list['ENE'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 90 and esaz < 135:
                sta_list['ESE'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 135 and esaz < 180:
                sta_list['SSE'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 180 and esaz < 225:
                sta_list['SSW'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 225 and esaz < 270:
                sta_list['WSW'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 270 and esaz < 315:
                sta_list['WNW'].append((sta, esaz, depth, distance_km, at))
            elif esaz >= 315 and esaz <= 360:
                sta_list['NNW'].append((sta, esaz, depth, distance_km, at))
        return sta_list

    def get_chan_data(self, wave_db, chan_to_use, esaz, at, filter_string, filter_timepad, stacode, size, clip_values):
        """Opens waveform database and returns
        trace objects based on sta_chan. Applies
        calibration, splice, filter, decimation
        and rotation of channels
        !!! NOTE: To avoid filter transient problems
                  we need to ensure we get a large enough
                  record section prior to the P-arrival time,
                  on the order of minutes (filter_timepad in secs). 
                  Once we filter, then subset the trace object 
                  to the window around the arrival time.
                  RLN (2012-02-02)
        """
        clip_max = float(eval(clip_values['max']))
        clip_min = float(eval(clip_values['min']))

        if self.verbosity > 0:
            logmt(1, '  - Get channel data (trace objects) for station (%s)' % stacode)
        '''
        !!! NOTE: Start time and endtime need to be calculated using ptime 
        Right now we are just getting size/2 either side of first arrival
        Instead use arrival.time - should be the same thing. RLN (2011-11-29)
        '''
        stachan_trace = {} # Holder dict for data

        '''
        !!! NOTE: Need to get a fresh view of 
                  the wfdisc as subsetted for P 
                  arrivals which are only done 
                  on the BH.* channels.
                  RLN (2011-11-29)
        '''
        try:
            wvdb = antdb.dbopen(wave_db, 'r')
        except Exception, e:
            logmt(3, 'Could not open waveform database %s. Exception: %s' % (wave_db, e))
            return False
        else:
            wvdb.lookup(table='wfdisc')
            wvdb.subset('sta =~ /%s/' % stacode)
            wvdb[3] = 0
            '''
            !!!! NOTE: Want maximium x3 the size 
                       to account for timeshifts
                       and add filter_timepad before
                       RLN (2012-02-02)
            st              p origin.time            et
            | <---------------->*<----><----><---->|
                filter_timepad
            '''
            st = at - filter_timepad
            et = at + (size * 3) 
            resample = 0
            vang = 0
            if self.verbosity > 0:
                logmt(1, '    - Retrieve %s data with time range (%s, %s)' % (stacode, st, et))
            wvdb.subset('sta=~/^%s$/ && chan=~/%s/' % (stacode, chan_to_use))
            trace = wvdb.load_css(st, et)
            trace.apply_calib()
            trace.splice() # Join all segments together
            trace.filter(filter_string)
            rotchan = ('R', 'T', 'Z')
            trace.rotate(esaz, vang, rotchan)
            trace.subset('chan =~ /R|T|Z/')
            if self.verbosity > 1:
                logmt(1, '    - Number of traces for %s: %s' % (stacode, trace.nrecs()))
            for j in range(trace.nrecs()):
                trace[3] = j
                sta, chan, ns, sr = trace.getv('sta', 'chan', 'nsamp', 'samprate')
                stachan = '%s_%s' % (sta, chan)
                ns_tra = 0
                ns_req = 0
                filtered_samples = list(trace.data())
                use_data = 0

                # {{{ Test for clipped data
                '''
                for i in filtered_samples:
                    if i > clip_max:
                        use_data = 1
                    elif i < clip_min:
                        use_data = -1
                '''
                # }}}

                # if use_data == 1:
                #     logmt(1, '    - %s: One or more samples > clip_values[max] (%s). IGNORE.' % (stachan, clip_values['max']))
                # elif use_data == -1:
                #     logmt(1, '    - %s: One or more samples < clip_values[min] (%s). IGNORE.' % (stachan, clip_values['min']))
                # else:
                #    '''
                #    !!! NOTE: Get three times the size to 
                #              allow for a positive time shift 
                #              RLN (2011-12-06)
                #    '''
                max_size = int(size * 3) # Only works if 1Hz (1 sample per second, LH.*)
                # max_size = int(size) # Only works if 1Hz (1 sample per second, LH.*)
                if self.verbosity > 1 and len(filtered_samples) > max_size:
                    logmt(1, '  - Sample size for %s: %s. Max: %s. Trimming.' % (chan, len(filtered_samples), max_size))
                stachan_trace[chan] = filtered_samples[filter_timepad:max_size]  
                if self.verbosity > 0:
                    logmt(1, '    - Trace extracted for %s samples:[%s], wanted:[%s]' % (chan, len(stachan_trace[chan]), max_size))
            trace.trdestroy()
            wvdb.free()
            wvdb.close()
        return stachan_trace

    def update_moment_tbl(self, strike, dip, rake):
        """Write out results to 
        database moment table
        """
        if self.verbosity > 0:
            logmt(1, 'Write out moment tensor for orid %s to database table %s.moment' % (self.orid, self.event_db))
            logmt(1, ' - MT for orid %s strike => %s' % (self.orid, strike))
            logmt(1, ' - MT for orid %s dip => %s' % (self.orid, dip))
            logmt(1, ' - MT for orid %s rake => %s' % (self.orid, rake))
        moment_dbptr = antdb.dbopen(self.event_db, 'r+')
        try:
            moment_dbptr.lookup(table='moment')
        except Exception, e:
            logmt(3, 'update_moment_tbl error: Error in lookup: %s' % e)
        else:
            orid_subset = antdb.dbsubset(moment_dbptr, 'orid == %s' % self.orid)
            if orid_subset.query('dbRECORD_COUNT') == 0:
                logmt(1, 'Adding new moment tensor to moment table with orid (%s)' % self.orid)
                try:
                    moment_dbptr.addv(
                        'orid', int(self.orid),
                        'str1', math.degrees(strike[0]),
                        'dip1', math.degrees(dip[0]),
                        'rake1', math.degrees(rake[0]),
                        'str2', math.degrees(strike[1]),
                        'dip2', math.degrees(dip[1]),
                        'rake2', math.degrees(rake[1])
                    )
                except Exception, e:
                    logmt(3, 'Adding record to moment table unknown error: %s' % e)
                else:
                    logmt(1, 'Successfully added record with orid (%s) to moment table' % self.orid)
            else:
                logmt(1, 'Updating moment tensor to moment table with orid (%s). Deleting current record and rewriting.' % self.orid)
                for i in range(moment_dbptr.query('dbRECORD_COUNT')):
                    moment_dbptr[3] = i
                    try:
                        moment_dbptr.putv(
                            'orid', int(self.orid),
                            'str1', math.degrees(strike[0]),
                            'dip1', math.degrees(dip[0]),
                            'rake1', math.degrees(rake[0]),
                            'str2', math.degrees(strike[1]),
                            'dip2', math.degrees(dip[1]),
                            'rake2', math.degrees(rake[1])
                        )
                    except Exception, e:
                        logmt(3, 'Update record in moment table unknown error: %s' % e)
                    else:
                        logmt(1, 'Successfully updated record with orid (%s) in moment table' % self.orid)
            orid_subset.free()
            moment_dbptr.free()
            moment_dbptr.close()
            return True

    def create_focal_mechanism(self, obspy_beachball, matrix_M=False, strike=False, dip=False, rake=False):
        """Write out focal mechanism
        to images directory
        """
        if self.verbosity > 0:
            logmt(1, 'Writing file to images dir (%s)' % self.mt_images_dir)

        '''
        !!! NOTE: Focal mechanism can be 
                  defined in ObsPy in two ways:
                  focal_mechanism = [strike[0], dip[0], rake[0]]
                      - or -
                  focal_mechanism = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
                  Allow both. RLN (2011-12-02)
        '''
        if len(matrix_M) < 1:
            focal_mechanism = [strike[0], dip[0], rake[0]]
        else:
            focal_mechanism = [ matrix_M[0, 0], 
                                matrix_M[1, 1],  
                                matrix_M[2, 2], 
                                matrix_M[0, 1], 
                                matrix_M[0, 2],
                                matrix_M[1, 2]
                              ]

        if self.verbosity > 0:
            logmt(1, 'Try to plot focal mechanism: %s' % focal_mechanism)

        beachball_vals = {}
        '''
        !!! NOTE: From the ObsPy website. 
                  This will need to be 
                  updated if the library 
                  updates. RLN (2011-12-06)
        '''
        beachball_defaults = { 
            'size': 200, 
            'linewidth': 2, 
            'facecolor': 'b', 
            'edgecolor': 'k', 
            'bgcolor': 'w', 
            'alpha': 1.0, 
            'xy': (0, 0),
            'width': 200, 
            'format': None, 
            'nofill': False, 
            'fig': None
        }
        if self.verbosity > 1:
            logmt(1, 'Beachball(): defaults %s' % beachball_defaults)

        # Test for defaults in the pf
        for k, v in beachball_defaults.iteritems():
            if not obspy_beachball[k]:
                if self.verbosity > 1:
                    logmt(1, 'write_results(): Beachball(): Setting default for %s' % k)
                beachball_vals[k] = beachball_defaults[k]
            else:
                if self.verbosity > 1:
                    logmt(1, 'write_results(): Beachball(): Using pf defined value for %s' % k)
                beachball_vals[k] = obspy_beachball[k]
            if self.verbosity > 1:
                logmt(1, 'write_results(): Beachball(): Arg: %s, Val: %s' % (k, beachball_vals[k]))

        my_outfile = '%s_focal_mechanism.%s' % (self.orid, beachball_vals['format'])
        my_outpath = '%s/%s' % (self.mt_images_dir, my_outfile)

        try:
            Beachball(focal_mechanism, 
                size = beachball_vals['size'],
                linewidth = beachball_vals['linewidth'], 
                facecolor = beachball_vals['facecolor'],
                edgecolor = beachball_vals['edgecolor'],
                bgcolor = beachball_vals['bgcolor'],
                alpha = beachball_vals['alpha'],
                xy = beachball_vals['xy'],
                width = beachball_vals['width'],
                outfile = my_outpath,
                format = beachball_vals['format'],
                nofill = beachball_vals['nofill'],
                fig = beachball_vals['fig']
            )
        except Exception,e:
            logmt(3, 'Error creating Beachball() %s: %s' % (Exception, e))
            return False
        else:
            if self.verbosity > 0:
                logmt(1, 'Successfully created focal mechanism image (%s)' % my_outpath)
            return my_outpath

    # def calculate_synthetics_to_plot(self, gg, azimuth_list, matrix_M, size):
    def calculate_synthetics_to_plot(self, gg, ev2sta, matrix_M, size):
        """Calculate the data series
        to plot against the station
        data in the plots.

        !!! NOTE: Translation to Dreger like
                  code makes debugging easier. 
                  RLN (2011-12-07)
        """
        mxx = matrix_M[0, 0]
        mxy = matrix_M[0, 1]
        mxz = matrix_M[0, 2]
        myy = matrix_M[1, 1]
        myz = matrix_M[1, 2]
        mzz = matrix_M[2, 2]

        if self.verbosity > 0:
            logmt(1, 'Calculating synthetics to plot (rotated) for orid: %s' % self.orid)
            logmt(1, 'Size of synthetics: %s, Override size: %s' % (len(gg), size))

        syn_plot_dict = defaultdict(lambda: defaultdict(defaultdict))
        rotated_comps = ['T', 'R', 'Z']

        for i, sta_tup in enumerate(ev2sta):
            sta = sta_tup[0]
            az = sta_tup[1]
            if self.verbosity > 0:
                logmt(1, ' - Determine plot synthetics for station (%s)' % sta)
            for rc in rotated_comps:
                if self.verbosity > 0:
                    logmt(1, '  - Work on rotated channel (%s)' % rc)
                for j in range(size):
                    '''
                    !!! FIX: In Dreger's code there is a timeshift
                             but that timeshift is set to 0, i.e. use
                             all datapoints from the Greens function.
                             RLN (2011-12-07)
                    '''
                    if rc == 'T':
                        syn_plot_dict[sta][rc][j] = (mxx*0.5*gg[i]['TSS'][j]*math.sin(2*az) 
                                                - myy*0.5*gg[i]['TSS'][j]*math.sin(2*az)
                                                - mxy*gg[i]['TSS'][j]*math.cos(2*az)
                                                - mxz*gg[i]['TDS'][j]*math.sin(az)
                                                + myz*gg[i]['TDS'][j]*math.cos(az))
                    elif rc == 'R':
                        syn_plot_dict[sta][rc][j] = (mxx*0.5*gg[i]['XDD'][j]
                                                - mxx*0.5*gg[i]['XSS'][j]*math.cos(2*az)
                                                + mxx*0.3333*gg[i]['REX'][j]
                                                + myy*0.5*gg[i]['XDD'][j]
                                                + myy*0.5*gg[i]['XSS'][j]*math.cos(2*az)
                                                + myy*0.3333*gg[i]['REX'][j]
                                                + mzz*0.3333*gg[i]['REX'][j]
                                                - mxy*gg[i]['XSS'][j]*math.sin(2*az)
                                                + mxz*gg[i]['XDS'][j]*math.cos(az)
                                                + myz*gg[i]['XDS'][j]*math.sin(az))
                    elif rc == 'Z':
                        syn_plot_dict[sta][rc][j] = (mxx*0.5*gg[i]['ZDD'][j]
                                                - mxx*0.5*gg[i]['ZSS'][j]*math.cos(2*az)
                                                + mxx*0.3333*gg[i]['ZEX'][j]
                                                + myy*0.5*gg[i]['ZDD'][j]
                                                + myy*0.5*gg[i]['ZSS'][j]*math.cos(2*az)
                                                + myy*0.3333*gg[i]['ZEX'][j]
                                                + mzz*0.3333*gg[i]['ZEX'][j]
                                                - mxy*gg[i]['ZSS'][j]*math.sin(2*az)
                                                + mxz*gg[i]['ZDS'][j]*math.cos(az)
                                                + myz*gg[i]['ZDS'][j]*math.sin(az))
        return syn_plot_dict

    def normalize_coefficient(self, data_as_list):
        """Determine the 
        normalization factor
        for all the data
        """
        normalizer = max([abs(x) for x in data_as_list])
        if self.verbosity > 1:
            logmt(1, "  - Normalization factor: %s" % normalizer)
        return normalizer

    def create_data_synthetics_plot(self, ss, ev2sta, mod_gg, size):
        """Create and save data vs. 
        synthetic waveform plots. 
        Equivalent to Dreger's function 
        mt_plot in mt_plot6iso2_linux2.c:

        mt_plot(ss,gg,nsta,Strike,Rake,Dip,
                St2,Rk2,Dp2,d_mt,Pdc,Pclvd,Piso,Mo, Mw, E, VR);

        There are as many TRZ 
        plots as stations
        """
        if self.verbosity > 0:
            logmt(1, 'Create plots of data vs. synthetics')

        # Init figure
        my_plot = plt.figure(figsize=(9, len(ev2sta)+0.5), dpi=100)
        my_plot.subplots_adjust(hspace=0.05, wspace=0.02,
                                bottom=0.02, top=0.96,
                                left=0.07, right=0.98)
        # Immutable rotated components
        # The order is important!
        rotated_components = ('T', 'R', 'Z')

        if self.verbosity > 1:
            logmt(1, "Synthetics (Green's functions) used in calculations:")
            pprint(mod_gg)

        rows = len(ev2sta)
        cols = len(rotated_components)
        axis_num = 0
        # Subplot headings
        my_plot.text(0.22, 0.985, 'Tangential', ha='center', va='top')
        my_plot.text(0.525, 0.985, 'Radial', ha='center', va='top')
        my_plot.text(0.83, 0.985, 'Vertical', ha='center', va='top')

        for i, tup in enumerate(ev2sta):
            s = tup[0]
            if self.verbosity > 0:
                logmt(1, 'Creating plots for station (%s)' % s)
            for j, rc in enumerate(rotated_components):
                axis_num += 1
                if self.verbosity > 0:
                    logmt(1, ' - Processing (%s)' % rc)
                    logmt(1, '  - ss[%s][%s], length:%s, timeshift: %s' % (i, rc, len(ss[i][rc]), ev2sta[i][3]))
                synthetic_vals = mod_gg[s][rc].values()
                if ev2sta[i][3] < 0:
                    real_start = 0
                    syn_start = int(ev2sta[i][3]) * -1
                    syn_end = len(mod_gg[s][rc])
                    synthetic_list = synthetic_vals[syn_start:syn_end]
                else:
                    real_start = ev2sta[i][3]
                    synthetic_list = synthetic_vals[0:len(mod_gg[s][rc])]
                real_end = real_start + size
                ax = plt.subplot(len(ev2sta), len(rotated_components), axis_num)

                try:
                    data_scale_factor = self.normalize_coefficient(ss[i][rc][real_start:real_end])
                except ValueError, e:
                    logmt(1, '  *** create_data_synthetics_plot(): Error: %s' % e)
                else:
                    new_data_list = [(v/data_scale_factor) for v in ss[i][rc][real_start:real_end]]

                try:
                    syn_scale_factor = self.normalize_coefficient(synthetic_list)
                except ValueError, e:
                    logmt(1, '  *** create_data_synthetics_plot(): Error: %s' % e)
                else:
                    new_synthetic_list = [(v/syn_scale_factor) for v in synthetic_list]
                '''
                if self.verbosity > 1:
                    print "\n\n  - NEW DATA LIST:"
                    pprint(new_data_list)
                    print "\n\n  - NEW SYNTHETICS LIST:"
                    pprint(new_synthetic_list)
                '''
                plt.plot(new_data_list, 'b', new_synthetic_list, 'g')
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel(s, rotation='horizontal')
        syn_plot_outfile = '%s/%s_synthetics_fit.png' % (self.mt_images_dir, self.orid)
        my_plot.savefig(syn_plot_outfile)
        if os.path.isfile(syn_plot_outfile) and self.verbosity > 0:
            logmt(1, 'Successfully created data vs. synthetics plot (%s)' % syn_plot_outfile)
        elif not os.path.isfile(syn_plot_outfile):
            logmt(3, 'Error creating data vs. synthetics plot (%s)' % syn_plot_outfile)
        return syn_plot_outfile

    def create_composite_plot(self, ttfont, img_anno, synthetics_img, focalmech_img):
        """
        Create a composite image 
        from focal mechanism and 
        synthetics plots and update 
        the database

        All PIL work
        """ 
        if self.verbosity > 0:
            logmt(1, 'Merging the image annotations, synthetics & focal mechanism plots together')

        # Create composite image
        size = (1400, 900)
        white = (255, 255, 255, 255)
        try:
            myfont = ImageFont.truetype(ttfont, 20)
        except IOError as e:
            logmt(3, 'Error importing font: %s' % e)

        final_file = '%s.png' % self.orid
        path_to_file = '%s/%s' % (self.mt_images_dir, final_file)
        syn_img = Image.open(synthetics_img, 'r')
        fm_img = Image.open(focalmech_img, 'r')
        fm_position = (size[0] - (fm_img.size)[0] - 50, size[1] - (fm_img.size)[1] - 50)

        composite = Image.new('RGBA', size, white)
        composite.paste(syn_img, (10, 10))
        composite.paste(fm_img, fm_position)

        draw = ImageDraw.Draw(composite)
        position = (1000, 20)
        incr = 0
        for anno in img_anno:
            new_position = (position[0], position[1] + (incr*25))
            complete_anno = '%s = %s' % (anno[0], anno[1])
            draw.text(new_position, complete_anno, font=myfont, fill='black')
            incr += 1

        try:
            composite.save(path_to_file, 'PNG')
        except IOError as e:
            logmt(3, 'Cannot save file (%s). Error: %s' % (final_file, e))
        else:
            # Update the database table
            mtimages_dbptr = antdb.dbopen(self.event_db, 'r+')
            try:
                mtimages_dbptr.lookup(table='moment_tensor_images')
            except Exception, e:
                logmt(3, "Cannot open table 'moment_tensor_images'. Do you have the schema extension correctly installed? Error: %s" % e)
            else:
                orid_subset = antdb.dbsubset(mtimages_dbptr, 'orid == %s' % self.orid)
                if orid_subset.query('dbRECORD_COUNT') == 0:
                    logmt(1, 'Adding new focal mechanism to moment_tensor_images table with orid (%s)' % self.orid)
                    try:
                        mtimages_dbptr.addv(
                            'sta', '109C',
                            'orid', int(self.orid),
                            'dir', self.mt_images_dir,
                            'dfile', final_file)
                    except Exception, e:
                        logmt(3, 'Adding record to moment_tensor_images table unknown error: %s' % e)
                    else:
                        logmt(1, 'Successfully added record with orid (%s) to moment_tensor_images table' % self.orid)
                else:
                    logmt(1, 'Updating focal mechanism to moment_tensor_images table with orid (%s). Deleting current record and rewriting.' % self.orid)
                    for i in range(mtimages_dbptr.query('dbRECORD_COUNT')):
                        mtimages_dbptr[3] = i
                        try:
                            mtimages_dbptr.putv(
                                'sta', '109C',
                                'orid', int(self.orid),
                                'dir', self.mt_images_dir,
                                'dfile', final_file)
                        except Exception, e:
                            logmt(3, 'Update record in moment_tensor_images table unknown error: %s' % e)
                        else:
                            logmt(1, 'Successfully updated record with orid (%s) to moment_tensor_images table' % self.orid)
                orid_subset.free()
                mtimages_dbptr.free()
                mtimages_dbptr.close()
        return True
    # }}}

def main():
    """Get the moment tensor solution
    for a particular origin
    """
    # Parse command line
    pf, orid, verbose, debug = configure()
    verbosity = determine_verbosity(verbose, debug)

    # Parse the parameter file
    (chan_to_use, trim_value, isoflag, model_name,
     model_type, statmax, clip_values, use_inc, event_db, 
     wave_db, green_db, resp_db, 
     mag_filters, mt_images_dir, 
     ttfont, obspy_beachball, 
     distance_weighting) = parse_pf(pf, verbosity)
    mt_images_dir = '%s/%s' % (mt_images_dir, orid)

    # Instantiate objects. RLN (2011-11-21)
    my_event = Event(orid, event_db, mt_images_dir, verbosity)
    my_mt = MomentTensor(distance_weighting, isoflag, trim_value, mt_images_dir, verbosity)
    evdbptr, evparams, filter_string, filter_timepad = my_event.extract_data(mag_filters)

    '''
    !!! NOTE: We have an inherent problem
              in only selecting the same
              number of samples as the
              Greens Function length. If
              we need to timeshift the 
              data to get a good fit, and
              the shift goes beyond the limit
              of the size of the stachan_traces
              we will error out. So make
              the stachan_traces have a length
              of 3 * size of Greens Functions.
              RLN (2011-12-06)
    '''

    ev2sta_azimuth_groups = my_event.get_stations_and_orientations(evdbptr)
    '''
    !!! NOTE: Go through each of the eight 
              groups until we get good 
              cross-correlation results 
              with the synthetics
              (that we have to generate for 
              each distance & depth)
              Ideally we need two per quad
              for a total of 8 stations
    '''
    ss = []
    gg = []
    stations_per_group = 1
    # stations_per_group = int(round(statmax/8))
    good_cross_corr_stations = {'NNE':0, 'ENE':0, 'ESE':0, 'SSE':0, 'SSW':0, 'WSW':0, 'WNW':0, 'NNW':0}
    '''
    !!! NOTE: This keeps track of what 
              we used and the metadata 
              associated with the station.
              It is important that it is
              a list - order is critical.
              RLN (2011-12-13)
    '''
    ev2sta = []

    if verbosity > 0:
        logmt(1, "Iterate over groups of stations in each azimuthal range. Requesting %s per group" % stations_per_group)
    for grp in ev2sta_azimuth_groups:
        if verbosity > 0:
            logmt(1, " - Working on group (%s) " % grp)
        sta_list = ev2sta_azimuth_groups[grp]
        if stations_per_group > len(sta_list):
            if verbosity > 0:
                logmt(1, " - # of stations requested (%s) > # of stations in %s azimuthal range (%s). Using max number in azimuthal range." % (stations_per_group, grp, len(sta_list)))
            stations_per_group = len(sta_list)
        for sta in sta_list:
            if good_cross_corr_stations[grp] < stations_per_group:
                stacode, esaz, depth, distance, arrival_time = sta
                depth = int(depth)
                distance = int(distance)
                if verbosity > 0:
                    logmt(1, "  - Generate Green's function for this depth (%skm) & distance (%skm)" % (depth, distance))
                '''
                !!! NOTE: Dynamic creation of Greens Functions
                          using Juan's new Python module. RLN (2011-12-02)
                          # {{{ Structure of ours vs Dreger:
                          TSS vs u1
                          TDS vs u2
                          XSS vs u3
                          XDS vs u4
                          XDD vs u5
                          ZSS vs u6
                          ZDS vs u7
                          ZDD vs u8
                          REX vs u9
                          ZEX vs u10
                          # }}}
                          Create for every stations distance
                          to the event.
                          RLN (2011-12-08)
                '''
                green = fkr.GreenFunctions(model_name)
                green.build(depth, distance, 1, 'v', filter_string)
                if verbosity > 1:
                    green.plot()
                synthetics = green.ALL
                samples = len(synthetics['TSS']) + filter_timepad
                real_data = my_event.get_chan_data(wave_db, chan_to_use, esaz, arrival_time, filter_string, filter_timepad, stacode, samples, clip_values)
                max_xcor, timeshift = my_mt.get_time_shift(real_data, synthetics, stacode, samples)
                shifted_real_data = {}
                for chan in real_data:
                    shifted_real_data[chan] = []
                    if timeshift < 0:
                        prepend_zeros = []
                        for i in range(abs(timeshift)):
                            prepend_zeros.append(0)
                        end = samples - filter_timepad
                        real_data[chan].insert(0, prepend_zeros)
                        shifted_real_data[chan] = real_data[chan][0:end]
                    else:
                        start = timeshift
                        end = samples + start - filter_timepad
                        shifted_real_data[chan] = real_data[chan][start:end]
                ss.append(shifted_real_data)
                gg.append(synthetics)
                ev2sta.append((stacode, esaz, distance, timeshift))
                good_cross_corr_stations[grp] += 1
            else:
                logmt(1, " - Found enough good matches between data and synthetics for group (%s)" % grp)
                break
    '''
    !!! NOTE: Dynamic creation of Greens Functions
              using Juan's new Python module. RLN (2011-12-02)
              # {{{ Info on structure
              Structure of ours vs Dreger:
              TSS vs u1
              TDS vs u2
              XSS vs u3
              XDS vs u4
              XDD vs u5
              ZSS vs u6
              ZDS vs u7
              ZDD vs u8
              REX vs u9
              ZEX vs u10
              # }}}
              Create for every stations distance
              to the event.
              RLN (2011-12-08)
    '''
    nl = len(gg[0]['TSS'])

    if verbosity > 0:
        print 'Stachan traces:'
        # pprint(ss)
        print len(ss[0]['Z'])
        print 'Event to station details:'
        pprint(ev2sta)
        print 'Greens Functions:'
        # pprint(gg)
        print len(gg[0]['XSS'])

    # !!! NOTE: Instantiate MomentTensor. RLN (2011-11-21)
    my_mt = MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
    # ss = my_mt.construct_data_matrix(stachan_traces)

    '''
    !!! FIX: Override the number of data points 
             like Dreger does to be 120 values.
             HACK!
             RLN (2011-12-07)
    '''
    nl = 512

    if len(ss) != 0:
        logmt(1, 'Data matrix S created --> %s stations used' % len(ss))

    '''
    !!! NOTE: New weighting factor based on timeshift:
              If timeshift < 50, include station in
              solution. Else reject.
              RLN (2011-12-02)

    !!! NOTE: INVERSION ROUTINE
              Dreger normalizes AtA (B?) and AIV (AIV) matrix. 
              We don't need to - just create default dictionary
              RLN (2011-08-24)
    '''

    AIV = defaultdict(dict)
    B = defaultdict(dict) 

    # AIV, B = my_mt.matrix_AIV_B(ss, gg, ev2sta_azimuths, timeshift, nl)
    # AIV, B = my_mt.matrix_AIV_B(ss, gg, ev2sta, nl)
    AIV, B = my_mt.matrix_AIV_B_SCIPY(ss, gg, ev2sta, nl)
    exit()
    # !!! NOTE: DETERMINE SOLUTION VECTOR. RLN (2011-08-24)
    M = my_mt.determine_solution_vector(AIV, B)

    gfscale = -1.0e+20
    gfscale = -1.0
    strike = []
    dip = []
    rake = []

    # !!! NOTE: DECOMPOSE MOMENT TENSOR INTO VARIOUS REPRESENTATIONS. RLN (2011-11-03)
    m0, Mw, strike, dip, rake, pcdc, pcclvd, pciso = my_mt.decompose_moment_tensor(M)
    # E, VR, VAR, svar, sdpower = my_mt.fitcheck(gg, ss, W, M, m0, timeshift, ev2sta_azimuths, nl)
    # E, VR, VAR, svar, sdpower = my_mt.fitcheck(gg, ss, W, M, m0, ev2sta, nl)
    E, VR, VAR, svar, sdpower = my_mt.fitcheck(gg, ss, M, m0, ev2sta, nl)

    qlt = my_mt.quality_check(VR)

    try:
        evdbptr.close()
    except Exception, e:
        logmt(1, 'Error closing event database. Already closed? Exception: %s' % e)

    '''
    !!! NOTE: In ObsPy there are two different 
              ways to create a focal mechanism.
              Allow both.  RLN (2011-12-02)
    '''
    my_event.update_moment_tbl(strike, dip, rake)

    focalmech_img = my_event.create_focal_mechanism(obspy_beachball, M, False, False, False)
    # focalmech_img = my_event.create_focal_mechanism(obspy_beachball, mt_images_dir, False, strike, dip, rake)

    # !!! NOTE: Use a list of tuples as order is important. RLN (2011-11-28)
    image_annotations = [ ('strike', '%.3f' % strike[0]), 
                          ('rake', '%.3f' % rake[0]), 
                          ('dip', '%.3f' % dip[0]), 
                          ('m0', '%.3f' % m0), 
                          ('Mw', '%.3f' % Mw), 
                          ('pcdc', '%.3f' % pcdc), 
                          ('pcclvd', '%.3f' % pcclvd), 
                          ('pciso', '%.3f' % pciso), 
                          ('VR', '%.3f' % VR), 
                          ('VAR', '%.3f' % VAR)
                          ]
    # !!! FIX: CREATE DATA vs. SYNTHETICS PLOTS - DONE. RLN (2011-11-03)
    mod_gg = my_event.calculate_synthetics_to_plot(gg, ev2sta, M, nl)
    # pprint(mod_gg)

    synthetics_img = my_event.create_data_synthetics_plot(ss, ev2sta, mod_gg, nl)

    # !!! NOTE: CREATE THE COMPOSITE (FINAL) IMAGE AND UPDATE DB. RLN (2011-11-28)
    my_event.create_composite_plot(ttfont, image_annotations, synthetics_img, focalmech_img)

    print "MOMENT TENSOR:"
    print M
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

    # !!! NOTE: From Dreger's fitcheck fprintf() commands. RLN (2011-08-24)
    for j in range(len(svar)):
        print 'Station (%s) = %s  %s' % (svar[j][0], svar[j][1], sdpower[j])
    print 'Var     = %s' % E
    print 'VR      = %s (UNWEIGHTED)' % VR
    print 'VR      = %s (WEIGHTED)' % VR

    print 'Var/Pdc = ', E/pcdc
    print 'Quality = %s' % qlt

    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)

"""
dbmoment.py

@authors  Gert-Jan van den Hazel <hazelvd@knmi.nl>
          Juan Reyes <jreyes1108@gmail.com>
          Rob Newman <robertlnewman@gmail.com>
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

# SIGNAL PROCESSING
from pylab import *
from scipy.signal.filter_design import butter, buttord
from scipy.signal import lfilter

# FOR CROSS CORRELATION
from scipy import fftpack

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
    model_type = stock.pfget_string(pfname, 'model_type')
    statmax = stock.pfget_int(pfname,'statmax')
    clip_values = stock.pfget_arr(pfname, 'clip_values')
    use_inc = stock.pfget_string(pfname, 'use_inc')
    event_db = stock.pfget_string(pfname,'event_db')
    try:
        wave_db = stock.pfget_string(self.pfname, 'wave_db')
    except:
        wave_db = False
    # !!! FIX: Currently not using a Greens db, just a file. RLN (2011-11-03)
    # green_db = stock.pfget_string(pfname, 'green_db')
    green_file = stock.pfget_string(pfname, 'green_file')
    green_pf = stock.pfget_string(pfname, 'green_pf')
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

    # return (chan_to_use, trim_value, isoflag, 
    #         model_type, statmax, use_inc, event_db, 
    #         wave_db, green_db, green_pf, resp_db, 
    #         mag_filters, mt_images_dir, obspy_beachball, 
    #         distance_weighting)
    return (chan_to_use, trim_value, isoflag, 
            model_type, statmax, clip_values, use_inc, event_db, 
            wave_db, green_file, green_pf, resp_db, 
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
        # print message
    else:
        print curtime, 'ERROR:',message,'--> exiting'
        # print 'ERROR:',message,'--> exiting'
        sys.exit(-1)
# }}}

class MomentTensor():
    """Class for building moment tensors
    and doing the inversion
    """
    # {{{

    def __init__(self, distance_weighting, isoflag, trim_value, verbosity=0):
        """Initialize"""
        self.distance_weighting = distance_weighting
        self.isoflag = isoflag
        self.trim_value = trim_value
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

    def cross_cor(self, a, b):
        """Calculate cross correlation 
        between data and accompanying 
        Green's function. Returns the 
        max_cor and the shift
        """
        if self.verbosity > 1:
            logmt(1, ' - Calculating cross correlation between station data and synthetic (Greens function) data')
        xcor = np.correlate(a, b, 'full')
        maxval = np.amax(xcor)
        maxshift = np.argmax(xcor)
        # Does this shift the data or the synthetic?
        shift = maxshift - (len(xcor)-1)/2

        if self.verbosity > 1:
            logmt(1, "  - Max val: %s. Timeshift (in samples): %s" % (maxval, shift))
            fig = plt.figure()
            ax = fig.add_subplot(311)
            # ax.title('Raw data')
            ax.plot(self.normalize(a), 'b-')
            ax.plot(self.normalize(b), 'g-')
            # Print shifted time series
            bx = fig.add_subplot(312)
            if shift < 0:
                [a.insert(0, 0) for x in range(abs(shift))]
            elif shift > 0:
                [b.insert(0, 0) for x in range(abs(shift))]
            bx.plot(self.normalize(a), 'r-')
            bx.plot(self.normalize(b), 'c-')
            # Print xcor
            cx = fig.add_subplot(313)
            # cx.title('Cross correlation results')
            cx.plot(xcor)
            plt.ylabel('xcor')
            plt.xlabel('point shift')
            plt.show()
        return maxval, shift

    def normalize(self, data_as_list):
        """Determine the 
        normalization factor
        for all the data
        """
        normalizer = max([abs(x) for x in data_as_list])
        if self.verbosity > 0:
            logmt(1, "  - Normalization factor: %s" % normalizer)
        return [x / normalizer for x in data_as_list]
        # return normalized_list

    def get_time_shift(self, data_dict, greens_dict, size):
        """Get time shift for each station
        using cross correlation of station 
        data points and the Green's Function
        data points. Make sure both matrices
        are the same size for the comparison.
        Return a list of tuples (stacode, timeshift).
        """
        if self.verbosity > 0:
            logmt(1, '  - Get the time shift for all components & return the mean shift shift')
        '''
        !!! NOTE: Loop over each stations data points.
                  Remember that the Green's Function 
                  matrix is just two dimensional - not 
                  replicating the number of stations
                  (3D matrix). RLN (2011-12-02)
        '''
        shift = 0
        xcor = 0
        if self.cross_cor(data_dict['T'][0:size], greens_dict['TSS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['T'][0:size], greens_dict['TSS'][0:size])
        if self.cross_cor(data_dict['T'][0:size], greens_dict['TDS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['T'][0:size], greens_dict['TDS'][0:size])
        if self.cross_cor(data_dict['R'][0:size], greens_dict['XSS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['R'][0:size], greens_dict['XSS'][0:size])
        if self.cross_cor(data_dict['R'][0:size], greens_dict['XDS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['R'][0:size], greens_dict['XDS'][0:size])
        if self.cross_cor(data_dict['R'][0:size], greens_dict['XDD'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['R'][0:size], greens_dict['XDD'][0:size])
        if self.cross_cor(data_dict['Z'][0:size], greens_dict['ZSS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['Z'][0:size], greens_dict['ZSS'][0:size])
        if self.cross_cor(data_dict['Z'][0:size], greens_dict['ZDS'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['Z'][0:size], greens_dict['ZDS'][0:size])
        if self.cross_cor(data_dict['Z'][0:size], greens_dict['ZDD'][0:size])[0] > xcor:
            xcor, shift = self.cross_cor(data_dict['Z'][0:size], greens_dict['ZDD'][0:size])
        if self.verbosity > 0:
            logmt(1, '  - Cross-correlation value: %s, timeshift value: %s' % (xcor, shift))
        return xcor, shift

    # def matrix_AIV_B(self, dict_s, dict_g, list_az, timeshift, size):
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
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            if self.verbosity > 0:
                logmt(1, '  - Working on inversion for station (%s)' % ev2sta[i][0])
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
                AJ[3][cnt1] = -math.sin(this_azimuth)*dict_g[i]['TDS'][i]
                AJ[3][cnt2] =  math.cos(this_azimuth)*dict_g[i]['XDS'][i]
                AJ[3][cnt3] =  math.cos(this_azimuth)*dict_g[i]['ZDS'][i]

                # Myz term
                AJ[4][cnt1] =  math.cos(this_azimuth)*dict_g[i]['TDS'][i]
                AJ[4][cnt2] =  math.sin(this_azimuth)*dict_g[i]['XDS'][i]
                AJ[4][cnt3] =  math.sin(this_azimuth)*dict_g[i]['ZDS'][i]

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
                            logmt(3, 'determine_solution_vector(): ERROR... 1')
            ipiv[icol] += 1
            if not irow == icol:
                for l in range(self.isoflag):
                    (dict_AIV[irow][l],dict_AIV[icol][l]) = swap(dict_AIV[irow][l],dict_AIV[icol][l])
                for l in range(1):
                    (dict_B[irow][l],dict_B[icol][l]) = swap(dict_B[irow][l],dict_B[icol][l])
            if dict_AIV[icol][icol] == 0.0:
                logmt(3, 'determine_solution_vector(): ERROR... 2')
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

    def fitcheck(self, dict_g, dict_s, list_W, matrix_M, m0, ev2sta, size):
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
    # }}}

class Event():
    """Class for extracting per event 
    data from the database
    and creating the plots
    """
    # {{{

    def __init__(self, orid, event_db, verbosity=0):
        """Initialize"""
        self.orid = orid
        self.event_db = event_db
        self.verbosity = verbosity

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
        except:
            logmt(3, 'Could not open database (%s)' % self.event_db)
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
                except:
                    logmt(3, 'Could not find field (%s) in join of origin & netmag tables for orid %s' % (f, self.orid))
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
            lower_mags = []
            upper_mags = []
            filters = []
            for line in mag_filters:
                splitline = line.split()
                if not len(splitline) == 3:
                    lower_magnitude, filter = splitline
                else:
                    lower_magnitude, upper_magnitude, filter = splitline
                lower_mags.append(lower_magnitude)
                upper_mags.append(upper_magnitude)
                filters.append(filter)
            min_mag = lower_mags[0] # !!! NOTE: Assumes order correct in pf. RLN (2011-07-28)
            max_mag = upper_mags[-1] # !!! NOTE: Assumes order correct in pf. RLN (2011-07-28)
            filter = filters[-1] # !!! NOTE: Assumes all filters are the same. RLN (2011-07-28)
            if float(evparams['magnitude']) > float(min_mag) and float(evparams['magnitude']) < float(max_mag):
                filter_string = filter.replace('_', ' ')
            else:
                logmt(3, 'Magnitude %s not within bounds (%s - %s) defined in the parameter file' % (evparams['magnitude'], float(min_mag), float(max_mag)))
            return evdb, evparams, filter_string

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
            sta, esaz, depth, at, ev_lat, ev_lon, site_lat, site_lon = dbptr.getv('sta', 'esaz', 'depth', 'arrival.time', 'lat', 'lon', 'site.lat', 'site.lon')
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

    def get_chan_data(self, wave_db, chan_to_use, esaz, at, filter_string, stacode, size, clip_values):
        """Opens waveform database and returns
        trace objects based on sta_chan. Applies
        calibration, splice, filter, decimation
        and rotation of channels
        """

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
        except:
            logmt(3, 'Could not open waveform database %s' % wave_db)
            return False
        else:
            wvdb.lookup(table='wfdisc')
            wvdb.subset('sta =~ /%s/' % stacode)
            wvdb[3] = 0
            '''
            !!!! NOTE: Want maximium x2 the size 
                       to account for timeshifts
                       and add 10 seconds before
                       RLN (2011-12-06)
            st   p arrival.time             et
            | <-------->*<--------><-------->|
            '''
            st = at - 10
            et = at + (size * 3) 
            resample = 0
            vang = 0
            if self.verbosity > 0:
                logmt(1, '    - Retrieve %s data with time range (%s, %s)' % (stacode, st, et))
            wvdb.subset('sta=~/^%s$/ && chan=~/%s/' % (stacode, chan_to_use))
            trace = wvdb.load_css(st, et)
            trace.apply_calib()
            trace.splice()
            # Comment out this filtering - need same filtering as GF
            # trace.filter(filter_string)
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
                samples = trace.data()

                # {{{ Bandpass filters
                ord_1 = 5
                wn_1 = 0.01 * math.pi * 2
                b_1, a_1 = butter(ord_1, wn_1, btype='low')

                ord_2 = 5
                wn_2 = 0.07 * math.pi * 2
                b_2, a_2 = butter(ord_2, wn_2, btype='low')

                samples = lfilter(b_1, a_1, samples)
                samples = (lfilter(b_2, a_2, samples)).tolist()
                # }}} Bandpass filters

                use_data = 0

                # {{{ Test for clipped data
                for i in samples:
                    if i > float(clip_values['max']):
                        use_data = 1
                    elif i < float(clip_values['min']):
                        use_data = -1
                # }}}

                if use_data == 1:
                    logmt(1, '    - %s: One or more samples > clip_values[max] (%s). IGNORE.' % (stachan, clip_values['max']))
                elif use_data == -1:
                    logmt(1, '    - %s: One or more samples < clip_values[min] (%s). IGNORE.' % (stachan, clip_values['min']))
                else:
                    '''
                    !!! NOTE: Get three times the size 
                              to allow for time shift of 
                              max length either side.
                              RLN (2011-12-06)
                    '''
                    max_size = (size * 3) + 10 # Only works if 1Hz (1 sample per second, LH.*)
                    if self.verbosity > 1 and len(samples) > max_size:
                        logmt(1, '  - Sample size for %s: %s. Max: %s. Trimming.' % (chan, len(samples), max_size))
                    stachan_trace[chan] = samples[0:max_size]  
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

    def create_focal_mechanism(self, obspy_beachball, mt_images_dir, matrix_M=False, strike=False, dip=False, rake=False):
        """Write out focal mechanism
        to images directory
        """
        if self.verbosity > 0:
            logmt(1, 'Writing file to images dir (%s)' % mt_images_dir)

        if not os.path.exists(mt_images_dir):
            if self.verbosity > 0:
                logmt(1, 'Images dir (%s) does not exist. Try to create...' % mt_images_dir)
            try:
                os.mkdir(mt_images_dir, 0775)
            except:
                logmt(3, 'Moment tensor images dir (%s) does not exist and cannot be created!' % mt_images_dir)

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
        my_outpath = '%s/%s' % (mt_images_dir, my_outfile)

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
        if self.verbosity > 0:
            logmt(1, "  - Normalization factor: %s" % normalizer)
        return normalizer

    def create_data_synthetics_plot(self, ss, ev2sta, mod_gg, size, mt_images_dir):
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
                    logmt(1, 'ss[%s][%s], length:%s, timeshift: %s' % (i, rc, len(ss[i][rc]), ev2sta[i][3]))
                start = ev2sta[i][3]
                end = start + size
                ax = plt.subplot(len(ev2sta), len(rotated_components), axis_num)
                print start
                print end
                print ss[i][rc][start:end]
                data_scale_factor = self.normalize_coefficient(ss[i][rc][start:end])
                syn_scale_factor = self.normalize_coefficient(mod_gg[s][rc].values())
                new_data_list = [(v/data_scale_factor) for v in ss[i][rc][start:end]]
                new_synthetic_list = [(v/syn_scale_factor) for v in mod_gg[s][rc].values()]
                if self.verbosity > 1:
                    print "\n\n  - NEW DATA LIST:"
                    pprint(new_data_list)
                    print "\n\n  - NEW SYNTHETICS LIST:"
                    pprint(new_synthetic_list)
                plt.plot(new_data_list, 'b', new_synthetic_list, 'g')
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel(s, rotation='horizontal')
        syn_plot_outfile = '%s/%s_synthetics_fit.png' % (mt_images_dir, self.orid)
        my_plot.savefig(syn_plot_outfile)
        if os.path.isfile(syn_plot_outfile) and self.verbosity > 0:
            logmt(1, 'Successfully created data vs. synthetics plot (%s)' % syn_plot_outfile)
        elif not os.path.isfile(syn_plot_outfile):
            logmt(3, 'Error creating data vs. synthetics plot (%s)' % syn_plot_outfile)
        return syn_plot_outfile

    def create_composite_plot(self, ttfont, img_anno, mt_images_dir, synthetics_img, focalmech_img):
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
        path_to_file = '%s/%s' % (mt_images_dir, final_file)
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
                            'dir', mt_images_dir,
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
                                'dir', mt_images_dir,
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
    (chan_to_use, trim_value, isoflag, 
     model_type, statmax, clip_values, use_inc, event_db, 
     wave_db, green_file, green_pf, resp_db, 
     mag_filters, mt_images_dir, 
     ttfont, obspy_beachball, 
     distance_weighting) = parse_pf(pf, verbosity)

    # !!! NOTE: Instantiate Event. RLN (2011-11-21)
    my_event = Event(orid, event_db, verbosity)
    my_mt = MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
    evdbptr, evparams, filter_string = my_event.extract_data(mag_filters)

    '''
    !!! NOTE: Is this relevant?
              In Dreger's code look
              for the comments:
              /*Note the vertical GF's are*/ 
              /*flipped in earqt1.f and TW's*/
              /* Blackbox.f DVH conv. z + down*/
              RLN (2011-12-07)
    '''
    '''
    if verbosity > 1:
        print "\n\n\nOLD GREENS FUNCTION FORMAT:"
        pprint(gg)
    new_zss = [val*-1 for val in gg['ZSS']]
    new_zds = [val*-1 for val in gg['ZDS']]
    new_zdd = [val*-1 for val in gg['ZDD']]
    gg['ZSS'] = new_zss
    gg['ZDS'] = new_zds
    gg['ZDD'] = new_zdd
    if verbosity > 1:
        print "\n\n\nNEWLY FORMATTED GREENS FUNCTION, VERTICAL COMPS * -1 (LOOK AT ZSS, ZDS & ZDD):"
        pprint(gg)
    '''

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

    !!! NOTE: We don't want to be limited by
              the number of stations in 
              determining a solution. Need a
              while loop that makes sure the 
              cross-correlation is good before
              accepting the data as part of 
              generating the MT solution.
              Would be nice to have a minimum
              of 8 stations with good cross-
              correlations results to use.
              RLN (2011-12-12)
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
    stations_per_group = int(round(statmax/8))
    print stations_per_group
    good_cross_corr_stations = {'NNE':0, 'ENE':0, 'ESE':0, 'SSE':0, 'SSW':0, 'WSW':0, 'WNW':0, 'NNW':0}
    '''
    !!! NOTE: This keeps track of what 
              we used and the metadata 
              associated with the station.
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
                    logmt(1, "  - Getting data for station %s (distance = %s, depth = %s)" % (stacode, distance, depth))
                # Not sure we still need to define 200, but leave for now
                real_data = my_event.get_chan_data(wave_db, chan_to_use, esaz, arrival_time, filter_string, stacode, 200, clip_values)
                if verbosity > 0:
                    logmt(1, "  - Generate Green's function for this depth (%s) & distance (%s)" % (depth, distance))
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
                green = fkr.GreenFunctions('SOCAL_MODEL')
                all_greens = green.generate(depth, distance)
                if verbosity > 1:
                    green.plot()
                synthetic_data = green.ALL
                # {{{ Bandpass filters on the GF
                ord_1 = 5
                wn_1 = 0.01 * math.pi * 2
                b_1, a_1 = butter(ord_1, wn_1, btype='low')

                ord_2 = 5
                wn_2 = 0.07 * math.pi * 2
                b_2, a_2 = butter(ord_2, wn_2, btype='low')
                # }}}
                for k in synthetic_data:
                    # Apply both filters to Greens Function
                    synthetic_data[k] = lfilter(b_1, a_1, synthetic_data[k])
                    synthetic_data[k] = lfilter(b_2, a_2, synthetic_data[k])
                    synthetic_data[k] = synthetic_data[k].tolist()

                # Do the cross correlation and get the time shift
                max_xcor, timeshift = my_mt.get_time_shift(real_data, synthetic_data, 120)

                # If a good match, add to the solution
                # if max_xcor > 0.8:
                ss.append(real_data)
                gg.append(synthetic_data)
                ev2sta.append((stacode, esaz, distance, timeshift))
                good_cross_corr_stations[grp] += 1
                # else:
                #     logmt(1, "  - Cross correlation (%s) did not result in a good match (>0.8). Trying another station..." % max_xcor)
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

    if verbosity > 1:
        print 'Stachan traces:'
        # pprint(stachan_traces)
        pprint(ss)
        print 'Event to station details:'
        pprint(ev2sta)
        print 'Greens Functions:'
        pprint(gg)

    # !!! NOTE: Instantiate MomentTensor. RLN (2011-11-21)
    my_mt = MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
    # ss = my_mt.construct_data_matrix(stachan_traces)

    '''
    !!! FIX: Override the number of data points 
             like Dreger does to be 120 values.
             HACK!
             RLN (2011-12-07)
    '''
    nl = 120

    if len(ss) != 0:
        logmt(1, 'Data matrix S created --> %s stations used' % len(ss))

    '''
    !!! NOTE: Weighting factor:
              number of stations to be used 
              in solution. If set to 1 then 
              station fully incorporated into 
              solution. If < 1, don't include. 
              For now, use ALL stations.
              RLN (2011-12-02)
    '''
    W = []
    for i in range(len(ss)):
        W.append(1.0)

    '''
    !!! NOTE: INVERSION ROUTINE
              Dreger normalizes AtA (B?) and AIV (AIV) matrix. 
              We don't need to - just create default dictionary
              RLN (2011-08-24)
    '''

    AIV = defaultdict(dict)
    B = defaultdict(dict) 

    # AIV, B = my_mt.matrix_AIV_B(ss, gg, ev2sta_azimuths, timeshift, nl)
    AIV, B = my_mt.matrix_AIV_B(ss, gg, ev2sta, nl)

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
    E, VR, VAR, svar, sdpower = my_mt.fitcheck(gg, ss, W, M, m0, ev2sta, nl)

    qlt = my_mt.quality_check(VR)

    try:
        evdbptr.close()
    except:
        logmt(1, 'Error closing event database. Already closed?')

    '''
    !!! NOTE: In ObsPy there are two different 
              ways to create a focal mechanism.
              Allow both.  RLN (2011-12-02)
    '''
    my_event.update_moment_tbl(strike, dip, rake)

    focalmech_img = my_event.create_focal_mechanism(obspy_beachball, mt_images_dir, M, False, False, False)
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
    pprint(mod_gg)

    synthetics_img = my_event.create_data_synthetics_plot(ss, ev2sta, mod_gg, nl, mt_images_dir)

    # !!! NOTE: CREATE THE COMPOSITE (FINAL) IMAGE AND UPDATE DB. RLN (2011-11-28)
    my_event.create_composite_plot(ttfont, image_annotations, mt_images_dir, synthetics_img, focalmech_img)

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
        print 'Station (%d) = %s  %s' % (j,svar[j],sdpower[j])
    print 'Var     = %s' % E
    print 'VR      = %s (UNWEIGHTED)' % VR
    print 'VR      = %s (WEIGHTED)' % VR

    print 'Var/Pdc = ', E/pcdc
    print 'Quality = %s' % qlt

    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)

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
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# OBSPY
try:
    from obspy.imaging.beachball import Beachball
except ImportError:
    print "Import Error: Do you have ObsPy installed correctly?"

# {{{
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

    def construct_data_matrix_hard_coded(self, stachan_traces):
        """Construct data matrix from files
        containing trace_objects per sta_chan
        Returns matrix s, which is a 3D matrix, 
        containing all data per channel
        """
        '''
        this_dict = defaultdict(lambda: defaultdict(defaultdict))
        if self.debug:
            logmt(1, 'Constructing data matrix from hard-coded files in db/data/data1-3')
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
        '''

    def cross_cor(self, a, b):
        """Calculate cross correlation 
        between data and accompanying 
        Green's function. Returns the 
        max_cor and the shift
        """
        if self.verbosity > 0:
            logmt(1, 'Calculate cross correlation between data and Greens function')
        xcor = np.correlate(a.values(), b.values(), 'full')
        pr = len(xcor)/2 
        maxval = 0
        maxcoef = 0
        for j in range(pr):
            if abs(xcor[j+pr]) > maxval:
                maxval = abs(xcor[j+pr])
                #maxcoef = j+1
                maxcoef = j
        print "maxval %s maxcoef %s" % (maxval,maxcoef)
        return (maxval, maxcoef)

    def get_time_shift(self, data_dict, greens_dict):
        """Get the time shift 
        and return a list
        """
        if self.verbosity > 0:
            logmt(1, 'Get the time shift and return a list')
            exit()
        timeshift = []
        for i in range(len(data_dict['T'])):
            shift = 0
            xcor  = 0
            if self.cross_cor(data_dict['T'][i], greens_dict[0][i])[0] > xcor:
                xcor, shift = self.cross_cor(data_dict['T'][i], greens_dict[0][i])
            if self.cross_cor(data_dict['T'][i], greens_dict[1][i])[0] > xcor:
                xcor, shift = self.cross_cor(data_dict['T'][i], greens_dict[1][i])
            if self.cross_cor(data_dict['R'][i], greens_dict[2][i])[0] > xcor:
                xcor, shift = self.cross_cor(data_dict['R'][i], greens_dict[2][i])
            if self.cross_cor(data_dict['R'][i], greens_dict[3][i])[0] > xcor:
                xcor,shift = self.cross_cor(data_dict['R'][i], greens_dict[3][i])
            if self.cross_cor(data_dict['R'][i], greens_dict[4][i])[0] > xcor:
                xcor,shift = self.cross_cor(data_dict['R'][i], greens_dict[4][i])
            if self.cross_cor(data_dict['Z'][i], greens_dict[5][i])[0] > xcor:
                xcor, shift = self.cross_cor(data_dict['Z'][i], greens_dict[5][i])
            if self.cross_cor(data_dict['Z'][i], greens_dict[6][i])[0] > xcor:
                xcor,shift = self.cross_cor(data_dict['Z'][i], greens_dict[6][i])
            if self.cross_cor(data_dict['Z'][i], greens_dict[7][i])[0] > xcor:
                xcor,shift = self.cross_cor(data_dict['Z'][i], greens_dict[7][i])
            timeshift.append(shift)
        return timeshift

    def matrix_AIV_B(self, dict_s, dict_g, list_az, this_timeshift):
        """Inversion routine
        Construct matrices AIV and B using 
        the data and Green's function matrices
        Return AIV and B for further processing 
        """
        if self.verbosity > 0:
            logmt(1, 'Construct matrices AIV and B using the data and Greens function matrices')
        AJ = defaultdict(dict) # RLN (2011-07-29): What is AJ? COPY OF DREGER
        trim = 0
        cnt1 = cnt2 = cnt3 = 0
        #trim = int(len(dict_s['T'][0]) * float(self.trim_value)) 
        if self.verbosity > 1:
            logmt(1, 'Timeshift: %s' % this_timeshift)
            logmt(1, 'Number of stations used in calculation: %s' % len(dict_s))

        # Loop over the number of stations in the dictionary
        for i in range(len(dict_s)):
            cnt1 = cnt2 = cnt3
            cnt2 += len(dict_s['T'][0]) - trim
            cnt3 += 2*len(dict_s['T'][0]) - 2*trim
            list_az[i][1] *= math.pi/180 # list_az is a tuple
            # Index over time
            for j in range(len(dict_s['T'][0])-trim):
                # Mxx term
                AJ[0][cnt1] =  math.sin(2*list_az[i][1])*dict_g[0][i][j]/2

                # Myy term
                AJ[1][cnt1] = -math.sin(2*list_az[i][1])*dict_g[0][i][j]/2

                # Mxy term
                AJ[2][cnt1] = -math.cos(2*list_az[i][1])*dict_g[0][i][j]
                AJ[2][cnt2] = -math.sin(2*list_az[i][1])*dict_g[2][i][j]
                AJ[2][cnt3] = -math.sin(2*list_az[i][1])*dict_g[5][i][j]

                # Mxz term
                AJ[3][cnt1] = -math.sin(list_az[i][1])*dict_g[1][i][j]
                AJ[3][cnt2] =  math.cos(list_az[i][1])*dict_g[3][i][j]
                AJ[3][cnt3] =  math.cos(list_az[i][1])*dict_g[6][i][j]

                # Myz term
                AJ[4][cnt1] =  math.cos(list_az[i][1])*dict_g[1][i][j]
                AJ[4][cnt2] =  math.sin(list_az[i][1])*dict_g[3][i][j]
                AJ[4][cnt3] =  math.sin(list_az[i][1])*dict_g[6][i][j]

                # Vary the other values depending on isoflag value
                if self.isoflag == 5:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[4][i][j])/2 - (math.cos(2*list_az[i][1])*dict_g[2][i][j])/2
                    AJ[0][cnt3] = (dict_g[7][i][j])/2 - (math.cos(2*list_az[i][1])*dict_g[5][i][j])/2
                    # Myy term
                    AJ[1][cnt2] = (dict_g[4][i][j])/2 + (math.cos(2*list_az[i][1])*dict_g[2][i][j])/2
                    AJ[1][cnt3] = (dict_g[7][i][j])/2 + (math.cos(2*list_az[i][1])*dict_g[5][i][j])/2
                if self.isoflag == 6:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[4][i][j])/6 - (math.cos(2*list_az[i][1])*dict_g[2][i][j])/2 + (dict_g[8][i][j])/3
                    AJ[0][cnt3] = (dict_g[7][i][j])/6 - (math.cos(2*list_az[i][1])*dict_g[5][i][j])/2 + (dict_g[9][i][j])/3
                    # Myy term
                    AJ[1][cnt2] = (dict_g[4][i][j])/6 + (math.cos(2*list_az[i][1])*dict_g[2][i][j])/2 + (dict_g[8][i][j])/3
                    AJ[1][cnt3] = (dict_g[7][i][j])/6 + (math.cos(2*list_az[i][1])*dict_g[5][i][j])/2 + (dict_g[9][i][j])/3
                    # RLN (2011-08-24): Where do these values come from and why?
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

        # Compute AtA
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
            # print "i is %d" % i
            cnt1 = cnt2 = cnt3
            cnt2 += len(dict_s['T'][0])-trim
            cnt3 += 2*(len(dict_s['T'][0])-trim)
            for j in range(len(dict_s['T'][0])-trim):
                #print "\nj is %s" % j
                #print "i is %s" % i
                #print "trim is %s" % trim
                #print "this_timeshift[i] is %s" % this_timeshift[i]
                #print "dict_s['T'][i] is %s" % len(dict_s['T'][i])
                #print "cnt1 is %s" % cnt1
                #print "tmp[cnt1] is %s" % tmp[cnt1]
                #print "cnt2 is %s" % cnt2
                #print "cnt3 is %s" % cnt3
                #tmp[cnt1] = dict_s['T'][i][j + this_timeshift[i]]
                #tmp[cnt2] = dict_s['R'][i][j + this_timeshift[i]]
                #tmp[cnt3] = dict_s['Z'][i][j + this_timeshift[i]]
                tmp[cnt1] = dict_s['T'][i][j]
                tmp[cnt2] = dict_s['R'][i][j]
                tmp[cnt3] = dict_s['Z'][i][j]
                #print "\n#### tmp[cnt1] is %s\n" % tmp[cnt1]
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1

        # Calculate Righthand Side
        for i in range(self.isoflag):
            for j in range(cnt3):
                B[i][0] += AJ[i][j] * tmp[j]
                #print "i is %s" % i
                #print "j is %s" % j
                #print "B[i][0] is %s" % B[i][0]
                #print "AJ[i][j] is %s" % AJ[i][j]
                #print "tmp[j] is %s" % tmp[j]
                #B[i][0] = B[i][0] + AJ[i][j] * tmp[j]

        # Some logging
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
        return(a,b)

    def dyadic(self, v, n1, n2, c):
        """Calculate the dyadic matrix 
        of eigenvectors v
        """
        if self.verbosity > 0:
            logmt(1, 'Compute the dyadic matrix of vector v')
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
            logmt(1, 'Moment tensor matrix M[3,3] created')
        return M

    def decompose_moment_tensor(self, matrix_M):
        """Decompose moment tensor into eigenvector/values.
        Calculate moment tensor parameters from eigenvalues and vectors. 
        Return M0, Mw, strike, slip, rake and precentages of present
        source characteristics
        From Dreger source: fmap_subs_linux.f
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
            logmt(1, 'Decomposition of moment tensor succesful')

        return m0, Mw, strike, dip, slip, pcdc, pcclvd, pciso

    def fitcheck(self, dict_g, dict_s, list_W, matrix_M, m0, timeshift, az):
        """Calculate the variance, variance reduction
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
            evdb.join('netmag')
            if evdb.nrecs() == 0:
                self.logmt(3, 'Could not join netmag table for orid (%s)' % self.orid)
            evparams = {}
            evdb[3] = 0 # Just get the first record
            for f in evdb.query('dbTABLE_FIELDS'):
                try: 
                    evparams[f] = evdb.getv(f)[0]
                except:
                    logmt(3, 'Could not find field (%s) in join of origin & netmag tables for orid %s' % (f, self.orid))
            evdb.join('assoc')
            evdb.join('arrival')
            evdb.sort('delta')
            evdb.subset('iphase=~/.*P.*|.*p.*/')
            if evdb.nrecs() == 0:
                logmt(3, 'No arrivals for selected origin %s' % self.orid)
            elif self.verbosity > 0:
                logmt(1, 'There are %d records that match this origin arrival and iphase' % evdb.nrecs())
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

    def get_chan_data(self, dbptr, wave_db, chan_to_use, filter_string, statmax, size, clip_values):
        """Opens waveform database and returns
        trace objects based on sta_chan. Applies
        calibration, splice, filter, decimation
        and rotation of channels
        """

        if self.verbosity > 0:
            logmt(1, 'Get channel data (trace objects) for statmax (%s) stations' % statmax)
        '''
        !!! NOTE: Start time and endtime need to be calculated using ptime 
        Right now we are just getting size/2 either side of first arrival
        Instead use arrival.time - should be the same thing. RLN (2011-11-29)
        '''
        stachan_traces = {} # Holder dict for data
        ev2sta_azimuths = [] # Holder list - needs to remain sorted

        if self.verbosity > 0:
            logmt(1, 'There are %s records for this event' % dbptr.nrecs())

        if dbptr.nrecs() < statmax*3:
            if self.verbosity > 0:
                logmt(1, 'Only %s station channel records available.' % dbptr.nrecs())
            statmax = dbptr.nrecs()

        # !!! NOTE: Need to get a fresh view of the wfdisc as subsetted for P arrivals
        #           which are only done on the BH.* channels. RLN (2011-11-29)
        try:
            wvdb = antdb.dbopen(wave_db, 'r')
            wvdb.lookup(table='wfdisc')
        except:
            logmt(3, 'Could not open waveform database %s' % wave_db)
        else:
            if self.verbosity > 0:
                logmt(1, 'Get each entry in the table for (%s) stations' % statmax)
            for i in range(statmax):
                dbptr[3] = i
                sta, esaz, at = dbptr.getv('sta', 'esaz', 'arrival.time')
                ev2sta_azimuths.append((sta, esaz))
                st = at - int(size/2)  
                et = at + int(size/2) 
                et += 10 # Give 10s of endtime buffer (10 samples @ 1Hz)
                resample = 0
                vang = 0
                if self.verbosity > 0:
                    logmt(1, ' - Retrieve %s data with time range (%s, %s)' % (sta, st, et))
                wvstadb = antdb.dbsubset(wvdb, 'sta=~/^%s$/ && chan=~/%s/' % (sta, chan_to_use))
                trace = wvstadb.load_css(st, et)
                trace.apply_calib()
                trace.splice()
                trace.filter(filter_string)
                rotchan = ('R', 'T', 'Z')
                trace.rotate(esaz, vang, rotchan)
                trace.subset('chan =~ /R|T|Z/')
                if self.verbosity > 1:
                    logmt(1, ' - Number of traces for %s: %s' % (sta, trace.nrecs()))
                for j in range(trace.nrecs()):
                    trace[3] = j
                    sta, chan, ns, sr = trace.getv('sta', 'chan', 'nsamp', 'samprate')
                    stachan = '%s_%s' % (sta, chan)
                    ns_tra = 0
                    ns_req = 0
                    samples = trace.data()
                    use_data = 0
                    # Test for clipped data
                    for i in samples:
                        if i > float(clip_values['max']):
                            use_data = 1
                        elif i < float(clip_values['min']):
                            use_data = -1
                    if use_data == 1:
                        logmt(1, ' - %s: One or more samples > clip_values[max] (%s)' % (stachan, clip_values['max']))
                    elif use_data == -1:
                        logmt(1, ' - %s: One or more samples < clip_values[min] (%s)' % (stachan, clip_values['min']))
                    else:
                        if self.verbosity > 1 and len(samples) > size:
                            logmt(1, '  - Sample size for %s: %s. Max: %s. Trimming.' % (stachan, len(samples), size))
                        stachan_traces[stachan] = samples[0:size]  
                        if self.verbosity > 0:
                            logmt(1, ' - Trace extracted for %s samples:[%s], wanted:[%s]' % (stachan,len(stachan_traces[stachan]),size))
                trace.trdestroy()
                wvstadb.free()
            wvdb.free()
            wvdb.close()
        return stachan_traces, ev2sta_azimuths

    def update_moment_tbl(self, strike, dip, rake):
        """Write out results to 
        database moment table
        """
        if self.verbosity > 0:
            logmt(1, 'Write out moment tensor for orid %s to database table %s.moment' % (self.orid, self.event_db))
            logmt(1, '## MT for orid %s strike => %s' % (self.orid, strike))
            logmt(1, '## MT for orid %s dip => %s' % (self.orid, dip))
            logmt(1, '## MT for orid %s rake => %s' % (self.orid, rake))
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

    def create_focal_mechanism(self, obspy_beachball, mt_images_dir, strike, dip, rake):
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

        focal_mechanism = [strike[0], dip[0], rake[0]]
        if self.verbosity > 0:
            logmt(1, 'Try to plot focal mechanism: %s' % focal_mechanism)

        beachball_vals = {}

        # From the ObsPy website. This will need to be updated if the library updates
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
                # outfile = 'test.png'
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

    def create_data_synthetics_plot(self, stachan_traces, synthetics, mt_images_dir):
        """Create and save data vs. 
        synthetic waveform plots. 
        Equivalent to Dreger's function 
        mt_plot in mt_plot6iso2_linux2.c:

        mt_plot(ss,gg,nsta,Strike,Rake,Dip,
                St2,Rk2,Dp2,d_mt,Pdc,Pclvd,Piso,Mo, Mw, E, VR);

        There should be as many TRZ 
        plots as stations
        """
        # Init figure
        my_plot = plt.figure(figsize=(9, 8.5), dpi=100)
        my_plot.subplots_adjust(hspace=0.05, wspace=0.02,
                                bottom=0.02, top=0.96,
                                left=0.07, right=0.98)
        # Immutable rotated components
        # The order is important!
        rotated_components = ('T', 'R', 'Z')

        if self.verbosity > 1:
            print "Synthetics (Green's functions) used in calculations:"
            print synthetics

        # Get a list of the stations from the data
        stations = []
        for k in sorted(stachan_traces.iterkeys()):
            sta, chan = k.split('_')
            if sta not in stations:
                stations.append(sta)

        rows = len(stations)
        cols = len(rotated_components)
        no_of_subs = len(stations) * len(rotated_components)
        axis_num = 0
        my_plot.text(0.22, 0.985, 'Tangential', ha='center', va='top')
        my_plot.text(0.525, 0.985, 'Radial', ha='center', va='top')
        my_plot.text(0.83, 0.985, 'Vertical', ha='center', va='top')

        for i, s in enumerate(stations):
            for j, rc in enumerate(rotated_components):
                axis_num += 1
                new_sta_chan = s+'_'+rc
                if new_sta_chan in stachan_traces:
                    # Only plotting the first synthetic matrix
                    simple_synthetic = [(v) for k, v in synthetics[rc][0].iteritems()]
                    ax = plt.subplot(len(rotated_components), len(stations), axis_num)
                    plt.plot(stachan_traces[new_sta_chan], 'b', simple_synthetic, 'r--')
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])
                    if j == 0:
                        ax.set_ylabel(s, rotation='horizontal')

                else:
                    print "Trace (%s) is not in stachan_traces matrix" % new_sta_chan
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

class GreensFuncs():
    """Methods for working with the
    Greens functions, whether hard
    coded or via a db"""
    # {{{

    def __init__(self, green_db, green_pf, verbosity=0):
        """Initialize"""
        self.green_db = green_db
        self.green_pf = green_pf
        self.verbosity = verbosity

    def retrieve_file(self, dbptr):
        """Open Green's functions db and return
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
                logmt(1, 'Cannot open up Greens functions parameter file %s' % self.green_pf)
                return False
        ddist = stock.pfget_string(self.green_pf, 'ddist')
        dazim = stock.pfget_string(self.green_pf, 'dazim')
        ddip = stock.pfget_string(self.green_pf, 'ddip')
        ddepth = stock.pfget_string(self.green_pf, 'ddepth')

        if not os.path.isfile(self.green_db):
            logmt(3, 'Database (%s) does not exist' % self.green_db)
        gdb = antdb.dbopen(green_db, 'r')
        gdb.lookup(table = 'moment_tensor_greensfuncs')
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
                logmt(1, 'No Greens function found for %s' % sta)
            else:
                subgdb[3] = 0
                file = subgdb.extfile()       
                tmp[sta] = file
        return tmp

    def parse_data_file(self, green_file):
        """Hard coded way of getting synthetics (Green's
        function) from a file, not db. Remove before production,
        unless we can't get Green's functions auto-generated!!!
        """
        this_dict = defaultdict(lambda: defaultdict(defaultdict))

        if self.verbosity > 0:
            logmt(1, 'Constructing Greens function matrix from hard-coded file (%s)' % green_file)
        # !!! NOTE: Container list for parsed file object values. RLN (2011-11-03)
        greens = []
        # !!! NOTE: Open up the Green function file as specified in the pf (for now). RLN (2011-11-03)
        green = open(green_file, 'r')
        '''
        nl = number of lines
        dt = sample rate
        t1 = time period
        '''
        nl, dt, t1 = (green.readline()).split() # !!! NOTE: Only in PGC Greens functions, there are three starting numbers. RLN (2011-11-03)
        for line in green:
            one_line = line.rstrip('\n')
            vals = one_line.split()
            for j in vals:
                greens.append(float(j))
        green.close()
        tarr = greens.pop() # !!! NOTE: Only in PGC Greens functions, there is a trailing float on the last line. RLN (2011-11-03)
        if self.verbosity > 0:
            logmt(1, 'nl: %s, dt: %s, t1: %s, tarr: %s' % (nl, dt, t1, tarr))

        # !!! NOTE: Just force into two dimensional matrix. RLN (2011-11-03)
        #           'blockettes': separate chunks (components) of the Greens functions
        #           In the PGC: blockettes = 12, number of elements (nl) = 101
        #           In Dreger:  blockettes = 8, number of elements (nl) = 200
        blockettes = len(greens) / int(nl)
        for i in range(3):
            for j in range(blockettes):
                for k in range(int(nl)):
                    this_dict[j][i][k] = float(greens[k + j*int(nl)])

        '''
        !!! NOTE: The section below is Gert-Jan's original code - which takes into account
                  if the Green's function is isotropic, which you define in dbmoment.pf.
                  He forces the matrix into a certain format if it is isotropic,
                  but this does not play nice with the PGC Green's functions which
                  already have the isotropic components in them, but in a DIFFERENT
                  part of the matrix. I am pretty sure this has knock on effects in
                  the code, but I am not sure where or how yet. RLN (2011-11-03)
        '''

        '''
        # !!! NOTE: Open up the Green function file
        green = open('db/data/green', 'r')
        for line in green:
            # !!! NOTE: Each value is a float of twelve (12) significant figures with 
            #           no spaces, so split line into units of twelve. RLN (2011-11-03)
            for j in range(len(line)/12):
                # !!! NOTE: Append one value at a time. RLN (2011-11-03)
                greens.append(line[j*12:j*12+12])
        green.close()

        for i in range(3):
            for k in range(8):
                for j in range(200):
                    this_dict[k][i][j] = float(greens[j + k*200])
                    if k in [5,6,7]:
                        this_dict[k][i][j] *= -1 # Note the vertical GF's are flipped in earqt1.f and TW's Blackbox.f DVH conv. z + down
            if self.isoflag == 5:
                for k in [8,9]:
                    for j in range(200):
                        this_dict[k][i][j] = 0.0
            if self.isoflag == 6:
                for k in [8,9]:
                    for j in range(200):
                        this_dict[k][i][j] = float(greens[j + k*200])
        '''
        return (nl, this_dict)
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

    # !!! NOTE: Instance of Event. RLN (2011-11-21)
    my_event = Event(orid, event_db, verbosity)
    evdbptr, evparams, filter_string = my_event.extract_data(mag_filters)

    # !!! NOTE: Instance of GreensFuncs. RLN (2011-11-21)
    green_db = '0' # Holder for now as we don't have an active Greens function db
    my_greens = GreensFuncs(green_db, green_pf, verbosity)
    greens_file = my_greens.retrieve_file(evdbptr)
    greens_file = 'greens/1419_008_PDS1.gn'
    nl, gg = my_greens.parse_data_file(greens_file)
    nl = int(nl)

    stachan_traces, ev2sta_azimuths = my_event.get_chan_data(evdbptr, wave_db, chan_to_use, filter_string, statmax, nl, clip_values)

    if verbosity > 1:
        print 'Stachan traces:'
        pprint(stachan_traces)
        print 'Event to station azimuths:'
        pprint(ev2sta_azimuths)
        print 'Greens Function:'
        pprint(gg)

    # !!! NOTE: Works up to here. RLN (2011-11-29)

    """
    Opens Green's database and returns a 
    dictionary which the 10 component ascii series.
    """
    # !!! FIX: Read the Green's function from database based on distance and depth. RLN (2011-11-03)
    # greens_funcs = my_mt.get_greens_functions(evdbptr)
    #   !!! NOTE: green_funcs dict is not currently 
    #   used anywhere, so why is it here? Possibly relates to 
    #   the commented out dict above? RLN (2011-07-28)
    # nl, gg = my_mt.get_greens_functions_hard_coded(g)

    if len(gg[1][1]) != nl:
        logmt(3, "Size of Green's functions (%s) does not match the reported value (%s)" % (len(gg[1][1]), nl))

    """
    Pull the data from the db into stachan_traces
    then...
    Name them ss and gg after Dregers code
    """
    # !!! FIX: Why only 24 values for T,R,Z from the Texas 
    #          event when pulling from the db?
    #          It needs to match the length of the Green's 
    #          function to ensure we are comparing apples to apples
    #          RLN (2011-11-03)

    # !!! NOTE: Now we instantiate MomentTensor. RLN (2011-11-21)
    my_mt = MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
    ss = my_mt.construct_data_matrix(stachan_traces)

    if len(ss) != 0:
        print 'Data matrix S created --> %s stations used' % len(ss)
    if len(gg[1]) != len(ss['T']):
        print "Number of Green's functions (%d) does not match the number of stations (%d)" % (len(gg[1]), len(ss['T']))
    else:
        print 'Matrix S (data) and G (synthetics) created'

    '''
    !!! NOTE: The assumption here is that the size of both
    ss and gg is the same as the number of stations that are used
    in determining the moment tensor. In the hard-coded example
    this is three. RLN (2011-08-24)
    '''

    # !!! NOTE: timeshift is a list? RLN (2011-11-03)
    timeshift = my_mt.get_time_shift(ss, gg)
    if len(timeshift) != len(ss['T']):
        print 'Error in correlating synthetic and measured data'
    else:
        print 'Timeshift between data and synthetics computed for %s stations' % len(timeshift)

    # !!! NOTE: Where does this come from without explanation? Dregers code? RLN (2011-07-28)
    # Examples just for the examples
    # - event to station azimuth, in this example used 3 stations
    az = [10, 40, 50]
    ev2sta_azimuths

    # !!! NOTE: What is W? Weighting factor - number of stations, long list of values.
    # If 1 then station fully incorporated into solution. If less than one, don't include. RLN (2011-07-28)

    # Allocate Memory for Station Weights
    W = []
    for i in range(len(ss['T'])):
        W.append(1.0)

    # !!! NOTE: INVERSION ROUTINE
    #     Dreger normalizes AtA (B?) and AIV (AIV) matrix. We don't need to - just create default dictionary
    #     RLN (2011-08-24)
    AIV = defaultdict(dict)
    B = defaultdict(dict) 

    AIV, B = my_mt.matrix_AIV_B(ss, gg, az, timeshift)

    # !!! NOTE: DETERMINE SOLUTION VECTOR. RLN (2011-08-24)
    M = my_mt.determine_solution_vector(AIV, B)

    gfscale = -1.0e+20
    gfscale = -1.0
    strike = []
    dip = []
    rake = []

    # !!! NOTE: DECOMPOSE MOMENT TENSOR INTO VARIOUS REPRESENTATIONS. RLN (2011-11-03)
    m0, Mw, strike, dip, rake, pcdc, pcclvd, pciso = my_mt.decompose_moment_tensor(M)
    E, VR, VAR, svar, sdpower = my_mt.fitcheck(gg, ss, W, M, m0, timeshift, az)

    qlt = my_mt.quality_check(VR)

    evdbptr.close()

    # !!! NOTE: UPDATED TABLES AND CREATE FOCAL MECHANISM IMAGE - DONE. RLN (2011-11-22)
    my_event.update_moment_tbl(strike, dip, rake)
    focalmech_img = my_event.create_focal_mechanism(obspy_beachball, mt_images_dir, strike, dip, rake)

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
    synthetics_img = my_event.create_data_synthetics_plot(stachan_traces, gg, mt_images_dir)

    # !!! NOTE: CREATE THE COMPOSITE (FINAL) IMAGE AND UPDATE DB. RLN (2011-11-28)
    my_event.create_composite_plot(ttfont, image_annotations, mt_images_dir, synthetics_img, focalmech_img)

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

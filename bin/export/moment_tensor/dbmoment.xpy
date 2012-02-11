"""
dbmoment.py

    Calculates the moment tensor for a given  origin-id  (orid) for a certain amount of stations 
    (defined in the parameter file). The accompanying Green's functions are either  read  from  a
    pre-constructed wfdisc table in a database or generated and stored in a database wfdisc table. 
    It relies significantly on the  Antelope  Python Interface (Datascope and Stock), NumPy, 
    Matplotlib, ObsPy and SciPy.

        dbmoment [-vV] [-p pfname] orid

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

#{{{ Import libraries

import os
import re
import sys
import math as math
from pprint import pprint
from time import gmtime, time
from datetime import datetime
from optparse import OptionParser
from collections import defaultdict

# ANTELOPE
try:
    import antelope.stock as stock
    import antelope.datascope as datascope
except Exception,e:
    sys.exit("Antelope Import Error: [%s] => [%s]" % (Exception,e))

# NUMPY
try:
    from PIL import Image, ImageDraw, ImageFont
except Exception,e:
    sys.exit("Import Error: [%s] => [%s] Do you have PIL installed correctly?" % (Exception,e))

# NUMPY
try:
    import numpy as np
except Exception,e:
    sys.exit("Import Error: [%s] => [%s] Do you have NumPy installed correctly?" % (Exception,e))

# MATPLOTLIB
try:
    import matplotlib as mpl
except Exception,e:
    sys.exit("Import Error: [%s] => [%s] Do you have MatplotLib installed correctly?" % (Exception,e))
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

# OBSPY
try:
    from obspy.imaging.beachball import Beachball
except Exception,e:
    sys.exit("Import Error: [%s] => [%s] Do you have ObsPy installed correctly?" % (Exception,e))

# Needed for fkrprog.py
from scipy.integrate import cumtrapz
from numpy.fft  import ifft
import matplotlib.pyplot as plt
from math import exp, log, sqrt, acos, asin, cos, sin
from cmath import log as clog
from cmath import sqrt as csqrt
from cmath import exp as cexp

#}}}

# {{{ Global functions
def configure():
    """ Function to configure parameters from command-line and the pf-file
    Check command line options Use Antelope's built-in PFPATH to
    determine the paths to search for the parameter file
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

    gf_lib = stock.pfget_string(pfname, 'gf_lib')
    inv_lib = stock.pfget_string(pfname, 'inv_lib')
    data_lib = stock.pfget_string(pfname, 'data_lib')

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

    try:
        green_db = stock.pfget_string(pfname, 'green_db')
    except:
        green_db = wave_db

    try:
        resp_db = stock.pfget_string(pfname, 'resp_db')
    except:
        resp_db = wave_db

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
            gf_lib,inv_lib, data_lib, wave_db, green_db, resp_db, 
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

def main():
#{{{
    """Get the moment tensor solution
    for a particular origin
    """
    # Parse command line
    pf, orid, verbose, debug = configure()
    verbosity = determine_verbosity(verbose, debug)
    # Parse the parameter file
    (chan_to_use, trim_value, isoflag, model_name,
     model_type, statmax, clip_values, use_inc, event_db, 
     gf_lib,inv_lib, data_lib, wave_db, green_db, resp_db, 
     mag_filters, mt_images_dir, 
     ttfont, obspy_beachball, 
     distance_weighting) = parse_pf(pf, verbosity)


    # FKRPROG - dynamic upload based on pf file value
    try:
        exec "import moment_tensor.%s as fkr" % gf_lib
    except Exception,e:
        sys.exit("Import Error: [%s] => [%s]" % (Exception,e))

    # Inversion - dynamic upload based on pf file value
    try:
        exec "import moment_tensor.%s as inversion" % inv_lib
    except Exception,e:
        sys.exit("Import Error: [%s] => [%s]" % (Exception,e))

    # Event Data - dynamic upload based on pf file value
    try:
        exec "import moment_tensor.%s as data" % data_lib
    except Exception,e:
        sys.exit("Import Error: [%s] => [%s]" % (Exception,e))


    # !!! NOTE: Instantiate Event. RLN (2011-11-21)
    my_event = data.Event(orid, event_db, verbosity)
    my_mt = inversion.MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
    evdbptr, evparams, filter_string, filter_timepad = my_event.extract_data(mag_filters)

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
                # Not sure we still need to define 200, but leave for now
                real_data = my_event.get_chan_data(wave_db, chan_to_use, esaz, arrival_time, filter_string, filter_timepad, stacode, 200, clip_values)
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
                # filtered_synthetic_data = filter_data(synthetics, 'dict', filter_string, verbosity)
                max_xcor, timeshift = my_mt.get_time_shift(real_data, synthetics, stacode, 120)
                ss.append(real_data)
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

    if verbosity > 1:
        print 'Stachan traces:'
        pprint(ss)
        print 'Event to station details:'
        pprint(ev2sta)
        print 'Greens Functions:'
        pprint(gg)

    # !!! NOTE: Instantiate MomentTensor. RLN (2011-11-21)
    my_mt = inversion.MomentTensor(distance_weighting, isoflag, trim_value, verbosity)
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
    !!! NOTE: New weighting factor based on timeshift:
              If timeshift < 50, include station in
              solution. Else reject.
              RLN (2011-12-02)
    '''
    # print ss
    # pprint(gg)
    cleaned_ss = []
    cleaned_gg = []
    cleaned_ev2sta = []
    for i in range(len(ss)):
        print "\n\n**** %s ***" % ev2sta[i][0]
        if abs(ev2sta[i][3]) > 50:
            print "**** Timeshift: %s. Ignoring! ***" % ev2sta[i][3]
        else:
            cleaned_ss.append(ss[i])
            cleaned_gg.append(gg[i])
            cleaned_ev2sta.append(ev2sta[i])

    # print cleaned_ss
    # print cleaned_gg
    # print cleaned_ev2sta
    '''
    !!! NOTE: INVERSION ROUTINE
              Dreger normalizes AtA (B?) and AIV (AIV) matrix. 
              We don't need to - just create default dictionary
              RLN (2011-08-24)
    '''

    AIV = defaultdict(dict)
    B = defaultdict(dict) 

    # AIV, B = my_mt.matrix_AIV_B(ss, gg, ev2sta_azimuths, timeshift, nl)
    AIV, B = my_mt.matrix_AIV_B(cleaned_ss, cleaned_gg, cleaned_ev2sta, nl)

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
    # E, VR, VAR, svar, sdpower = my_mt.fitcheck(cleaned_gg, cleaned_ss, W, M, m0, cleaned_ev2sta, nl)
    E, VR, VAR, svar, sdpower = my_mt.fitcheck(cleaned_gg, cleaned_ss, M, m0, cleaned_ev2sta, nl)

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
    mod_gg = my_event.calculate_synthetics_to_plot(cleaned_gg, cleaned_ev2sta, M, nl)
    # pprint(mod_gg)

    synthetics_img = my_event.create_data_synthetics_plot(cleaned_ss, cleaned_ev2sta, mod_gg, nl, mt_images_dir)

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
        print 'Station (%s) = %s  %s' % (svar[j][0], svar[j][1], sdpower[j])
    print 'Var     = %s' % E
    print 'VR      = %s (UNWEIGHTED)' % VR
    print 'VR      = %s (WEIGHTED)' % VR

    print 'Var/Pdc = ', E/pcdc
    print 'Quality = %s' % qlt

    return 0
#}}}

if __name__ == '__main__':
    sys.exit(main())

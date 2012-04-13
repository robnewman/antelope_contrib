"""
dbmoment.py

    Calculates the moment tensor for a given  origin-id  (orid) for a certain amount of stations 
    (defined in the parameter file). The accompanying Green's functions are either  read  from  a
    pre-constructed wfdisc table in a database or generated and stored in a database wfdisc table. 
    It relies significantly on the  Antelope  Python Interface (Datascope and Stock), NumPy, 
    Matplotlib, ObsPy and SciPy.
    UPDATE: Antelope 5.2-64 contains the library PyLab. 

        dbmoment [-vV] [-p pfname] orid

    @authors  
            Gert-Jan van den Hazel <hazelvd@knmi.nl>
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
from math import exp, log, sqrt, acos, asin, cos, sin
from cmath import log as clog
from cmath import sqrt as csqrt
from cmath import exp as cexp

# ANTELOPE
try:
    import antelope.stock as stock
    import antelope.datascope as datascope
except Exception,e:
    sys.exit("Import Error: [%s] Do you have ANTELOPE installed correctly?" % e)

# PYLAB
try:
    import pylab as pylab
    from pylab  import ifft, trapz, linalg, array, zeros, matplotlib
    from matplotlib  import pyplot
except Exception,e:
    sys.exit("Import Error: [%s] Do you have PYLAB installed correctly?" % e)

# NUMPY
try:
    from PIL import Image, ImageDraw, ImageFont
except Exception,e:
    sys.exit("Import Error: [%s] Do you have PIL installed correctly?" % e)

# OBSPY
try:
    from obspy.imaging.beachball import Beachball
except Exception,e:
    sys.exit("Import Error: [%s] Do you have ObsPy installed correctly?" % e)

#}}}

def log(message):
    """Global function to handle 
    log output. Prints messages 
    with a timestamp.  
    """
#{{{

    if not options.verbose: return
    curtime = stock.epoch2str( time(), '%m/%d/%Y %H:%M:%S')

    print '%s dbmoment: %s' % (curtime,message)

#}}}

def debug(message):
    """Global function to handle 
    debug output. Just a wrapper
    for the log function.
    """
#{{{

    if options.debug: log(message)

#}}}

class DbMoment():
    """
    Main class for calculating moment tensors of events
    """
#{{{

    def __init__(self,pf):
#{{{

        log( "---------------------------------------" )
        log( "DBMOMENT.PY: DbMoment.__init__(%s)" % (pf) )

        self.pf = pf

        try:
            self._parse_pf()
        except Exception,e:
            sys.exit('ERROR: problem during parsing of pf file(%s) class.[%s]' % (self.pf,e) )

#}}}

    def _parse_pf(self):
        """Parse the parameter file
        and assign to vars that will be
        used throughout the script
        """
#{{{

        log( "---------------------------------------" )
        log( "DBMOMENT.PY: parse parameter file(%s)" % (self.pf) )


        self.gf_lib = stock.pfget_string(self.pf, 'gf_lib')
        log('\tgf_lib => %s' % self.gf_lib)

        self.inv_lib = stock.pfget_string(self.pf, 'inv_lib')
        log('\tinv_lib => %s' % self.inv_lib)

        self.data_lib = stock.pfget_string(self.pf, 'data_lib')
        log('\tdata_lib => %s' % self.data_lib)

        try:
            self.chan_to_use = stock.pfget_string(self.pf, 'chan_to_use')
        except:
            self.chan_to_use = 'LH.*'
        log('\tchan_to_use => %s' % self.chan_to_use)

        self.trim_value = stock.pfget_string(self.pf, 'trim_value')
        log('\ttrim_value => %s' % self.trim_value)

        if stock.pfget_int(self.pf, 'isoflag') == 1: 
            self.isoflag = 6
        else:
            self.isoflag = 5
        log('\tisoflag => %s' % self.isoflag)

        self.model_name = stock.pfget_string(self.pf, 'model_name')
        log('\tmodel_name => %s' % self.model_name)

        self.model_type = stock.pfget_string(self.pf, 'model_type')
        log('\tmodel_type => %s' % self.model_type)

        self.statmax = stock.pfget_int(self.pf,'statmax')
        log('\tstatmax => %s' % self.statmax)

        self.clip_values = stock.pfget_arr(self.pf, 'clip_values')
        log('\tclip_values => %s' % self.clip_values)

        self.use_inc = stock.pfget_string(self.pf, 'use_inc')
        log('\tuse_inc => %s' % self.use_inc)

        self.event_db = stock.pfget_string(self.pf,'event_db')
        log('\tevent_db => %s' % self.event_db)

        try:
            self.wave_db = stock.pfget_string(self.self.pf, 'wave_db')
        except:
            self.wave_db = self.event_db
        log('\twave_db => %s' % self.wave_db)

        try:
            self.resp_db = stock.pfget_string(self.pf, 'resp_db')
        except:
            self.resp_db = self.event_db
        log('\tresp_db => %s' % self.resp_db)

        self.mag_filters = stock.pfget_arr(self.pf, 'mag_filters')
        [log('\tmag_filters => %s' % x) for x in self.mag_filters]

        self.mt_images_dir = stock.pfget_string(self.pf, 'mt_images_dir')
        log('\tmt_images_dir => %s' % self.mt_images_dir)

        self.ttfont = stock.pfget_string(self.pf, 'ttfont')
        log('\tttfont => %s' % self.ttfont)

        self.obspy_beachball = stock.pfget_arr(self.pf, 'obspy_beachball')
        log('\tobspy_beachball => %s' % self.obspy_beachball)

        if stock.pfget_string(self.pf, 'distance_weighting') == 'on':
            self.distance_weighting = True
        else:
            self.distance_weighting = False
        log('\tdistance_weighting => %s' % self.distance_weighting)

#}}}

    def mt(self,orid,event_table=False,select='',reject=''):
        """ Get the moment tensor solution for a 
        particular origin. If event is True
        then we get the origin from the 
        event table.
        """
#{{{
        log( "---------------------------------------" )
        log( "DBMOMENT.PY: mt(%s,select=%s,reject=%s)" % (orid,select,reject) )

        # FKRPROG - dynamic upload based on pf file value
        try:
            log( "DBMOMENT.PY: Load module for GreenFunctions(%s)" % self.gf_lib )
            exec "import moment_tensor.%s as gf_lib" % self.gf_lib
        except Exception,e:
            sys.exit("Import Error: [%s]" % e)

        # Inversion - dynamic upload based on pf file value
        try:
            log( "DBMOMENT.PY: Load module for Inversion(%s)" % self.inv_lib )
            exec "import moment_tensor.%s as inv_lib" % self.inv_lib
        except Exception,e:
            sys.exit("Import Error: [%s]" % e)

        # Event Data - dynamic upload based on pf file value
        try:
            log( "DBMOMENT.PY: Load module for data_extraction(%s)" % self.data_lib )
            exec "import moment_tensor.%s as data_lib" % self.data_lib
        except Exception,e:
            sys.exit("Import Error: [%s]" % e)

        # Instantiate Classes
        my_event = data_lib.Event(self.chan_to_use,self.clip_values)

        # Instantiate Classes
        my_inv = inv_lib.MomentTensor(self.distance_weighting, self.isoflag, self.trim_value)

        # Instantiate Classes
        green = gf_lib.GreenFunctions(self.model_name)


        # Get waveforms from database
        log( "DBMOMENT.PY: Get event [%s]" % orid )
        if not my_event.subset_event(orid, self.event_db, event_table, self.mag_filters,select,reject):
            sys.exit('Problem creating subset for event [%s].' % orid)

        '''
        !!! NOTE: 
        Go through each of the eight groups until we get good 
        cross-correlation results with the synthetics
        (that we have to generate for each distance & depth)
        Ideally we need two per quad for a total of 8 stations.

        '''
        ss = []
        gg = []
        ev2sta = []
        stations_per_group = int(round(self.statmax/8))
        sta_list = my_event.sta_list

        log("\tIterate over groups of stations in each azimuthal range. Requesting %s per group" % stations_per_group)

        #
        # Start loop over each quadrant
        #
        for grp in sta_list:

            log("\t\tWorking on group (%s) " % grp)
            log("\t\t(%s) [%s] " % (grp,sta_list[grp]))
            good_cross_corr_stations = 0

            if not len(sta_list[grp]):
                print "\nERROR: No stations for quadrant[%s].\n" % grp
                continue

            #
            # Start loop over each station in the quadrant
            #
            log("\t\tStart loop over quadrant [%s]." % grp)
            for sta in sta_list[grp]:
            #{{{
                log("\t\t\tWorking on sta (%s) already(%s) " % (sta,good_cross_corr_stations))
                if good_cross_corr_stations >= stations_per_group:
                    log("\t\t\tFound enough good matches for quadrant (%s)" % grp)
                    break
                else:
                    stacode, esaz, depth, distance, arrival_time = sta

                    # Get new Green's Functions for this distance, depth
                    # Clean depth and distance to create clean green's functions
                    depth = int(depth)
                    distance = int(distance)
                    log("\t\t\tGenerate GF for this depth (%s km) & distance (%s km)" % (depth, distance))
                    green.build(depth, distance, 1, 'v', my_event.filter_string)
                    synthetics = green.ALL
                    '''
                    !!! NOTE: 
                    Dynamic creation of Greens Functions
                    #  Structure of ours vs Dreger:
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

                    Create for every stations distance to the event.  RLN (2011-12-08)
                    '''


                    # Plot greens functions fundamental traces.  
                    if options.debug: 
                        pprint(synthetics)
                        green.plot()
                        log('Greens Functions:')

                    '''
                    !!! NOTE: 
                    Since we hardcoded 1 sps GF we can use the value of the total
                    samples in the trace to calculate the time range of the 
                    requested data object. JCR 
                    '''
                    # Get waveform data for station
                    log("\t\t\tGet data for (%s)" % stacode)
                    real_data = my_event.get_chan_data(self.wave_db, esaz, arrival_time, stacode, len(synthetics['TSS']) )

                    if options.debug: 
                        debug("\t\t\tPlot traces for (%s)" % stacode)
                        now = 0
                        for trace in real_data:
                            try: 
                                now += 1
                                pyplot.subplot(len(real_data),1,now)
                                pyplot.plot(real_data[trace])
                                pyplot.legend([trace])
                            except Exception,e:
                                sys.exit('ERROR: plot real data: %s . Exception: %s' % (trace, e))


                        pyplot.suptitle("Real Data: Trace %s" % stacode)
                        pyplot.show()



                    # Get correlation from data to greenfuncs
                    delta = int(arrival_time - my_event.event_time)
                    max_xcor, timeshift = my_inv.get_time_shift(real_data, synthetics, delta)
                    
                    if max_xcor < 25: 
                        log("\t\t\tREJECT: %s max-correlation:[%s] and timeshift:[%s]" % (stacode,max_xcor,timeshift))
                        continue

                    log("\t\t\tUSE: %s max-correlation:[%s] and timeshift:[%s]" % (stacode,max_xcor,timeshift))

                    ss.append(real_data)
                    gg.append(synthetics)
                    ev2sta.append((stacode, esaz, distance, timeshift))
                    good_cross_corr_stations += 1
                #}}}

            debug('\tEvent to station details:[%s]' % ev2sta)

        '''
        !!! FIX: Override the number of data points 
                like Dreger does to be 120 values.
                HACK!
                RLN (2011-12-07)
        '''
        # limit the samples
        #nl = 120
        nl = int(len(gg[0]['TSS'])/2)

        '''
        !!! NOTE: INVERSION ROUTINE
                Dreger normalizes AtA (B?) and AIV (AIV) matrix. 
                We don't need to - just create default dictionary
                RLN (2011-08-24)
        '''
        AIV = defaultdict(dict)
        B = defaultdict(dict) 

        AIV, B = my_inv.matrix_AIV_B(ss, gg, ev2sta, nl)

        M = my_inv.determine_solution_vector(AIV, B)

        gfscale = -1.0e+20
        gfscale = -1.0
        strike = []
        dip = []
        rake = []

        # !!! NOTE: DECOMPOSE MOMENT TENSOR INTO VARIOUS REPRESENTATIONS. RLN (2011-11-03)
        m0, Mw, strike, dip, rake, pcdc, pcclvd, pciso = my_inv.decompose_moment_tensor(M)
        # E, VR, VAR, svar, sdpower = my_inv.fitcheck(gg, ss, W, M, m0, timeshift, ev2sta_azimuths, nl)
        # E, VR, VAR, svar, sdpower = my_inv.fitcheck(cleaned_gg, cleaned_ss, W, M, m0, cleaned_ev2sta, nl)
        #E, VR, VAR, svar, sdpower = my_inv.fitcheck(cleaned_gg, cleaned_ss, M, m0, cleaned_ev2sta, nl)
        E, VR, VAR, svar, sdpower = my_inv.fitcheck(gg, ss, M, m0, ev2sta, 120)

        qlt = my_inv.quality_check(VR)

        #try:
        #    evdbptr.close()
        #except Exception, e:
        #    log('Error closing event database. Already closed? Exception: %s' % e)

        '''
        !!! NOTE: In ObsPy there are two different 
                ways to create a focal mechanism.
                Allow both.  RLN (2011-12-02)
        '''
        my_event.update_moment_tbl(strike, dip, rake)

        print 'my_event.create_focal_mechanism()'
        #focalmech_img = my_event.create_focal_mechanism(self.obspy_beachball, self.mt_images_dir, M, False, False, False)
        focalmech_img = my_event.create_focal_mechanism(self.obspy_beachball, self.mt_images_dir, M, strike, dip, rake)

        print 'build image_annotations'
        # !!! NOTE: Use a list of tuples as order is important. RLN (2011-11-28)
        image_annotations = [ ('orid', '%s' % orid), 
                            ('time', '%s' % stock.epoch2str(my_event.event_time,'%m/%d/%Y %H:%M:%S')), 
                            ('M00', '%0d' % M[0,0]), 
                            ('M11', '%0d' % M[1,1]), 
                            ('M22', '%0d' % M[2,2]), 
                            ('M01', '%0d' % M[0,1]), 
                            ('M02', '%0d' % M[0,2]), 
                            ('M12', '%0d' % M[1,2]), 
                            ('1 Strike ', '%0d' % math.degrees(strike[0])), 
                            ('1 Dip', '%0d' % math.degrees(dip[0])), 
                            ('1 Rake', '%0d' % math.degrees(rake[0])), 
                            ('2 Strike', '%0d' % math.degrees(strike[1])), 
                            ('2 Dip', '%0d' % math.degrees(dip[1])), 
                            ('2 Rake', '%0d' % math.degrees(rake[1])), 
                            ('Mo', '%0.3E' % m0), 
                            ('Mw', '%0.3f' % Mw), 
                            ('% DC', '%0.3f' % pcdc), 
                            ('% CLVD', '%0.3f' % pcclvd), 
                            ('% ISO', '%0.3f' % pciso), 
                            ('VR', '%0.3E' % VR), 
                            ('VAR', '%0.3E' % VAR)
                            ]
        print 'my_event.calculate_synthetics_to_plot()'
        # !!! FIX: CREATE DATA vs. SYNTHETICS PLOTS - DONE. RLN (2011-11-03)
        #mod_gg = my_event.calculate_synthetics_to_plot(cleaned_gg, cleaned_ev2sta, M, nl)
        mod_gg = my_event.calculate_synthetics_to_plot(gg, ev2sta, M, nl)
        # pprint(mod_gg)

        #synthetics_img = my_event.create_data_synthetics_plot(cleaned_ss, cleaned_ev2sta, mod_gg, nl, mt_images_dir)
        synthetics_img = my_event.create_data_synthetics_plot(ss, ev2sta, mod_gg, nl, self.mt_images_dir)

        # !!! NOTE: CREATE THE COMPOSITE (FINAL) IMAGE AND UPDATE DB. RLN (2011-11-28)
        my_event.create_composite_plot(self.ttfont, image_annotations, self.mt_images_dir, synthetics_img, focalmech_img)

        #log("MOMENT TENSOR:")
        #log(M)
        #log('M0      = %s' % m0)
        #log('Mw      = %s' % Mw )
        #log('Strike  = %s' % math.degrees(strike[0]))
        #log('Rake    = %s' % math.degrees(rake[0]))
        #log('Dip     = %s' % math.degrees(dip[0]))
        #log('Strike2 = %s' % math.degrees(strike[1]))
        #log('Rake2   = %s' % math.degrees(rake[1]))
        #log('Dip2    = %s' % math.degrees(dip[1]))
        #log('Pdc     = %s' % pcdc)
        #log('Pclvd   = %s' % pcclvd)
        #log('Piso    = %s' % pciso)

        # !!! NOTE: From Dreger's fitcheck fprintf() commands. RLN (2011-08-24)
        for j in range(len(svar)):
            log('Station (%s) = %s  %s' % (svar[j][0], svar[j][1], sdpower[j]))
        log('Var     = %s' % E)
        log('VR      = %s (UNWEIGHTED)' % VR)
        log('VR      = %s (WEIGHTED)' % VR)

        log('Var/Pdc = %s' % (E/pcdc))
        log('Quality = %s' % qlt)

        return 0
#}}}
#}}}

if __name__ == '__main__':
    """ Configure parameters from command-line and the pf-file
    Check command line options Use Antelope's built-in PFPATH to
    determine the paths to search for the parameter file
    """
#{{{
    usage = "Usage: dbmoment [-vVde] [-p pfname] orid"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="verbose output")
    parser.add_option("-V", "--veryverbose", action="store_true", dest="veryverbose", default=False, help="debug application")
    parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="debug application")
    parser.add_option("-p", "--pf", action="store", dest="pf", type="string", help="parameter file path")
    parser.add_option("-e", "--event", action="store", default=False, dest="event", type="string", help="use event table")
    parser.add_option("-s", "--select", action="store", default='', dest="select", type="string", help="only select these stations")
    parser.add_option("-r", "--reject", action="store", default='', dest="reject", type="string", help="reject these stations")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit( usage );
    else:
        orid = args[0]

    if not options.pf: options.pf = 'dbmoment'

    if options.veryverbose: options.debug = True
    if options.debug: options.verbose = True


    try: 
        options.pf = stock.pffiles(options.pf)[-1]
    except Exception,e:
        sys.exit('ERROR: problem loading pf(%s) class.[%s => %s]' % (options.pf,Exception,e) )

    log("Parameter file to use [%s]" % options.pf)

    if not os.path.isfile(options.pf): sys.exit('ERROR: Cannot find pf(%s)' % options.pf )

    try:
        log("Loading module [ DbMoment ]" )
        dbmnt = DbMoment(options.pf)
    except Exception,e:
        sys.exit('ERROR: problem loading main DbMoment(%s) class.[%s => %s]' % (options.pf,Exception,e) )

    try:
        log("Start calculation [ %s ]" % orid )
        dbmnt.mt(orid,options.event,options.select,options.reject)
    except Exception,e:
        sys.exit('ERROR: problem during calculation of Moment Tensor: mt(%s,%s) class.[%s => %s]' % (orid,options.event,Exception,e) )

    sys.exit()

else:
    sys.exit('ERROR: Cannot import class DbMoment into your code!!!')
#}}}

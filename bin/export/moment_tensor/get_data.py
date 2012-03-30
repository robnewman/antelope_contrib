from __main__ import *      # _et all the libraries from parent

class Event():
    # {{{
    """Class for extracting per event data 
    from the database and creating the plots.
    """

    def __init__(self, chan_to_use, clip_values, verbose=False, debug=False):
#{{{
        """Initialize"""
        self.verbose = verbose
        self.debug = debug
        self.clip_values = clip_values
        self.chan_to_use = chan_to_use

        if self.debug: self.verbose = True
#}}}

    def _log(self,message):
    #{{{
        """Global function to handle log output. Prints messages with a timestamp.  """

        curtime = stock.epoch2str( time(), '%m/%d/%Y %H:%M:%S')
        
        print '%s dbmoment: \tEvent() => %s' % (curtime,message)
    #}}}

    def subset_event(self, orid, event_db, event_table, filters,select='',reject=''):
    #{{{
        """Open event database and get 
        event data for given origin-id
        Returns event parameters or 
        exit when an error occurs
        """

        if self.debug: self._log('Get canned view from database')

        self.orid = orid
        self.event_db = event_db

        if not os.path.isfile(event_db):
            sys.exit('Database (%s) does not exist' % event_db)

        # Open the database
        try:
            evdb = datascope.dbopen(event_db, 'r')
        except Exception, e:
            sys.exit('Could not open database (%s). Exception: %s' % (event_db, e))

        if event_table:
            # For events from event table
            if self.debug: self._log('open event table and join with origin ')
            try:
                evdb.lookup(table='event')
                evdb.join('origin')
            except Exception, e:
                sys.exit('Could not open(event)/join(origin) for (%s). Exception: %s' % (event_db, e))

            if self.debug: self._log('subset(event==%s' % orid)
            try:
                evdb.subset('evid == %s' % orid)
            except Exception, e:
                sys.exit('Could not subset(evid == %s). Exception: %s' % (orid, e))

            if evdb.nrecs() == 0:
                sys.exit('evid (%s) does not exist after join' % orid)

        else:
            # For events from origin table
            if self.debug: self._log('open origin')
            try:
                if self.debug: self._log('open origin table')
                evdb.lookup(table='origin')
            except Exception, e:
                sys.exit('Could not open(origin) for (%s). Exception: %s' % (event_db, e))

            if self.debug: self._log('subset(orid==%s' % orid)
            try:
                evdb.subset('orid == %s' % orid)
            except Exception, e:
                sys.exit('Could not subset(orid == %s). Exception: %s' % (orid, e))

            if evdb.nrecs() == 0:
                sys.exit('Orid (%s) does not exist in origin table' % orid)

        # Join with netmag table
        if self.debug: self._log('join(netmag,outer=True)')
        try: 
            evdb.join('netmag', outer=True)
        except Exception, e:
            sys.exit('Could not join(netmag,outer=True). Exception: %s' % e)

        if evdb.nrecs() == 0:
            sys.exit('Could not join netmag table for orid (%s)' % orid)

        if evdb.nrecs() > 1:
            self._log('More than one entry (%s) in the netmag table for orid (%s)' % (evdb.nrecs(), orid))
            evdb[3] = 0 # Force to first record
            magid = evdb.getv('magid') 
            self._log('Subset to get unique entry. Using magid (%s).' % magid)
            evdb.subset('magid == %d' % magid)

        evdb[3] = 0 # Should only be one record now
        self.evparams = {} # Prep variable

        # Get all event data into evparams variable
        for f in evdb.query('dbTABLE_FIELDS'):
            try: 
                self.evparams[f] = evdb.getv(f)[0]
                if self.verbose: self._log('%s => %s' % (f,self.evparams[f]))
            except Exception, e:
                sys.exit('Could not find (%s) after origin netmag join for orid %s. Exception: %s' % (f,orid,e))

        #save event time into special variable
        self.event_time = self.evparams['time']

        if self.debug: self._log('join(assoc,arrival,site), sort(delta), subset(iphase=~P)')
        try:
            evdb.join('assoc')
            evdb.join('arrival')
            evdb.join('site')
            evdb.sort('delta')
            evdb.subset('iphase=~/.*P.*|.*p.*/')
        except Exception, e:
            sys.exit('in join (assoc,arrival,site) sort(delta), subset(iphase=~P). Exception: %s' % (f,orid,e))

        # select only stations listed on the -s flag at command-line
        if select:
            if self.debug: self._log('subset(sta =~ ==%s' % select)
            try:
                evdb.subset('sta =~ /%s/' % select)
            except Exception, e:
                sys.exit('Could not subset(sta =~/%s/). Exception: %s' % (select, e))

        # reject all stations listed on the -r flag at command-line
        if reject:
            if self.debug: self._log('subset(sta !~ ==%s' % reject)
            try:
                evdb.subset('sta !~ /%s/' % reject)
            except Exception, e:
                sys.exit('Could not subset(sta !~/%s/). Exception: %s' % (reject, e))

        '''
        !!! NOTE: This is regional moment
                    tensor. Will fail if stations
                    are closer than 1 deg (clipping)
                    or further away than 10 deg
                    (teleseismic) to the event.
                    RLN (2011-12-13)
        '''
        if self.debug: self._log('subset(delta>=1 && delta < 10)')
        try:
            evdb.subset('delta >= 1 && delta < 10') # REGIONAL moment tensor,
        except Exception, e:
            sys.exit('subset(delta>=1 && delta < 10). Exception: %s' % (f,orid,e))


        if evdb.nrecs() == 0:
            sys.exit('No arrivals for selected origin %s' % orid)

        if self.verbose: self._log('%d records that match this origin arrival and iphase' % evdb.nrecs())

        # Get magnitude filter options from pf values
        if self.debug: self._log('Parse magnitude filters')
        mag_filter_list = []

        for line in filters:
            if self.debug: self._log('Raw filter line [%s]' % line)
            splitline = line.split()
            # No upper limit for the last element
            if len(splitline) == 2:
                lower_magnitude, filter = splitline
                upper_magnitude = 99
            elif len(splitline) == 3:
                lower_magnitude, upper_magnitude, filter = splitline
            else:
                sys.exit('Filter format in pf file not recognized [%s]' % filters)

            mag_filter_list.append({'lowermag':lower_magnitude, 'uppermag':upper_magnitude, 'filter':filter})

            for mfl in mag_filter_list:
                if float(self.evparams['magnitude']) >= float(mfl['lowermag']) and float(self.evparams['magnitude']) < float(mfl['uppermag']):
                    # Parse the filter string
                    if self.debug: self._log('Parsed filter line [%s]' % mfl['filter'])
                    self.filter_string = (mfl['filter']).replace('_', ' ')

        if not self.filter_string:
            sys.exit('Magnitude %s not within magnitude bounds defined in the parameter file:' % self.evparams['magnitude'])

        """Build a list of all stations that
        recorded the event. Split into eight sections:

                 0
             NNW | NNE
           WNW \ | / ENE
        270 ----------- 90
           WSW / | \ ESE
             SSW | SSE  
                180

        """

        self.sta_list = {'NNE':[], 'ENE':[], 'ESE':[], 'SSE':[], 'SSW':[], 'WSW':[], 'WNW':[], 'NNW':[]}
        for i in range(evdb.nrecs()):
            evdb[3] = i
            sta, esaz, depth, at, ev_lat, ev_lon, site_lat, site_lon = evdb.getv('sta', 'esaz', 'depth', 'arrival.time', 'lat', 'lon', 'site.lat', 'site.lon')
            distance_deg = evdb.ex_eval('distance(%s, %s, %s, %s)' % (ev_lat, ev_lon, site_lat, site_lon))
            distance_km = evdb.ex_eval('deg2km(%s)' % (distance_deg))
            if esaz >= 0 and esaz < 45:
                self.sta_list['NNE'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station NNE [%s]' % sta)
            elif esaz >= 45 and esaz < 90:
                self.sta_list['ENE'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station ENE [%s]' % sta)
            elif esaz >= 90 and esaz < 135:
                self.sta_list['ESE'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station ESE [%s]' % sta)
            elif esaz >= 135 and esaz < 180:
                self.sta_list['SSE'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station SSE [%s]' % sta)
            elif esaz >= 180 and esaz < 225:
                self.sta_list['SSW'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station SSW [%s]' % sta)
            elif esaz >= 225 and esaz < 270:
                self.sta_list['WSW'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station WSW [%s]' % sta)
            elif esaz >= 270 and esaz < 315:
                self.sta_list['WNW'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station WNW [%s]' % sta)
            elif esaz >= 315 and esaz <= 360:
                self.sta_list['NNW'].append((sta, esaz, depth, distance_km, at))
                if self.debug: self._log('Add station NNW [%s]' % sta)

        evdb.close()

        if len(self.sta_list): 
            return True
        else: 
            return False
    #}}}

    def get_chan_data(self, wave_db,  esaz, at,  stacode, size):
#{{{
        """Opens waveform database and returns trace objects based on sta_chan. Applies
        calibration, splice, filter, decimation and rotation of channels
        !!! NOTE: 
        To avoid filter transient problems we need to ensure we get a large enough
        record section prior to the P-arrival time, on the order of minutes (size in secs). 
        Once we filter, then subset the trace object to the window around the arrival time.
        RLN (2012-02-02)

        """
        clip_max = float(eval(self.clip_values['max']))
        clip_min = float(eval(self.clip_values['min']))

        if self.verbose: self._log('Get channel data (trace objects) for station (%s)' % stacode)

        #'''
        #!!! NOTE: Start time and endtime need to be calculated using ptime 
        #Right now we are just getting size/2 either side of first arrival
        #Instead use arrival.time - should be the same thing. RLN (2011-11-29)
        #'''

        stachan_trace = {} # Holder dict for data

        #'''
        #!!! NOTE: Need to get a fresh view of 
        #          the wfdisc as subsetted for P 
        #          arrivals which are only done 
        #          on the BH.* channels.
        #          RLN (2011-11-29)
        #'''

        # opend database
        try:
            wvdb = datascope.dbopen(wave_db, 'r')
        except Exception, e:
            sys.exit('Could not open waveform database %s. Exception: %s' % (wave_db, e))

        # subset for wfdisc and station
        try: 
            wvdb.lookup(table='wfdisc')
            wvdb.subset('sta =~ /%s/' % stacode)
        except Exception, e:
            sys.exit('Could not subset wfdisc for %s [%s]' % (stacode, e))

        wvdb[3] = 0

        #'''
        #!!!! NOTE: Want maximium x3 the size 
        #            to account for timeshifts
        #            and add size before
        #            RLN (2012-02-02)
        #st              p arr.time            et
        #| <---------------->*<----><----><---->|
        #    size
        #'''

        """ 

        ***** OVERRIDE ***** 

        Lets try using the event time.  

        """
        at = self.evparams['time']

        st = at - size
        et = at + size
        resample = 0
        vang = 0

        try:
            if self.verbose: self._log('Retrieve %s data with time range (%s, %s)' % (stacode, st, et))
            wvdb.subset('sta=~/^%s$/ && chan=~/%s/' % (stacode, self.chan_to_use))
            trace = wvdb.load_css(st, et)
            trace.apply_calib()
            trace.splice() # Join all segments together
            trace.filter(self.filter_string)
            rotchan = ('R', 'T', 'Z')
            trace.rotate(esaz, vang, rotchan)
            trace.subset('chan =~ /R|T|Z/')
        except Exception, e:
            sys.exit('Could not prepare data for %s:%s [%s]' % (stacode,self.chan_to_use, e))

        if self.verbose: self._log('Number of traces for %s: %s' % (stacode, trace.nrecs()))

        for j in range(trace.nrecs()):
            try:
                trace[3] = j
                sta, chan, ns, sr = trace.getv('sta', 'chan', 'nsamp', 'samprate')
                stachan = '%s_%s' % (sta, chan)
                ns_tra = 0
                ns_req = 0
                filtered_samples = list(trace.data())
                use_data = 0
            except Exception, e:
                sys.exit('Could not extract data for %s:%s [%s]' % (stacode,self.chan_to_use, e))

            #  Test for clipped data
            #for i in filtered_samples:
            #    if i > clip_max:
            #        use_data = 1
            #    elif i < clip_min:
            #        use_data = -1
            # 

            # if use_data == 1:
            #     self._log('    - %s: One or more samples > clip_values[max] (%s). IGNORE.' % (stachan, clip_values['max']))
            # elif use_data == -1:
            #     self._log('    - %s: One or more samples < clip_values[min] (%s). IGNORE.' % (stachan, clip_values['min']))
            # else:
            #    '''
            #    !!! NOTE: Get three times the size to 
            #              allow for a positive time shift 
            #              RLN (2011-12-06)
            #    '''


            #max_size = int(size * 3) # Only works if 1Hz (1 sample per second, LH.*)

            # max_size = int(size) # Only works if 1Hz (1 sample per second, LH.*)

            #if self.verbose:
            #    self._loge('Sample size for %s: %s. Max: %s. Trimming.' % (chan, len(filtered_samples), max_size))

            #stachan_trace[chan] = filtered_samples[filter_timepad:max_size]  

            stachan_trace[chan] = filtered_samples[size:-1]  
            #if self.verbose:
            #    self._log('    - Trace extracted for %s samples:[%s], wanted:[%s]' % (chan, len(stachan_trace[chan]), max_size))

        trace.trdestroy()
        wvdb.free()
        wvdb.close()

        if self.verbose: self._log('Got data for [%s]: %s elements' % (stacode, len(stachan_trace)))

        return stachan_trace
#}}}

    def update_moment_tbl(self, strike, dip, rake):
#{{{
        """Write out results to 
        database moment table
        """
        if self.verbose:
            self._log('Write out moment tensor for orid %s to database table %s.moment' % (self.orid, self.event_db))
            self._log(' - MT for orid %s strike => %s' % (self.orid, strike))
            self._log(' - MT for orid %s dip => %s' % (self.orid, dip))
            self._log(' - MT for orid %s rake => %s' % (self.orid, rake))
        moment_dbptr = datascope.dbopen(self.event_db, 'r+')
        try:
            moment_dbptr.lookup(table='moment')
        except Exception, e:
            sys.exit('update_moment_tbl error: Error in lookup: %s' % e)
        else:
            orid_subset = datascope.dbsubset(moment_dbptr, 'orid == %s' % self.orid)
            if orid_subset.query('dbRECORD_COUNT') == 0:
                self._log('Adding new moment tensor to moment table with orid (%s)' % self.orid)
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
                    sys.exit('Adding record to moment table unknown error: %s' % e)
                else:
                    self._log('Successfully added record with orid (%s) to moment table' % self.orid)
            else:
                self._log('Updating moment tensor to moment table with orid (%s). Deleting current record and rewriting.' % self.orid)
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
                        sys.exit('Update record in moment table unknown error: %s' % e)
                    else:
                        self._log('Successfully updated record with orid (%s) in moment table' % self.orid)
            orid_subset.free()
            moment_dbptr.free()
            moment_dbptr.close()
            return True
#}}}

    def create_focal_mechanism(self, obspy_beachball, mt_images_dir, matrix_M=False, strike=False, dip=False, rake=False):
#{{{
        """Write out focal mechanism
        to images directory
        """
        if self.verbose:
            self._log('Writing file to images dir (%s)' % mt_images_dir)

        if not os.path.exists(mt_images_dir):
            if self.verbose:
                self._log('Images dir (%s) does not exist. Try to create...' % mt_images_dir)
            try:
                os.mkdir(mt_images_dir, 0775)
            except Exception, e:
                sys.exit('Moment tensor images dir (%s) does not exist and cannot be created! Exception: %s' % (mt_images_dir, e))

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

        if self.verbose:
            self._log('Try to plot focal mechanism: %s' % focal_mechanism)

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
        if self.verbose:
            self._log('Beachball(): defaults %s' % beachball_defaults)

        # Test for defaults in the pf
        for k, v in beachball_defaults.iteritems():
            if not obspy_beachball[k]:
                if self.verbose:
                    self._log('write_results(): Beachball(): Setting default for %s' % k)
                beachball_vals[k] = beachball_defaults[k]
            else:
                if self.verbose:
                    self._log('write_results(): Beachball(): Using pf defined value for %s' % k)
                beachball_vals[k] = obspy_beachball[k]
            if self.verbose:
                self._log('write_results(): Beachball(): Arg: %s, Val: %s' % (k, beachball_vals[k]))

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
            sys.exit('Error creating Beachball() %s: %s' % (Exception, e))
            return False
        else:
            if self.verbose:
                self._log('Successfully created focal mechanism image (%s)' % my_outpath)
            return my_outpath
#}}}

    # def calculate_synthetics_to_plot(self, gg, azimuth_list, matrix_M, size):
    def calculate_synthetics_to_plot(self, gg, ev2sta, matrix_M, size):
#{{{
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

        if self.verbose:
            self._log('Calculating synthetics to plot (rotated) for orid: %s' % self.orid)
            self._log('Size of synthetics: %s, Override size: %s' % (len(gg), size))

        syn_plot_dict = defaultdict(lambda: defaultdict(defaultdict))
        rotated_comps = ['T', 'R', 'Z']

        for i, sta_tup in enumerate(ev2sta):
            sta = sta_tup[0]
            az = sta_tup[1]
            if self.verbose:
                self._log(' - Determine plot synthetics for station (%s)' % sta)
            for rc in rotated_comps:
                if self.verbose:
                    self._log('  - Work on rotated channel (%s)' % rc)
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
#}}}

    def normalize_coefficient(self, data_as_list):
#{{{
        """Determine the 
        normalization factor
        for all the data
        """
        normalizer = max([abs(x) for x in data_as_list])
        if self.verbose:
            self._log("  - Normalization factor: %s" % normalizer)
        return normalizer
#}}}

    def create_data_synthetics_plot(self, ss, ev2sta, mod_gg, size, mt_images_dir):
#{{{
        """Create and save data vs. 
        synthetic waveform plots. 
        Equivalent to Dreger's function 
        mt_plot in mt_plot6iso2_linux2.c:

        mt_plot(ss,gg,nsta,Strike,Rake,Dip,
                St2,Rk2,Dp2,d_mt,Pdc,Pclvd,Piso,Mo, Mw, E, VR);

        There are as many TRZ 
        plots as stations
        """
        if self.verbose:
            self._log('Create plots of data vs. synthetics')

        # Init figure
        my_plot = plt.figure(figsize=(9, len(ev2sta)+0.5), dpi=100)
        my_plot.subplots_adjust(hspace=0.05, wspace=0.02,
                                bottom=0.02, top=0.96,
                                left=0.07, right=0.98)
        # Immutable rotated components
        # The order is important!
        rotated_components = ('T', 'R', 'Z')

        if self.verbose:
            self._log("Synthetics (Green's functions) used in calculations:")
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
            if self.verbose:
                self._log('Creating plots for station (%s)' % s)
            for j, rc in enumerate(rotated_components):
                axis_num += 1
                if self.verbose:
                    self._log(' - Processing (%s)' % rc)
                    self._log('  - ss[%s][%s], length:%s, timeshift: %s' % (i, rc, len(ss[i][rc]), ev2sta[i][3]))
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
                    self._log('  *** create_data_synthetics_plot(): Error: %s' % e)
                else:
                    new_data_list = [(v/data_scale_factor) for v in ss[i][rc][real_start:real_end]]

                try:
                    syn_scale_factor = self.normalize_coefficient(synthetic_list)
                except ValueError, e:
                    self._log('  *** create_data_synthetics_plot(): Error: %s' % e)
                else:
                    new_synthetic_list = [(v/syn_scale_factor) for v in synthetic_list]
                '''
                if self.verbose:
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
        syn_plot_outfile = '%s/%s_synthetics_fit.png' % (mt_images_dir, self.orid)
        my_plot.savefig(syn_plot_outfile)
        if os.path.isfile(syn_plot_outfile) and self.verbose:
            self._log('Successfully created data vs. synthetics plot (%s)' % syn_plot_outfile)
        elif not os.path.isfile(syn_plot_outfile):
            sys.exit('Error creating data vs. synthetics plot (%s)' % syn_plot_outfile)
        return syn_plot_outfile
#}}}

    def create_composite_plot(self, ttfont, img_anno, mt_images_dir, synthetics_img, focalmech_img):
#{{{
        """
        Create a composite image 
        from focal mechanism and 
        synthetics plots and update 
        the database

        All PIL work
        """ 
        if self.verbose:
            self._log('Merging the image annotations, synthetics & focal mechanism plots together')

        # Create composite image
        #size = (1400, 900)
        white = (255, 255, 255, 255)
        try:
            myfont = ImageFont.truetype(ttfont, 20)
        except IOError as e:
            sys.exit('Error importing font: %s' % e)

        final_file = '%s.png' % self.orid
        path_to_file = '%s/%s' % (mt_images_dir, final_file)
        syn_img = Image.open(synthetics_img, 'r')
        fm_img = Image.open(focalmech_img, 'r')

        size = (1400, syn_img.size[1])
        composite = Image.new('RGBA', size, white)
        composite.paste(syn_img, (10, 10))
        fm_position = (size[0] - (fm_img.size)[0] - 50, syn_img.size[1] - (fm_img.size)[1] - 50)
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
            sys.exit('Cannot save file (%s). Error: %s' % (final_file, e))
        else:
            # Update the database table
            try:
                mtimages_dbptr = datascope.dbopen(self.event_db, 'r+')
                mtimages_dbptr.lookup(table='moment_tensor_images')
            except Exception, e:
                sys.exit("Cannot open table 'moment_tensor_images'. Do you have the schema extension correctly installed? Error: %s" % e)

            try:
                mtimages_dbptr.subset('orid == %s' % self.orid)
                records = mtimages_dbptr.query('dbRECORD_COUNT')
                self._log('Subset moment_tensor_images table for orid == %s => [%s]' % (self.orid,records))
            except Exception,e:
                self._log('Subset moment_tensor_images table for orid == %s :  %s=>[%s]' % (self.orid,Exception,e))
                records = 0

            if not records:
                self._log('Adding new focal mechanism to moment_tensor_images table with orid (%s)[%s,%s,%s,%s]' % (self.orid,'test',int(self.orid),mt_images_dir,final_file))
                try:
                    # print 'sta:', 'test ', '\norid:', int(self.orid), '\ndir:', mt_images_dir, '\ndfile:', final_file
                    # print mtimages_dbptr
                    mtimages_dbptr.addv(
                        'sta', 'test',
                        'orid', int(self.orid),
                        'dir', mt_images_dir,
                        'dfile', final_file )
                    # print mtimages_dbptr
                    self._log('Successfully added record with orid (%s) to moment_tensor_images table' % self.orid)

                except Exception, e:
                    sys.exit('Adding record to moment_tensor_images table unknown error: %s=>[%s]' % (Exception,e))

            else:
                self._log('Updating focal mechanism to moment_tensor_images table with orid (%s). Deleting current record and rewriting.' % self.orid)
                for i in range(records):
                    mtimages_dbptr[3] = i
                    try:
                        mtimages_dbptr.putv(
                            'sta', 'test',
                            'orid', int(self.orid),
                            'dir', mt_images_dir,
                            'dfile', final_file)
                        self._log('Successfully updated record with orid (%s) to moment_tensor_images table' % self.orid)
                    except Exception, e:
                        sys.exit('Update record in moment_tensor_images table unknown error: %s' % e)

            mtimages_dbptr.close()

        return True
#}}}

    def __getitem__(self,i):
#{{{
        try:
            return self.i
        except:
            sys.exit('Event() => Wrong name of called element. [%s]\n'% i)
#}}}

    def __call__(self,i):
#{{{
        try:
            return self.i
        except:
            sys.exit('Event() => Wrong name of called element. [%s]\n'% i)
#}}}

    def __getattr__(self,i):
#{{{
        try:
            return self.i
        except:
            sys.exit('Event() => Wrong name of called element. [%s]\n'% i)
#}}}

    # }}}


"""
If we call this script directly, then output help.
"""
if __name__ == "__main__":
#{{{
    print "Moment Tensor Data Extaction Library:"
    print "\t1) Import the library.  "
    print "\t\t"
    print "\t2) "
    print "\t\t"
    print "\t3) "
    print "\t\t"
    print "\t3) "
    print "\t\t"

    print "\n\n"
    print "No BRTT support."
    print "Juan Reyes <reyes@ucsd.edu>"
    print "Rob Newman <rlnewman@ucsd.edu>"

#}}}

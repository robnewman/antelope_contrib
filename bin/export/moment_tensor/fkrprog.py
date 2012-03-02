from __main__ import *      # Get all the libraries from parent


class GreenFunctions():
#{{{
    """ 
      FREQUENCY WAVENUMBER-INTEGRATION PROGRAM.

    ORIGINAL FORTRAN CODE WRITTEN BY CHANDAN K. SAIKIA
    TRANSLATED TO PYTHON BY JUAN C. REYES 1/2012 <REYES@UCSD.EDU>

    NOTES:

    Call generate to build your 10 element matrix. The produced elements are the same
    as in Jost and Hermann (1989); transverse component vertical strike-slip(TSS) and vertical dip-slip (TDS) faluts,
    radial component vertical strike-slip (XSS), vertical dip-slip (XDS) and 45-degree dip-slip (XDD) faults, and 
    vertical component vertical strike-slip (ZSS), vertical dip-slip (ZDS), 45-degree (ZDD) faults, and radial (REX) 
    and vertical (ZEX) explosions.

    USAGE:

    Initialize the library with:
        GF =  fkr.GreenFunctions(pf_file)

    The pf_file is the parameter file that contains the model that we want to use.
    Generate GFs for depth of 8km and distance of 1km.
        GF.generate_python(depth=8,distance=10)

    Plot the GFs for the first 150 samples. 
        GF.plot(0,150)

    """
    def __init__(self,model,verbose=False):
#{{{
        self.verbose = verbose
        self.pf = model

        self.IVEL = 'v'
        self.IB  = 100
        self.IFREQ = 0
        self.DATA = defaultdict(lambda: defaultdict(dict))
        self.ALLOMEGA = defaultdict(lambda: defaultdict(dict))
        self.ELEMENTS = defaultdict(dict)


        if self.verbose: self._log("  init(%s)" % model)

        """ Read configuration from parameter file self.pf """
        self._read_pf()

        """ Return a normalized absolutized version of the path. """
        self.gf_db_path = os.path.abspath(self.gf_db_path)

        """ Verify database """
        db = self._open_db()

        db.close()

        if self.verbose: 
            self._log("  Using database [%s]" % os.path.normpath(self.gf_db_path + '/' + self.gf_db + '.wfdisc'))

#}}}


    def _log(self,message):
    #{{{
        """Global function to handle log output. Prints messages with a timestamp.  """

        curtime = stock.epoch2str( time(), '%m/%d/%Y %H:%M:%S')
            
        print '%s dbmoment: \tGreenFunction() => %s' % (curtime,message)
    #}}}



    def _open_db(self):
#{{{
        """  
        Open the database that will save our functions.

            From PF file:
            self.gf_db      => database name
            self.gf_db_path => database path *

            * Already nomalized by __init__()

        """

        """ Recursive directory creation function. """
        try:
            if not os.path.exists(self.gf_db_path): os.makedirs(self.gf_db_path)
        except Exception,e:
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot create direcotry (%s) %s => %s \n\n' % (self.gf_db_path,Exception,e))

        """ Verify the directory. """
        if not os.path.exists(self.gf_db_path):
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot find direcotry (%s)\n\n' % self.gf_db_path)

        if self.verbose: self._log("  Open database: %s/%s\n" % (self.gf_db_path,self.gf_db))

        """ Open database and keep pointer in db object"""
        try: 
            db = datascope.dbopen( os.path.normpath(self.gf_db_path + '/' + self.gf_db), "r+" )
            db.lookup( table='wfdisc' )

        except Exception, e:
            raise SystemExit('\n\nERROR: dbopen() %s => %s\n\n' % (Exception,e))

        return db
#}}}

    def build(self,depth,distance,sps=1,type='v',filter=None):
#{{{
        """
        Main function to retrieve ( and produce if missing )
        the requested elements.

        sps: The requested samplerate for the data
        type: Output velocity of displacepment
        filter: Apply this filter to the data 

        """

        self.DEPTH = depth
        self.DISTANCE = distance
        self.IVEL = type
        if self.verbose: self._log("  build(%s,%s,%s,%s,%s)" % (depth,distance,sps,type,filter))

        if self.verbose: self._log("  Get from database.")
        if not self._get_from_db(depth,distance,sps,type,filter): 
            if self.verbose: self._log("  Not in database. Generate new.")
            self._generate(depth,distance)
            if self.verbose: self._log("  Get from database.")
            self._get_from_db(depth,distance,sps,type,filter)

#}}}

    def _generate_python(self,depth,distance,type='v'):
#{{{
        self.IVEL = type
        self.DEPTH = depth
        self.DISTANCE = distance

        tempD =  self.D
        tempA =  self.A
        tempB =  self.B
        tempRHO =  self.RHO
        tempQA =  self.QA
        tempQB =  self.QB

        # Init vars for model layers
        self.D   = defaultdict(int)
        self.A   = defaultdict(int)
        self.B   = defaultdict(int)
        self.RHO = defaultdict(int)
        self.QA  = defaultdict(int)
        self.QB  = defaultdict(int)

        d = 0
        delta = 0
        for y in range(len(tempD)):
            if self.DEPTH == (d + tempD[y]):
                delta  = 1

                self.LMAX = y + 2
                for x in range(len(tempD)):
                    if x < y:
                        self.D[x]   =  tempD[x]
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                    elif x == y:
                        # new layer
                        self.D[x]   =  tempD[x] - delta
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                        # old layer
                        self.D[x+1]   =  delta
                        self.A[x+1]   =  tempA[x+1]
                        self.B[x+1]   =  tempB[x+1]
                        self.RHO[x+1] =  tempRHO[x+1]
                        self.QA[x+1]  =  tempQA[x+1]
                        self.QB[x+1]  =  tempQB[x+1]
                    else:
                        self.D[x+1]   =  tempD[x]
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                break

            elif self.DEPTH < (d + tempD[y]):
                self.LMAX = y + 1
                delta = d + tempD[y] - self.DEPTH

                for x in range(len(tempD)):
                    if x < y:
                        self.D[x]   =  tempD[x]
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                    elif x == y:
                        # new layer
                        self.D[x]   =  tempD[x] - delta
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                        # old layer
                        self.D[x+1]   =  delta
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                    else:
                        self.D[x+1]   =  tempD[x]
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                break

            else:
                d += tempD[y]


        if self.verbose:
            for I in range(len(self.D)):
                self._log("%s: %s\t%s\t%s\t%s\t%s\t%s" % (self.LMAX,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I]))

        
        if self.verbose: self._log("LMAX => %s" % self.LMAX)

        self.MMAX = len(self.D) - 1
        if self.verbose: self._log("MMAX => %s" % self.MMAX)

        if self.verbose: self._log("DEPTH => %s" % self.DEPTH)
        if self.verbose: self._log("RANGE => %s" % self.DISTANCE)

        self._setup_vars()


        #self._excit()
        self._wvint9()

#}}}

    def _generate(self,depth,distance,type='v'):
#{{{

        if self.verbose: self._log("generate(%s,%s,%s)" %(depth,distance,type))
        try:
            os.remove('./.greens_funcs')
            os.remove('./GREEN.1')
            os.remove('./junk')
        except:
            pass

        #
        # Build model
        # output new model to ./greens_funcs 
        #
#{{{

        self.IVEL = type
        self.DEPTH = depth
        self.DISTANCE = distance
        self.NYQ = self.N / 2 + 1
        self.DF = 1 / (self.N * (1/self.DT))
        self.FU = (self.NYQ - 1) * self.DF
        self.TLEN = self.N * (1/self.DT)
        #self.ALPHA = -log(self.DECAY) / self.TLEN
        if self.verbose: self._log("Build new model.")

        """
        Build new layers on earth model.
        """
        tempD =  self.D
        tempA =  self.A
        tempB =  self.B
        tempRHO =  self.RHO
        tempQA =  self.QA
        tempQB =  self.QB

        # Init vars for model layers
        self.D   = defaultdict(int)
        self.A   = defaultdict(int)
        self.B   = defaultdict(int)
        self.RHO = defaultdict(int)
        self.QA  = defaultdict(int)
        self.QB  = defaultdict(int)

        if self.verbose: self._log("Setup vars.")
        d = 0
        delta = 0
        for y in range(len(tempD)):
            if self.DEPTH == (d + tempD[y]):
                delta  = 1

                self.LMAX = y + 2
                for x in range(len(tempD)):
                    if x < y:
                        self.D[x]   =  tempD[x]
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                    elif x == y:
                        # new layer
                        self.D[x]   =  tempD[x] - delta
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                        # old layer
                        self.D[x+1]   =  delta
                        self.A[x+1]   =  tempA[x+1]
                        self.B[x+1]   =  tempB[x+1]
                        self.RHO[x+1] =  tempRHO[x+1]
                        self.QA[x+1]  =  tempQA[x+1]
                        self.QB[x+1]  =  tempQB[x+1]
                    else:
                        self.D[x+1]   =  tempD[x]
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                break

            elif self.DEPTH < (d + tempD[y]):
                self.LMAX = y + 1
                delta = d + tempD[y] - self.DEPTH

                for x in range(len(tempD)):
                    if x < y:
                        self.D[x]   =  tempD[x]
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                    elif x == y:
                        # new layer
                        self.D[x]   =  tempD[x] - delta
                        self.A[x]   =  tempA[x]
                        self.B[x]   =  tempB[x]
                        self.RHO[x] =  tempRHO[x]
                        self.QA[x]  =  tempQA[x]
                        self.QB[x]  =  tempQB[x]
                        # old layer
                        self.D[x+1]   =  delta
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                    else:
                        self.D[x+1]   =  tempD[x]
                        self.A[x+1]   =  tempA[x]
                        self.B[x+1]   =  tempB[x]
                        self.RHO[x+1] =  tempRHO[x]
                        self.QA[x+1]  =  tempQA[x]
                        self.QB[x+1]  =  tempQB[x]
                break

            else:
                d += tempD[y]

        """
        self._log(new layers and other params.)
        """

        #if self.verbose:
        #    for I in range(len(self.D)):
        #        self._log("%s: %s\t%s\t%s\t%s\t%s\t%s" % (self.LMAX,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I]))

        
        if self.verbose: self._log("LMAX => %s" % self.LMAX)

        self.MMAX = len(self.D) - 1
        if self.verbose: self._log("MMAX => %s" % self.MMAX)

        if self.verbose: self._log("DEPTH => %s" % self.DEPTH)
        if self.verbose: self._log("RANGE => %s" % self.DISTANCE)


        fortran_gg = defaultdict(lambda: defaultdict(dict))
        self.GG = defaultdict(dict)
        self.ALLOMEGA = []
        self.LMAX += 1
        self.MMAX += 1
        template = self._file_format()
        model = ''
        for I in range(len(self.D)):
            model += " %1.4E %1.4E %1.4E %1.4E   600.00    300.00\n" % (self.D[I],self.A[I],self.B[I],self.RHO[I])

        model =  template % (self.DEPTH,(1/self.DT),self.MMAX,model,self.LMAX,self.DISTANCE)

        """ 
        Open temp file to put model.
        """
        try:
            os.remove('./.greens_funcs')
        except Exception, e:
            pass

        try: 
            f = open('./.greens_funcs', 'w')
            f.write(model)
            f.close()
        except Exception,e:
            raise SystemExit('\n\nERROR: Cannot open temp file .greens_funcs %s %s\n'% (Exception,e))
#}}}

        self._log("\n%s" % model)
        if self.verbose: self._log('Running: fortran_fkrprog < ./greens_func')
        p = os.popen('fortran_fkrprog < ./.greens_funcs',"r")
        while 1:
            line = p.readline()
            self._log(line)
            if not line: break

        p.close()

        if self.verbose: self._log('Running: wvint9')
        p = os.popen('wvint9',"r")
        while 1:
            line = p.readline()
            line = line.strip('\n')
            line = line.replace(' ','')
                
            if not line: break

            if re.match("GG",line):
                m = re.match("GG\[(\d+)\]\[(\d+)\]=(.+)",line)
                if not m: continue
                l = int(m.group(1))-1
                n = int(m.group(2))-1
                self.DATA[l][n]= eval(m.group(3))
                self._log("DATA[%s][%s] = %s" % (l,n,m.group(3)))
            else: 
                if self.verbose: self._log("ERROR ** no match ** : %s" % line)

        p.close()

        try:
            os.remove('./.greens_funcs')
            os.remove('./GREEN.1')
            os.remove('./junk')
        except:
            pass

        self._reformat()

#}}}

    def _get_from_db(self,depth,distance,sps=1,type='v',filter=None):
#{{{

        """ 
        Open the database and extract the traces for requested depth and distance. 

        The functions are archived in ascii files referenced by Datascope using a simple 
        wfdisc table. The value for the station is our DISTANCE. The value for the channel  
        is our DEPTH and the element is specified in the location code. 
        i.e.   
            depth: 8
            distance: 10
            element: TDS
            => 10_8_TDS ( format: sta_chan_loc )


        All data will be extracted from the database and archived in memory internal to the class. We 
        will use the dictionary self.ELEMENT for this. Each key in the dictionary will be dedicated
        to a different component of the fundamental Green's Functions and will include objects for metadata. 
        """

        """ Make local copy of database pointer. """
        db = self._open_db()

        try:
            records = db.query(datascope.dbRECORD_COUNT)
        except Exception,e:
            records = 0 

        if not records:
            return False

        """ Subset for distance and depth. """

        try: 
            db.subset('sta =~ /%s/ && chan =~ /%s_.*/' % (distance,depth))
            records = db.query(datascope.dbRECORD_COUNT)
        except: 
            records = 0

        if records: 
            """ Get list of elements and get values for time, endtime, nsamp, samprate and LOC_CODE. """
            for i in range(records):

                db.record = i

                try:
                    (db_sta,db_chan,db_nsamp,db_time,db_endtime,db_samprate) = \
                            db.getv('sta','chan','nsamp','time','endtime','samprate')
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): Problems while extracting from database %s: %s\n\n' % (Exception,e))

                if self.verbose:
                    self._log('GreenFunctions(): getv()=> (%s,%s,%s,%s,%s,%s)'%(db_sta,db_chan,db_nsamp,db_time,db_endtime,db_samprate))

                """ Get steps needed to get to new samplerate. """
                decimate_files = self._decimate_file(db_samprate,sps)
                if self.verbose: self._log('GreenFunctions(): decimate factor [%s] file [%s]' % ( db_samprate/sps,decimate_files))

                """ Extract the elemnt name from the channel text. """
                try:
                    m = re.match(".*_(...)",db_chan)
                    comp = m.group(1)
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): Problems in regex [.*_(...)] on [%s] %s: %s\n\n' % (db_chan,Exception,e))

                if not comp: 
                    raise SystemExit('\n\nERROR: GreenFunctions(): Cannot find component name in wfdisc entry: %s_%s\n\n' % (db_sta,db_chan))

                """ Build/clean object for selected channel. """
                self.ELEMENTS[comp] = {'data':[], 'samplerate':db_samprate, 'filter':filter} 


                # With TRACEOBJECTS
                try:
                    if self.verbose: self._log('GreenFunctions(): trloadchan(%s,%s,%s,%s)' % (db_time,db_endtime,db_sta,db_chan))
                    tr = datascope.trloadchan(db,db_time,db_endtime,db_sta,db_chan)
                    if self.verbose: self._log('GreenFunctions(): tr.slice()')
                    tr.splice()
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): trace object %s: %s\n\n' % (Exception,e))

                try:
                    # integrate data
                    if type == 'D' or type == 'd':
                        if self.verbose: self._log('GreenFunctions(): integrate for velocity()')
                        tr.filter('INT')
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): integrate %s: %s\n\n' % (Exception,e))

                try:
                    # filter data
                    if filter:
                        if self.verbose: self._log('GreenFunctions(): filter(%s)' % filter)
                        tr.filter(filter)
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): filter %s: %s\n\n' % (Exception,e))


                #
                # TEST DECIMATION
                #
                #tr[3] = 0
                #x =  tr.data()
                #plt.plot(x[1000:1200])


                try:
                    # decimate data
                    for f in decimate_files:
                        full_path = os.environ['ANTELOPE'] + '/data/responses/' + f
                        if self.verbose: self._log('GreenFunctions(): filter(DECIMATE %s)' % full_path)
                        tr.filter('DECIMATE %s' % full_path)
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): decimate %s: %s\n\n' % (Exception,e))

                #tr[3] = 0
                #y =  tr.data()
                #plt.plot(y[500:700])
                #plt.legend(['original','decimated'])

                #plt.suptitle("Green Functions: ")
                #plt.show()

                try:
                    tr[3] = 0
                    self.ELEMENTS[comp]['data'] = tr.data()
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): extract data %s: %s\n\n' % (Exception,e))

                try:
                    tr.trdestroy()
                except:
                    pass

                self.ELEMENTS[comp]['data'] = self.ELEMENTS[comp]['data'][int(len(self.ELEMENTS[comp]['data'])/2):-1] 

                if not self.ELEMENTS[comp]['data']:
                    raise SystemExit('\n\nERROR: GreenFunctions(): Cannot build component [%s]: %s_%s\n\n' % (m.group(1),db_sta,db_chan))

        else:
            self._log("\n\nERROR: GreenFunctions(): ERROR: No records for subset [sta =~ /%s/ && chan =~ /%s/].\n\n" % (distance,depth))
            return False

        try:
            db.close()
        except:
            pass

        return True
#}}}

    def _save_to_db(self):
#{{{
        """ 
        Open the database and save the new traces.

        The functions are archived in ascii files referenced by Datascope using a simple 
        wfdisc table. The value for the station is our DISTANCE. The value for the channel  
        is our DEPTH and the element is specified in the location code. 
        i.e.   
            depth: 8
            distance: 10
            element: TDS
            => 10_8_TDS ( format: sta_chan_loc )

        """

        """ Open the database. """
        try:
            db = self._open_db()
        except Exception,e:
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot open database (%s) %s => %s \n\n' % (db,Exception,e))

        """ Return a normalized absolutized version of the path to save the files. """
        path = os.path.abspath(self.gf_db_path + '/' + self.gf_db + '/' + str(self.DEPTH)  )


        """ Recursive directory creation function. """
        try:
            if not os.path.exists(path): os.makedirs(path)
        except Exception,e:
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot create direcotry (%s) %s => %s \n\n' % (path,Exception,e))


        """ Save all data to file. """
        dfile = "%s_%s_%s.gf" % (self.DISTANCE,self.DEPTH,self.gf_db)
        try:
            os.remove("%s/%s"%(path,dfile))
        except Exception, e:
            pass

        try: 
            f = open("%s/%s"%(path,dfile), 'w')
        except Exception,e:
            raise SystemExit('\n\nERROR: Cannot open file %s %s %s\n'% (dfile,Exception,e))

        """ Return a relative path to add to wfdisc. """
        dbpath = os.path.relpath(path, self.gf_db_path)

        for element in sorted(self.ELEMENTS.keys()):

            #self._log('%s' % self.ELEMENTS[element])
            if self.ELEMENTS[element]['data']:
                sta = '%s' % self.DISTANCE
                chan_loc = '%s_%s' % (self.DEPTH,element)
                nsamp = len(self.ELEMENTS[element]['data'])
                samprate = self.ELEMENTS[element]['samplerate']
                endtime = (nsamp*samprate)+1.00
                f.write('%s\t%s\t%s\n'%(element,nsamp,samprate))
                start = f.tell()

                try:
                    [f.write('%s\n' % x )for x in self.ELEMENTS[element]['data'] ]
                except Exception,e:
                    raise SystemExit('\n\nERROR: Cannot add to file [%s] => %s\n'% (Exception,e))

                if self.verbose:
                    self._log('Add new record to wfdisc table: (%s,%s,%s,%s,%s,%s)' % (sta,chan_loc,nsamp,samprate,endtime,start))
                    #self._log('\tsta=',sta,'chan=',chan_loc,'time=',1.00,'endtime=',endtime)
                    #self._log('\tnsamp=',nsamp,'samprate=',samprate,'calib=',1,'datatype=','as')
                    #self._log('\tsegtype=','c','dir=',dbpath,'dfile=',dfile,'foff=',start)

                try:
                    db.addv('sta',sta,'chan',chan_loc,'time',1.00,'endtime',endtime,'nsamp',nsamp,'samprate',samprate,'calib',1,'datatype','as','segtype','V','dir',dbpath,'dfile',dfile,'foff',start)
                except Exception,e:
                    raise SystemExit('\n\nERROR: Cannot add new line [%s] %s %s\n'% (element,Exception,e))

            else: 
                raise SystemExit('\n\nERROR: Empty element to dump into file %s %s\n'% (element,file))

        try: 
            f.close()
        except Exception,e:
            raise SystemExit('\n\nERROR: Cannot close file %s %s %s\n'% (file,Exception,e))



        try:
            records = db.query(datascope.dbRECORD_COUNT)

        except Exception,e:
            raise SystemExit('\n\nERROR: GreenFunctions(): Problems with database %s: %s\n\n' % (Exception,e))

        try: 
            db.close()
        except Exception,e:
            raise SystemExit('\n\nERROR: Cannot close database %s %s %s\n'% (db,Exception,e))

        #if not records:
        #    raise SystemExit("\n\nERROR: GreenFunctions(): ERROR: No records in databse!.\n\n")

        #""" Subset for distance and depth. """

        #try: 
        #    db.subset('sta =~ /%s/ && chan =~ /%s_.*/' % (distance,depth)
        #    records = db.query(datascope.dbRECORD_COUNT)
        #except: 
        #    records = 0

        #if records: 
        #    """ Get list of elements and get values for time, endtime, nsamp, samprate and LOC_CODE. """
        #    for i in range(records):

        #        db.record = i

        #        try:
        #            (db_sta,db_chan,db_nsamp,db_time,db_endtime,db_samprate) = \
        #                    db.getv('sta','chan','nsamp','time','endtime','samprate')
        #        except Exception,e:
        #            raise SystemExit('\n\nERROR: GreenFunctions(): Problems while extracting from database %s: %s\n\n' % (Exception,e))

        #        """ 
        #        Save each trace in memory. 
        #        """

        #        """ Extract the elemnt name from the channel text. """
        #        m = re.match(".?_(...)",db_chan)
        #        comp = m.group(1)

        #        if not comp: 
        #            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot find component name in wfdisc entry: %s_%s\n\n' % (db_sta,db_chan))

        #        if not self.ELEMENTS[comp]:
        #            self.ELEMENTS[comp] = {'data':[], 'samplerate':1, 'filter':filter} 

        #        try:
        #            self.ELEMENTS[comp]['data'] = db.sample(db_time,db_endtime,db_sta,db_chan,False,filter)
        #        except Exception,e:
        #            self._log('\n %s = db.sample(%s,%s,%s,%s) => %s: %s\n' % (m.group(1),db_time,db_endtime,db_sta,db_chan,Exception,e) )
        #            raise SystemExit('\n\nERROR: GreenFunctions(): Problems while extracting trace object %s: %s\n\n' % (Exception,e))


        #        if not self.ELEMENTS[comp]['data']:
        #            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot build component [%s]: %s_%s\n\n' % (m.group(1),db_sta,db_chan))

        #else:
        #    self._log("\n\nERROR: GreenFunctions(): ERROR: No records for subset [sta =~ /%s/ && chan =~ /%s/].\n\n" % (distance,depth))
        #    return False

        #return True
#}}}

    def _decimate_file(self,have=1,want=1):
#{{{

        """
        Function to return needed response file for proper decimation of the sample 
        during the filtering process in the extraction routine. Reead manpage for wffilbrtt.
        ...
            DECIMATE fir_file
                This implements a decimation filter. The decimation is performed
                after applying one or more FIR anti-alias filters to  the  data.
                The  FIR  filters  are  specified  in  the  standard  instrument
                response  file,   fir_file   (the   format   is   specified   in
                response(5)).   The  response  function  specified  in  fir_file
                should only contain FIR stages with  their  decimation  factors.
                The  output  filtered  data  will  have the new decimated sample
                rate. Therefore any software that assumes sample  rates  do  not
                change  with  filtering will not work properly with this filter.
                The output filtered results are only computed for  time  samples
                for  which  the FIR filter sequences do not encounter data gaps.
                Therefore gaps in the data will grow larger when using this fil-
                ter to account for the edge effects of the FIR filtering.
        ...
        

        I'm selecting a series of files with simple decimation factors to build our
        new time-series with the requested samplerate. 

            % egrep "decimation factor" * |grep trident
                (standard input):93:trident_1000sps_fir1:# decimation factor     5 
                (standard input):94:trident_1000sps_fir2:# decimation factor     3 
                (standard input):95:trident_1000sps_fir3:# decimation factor     2 
                (standard input):96:trident_100sps_fir1:# decimation factor     15 
                (standard input):97:trident_100sps_fir2:# decimation factor     10 
                

        """
        factor = int(have/want)

        if factor == 1: 
            return [] 
        elif factor == 2:
            return ['trident_1000sps_fir3']
        elif factor == 3:
            return ['trident_1000sps_fir2']
        elif factor == 4:
            return ['trident_1000sps_fir3','trident_1000sps_fir3']
        elif factor == 5:
            return ['trident_1000sps_fir1']
        elif factor == 6:
            return ['trident_1000sps_fir3','trident_1000sps_fir2']
        elif factor == 10:
            return ['trident_100sps_fir2']
        elif factor == 15:
            return ['trident_100sps_fir1']
        else:
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot decimate by (%s) only [2,3,4,5,6,10,15]\n\n' % factor)

#}}}

    def _file_format(self):
    #{{{
        return ".F.\n\
    0   64\n\
GREEN.1\n\
    6.0      %05.2f      1  512 1024    %03.1f      %2d    1\n\
    1    1    1    1    1    1    1    1    1    1    0\n\
%s\
    %d\n\
  0.4000000E+03  1.500000E+00         0\n\
    1  10000.0     30.0      2.9       2.5\n\
  %03.1f        0.0      10.0\n"
    #}}}

    def _read_pf(self):
    #{{{
        # NEED TO FIX...
        #    self.MMAX   = int(self.MMAX) total number of layers
        # LMAX
        # Layer number below the source
        #    self.LMAX = 3

        #
        # Read data from file
        #
        if self.verbose: self._log("Read PF file: %s.pf" % self.pf)

        try:
            self.gf_db = stock.pfget_string(self.pf,'database')
            self.gf_db_path = stock.pfget_string(self.pf,'db_path')
        except Exception, e:
            raise SystemExit('\n\nWrong Format of PF file[%s]. %s %s\n'% (self.pf,Exception,e))

        try:
            self.DECAY  = stock.pfget_double(self.pf,'decay')
            self.N1     = stock.pfget_int(self.pf,'start_frequency') 
            self.N2     = stock.pfget_int(self.pf,'end_frequency') 
            self.N      = (self.N2-self.N1+1)*2 # 2 times the total number of freqs
            self.DT     = stock.pfget_double(self.pf,'samplerate') 
        except Exception,e: 
            raise SystemExit('\n\nWrong Format of PF file[%s]. %s %s\n'% (self.pf,Exception,e))

        if self.verbose: self._log("%s: DECAY=%s N1=%s N2=%s N=%s DT=%s " % (self.pf,self.DECAY,self.N1,self.N2,self.N,self.DT))


        # ISRC, JBDRY 
        self.ISRC = [ 1 for x in range(10) ]
        self.JBDRY = 0
        if self.verbose: self._log("%s: ISRC=%s JBDRY=%s" % (self.pf,self.ISRC,self.JBDRY))


        # Init vars for model layers
        self.D   = defaultdict(int)
        self.A   = defaultdict(int)
        self.B   = defaultdict(int)
        self.RHO = defaultdict(int)
        self.QA  = defaultdict(int)
        self.QB  = defaultdict(int)

        try:
            temp  = stock.pfget_string(self.pf,'model').splitlines()
            t  = defaultdict(int)
            for x in range(len(temp)):
                t = temp[x].split()
                if not t: continue
                if self.verbose: self._log("%s: RAW MODEL LAYER: %s" % (self.pf,t))
                self.D[x]   = float(t[0])
                self.A[x]   = float(t[1])
                self.B[x]   = float(t[2])
                self.RHO[x] = float(t[3])
                self.QA[x]  = 1/float(t[4])
                self.QB[x]  = 1/float(t[5])

        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n RAW: %s'% (self.pf,Exception,e,temp))


        if self.verbose:
            for I in range(len(self.D)):
                self._log("%s: %s\t%s\t%s\t%s\t%s\t%s" % (self.pf,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I]))


        #  NFC...
        self.XLENG  = 400 
        self.XFAC  = 1.5


        self.CMAX = stock.pfget_double(self.pf,'cmax')
        self.C1   = stock.pfget_double(self.pf,'c1')
        self.C2   = stock.pfget_double(self.pf,'c2')
        self.CMIN = stock.pfget_double(self.pf,'cmin')
        if self.verbose: self._log("%s: CMAX=%s C1=%s C2=%s CMIN=%s" % (self.pf,self.CMAX,self.C1,self.C2,self.CMIN))

        # need to fix this while running
        self.T0    = defaultdict(int) 


        ##
        ## FIX VALUES TO MATCH PYTHON INDEXING
        ##
        self.LMAX = len(self.D)-1 # layer number below the source

        self.MMAX = len(self.D)-1 # Total number of layers

        return

    #}}}

    def _setup_vars(self):
    #{{{

        # Fix values for C1, C2, CMAX and CMIN
        if self.verbose: self._log("CMAX=%s" % self.CMAX)
        try: self.CMAX = 1/self.CMAX
        except Exception,e: raise SystemExit('\n\nCMAX. %s(%s) \n'% (Exception,e))
        try: self.C1 = 1/self.C1
        except Exception,e: raise SystemExit('\n\nC1. %s(%s) \n'% (Exception,e))
        try: self.C2 = 1/self.C2
        except Exception,e: raise SystemExit('\n\nC2. %s(%s) \n'% (Exception,e))
        try: self.CMIN = 1/self.CMIN
        except Exception,e: raise SystemExit('\n\nCMIN. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log("CMAX=%s CMIN=%s C1=%s C2=%s" % (self.CMAX,self.CMIN,self.C1,self.C2))

        # FACT(I)
        self.FACT = [0.0 for x in range(self.IB)]
        if self.CMAX < 0.0:
            for I in range(self.IB):
                self.FACT[I] = 1.0

        if self.verbose: self._log("FACT=%s" % self.FACT)

        #
        # CHECK IF DEPTH DOES CORRESPOND TO DEPTH INTERFACE
        #
        self.DPH = 0.0
        if self.verbose: self._log("Adding layers on top of source[%s]" % self.LMAX)

        if self.verbose: self._log("\tLMAX= %s" % self.LMAX)
        for I in range(self.LMAX):
            self.DPH += self.D[I]
            if self.verbose: self._log("\tDPH= %s" % self.DPH)

        if abs(self.DPH-self.DEPTH) <= 0.1: self.DEPTH = self.DPH

        #if abs(self.DPH-self.DEPTH) > 0.1:
        #    raise SystemExit('\n\nCHECK TWO PARAMETERS <+> LMAX [%s] AND DEPTH [%s]\n'% (self.LMAX,self.DEPTH))

        #
        # CALCULATE ALPHA FROM DECAY
        #
        self.TLEN = self.N * self.DT
        if self.verbose: self._log('TLEN = %s' % self.TLEN)

        try: self.DECAY = 1 / self.DECAY
        except Exception,e: raise SystemExit('\n\nDECAY. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log('DECAY = %s' % self.DECAY)

        try:
            self.ALPHA = -log(self.DECAY)
            self.ALPHA = self.ALPHA / self.TLEN
        except Exception,e: raise SystemExit('\n\nALPHA. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log('ALPHA = %s' % self.ALPHA)

        #
        # CHECK BOUCHON'S CRITERIA
        #
        try:
            self.RMAX  = self.DISTANCE
            self.VAMIN = min([x for x in self.A.values()])
            self.VAMAX = max([x for x in self.A.values()])
            self.VBMIN = min([x for x in self.B.values()])
            self.VBMAX = max([x for x in self.B.values()])
        except Exception,e: raise SystemExit('\n\nRMAX VMIN VMAX VBMIN VBMAX. %s(%s) \n'% (Exception,e))

        #if self.verbose: self._log('A= %s' % self.A)
        #if self.verbose: self._log('B= %s' % self.B)
        #if self.verbose: self._log('RANGE= %s' % self.DISTANCE)
        if self.verbose: self._log('RMAX = %s' % self.RMAX)
        if self.verbose: self._log('VAMAX = %s' % self.VAMAX)
        if self.verbose: self._log('VBMAX = %s' % self.VBMAX)
        if self.verbose: self._log('VAMIN = %s' % self.VAMIN)
        if self.verbose: self._log('VBMIN = %s' % self.VBMIN)


        if self.verbose: self._log('RANGE= %s' % self.DISTANCE)
        self.YMAX = 2.0 * self.RMAX
        if self.verbose: self._log('RMAX= %s' % self.RMAX)
        if self.verbose: self._log('YMAX = %s' % self.YMAX)

        if self.XLENG < self.YMAX: self.XLENG = self.YMAX
        if self.verbose: self._log('XLENG = %s' % self.XLENG)

        self.RMAX = self.RMAX + sqrt( self.VAMAX * self.VAMAX * self.TLEN * self.TLEN - self.DEPTH * self.DEPTH )
        if self.verbose: self._log('RMAX = %s' % self.RMAX)

        if self.XLENG < self.RMAX:
                self.XLENG=self.RMAX
                self.MM = int(self.XLENG / 100.0)
                self.XLL = int(self.MM) * 100.0
                if self.XLL != self.XLENG: 
                    self.XLENG = self.XLL + 100.0
        if self.verbose: 
            self._log('MM = %s' % self.MM)
            self._log('XLL = %s' % self.XLL)
            self._log('XLENG = %s' % self.XLENG)

        try: self.DF = 1 / (self.N * self.DT)
        except Exception,e: raise SystemExit('\n\nDF. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log('DF = %s' % self.DF)

        try: self.NYQ = self.N / 2 + 1
        except Exception,e: raise SystemExit('\n\nNYQ. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log('NYQ = %s' % self.NYQ)

        self.FL =0.0

        if self.verbose: self._log('FL = %s' % self.FL)

        try: self.FU = (self.NYQ - 1) * self.DF
        except Exception,e: raise SystemExit('\n\nFU. %s(%s) \n'% (Exception,e))

        if self.verbose: self._log('FU = %s' % self.FU)


        if self.IFREQ > self.N1: self.N1 = self.IFREQ + 1


        if self.verbose: self._log('N1 = %s' % self.N1)
        if self.verbose: self._log('N2 = %s' % self.N2)
        if self.verbose: self._log('DF = %s' % self.DF)

        return
        #for x in range(self.N2):
        for x in range(self.N1,self.N2+1):
            if self.verbose: self._log('I=%s' % x)

            #self.I = x
            self.FREQ = (x-1) * self.DF
            if self.FREQ < self.DF: self.FREQ = 0.01 * self.DF
            if self.verbose: self._log('FREQ => %s' % self.FREQ)
            self.OMEGA = 2.0 * acos(-1.0) * self.FREQ
            if self.verbose: self._log('OMEGA => %s' % self.OMEGA)

            #self._log('X = %s' % x)
            #self._log('DF = %s' % self.DF)
            #self._log('FREQ => %s' % self.FREQ)
            #self._log('OMEGA => %s' % self.OMEGA)

            self.WVCM = self.OMEGA * self.CMAX
            self.WVC1 = self.OMEGA * self.C1
            self.WVC2 = self.OMEGA * self.C2
            self.WVCN = self.OMEGA * self.CMIN

            if self.verbose: self._log('WVCM => %s' % self.WVCM)
            if self.verbose: self._log('WVC1 => %s' % self.WVC1)
            if self.verbose: self._log('WVC2 => %s' % self.WVC2)
            if self.verbose: self._log('WVCN => %s' % self.WVCN)

            #self._excit()
            self.ALLOMEGA[x-1] = [self.OMEGA,self.NK]
            if self.verbose: self._log("\t\t\t\t%s  F =  %12.4E NBLOCK =  %5s IBLOCK= %3s" % (x,self.FREQ,self.NBLOCK,self.IBLOCK))
            #self._log("\t\t\t\t%s  DATA :  %s" % (x,self.DATA[x]))
            if self.verbose:
                if self.verbose: self._log("\tDIST =>    1")
            for z in range(10):
                #self._log("\t\t%-10s%-16.9E%-4s%-16.9E" % self.GG[y][z])
                self.DATA[x-1][z] = self.GG[1][z]
                if self.verbose: self._log("\t\t%16.9E\t%16.9E" % (self.GG[1][z].real,self.GG[1][z].imag))


        self.XX = -999.0
        self.T0X = 0.0

        #self._log("%s %s" % (self.XX,self.T0X))


    #}}}

    def _excit(self):
#{{{
        # SUBROUTINE EXCIT(OMEGA,ISRC,LMAX,RANGE,ND,NBLK,IBLOCK)
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # VECTORIZED FREQUENCY-WAVENUMBER CODE +++++ THIS SUBROUTINE CALCULATES
        # THE KERNEL FOR BOTH P-SV AND SH SYSTEM I.E., EIGHT BASIC FAULTS AND
        # FOR THE EXPLOSIVE SOURCE.  SOURCE SHOULD LOCATE AT AN INTERFACE. I HAVE
        # USED HASKELL'S 4X4 MATRIX AND WANG AND HERRMANN'S 6X6 COMPOUND MATRICES.
        # MATRIX MULTIPLICATION IS REDUCED SIGNIFICANTLY BY USING WATSON (1970).
        #
        #               OMEGA    = ANGULAR FREQUENCY
        #               ISRC(I)  = I VARIES BETWEEN 1 TO 10
        #               LMAX     = INTERFACE WHERE THE SOURCE IS BUIRED
        #               RANGE    = DISTANCES
        #               ND       = NUMBER OF TOTAL DISTANCES
        #               NBLK   = TOTAL NUMBER OF BLOCKS FOR THIS FREQUENCY
        #               IBLOCK   = COUNT IN NBLK TH BLOCK
        #
        # WRITTEN  BY CHANDAN K. SAIKIA
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        #{{{

#{{{ BUILD VARS 
        self.ICNT = defaultdict(lambda: defaultdict(dict))
        self.SMM = defaultdict(lambda: defaultdict(dict))
        self.RA = defaultdict(lambda: defaultdict(dict))
        self.RB = defaultdict(lambda: defaultdict(dict))
        self.GAM = defaultdict(lambda: defaultdict(dict))
        self.GAM2 = defaultdict(lambda: defaultdict(dict))
        self.GAM3 = defaultdict(lambda: defaultdict(dict))
        self.GAM4 = defaultdict(lambda: defaultdict(dict))
        self.GAMM1 = defaultdict(lambda: defaultdict(dict))
        self.GAMM2 = defaultdict(lambda: defaultdict(dict))
        self.GAMM3 = defaultdict(lambda: defaultdict(dict))
        self.GAMX1 = defaultdict(lambda: defaultdict(dict))
        self.GAMX2 = defaultdict(lambda: defaultdict(dict))
        self.EXE = defaultdict(lambda: defaultdict(dict))
        self.EXEL = defaultdict(lambda: defaultdict(dict))
        self.EXL = defaultdict(lambda: defaultdict(dict))
        self.EXLL = defaultdict(lambda: defaultdict(dict))
        self.COSQL = defaultdict(lambda: defaultdict(dict))
        self.YL = defaultdict(lambda: defaultdict(dict))
        self.YLI = defaultdict(lambda: defaultdict(dict))
        self.ZL = defaultdict(lambda: defaultdict(dict))
        self.ZX = defaultdict(lambda: defaultdict(dict))
        self.K = defaultdict(lambda: defaultdict(dict))
        self.SNGL = defaultdict(lambda: defaultdict(dict))
        self.C = defaultdict(lambda: defaultdict(dict))

        self.GG = defaultdict(lambda: defaultdict(complex)) 
        self.J0 = defaultdict(complex) 
        self.J1 = defaultdict(complex) 
        self.J2 = defaultdict(complex) 
        self.SUME = defaultdict(complex) 
        self.SUMC = defaultdict(complex) 
        self.WVNO  = defaultdict(complex)
        self.WVNO2 = defaultdict(complex)
        self.WVNO3 = defaultdict(complex)
        self.S1 =    defaultdict(complex)
        self.S2 =    defaultdict(complex)
        self.S3 =    defaultdict(complex)
        self.S4 =    defaultdict(complex)
        self.S5 =    defaultdict(complex)
        self.S6 =    defaultdict(complex)
        self.S7 =    defaultdict(complex)
        self.S8 =    defaultdict(complex)
        self.S9 =    defaultdict(complex)
        self.S10=    defaultdict(complex)
        self.G1 =    defaultdict(complex)
        self.G2 =    defaultdict(complex)
        self.G3 =    defaultdict(complex)
        self.G4 =    defaultdict(complex)
        self.G5 =    defaultdict(complex)
        self.G6 =    defaultdict(complex)
        self.G7 =    defaultdict(complex)
        self.G8 =    defaultdict(complex)
        self.G9 =    defaultdict(complex)
        self.G10=    defaultdict(complex)
        self.FACT  = defaultdict(complex)
        self.FACX  = defaultdict(complex)
        self.S12   = defaultdict(complex)
        self.S21   = defaultdict(complex)
        self.S32   = defaultdict(complex)
        self.S14   = defaultdict(complex)
        self.S23   = defaultdict(complex)
        self.S34   = defaultdict(complex)
        self.S32E  = defaultdict(complex)
        self.S34E  = defaultdict(complex)
        self.R1 =    defaultdict(complex)
        self.R2 =    defaultdict(complex)
        self.R3 =    defaultdict(complex)
        self.R4 =    defaultdict(complex)
        self.R5 =    defaultdict(complex)
        self.R6 =    defaultdict(complex)
        self.R7 =    defaultdict(complex)
        self.X1 =    defaultdict(complex)
        self.X2 =    defaultdict(complex)
        self.X3 =    defaultdict(complex)
        self.X4 =    defaultdict(complex)
        self.X5 =    defaultdict(complex)
        self.X6 =    defaultdict(complex)
        self.X7 =    defaultdict(complex)
        self.Y11 =   defaultdict(complex)
        self.Y12 =   defaultdict(complex)
        self.Y21 =   defaultdict(complex)
        self.Y31 =   defaultdict(complex)
        self.Y41 =   defaultdict(complex)
        self.Y42 =   defaultdict(complex)
        self.Y32 =   defaultdict(complex)
        self.Y22 =   defaultdict(complex)
        self.ATNA =  defaultdict(complex)
        self.ATNB =  defaultdict(complex)
        self.XKA  =  defaultdict(complex)
        self.XKB  =  defaultdict(complex)
        self.AA =    defaultdict(int)
        self.BB =    defaultdict(int)
        self.CC =    defaultdict(int)



#}}}

        self.IASYMP = True
        self.ID = 0
        self.ISYST = 0
        for I in range(8):
            if self.ISRC[I] == 1: self.ID = 1
        if self.verbose: self._log('ISRC => %s' % self.ISRC)

        if self.ID == 1: self.ISYST = 1
        if self.verbose: self._log('ISYST => %s' % self.ISYST)
        self.ID = 0
        if self.ISRC[8] == 1: self.ID=1
        if self.ISRC[9] == 1: self.ID=1
        if self.ID == 0 and self.ISYST == 0: return
        if self.ID == 1 and self.ISYST == 0: self.ISYST = 2
        if self.ID == 1 and self.ISYST == 1: self.ISYST = 3

        if self.verbose: self._log('ISYST => %s' % self.ISYST)

        self.DK = 2.0 * acos(-1.0) / self.XLENG
        if self.verbose: self._log('DK => %s' % self.DK)

        self.WVBM = self.OMEGA/self.VBMIN
        self.WVMM = (5.0/self.DEPTH)+self.XFAC*self.WVBM
        self.NK   = self.WVMM/self.DK
        self.NK   = int(self.NK+2)
        if self.verbose: self._log('WVBM => %s' % self.WVBM)
        if self.verbose: self._log('WVMM => %s' % self.WVMM)
        if self.verbose: self._log('NK => %s' % self.NK)

        self.V = 6.0 / self.DEPTH

        self.WVNO[1]  = complex( self.V, 0 )
        self.WVNO2[1] = self.WVNO[1] * self.WVNO[1]
        self.WVNO3[1] = self.WVNO2[1] * self.WVNO[1]
        self.AKL      = ( self.NK - 3 ) * self.DK + 0.218 * self.DK

        if self.verbose: self._log('WVNO[1] => %s' % self.WVNO[1])
        if self.verbose: self._log('WVNO2[1] => %s' % self.WVNO2[1])
        if self.verbose: self._log('WVNO3[1] => %s' % self.WVNO3[1])
        if self.verbose: self._log('AKL => %s' % self.AKL)


        if self.AKL > self.V: self.IASYMP = False
        if self.verbose: self._log('IASYMP => %s' % self.IASYMP)

        self.V = 2.5 / self.DEPTH
        if self.verbose: self._log('V => %s' % self.V)
        self.WVNO[0]  = complex( self.V , 0 )
        self.WVNO2[0] = self.WVNO[0] * self.WVNO[0]
        if self.verbose: self._log('WVNO[1] => %s' % self.WVNO[0])
        if self.verbose: self._log('WVNO2[1] => %s' % self.WVNO2[0])

        if self.CMAX > 0.0: self.MM = int(self.WVCN / self.DK) + 2
        if self.NK > self.MM: self.NK = self.MM
        if self.verbose: self._log('MM => %s' % self.MM)
        if self.verbose: self._log('NK => %s' % self.NK)
        if self.verbose: self._log('OMEGA => %s' % self.OMEGA)
        if self.verbose: self._log('OMEGA:%s NK:%s' % (self.OMEGA,self.NK))
        if self.verbose: self._log('\t\t\tOMEGA:%s NK:%s' % (self.OMEGA,self.NK))

        #
        # DETERMINE INDEX FOR INTEGRATION  I.E., WHEN KR IS GT. 3.0
        #
#{{{
        
        for J in range(1,50000):
            self.WVV = ( J - 1 ) * self.DK + 0.218 * self.DK
            self.AKR = self.WVV * self.DISTANCE
            if self.verbose: self._log('')
            if self.verbose: self._log('J => %s' % J)
            if self.verbose: self._log('RANGE => %s' % self.DISTANCE)
            if self.verbose: self._log('WVV => %s' % self.WVV)
            if self.verbose: self._log('AKR => %s' % self.AKR)
            if self.AKR <= 3.0: self.ICNT[1] = J + 2
            if self.AKR > 3.0: break
        if self.verbose: self._log('ICNT => %s' % self.ICNT)

        self.OM     = complex(self.OMEGA, -1 * self.ALPHA)
        self.FOURPI = 12.5663706 * self.OM * self.OM
        self.PI     = 3.141592654
        self.OM1    = 6.283185307
        self.OML    = 0.06283185307

        if self.verbose: self._log('OM => %s' % self.OM)
        if self.verbose: self._log('FOURPI => %s' % self.FOURPI)
        if self.verbose: self._log('PI => %s' % self.PI)
        if self.verbose: self._log('OM1 => %s' % self.OM1)
        if self.verbose: self._log('OML => %s' % self.OML)
#}}}

        #
        # COMPUTING THE ATTENUATION OF THE MDEIUM USING QA, QB
        #

        self.AT =complex() 
        #self._log('OM => %s' % self.OM)
        #self._log('abs(OM) => %s' % abs(self.OM))
        #self._log('OML => %s' % self.OML)
        #self._log('OM => %s' % self.OM)
        #self._log('OM1 => %s' % self.OM1)
        if (abs(self.OM) > self.OML ): self.AT=log(self.OM.real/self.OM1.real)/self.PI.real
        if self.verbose: self._log('AT => %s' % self.AT)
        if (abs(self.OM) < self.OML ):
            self.FA = sqrt(self.OML**2 + self.ALPHA**2)/self.OML
            self.FAC = float( self.FA )
            self.AT = clog(complex( self.OML, -self.ALPHA)/(self.OM1*self.FAC)) / self.PI

        if self.verbose: self._log('AT => %s' % self.AT)
        if self.verbose: self._log('FA => %s' % self.FA)
        if self.verbose: self._log('FAC => %s' % self.FAC)

        for I in range(self.MMAX+1):
            self.ATNA[I] = complex(1)
            self.ATNB[I] = complex(1)
            if self.verbose: self._log('ATNA[I] => %s' % self.ATNA[I])
            if self.verbose: self._log('ATNB[I] => %s' % self.ATNB[I])
            if self.QA[I] > 0.0:
                self.PAL = complex(self.QA[I],0) * self.AT
                self.BA = float(0.5 * self.QA[I])
                self.QAL = complex(0, self.BA)
                self.ATNA[I] = complex(1) + self.PAL + self.QAL
                if self.verbose: self._log('PAL => %s' % self.PAL )
                if self.verbose: self._log('BA => %s' % self.BA )
                if self.verbose: self._log('QAL => %s' % self.QAL )
                if self.verbose: self._log('ATNA[I] => %s' % self.ATNA[I])

            if self.QB[I] > 0.0:
                self.PAL = complex(self.QB[I],0) * self.AT
                self.BA = float(0.5 * self.QB[I])
                self.QAL = complex(0, self.BA)
                self.ATNB[I] = complex(1) + self.PAL + self.QAL
                if self.verbose: self._log('PAL => %s' % self.PAL )
                if self.verbose: self._log('BA => %s' % self.BA )
                if self.verbose: self._log('QAL => %s' % self.QAL )
                if self.verbose: self._log('ATNB[I] => %s' % self.ATNB[I])

            if self.verbose: self._log('A[I] => %s' % self.A[I])
            if self.verbose: self._log('ATNA[I] => %s' % self.ATNA[I])
            self.XKA[I] = self.OM / (self.A[I] * self.ATNA[I] )
            self.XKB[I] = self.OM / (self.B[I] * self.ATNB[I] )
            if self.verbose: self._log('XKA[I] => %s' % self.XKA[I])
            if self.verbose: self._log('XKB[I] => %s' % self.XKB[I])

        if self.verbose: self._log('PAL => %s' % self.PAL)
        if self.verbose: self._log('BA => %s' % self.BA)
        if self.verbose: self._log('QAL => %s' % self.QAL)
        if self.verbose: self._log('ATNA => %s' % self.ATNA)
        if self.verbose: self._log('ATNB => %s' % self.ATNB)
        if self.verbose: self._log('XKA => %s' % self.XKA)
        if self.verbose: self._log('XKB => %s' % self.XKB)


        if self.verbose: self._log('LMAX => %s' % self.LMAX)
        if self.verbose: self._log('RHO => %s' % self.RHO)
        if self.verbose: self._log('RHO[0] => %s' % self.RHO[0])
        if self.verbose: self._log('RHO[3] => %s' % self.RHO[3])
        if self.verbose: self._log('RHO[LMAX] => %s' % self.RHO[self.LMAX])
        if self.verbose: self._log('complex(RHO[LMAX]) => %s' % complex(self.RHO[self.LMAX]))
        self.FRHO = complex( self.RHO[self.LMAX] ) * self.FOURPI
        self.XKA2 = self.XKA[self.LMAX] * self.XKA[self.LMAX]
        self.XKB2 = self.XKB[self.LMAX] * self.XKB[self.LMAX]
        self.XKK = self.ATNB[self.LMAX] * self.ATNB[self.LMAX]
        self.FACXX = 1.0 / ( 12.5663706 * self.B[self.LMAX] ** 2 )

        if self.verbose: self._log('FRHO => %s' % self.FRHO)
        if self.verbose: self._log('XKA2 => %s' % self.XKA2)
        if self.verbose: self._log('XKB2 => %s' % self.XKB2)
        if self.verbose: self._log('XKK => %s' % self.XKK)
        if self.verbose: self._log('FACXX => %s' % self.FACXX)

        if self.NK < self.IB:
            self.IBLOCK = self.NK
            self.NBLOCK = 1
            self.LEFT = 0
        else:
            self.IBLOCK = self.IB
            self.NBLOCK = self.NK / self.IB
            self.LEFT = self.NK - self.NBLOCK * self.IB
            if self.LEFT: self.NBLOCK += 1

        if self.verbose: self._log('IBLOCK => %s' % self.IBLOCK)
        if self.verbose: self._log('NBLOCK => %s' % self.NBLOCK)
        if self.verbose: self._log('LEFT => %s' % self.LEFT)
        if self.verbose: self._log('\n\n')

        self.NBLK = 1
        self.LR = 0
        #}}}

        # Loop includes half-space also
        self.WVNO3[0] = self.WVNO2[0] * self.WVNO[0]
        self.NN1 = 2

        # 500 loop.
        if self.verbose: self._log("Loop 500: while: %s <= %s" % (self.NBLK,self.NBLOCK))
        while self.NBLK <= self.NBLOCK:
            if self.verbose: self._log("\tLoop 500: while: self.NBLK[%s] ?< self.NBOCLK[%s]" % (self.NBLK,self.NBLOCK))
        #{{{ A

            if self.verbose: self._log("Loop 500 section AR => self.ISYST:%s" % self.ISYST)
            #{{{ AR

            if self.NBLK == self.NBLOCK and self.LEFT != 0: self.IBLOCK = self.LEFT
            if self.NBLK > 1: self.NN1 = 0
            if self.verbose: self._log("\tLoop 500: NBKL = %s" % self.NBLK)
            if self.verbose: self._log("\tLoop 500: IBLOCK = %s" % self.IBLOCK)
            if self.verbose: self._log("\tLoop 500: N1 = %s" % self.NN1)

            if self.verbose: self._log("\tLoop 500: self.NN1:[%s] self.IBLOCK[%s]" % (self.NN1,self.IBLOCK))
            #self._log("\tLoop 500: self.NN1 = %s" % self.NN1)
            #self._log("\tLoop 500: self.IBLOCK = %s" % self.IBLOCK)
            for J in range(self.NN1,self.IBLOCK):
            #{{{ D
                #self._log("\tLoop 500: J = %s" % J)
                self.LR += 1
                self.WV = (self.LR - 1) * self.DK + 0.218 * self.DK
                self.WVNO[J] = complex(self.WV)
                self.WVNO2[J] = self.WVNO[J]**2
                self.WVNO3[J] = self.WVNO2[J] * self.WVNO[J]
                if self.verbose: self._log("\tLoop 500: J = %s" % J)
                if self.verbose: self._log("\tLoop 500: LR = %s" % self.LR)
                if self.verbose: self._log("\tLoop 500: WVNO[J] = %s" % self.WVNO[J])
                if self.verbose: self._log("\tLoop 500: WVNO2[J] = %s" % self.WVNO2[J])
                if self.verbose: self._log("\tLoop 500: WVNO3[J] = %s" % self.WVNO3[J])
                if self.verbose: self._log("\tLoop 500: WVNO2 = %s" % self.WVNO2)
                #if J == 4: exit()
                #self._log("\tLoop 500: WVNO2[J] = %s" % self.WVNO2[J])
            #}}} D

                self.R1 =    defaultdict(complex)
                self.R2 =    defaultdict(complex)
                self.R3 =    defaultdict(complex)
                self.R4 =    defaultdict(complex)
                self.R5 =    defaultdict(complex)
                self.R6 =    defaultdict(complex)
                self.R7 =    defaultdict(complex)
                self.X1 =    defaultdict(complex)
                self.X2 =    defaultdict(complex)
                self.X3 =    defaultdict(complex)
                self.X4 =    defaultdict(complex)
                self.X5 =    defaultdict(complex)
                self.X6 =    defaultdict(complex)
                self.X7 =    defaultdict(complex)

            #self._log("self.WVNO2:%s" % self.WVNO2)
            if self.verbose: self._log("self.WVNO2:%s" % self.WVNO2)
            for I in range(self.MMAX+1):
                for J in range(self.IBLOCK):
                #{{{ F
                    self.RA[I][J] = csqrt(self.WVNO2[J] - self.XKA[I]**2)
                    self.RB[I][J] = csqrt(self.WVNO2[J] - self.XKB[I]**2)
                    self.GM = float(self.B[I]) * self.WVNO[J] / self.OM
                    self.GAM[I][J] = self.GM * self.ATNB[I]
                    self.GAM[I][J] = 2.0 * self.GAM[I][J]**2
                    self.GAM2[I][J] = self.GAM[I][J]**2
                    self.GAM3[I][J] = self.GAM2[I][J] * self.GAM[I][J]
                    self.GAM4[I][J] =  self.GAM3[I][J] * self.GAM[I][J]
                    self.GAMM1[I][J] =  self.GAM[I][J] - complex(1)
                    self.GAMM2[I][J] =  self.GAMM1[I][J] * self.GAMM1[I][J]
                    self.GAMM3[I][J] =  self.GAMM1[I][J] * self.GAMM2[I][J]
                    self.GAMX1[I][J] =  self.GAM[I][J] * self.GAMM1[I][J]
                    self.GAMX2[I][J] =  self.GAM2[I][J] * self.GAMM2[I][J]
                    if self.verbose: self._log("")
                    if self.verbose: self._log("Loop 500: RA[%s][%s]=> %s" % (I,J,self.RA[I][J]))
                    if self.verbose: self._log("Loop 500: RB[%s][%s]=> %s" % (I,J,self.RB[I][J]))
                    if self.verbose: self._log("Loop 500: GAM[%s][%s]=> %s" % (I,J,self.GAM[I][J]))
                    if self.verbose: self._log("Loop 500: GAM2[%s][%s]=> %s" % (I,J,self.GAM2[I][J]))
                    if self.verbose: self._log("Loop 500: GAM3[%s][%s]=> %s" % (I,J,self.GAM3[I][J]))
                    if self.verbose: self._log("Loop 500: GAM4[%s][%s]=> %s" % (I,J,self.GAM4[I][J]))
                    if self.verbose: self._log("Loop 500: GAMM1[%s][%s]=> %s" % (I,J,self.GAMM1[I][J]))
                    if self.verbose: self._log("Loop 500: GAMM2[%s][%s]=> %s" % (I,J,self.GAMM2[I][J]))
                    if self.verbose: self._log("Loop 500: GAMM3[%s][%s]=> %s" % (I,J,self.GAMM3[I][J]))
                    if self.verbose: self._log("Loop 500: GAMX1[%s][%s]=> %s" % (I,J,self.GAMX1[I][J]))
                    if self.verbose: self._log("Loop 500: GAMX2[%s][%s]=> %s" % (I,J,self.GAMX2[I][J]))
                    if self.verbose: self._log("")
                #}}} F

            #if self.verbose: self._log("self.WVNO2:%s" % self.WVNO2)
            #exit()

            self.EX = 0.0
            for J in range(self.IBLOCK):
            #{{{ G
                self.EXE[J] = 0.0
                self.EXEL[J] = 0.0
                self.EXL[J] = 0.0
                self.EXLL[J] = 0.0
            #}}} G

            #}}} AR

            if self.verbose: self._log("Loop 500 section BR => self.ISYST:%s" % self.ISYST)
            if self.ISYST != 2:
            #{{{ BR
                # Computation for P-SV system starts
                #{{{ H
                for J in range(self.IBLOCK):
                    self.S14[J] = -complex(2) * self.WVNO[J] / self.FOURPI
                    self.S21[J] =  complex(2) * self.XKB2 / self.FRHO
                    self.S32[J] = self.WVNO[J] * 4.0 * self.XKA2 / self.FRHO
                    self.TT = (2.0 * self.B[self.LMAX] / self.A[self.LMAX]) ** 2 - 3.
                    #self._log("\t\t\t B(LMAX)=>%s" % self.B[self.LMAX])
                    #self._log("\t\t\t A(LMAX)=>%s" % self.A[self.LMAX])
                    #self._log("\t\t\t TT=>%s" % self.TT)
                    self.S34[J] = complex(2) * self.WVNO[J] * complex(self.TT) / self.FOURPI 
                    #self._log("\t\t\t WVNO[J]=>%s" % self.WVNO[J])
                    #self._log("\t\t\t TT=>%s" % self.TT)
                    #self._log("\t\t\t complex(self.TT)=>%s" % complex(self.TT))
                    #self._log("\t\t\t FOURPI=>%s" % self.FOURPI)
                    if self.verbose: self._log("")
                    if self.verbose: self._log("Loop 500: S14=> %s" % self.S14)
                    if self.verbose: self._log("Loop 500: S21=> %s" % self.S21)
                    if self.verbose: self._log("Loop 500: S32=> %s" % self.S32)
                    if self.verbose: self._log("Loop 500: S34=> %s" % self.S34)
                    if self.verbose: self._log("")
                #}}} H

                #     SET UP HALFSPACE BOUNDARY CONDITIONS
                #     JBDRY=-1 RIGID;  0= ELASTIC;  +1= FREE SURFACE
                #{{{ I 
                M = self.MMAX
                if self.verbose: self._log("M:[%s] " % M)
                if self.verbose: self._log("RHO:%s" % self.RHO)
                if self.JBDRY == 0:
                    for J in range(self.IBLOCK):
                        if self.verbose: self._log("M:%s J:%s" % (M,J))
                        #if self.verbose: self._log(" self.GAM2[M][J] :%s" % self.GAM2[M][J])
                        #if self.verbose: self._log(" self.RA[M][J] :%s" % self.RA[M][J])
                        #if self.verbose: self._log(" self.RB[M][J] :%s" % self.RB[M][J])
                        #if self.verbose: self._log(" self.WVNO2[J] :%s" % self.WVNO2[J])
                        #if self.verbose: self._log(" self.GAMM2[M][J] :%s" % self.GAMM2[M][J])
                        self.RH = self.RHO[M]**2
                        if self.verbose: self._log(" self.RH :%s" % self.RH)
                        self.R1[J] = self.RH * (-self.GAM2[M][J] * self.RA[M][J] * self.RB[M][J] + self.WVNO2[J] * self.GAMM2[M][J])
                        self.R2[J] = -1 * complex(self.RHO[M]) * self.WVNO2[J] * self.RA[M][J]
                        self.R3[J] = -1 * self.GAM[M][J] * self.RA[M][J] * self.RB[M][J] + self.WVNO2[J] * self.GAMM1[M][J]
                        self.R3[J] = complex(self.RHO[M]) * self.R3[J]
                        self.R4[J] = complex(self.RHO[M]) * self.WVNO2[J] * self.RB[M][J]
                        self.R5[J] = self.WVNO2[J] * ( self.WVNO2[J] - self.RA[M][J] * self.RB[M][J])
                        if self.verbose: self._log("Loop 500: R1[J]=> %s" % self.R1[J])
                        if self.verbose: self._log("Loop 500: R2[J]=> %s" % self.R2[J])
                        if self.verbose: self._log("Loop 500: R3[J]=> %s" % self.R3[J])
                        if self.verbose: self._log("Loop 500: R4[J]=> %s" % self.R4[J])
                        if self.verbose: self._log("Loop 500: R5[J]=> %s" % self.R5[J])
                        if self.verbose: self._log("\n")
                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R1[J] = complex(1)
                        self.R2[J] = complex()
                        self.R3[J] = complex()
                        self.R4[J] = complex()
                        self.R5[J] = complex()
                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R1[J] = complex() 
                        self.R2[J] = complex() 
                        self.R3[J] = complex() 
                        self.R4[J] = complex() 
                        self.R5[J] = complex(1)

                if self.verbose: self._log("Loop 500: RH=> %s" % self.RH)
                if self.verbose: self._log("Loop 500: R1=> %s" % self.R1)
                if self.verbose: self._log("Loop 500: R2=> %s" % self.R2)
                if self.verbose: self._log("Loop 500: R3=> %s" % self.R3)
                if self.verbose: self._log("Loop 500: R4=> %s" % self.R4)
                if self.verbose: self._log("Loop 500: R5=> %s" % self.R5)
                if self.verbose: self._log("")
                #}}} I 

                # MATRIX-MULTIPLICATION - LAYERWISE ++++> FROM BOTTOM UPWARD TO LMAX
                #{{{ J 
                if self.verbose: self._log("Loop 500: LMAX = >%s" % self.LMAX)
                if self.verbose: self._log("Loop 500: MMAX = >%s" % self.MMAX)
                #for I in range(self.LMAX, self.MMAX):
                for I in range(self.LMAX, self.MMAX):
                    # MATRIX MULTIPLICATION FROM BOTTOM LAYER UPWARD
                    M = self.MMAX -1 + self.LMAX - I
                    self.DPTH = self.D[M]
                    if self.verbose: self._log("Loop 500: I=> %s" % I)
                    if self.verbose: self._log("Loop 500: M=> %s" % M)
                    if self.verbose: self._log("Loop 500: DPTH=> %s" % self.DPTH)
                    if self.verbose: self._log("Loop 500: D=> %s" % self.D)
                    if self.verbose: self._log("Loop 500: IBLOCK=> %s" % self.IBLOCK)
                    if self.verbose: self._log("" )

                    if self.verbose: self._log("IBLOCK=%s" % self.IBLOCK)
                    #self._log("self.WVNO2 = %s " % self.WVNO2)
                    for J in range(self.IBLOCK):
                        if self.verbose: self._log("J=%s" % J)
                        #self._log("self.WVNO2=%s" % self.WVNO2[J])
                        #{{{ K
                        self.P = self.RA[M][J] * complex(self.DPTH)
                        self.Q = self.RB[M][J] * complex(self.DPTH)
                        #if self.verbose: self._log("Loop 500:   RA[%s][%s]=> %s" % (M,J,self.RA[M][J]))
                        #if self.verbose: self._log("Loop 500:   RB[%s][%s]=> %s" % (M,J,self.RB[M][J]))
                        #if self.verbose: self._log("Loop 500:   complex(DEPTH)=> %s" % complex(self.DEPTH))
                        if self.verbose: self._log("Loop 500:   P=> %s" % self.P)
                        if self.verbose: self._log("Loop 500:   Q=> %s" % self.Q)

                        self.A0 = 0.0
                        self.PR = self.P.real
                        self.PI = self.P.imag
                        self.QR = self.Q.real
                        self.QI = self.Q.imag

                        if self.verbose: self._log("Loop 500:   A0=> %s" % self.A0)
                        if self.verbose: self._log("Loop 500:   PR=> %s" % self.PR)
                        if self.verbose: self._log("Loop 500:   PI=> %s" % self.PI)
                        if self.verbose: self._log("Loop 500:   QR=> %s" % self.QR)
                        if self.verbose: self._log("Loop 500:   QI=> %s" % self.QI)

                        self.EPP = complex(cos(self.PI),sin(self.PI)) / 2.0
                        self.EPM = self.EPP.conjugate()
                        self.EQP = complex( cos(self.QI),sin(self.QI) ) / 2.0
                        self.EQM = self.EQP.conjugate()
                        self.EXB = self.QR
                        self.FAC = 0.0

                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   EPP=> %s" % self.EPP)
                        if self.verbose: self._log("Loop 500:   EPM=> %s" % self.EPM)
                        if self.verbose: self._log("Loop 500:   EQP=> %s" % self.EQP)
                        if self.verbose: self._log("Loop 500:   EQM=> %s" % self.EQM)
                        if self.verbose: self._log("Loop 500:   EXB=> %s" % self.EXB)

                        if self.PR < 15.0: self.FAC = exp(-(2.0) * self.PR).real
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)


                        self.COSP = self.EPP + self.FAC * self.EPM
                        self.SINP = self.EPP - self.FAC * self.EPM
                        self.W = self.SINP / self.RA[M][J]
                        self.X = self.RA[M][J] * self.SINP
                        self.FAC = 0.0

                        if self.verbose: self._log("Loop 500:   COSP=> %s" % self.COSP)
                        if self.verbose: self._log("Loop 500:   SINP=> %s" % self.SINP)
                        if self.verbose: self._log("Loop 500:   W=> %s" % self.W)
                        if self.verbose: self._log("Loop 500:   X=> %s" % self.X)

                        if self.QR < 15.0: self.FAC = exp( -(2.0) * self.QR).real
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)

                        if self.verbose: self._log("")

                        self.COSQL[M][J] = self.EQP + self.FAC * self.EQM
                        self.SINQ = self.EQP - self.FAC * self.EQM
                        self.YL[M][J] = self.SINQ / self.RB[M][J]
                        if self.verbose: self._log("Loop 500:   M=> %s" % M)
                        if self.verbose: self._log("Loop 500:   J=> %s" % J)
                        if self.verbose: self._log("Loop 500:   COSQL[M][J]=> %s" % self.COSQL[M][J])
                        if self.verbose: self._log("Loop 500:   SINQ=> %s" % self.SINQ)
                        if self.verbose: self._log("Loop 500:   YL[M][J]=> %s" % self.YL[M][J])
                        if self.verbose: self._log("")
                        #self._log("Loop 500:   self.YL[%s][%s]=> %s" % (M,J,self.YL[M][J]))
                        self.ZL[M][J] = self.RB[M][J] * self.SINQ
                        self.EXA = self.PR + self.QR
                        self.CPCQ = self.COSP * self.COSQL[M][J]
                        self.CPY = self.COSP * self.YL[M][J]
                        self.CPZ = self.COSP * self.ZL[M][J]
                        self.CQW = self.COSQL[M][J] * self.W
                        self.CQX = self.COSQL[M][J] * self.X
                        self.XY = self.X * self.YL[M][J]
                        self.XZ = self.X * self.ZL[M][J]
                        self.WY = self.W * self.YL[M][J]
                        self.WZ = self.W * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500:   ZL[M][J]=> %s" % self.ZL[M][J])
                        if self.verbose: self._log("Loop 500:   EXA=> %s" % self.EXA)
                        if self.verbose: self._log("Loop 500:   CPCQ=> %s" % self.CPCQ)
                        if self.verbose: self._log("Loop 500:   CPY=> %s" % self.CPY)
                        if self.verbose: self._log("Loop 500:   CPZ=> %s" % self.CPZ)
                        if self.verbose: self._log("Loop 500:   CQW=> %s" % self.CQW)
                        if self.verbose: self._log("Loop 500:   CQX=> %s" % self.CQX)
                        if self.verbose: self._log("Loop 500:   XY=> %s" % self.XY)
                        if self.verbose: self._log("Loop 500:   XZ=> %s" % self.XZ)
                        if self.verbose: self._log("Loop 500:   WY=> %s" % self.WY)
                        if self.verbose: self._log("Loop 500:   WZ=> %s" % self.WZ)
                        if self.verbose: self._log("")
                        self.FAC = 0.0
                        self.QMP = self.QR - self.PR
                        if self.verbose: self._log("Loop 500:   QMP=> %s" % self.QMP)

                        if self.QMP > -40.0: self.FAC = exp(self.QMP)
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)

                        self.COSQ = self.COSQL[M][J] * self.FAC
                        if self.verbose: self._log("Loop 500:   COSQ=> %s" % self.COSQ)
                        self.Y = self.FAC * self.YL[M][J]
                        self.Z = self.FAC * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500:   Y=> %s" % self.Y)
                        if self.verbose: self._log("Loop 500:   Z=> %s" % self.Z)
                        self.FAC = 0.0

                        if self.EXA < 60.0: self.A0 = exp(-self.EXA).real
                        if self.verbose: self._log("Loop 500:   A0=> %s" % self.A0)

                        self.RHO2  = self.RHO[M] * self.RHO[M]
                        if self.verbose: self._log("Loop 500:   RHO2=> %s" % self.RHO2)
                        if self.verbose: self._log("Loop 500:   WVNO2=%s" % self.WVNO2[J])
                        if self.verbose: self._log("Loop 500:   WVNO2=%s" % self.WVNO2)
                        if self.verbose: self._log("Loop 500:   J=%s" % J)
                        self.A0C  = 2.0 * (self.A0 - self.CPCQ)
                        if self.verbose: self._log("Loop 500:   A0C=> %s" % self.A0C)
                        self.XZ2  = self.XZ / self.WVNO2[J]
                        if self.verbose: self._log("Loop 500:   XZ2=> %s" % self.XZ2)
                        self.WY2  = self.WY * self.WVNO2[J]
                        if self.verbose: self._log("Loop 500:   WY2=> %s" % self.WY2)
                        self.TEMP = self.A0C * self.WVNO2[J] + self.XZ + self.WY2 * self.WVNO2[J]
                        if self.verbose: self._log("Loop 500:   TEMP=> %s" % self.TEMP)

                        self.CA15 = -self.TEMP / self.RHO2
                        if self.verbose: self._log("Loop 500:   CA15=> %s" % self.CA15)

                        self.TEMP = 0.5 * self.A0C * ( self.GAM[M][J] + self.GAMM1[M][J] ) + self.GAM[M][J] * self.XZ2 + self.GAMM1[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500:   TEMP=> %s" % self.TEMP)

                        if self.verbose: self._log("")


                        self.CA13= -self.TEMP / self.RHO[M]
                        if self.verbose: self._log("Loop 500:   CA13=> %s" % self.CA13)
                        if self.verbose: self._log("Loop 500:   GAMX1[M][J]=> %s" % self.GAMX1[M][J])
                        if self.verbose: self._log("Loop 500:   GAM2[M][J]=> %s" % self.GAM2[M][J])
                        if self.verbose: self._log("Loop 500:   GAMM2[M][J]=> %s" % self.GAMM2[M][J])

                        self.TEMP = self.A0C * self.GAMX1[M][J] + self.GAM2[M][J] * self.XZ2 + self.GAMM2[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500:   TEMP=> %s" % self.TEMP)

                        self.CA33 = self.A0 + self.TEMP*2
                        self.CA11 = self.CPCQ - self.TEMP
                        if self.verbose: self._log("Loop 500:   CA33=> %s" % self.CA33)
                        if self.verbose: self._log("Loop 500:   CA11=> %s" % self.CA11)
                        #if self.verbose: self._log("Loop 500:   CPCQ=> %s" % self.CPCQ)

                        self.TEMP = 0.5 * self.A0C * self.GAMX1[M][J] * ( self.GAM[M][J] + self.GAMM1[M][J]) + self.GAM3[M][J] * self.XZ2 + self.GAMM3[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500:   TEMP=> %s" % self.TEMP)

                        if self.verbose: self._log("")
                        self.CA31 = 2.0 * self.TEMP * self.RHO[M]
                        self.TEMP = self.A0C * self.GAMX2[M][J] + self.GAM4[M][J] * self.XZ2 +self.GAMM3[M][J] * self.GAMM1[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500:   CA31=> %s" % self.CA31)
                        if self.verbose: self._log("Loop 500:   TEMP=> %s" % self.TEMP)

                        self.CA51 = -self.RHO2 * self.TEMP / self.WVNO2[J]
                        self.CA14 = (-self.WVNO2[J] * self.CQW + self.CPZ) / self.RHO[M]
                        self.CA21 = (-self.GAMM2[M][J] * self.CQW + self.GAM2[M][J] * self.CPZ / self.WVNO2[J]) * self.RHO[M]
                        self.CA23 = -( self.WVNO2[J] * self.GAMM1[M][J] * self.CQW - self.GAM[M][J] * self.CPZ ) / self.WVNO2[J]
                        self.CA12 = (-self.CQX + self.WVNO2[J] * self.CPY ) / self.RHO[M]
                        self.CA32 = self.WVNO2[J] * ( self.GAM[M][J] * self.CQX / self.WVNO2[J] - self.GAMM1[M][J] * self.CPY ) * 2.0
                        self.CA41 = (-self.GAM2[M][J] * self.CQX / self.WVNO2[J] + self.GAMM2[M][J] * self.CPY) * self.RHO[M]
                        self.CA22 = self.CPCQ
                        self.CA24 = -self.WZ
                        self.CA42 = -self.XY
                        self.CA25 = self.CA14
                        self.CA55 = self.CA11
                        self.CA54 = self.CA21
                        self.CA53 = -self.CA31 / (2.0 * self.WVNO2[J])
                        self.CA52 = self.CA41
                        self.CA43 =  -self.CA32 / (2.0 * self.WVNO2[J])
                        self.CA45 = self.CA12
                        self.CA44 = self.CA22
                        self.CA34 = -2.0 * self.WVNO2[J] * self.CA23
                        self.CA35 = -2.0 * self.WVNO2[J] * self.CA13
                        if self.verbose: self._log("Loop 500:   CA51=> %s" % self.CA51)
                        if self.verbose: self._log("Loop 500:   CA14=> %s" % self.CA14)
                        if self.verbose: self._log("Loop 500:   CA21=> %s" % self.CA21)
                        if self.verbose: self._log("Loop 500:   CA23=> %s" % self.CA23)
                        if self.verbose: self._log("Loop 500:   CA12=> %s" % self.CA12)
                        if self.verbose: self._log("Loop 500:   CA32=> %s" % self.CA32)
                        if self.verbose: self._log("Loop 500:   CA41=> %s" % self.CA41)
                        if self.verbose: self._log("Loop 500:   CA22=> %s" % self.CA22)
                        if self.verbose: self._log("Loop 500:   CA24=> %s" % self.CA24)
                        if self.verbose: self._log("Loop 500:   CA42=> %s" % self.CA42)
                        if self.verbose: self._log("Loop 500:   CA25=> %s" % self.CA25)
                        if self.verbose: self._log("Loop 500:   CA55=> %s" % self.CA55)
                        if self.verbose: self._log("Loop 500:   CA54=> %s" % self.CA54)
                        if self.verbose: self._log("Loop 500:   CA53=> %s" % self.CA53)
                        if self.verbose: self._log("Loop 500:   CA52=> %s" % self.CA52)
                        if self.verbose: self._log("Loop 500:   CA43=> %s" % self.CA43)
                        if self.verbose: self._log("Loop 500:   CA45=> %s" % self.CA45)
                        if self.verbose: self._log("Loop 500:   CA44=> %s" % self.CA44)
                        if self.verbose: self._log("Loop 500:   CA34=> %s" % self.CA34)
                        if self.verbose: self._log("Loop 500:   CA35=> %s" % self.CA35)



                        self.EXE[J] = self.EXE[J] + self.EXA

                        self.ER1 = self.R1[J] * self.CA11 + self.R2[J] * self.CA21 + self.R3[J] * self.CA31 + self.R4[J] * self.CA41 + self.R5[J] * self.CA51
                        self.ER2 = self.R1[J] * self.CA12 + self.R2[J] * self.CA22 + self.R3[J] * self.CA32 + self.R4[J] * self.CA42 + self.R5[J] * self.CA52
                        self.ER3 = self.R1[J] * self.CA13 + self.R2[J] * self.CA23 + self.R3[J] * self.CA33 + self.R4[J] * self.CA43 + self.R5[J] * self.CA53
                        self.ER4 = self.R1[J] * self.CA14 + self.R2[J] * self.CA24 + self.R3[J] * self.CA34 + self.R4[J] * self.CA44 + self.R5[J] * self.CA54
                        self.ER5 = self.R1[J] * self.CA15 + self.R2[J] * self.CA25 + self.R3[J] * self.CA35 + self.R4[J] * self.CA45 + self.R5[J] * self.CA55




                        if self.verbose: self._log("Loop 500:   ER1=> %s" % self.ER1)
                        if self.verbose: self._log("Loop 500:   ER2=> %s" % self.ER2)
                        if self.verbose: self._log("Loop 500:   ER3=> %s" % self.ER3)
                        if self.verbose: self._log("Loop 500:   ER4=> %s" % self.ER4)
                        if self.verbose: self._log("Loop 500:   ER5=> %s" % self.ER5)
                        if self.verbose: self._log("")

                        self.TESTT = 0.0
                        if abs(self.ER1.real) > self.TESTT: self.TESTT = abs(self.ER1.real)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER2.real) > self.TESTT: self.TESTT = abs(self.ER2.real)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER3.real) > self.TESTT: self.TESTT = abs(self.ER3.real)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER4.real) > self.TESTT: self.TESTT = abs(self.ER4.real)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER5.real) > self.TESTT: self.TESTT = abs(self.ER5.real)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER1.imag) > self.TESTT: self.TESTT = abs(self.ER1.imag)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER2.imag) > self.TESTT: self.TESTT = abs(self.ER2.imag)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER3.imag) > self.TESTT: self.TESTT = abs(self.ER3.imag)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER4.imag) > self.TESTT: self.TESTT = abs(self.ER4.imag)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)
                        if abs(self.ER5.imag) > self.TESTT: self.TESTT = abs(self.ER5.imag)
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)

                        if self.TESTT < 1.0E-30: self.TESTT = 1.0
                        if self.verbose: self._log("Loop 500:   TESTT=> %s" % self.TESTT)


                        self.TEST = 0.0

                        self.FAC = abs(self.ER1) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER2) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER3) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER4) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER5) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)
                        self.TEST = self.TEST * self.TESTT
                        if self.TEST < 1.0E-30: self.TEST = 1.0
                        if self.verbose: self._log("Loop 500:   FAC=> %s" % self.FAC)
                        if self.verbose: self._log("Loop 500:   TEST=> %s" % self.TEST)

                        if self.verbose: self._log("self.TEST:%s" % self.TEST)

                        self.XNORM = 1.0 / self.TEST
                        if self.verbose: self._log("Loop 500:   XNORM=> %s" % self.XNORM)
                        self.EXA  = -log(self.XNORM)
                        if self.verbose: self._log("Loop 500:   EXA=> %s" % self.EXA)
                        self.R1[J] = self.ER1 * self.XNORM
                        self.R2[J] = self.ER2 * self.XNORM
                        self.R3[J] = self.ER3 * self.XNORM
                        self.R4[J] = self.ER4 * self.XNORM
                        self.R5[J] = self.ER5 * self.XNORM

                        self.EXE[J] = self.EXE[J] + self.EXA
                        self.EXEL[J] = self.EXEL[J] + self.EXB
                        if self.verbose: self._log("Loop 500:   EXE[%s]=> %s" % (J,self.EXE[J]))
                        if self.verbose: self._log("Loop 500:   EXEL[%s]=> %s" % (J,self.EXEL[J]))
                        #}}} K
                #}}} J 

                # STORE PREVIOUS MATRIX INTO A SEPARATE BUFFER +++++++++++>
                #{{{ L
                for J in range(self.IBLOCK):
                    self.X1[J] = self.R1[J]
                    self.X2[J] = self.R2[J]
                    self.X3[J] = self.R3[J]
                    self.X4[J] = self.R4[J]
                    self.X5[J] = self.R5[J]
                    self.EXL[J] = self.EXE[J]
                    self.EXLL[J] = self.EXEL[J]
                    for K in range(4):
                        for L in range(4):
                            self.ZX[J][K][L] = complex(0.0,0.0)
                            if K == L: self.ZX[J][K][L] = complex(1.0,0.0)
                if self.verbose: self._log("Loop 500: X1=> %s" % self.X1)
                if self.verbose: self._log("Loop 500: X2=> %s" % self.X2)
                if self.verbose: self._log("Loop 500: X3=> %s" % self.X3)
                if self.verbose: self._log("Loop 500: X4=> %s" % self.X4)
                if self.verbose: self._log("Loop 500: X5=> %s" % self.X5)
                if self.verbose: self._log("Loop 500: EXL=> %s" % self.EXL)
                if self.verbose: self._log("Loop 500: EXLL=> %s" % self.EXLL)
                if self.verbose: self._log("Loop 500: ZX=> %s" % self.ZX)
                #}}} L

                # MATRIX MULTIPLICATION FROM THE SOURCE INTERFACE UPWARD TO THE SURFACE
                #{{{ O
                for I in range(self.LMAX):
                    M = self.LMAX - I -1 
                    self.DPTH = self.D[M]
                    if self.verbose: self._log("Loop 500 14: I=> %s" % I)
                    if self.verbose: self._log("Loop 500 14: LMAX=> %s" % self.LMAX)
                    if self.verbose: self._log("Loop 500 14: M=> %s" % M)
                    if self.verbose: self._log("Loop 500 14: DPTH=> %s" % self.DPTH)

                    #{{{ P
                    for J in range(self.IBLOCK):
                        self.P = self.RA[M][J] * self.DPTH
                        self.Q = self.RB[M][J] * self.DPTH

                        if self.verbose: self._log("Loop 500 15:   J=> %s" % J)
                        if self.verbose: self._log("Loop 500 15:   P=> %s" % self.P)
                        if self.verbose: self._log("Loop 500 15:   P=> %s" % self.Q)

                        A0 = 0.0
                        self.PR = self.P.real
                        if self.verbose: self._log("Loop 500 15: PR=> %s" % self.PR)
                        self.PI = self.P.imag
                        if self.verbose: self._log("Loop 500 15: PI=> %s" % self.PI)
                        self.QR = self.Q.real
                        if self.verbose: self._log("Loop 500 15: QR=> %s" % self.QR)
                        self.QI = self.Q.imag
                        if self.verbose: self._log("Loop 500 15: QI=> %s" % self.QI)
                        self.EPP = complex(cos(self.PI),sin(self.PI)) / 2.0
                        if self.verbose: self._log("Loop 500 15: EPP=> %s" % self.EPP)
                        self.EPM = self.EPP.conjugate()
                        if self.verbose: self._log("Loop 500 15: EPM=> %s" % self.EPM)
                        self.EQP = complex(cos(self.QI),sin(self.QI)) / 2.0
                        if self.verbose: self._log("Loop 500 15: EQP=> %s" % self.EQP)
                        self.EQM = self.EQP.conjugate()
                        if self.verbose: self._log("Loop 500 15: EQM=> %s" % self.EQM)
                        self.EX = self.PR
                        if self.verbose: self._log("Loop 500 15: EX=> %s" % self.EX)
                        self.EXB = self.QR
                        if self.verbose: self._log("Loop 500 15: EXB=> %s" % self.EXB)
                        self.FAC = 0.0
                        if self.verbose: self._log("")

                        if self.PR < 15.0: self.FAC = exp(-2*self.PR).real 
                        if self.verbose: self._log("Loop 500 15: FAC=> %s" % self.FAC)

                        self.COSP = self.EPP + self.FAC * self.EPM
                        if self.verbose: self._log("Loop 500 15: COSP=> %s" % self.COSP)
                        self.SINP = self.EPP - self.FAC * self.EPM
                        if self.verbose: self._log("Loop 500 15: SINP=> %s" % self.SINP)
                        self.W = self.SINP / self.RA[M][J]
                        if self.verbose: self._log("Loop 500 15: W=> %s" % self.W)
                        self.X = self.RA[M][J] * self.SINP
                        if self.verbose: self._log("Loop 500 15: X=> %s" % self.X)
                        self.FAC = 0.0
                        if self.QR < 15.0: self.FAC = exp(-2.0 * self.QR)
                        if self.verbose: self._log("Loop 500 15: FAC=> %s" % self.FAC)
                        self.COSQL[M][J] = self.EQP + self.FAC * self.EQM
                        if self.verbose: self._log("Loop 500 15: COSQL[M][J]=> %s" % self.COSQL[M][J])
                        self.SINQ = self.EQP - self.FAC * self.EQM
                        if self.verbose: self._log("Loop 500 15: SINQ=> %s" % self.SINQ)
                        self.YL[M][J] = self.SINQ / self.RB[M][J]
                        if self.verbose: self._log("Loop 500 15: YL[M][J]=> %s" % self.YL[M][J])
                        self.ZL[M][J] = self.RB[M][J] * self.SINQ
                        if self.verbose: self._log("Loop 500 15: ZL[M][J]=> %s" % self.ZL[M][J])
                        self.EXA = self.PR + self.QR
                        if self.verbose: self._log("Loop 500 15: EXA=> %s" % self.EXA)
                        self.CPCQ = self.COSP * self.COSQL[M][J]
                        if self.verbose: self._log("Loop 500 15: CPCQ=> %s" % self.CPCQ)
                        self.CPY = self.COSP * self.YL[M][J]
                        if self.verbose: self._log("Loop 500 15: CPY=> %s" % self.CPY)
                        self.CPZ = self.COSP * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500 15: CPZ=> %s" % self.CPZ)
                        self.CQW = self.COSQL[M][J] * self.W
                        if self.verbose: self._log("Loop 500 15: CQW=> %s" % self.CQW)
                        self.CQX = self.COSQL[M][J] * self.X
                        if self.verbose: self._log("Loop 500 15: CQX=> %s" % self.CQX)
                        self.XY = self.X * self.YL[M][J]
                        if self.verbose: self._log("Loop 500 15: XY=> %s" % self.XY)
                        self.XZ = self.X * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500 15: XZ=> %s" % self.XZ)
                        self.WY = self.W * self.YL[M][J]
                        if self.verbose: self._log("Loop 500 15: WY=> %s" % self.WY)
                        self.WZ = self.W * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500 15: WZ=> %s" % self.WZ)
                        if self.verbose: self._log("")

                        self.FAC = 0.0
                        self.QMP = self.QR - self.PR
                        if self.verbose: self._log("Loop 500 15: QMP=> %s" % self.QMP)
                        if self.QMP > -40.0: self.FAC = exp(self.QMP).real
                        if self.verbose: self._log("Loop 500 15: FAC=> %s" % self.FAC)
                        self.COSQ = self.COSQL[M][J] * self.FAC
                        if self.verbose: self._log("Loop 500 15: COSQ=> %s" % self.COSQ)
                        self.Y = self.FAC * self.YL[M][J]
                        if self.verbose: self._log("Loop 500 15: Y=> %s" % self.Y)
                        self.Z = self.FAC * self.ZL[M][J]
                        if self.verbose: self._log("Loop 500 15: Z=> %s" % self.Z)

                        if self.verbose: self._log("")

                        self.FAC=0.0
                        if self.EXA < 60.0: self.A0 = exp(-self.EXA)
                        self.RHO2 = self.RHO[M] ** 2 
                        self.A0C  = 2.0 * (self.A0 - self.CPCQ)
                        self.XZ2  = self.XZ / self.WVNO2[J]
                        self.WY2  = self.WY * self.WVNO2[J]
                        if self.verbose: self._log("Loop 500 15: A0=> %s" % self.A0)
                        if self.verbose: self._log("Loop 500 15: RHO2=> %s" % self.RHO2)
                        if self.verbose: self._log("Loop 500 15: A0C=> %s" % self.A0C)
                        if self.verbose: self._log("Loop 500 15: XZ2=> %s" % self.XZ2)
                        if self.verbose: self._log("Loop 500 15: WY2=> %s" % self.WY2)
                        self.TEMP = self.A0C * self.WVNO2[J] + self.XZ + self.WY2 * self.WVNO2[J]
                        if self.verbose: self._log("Loop 500 15: TEMP=> %s" % self.TEMP)
                        self.CA15 = -self.TEMP / self.RHO2
                        if self.verbose: self._log("Loop 500 15: CA15=> %s" % self.CA15)
                        self.TEMP = 0.5 * self.A0C * (self.GAM[M][J] + self.GAMM1[M][J]) + self.GAM[M][J] * self.XZ2 + self.GAMM1[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500 15: TEMP=> %s" % self.TEMP)
                        self.CA13 = -self.TEMP / self.RHO[M]
                        if self.verbose: self._log("Loop 500 15: CA13=> %s" % self.CA13)
                        self.TEMP = self.A0C * self.GAMX1[M][J] + self.GAM2[M][J] * self.XZ2 + self.GAMM2[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500 15: TEMP=> %s" % self.TEMP)
                        self.CA33 = self.A0 + self.TEMP + self.TEMP
                        if self.verbose: self._log("Loop 500: CA33=> %s" % self.CA33)
                        self.CA11 = self.CPCQ - self.TEMP
                        if self.verbose: self._log("Loop 500 15: CA11=> %s" % self.CA11)
                        self.TEMP = 0.5 * self.A0C * self.GAMX1[M][J] * (self.GAM[M][J] + self.GAMM1[M][J])+ self.GAM3[M][J] * self.XZ2 + self.GAMM3[M][J] * self.WY2
                        if self.verbose: self._log("Loop 500 15: TEMP=> %s" % self.TEMP)
                        self.CA31 = 2.0 * self.TEMP * self.RHO[M]
                        if self.verbose: self._log("Loop 500 15: CA31=> %s" % self.CA31)
                        self.TEMP = self.A0C * self.GAMX2[M][J] + self.GAM4[M][J] * self.XZ2 + self.GAMM3[M][J] * self.GAMM1[M][J] * self.WY2
                        #if self.verbose: self._log("Loop 500 15: GAMX2=> %s" % self.GAMX2[M][J])
                        #if self.verbose: self._log("Loop 500 15: GAM4=> %s" % self.GAM4[M][J])
                        #if self.verbose: self._log("Loop 500 15: XZ2=> %s" % self.XZ2)
                        #if self.verbose: self._log("Loop 500 15: GAMM3=> %s" % self.GAMM3[M][J])
                        #if self.verbose: self._log("Loop 500 15: GAMM1=> %s" % self.GAMM1[M][J])
                        #if self.verbose: self._log("Loop 500 15: WY2=> %s" % self.WY2)
                        if self.verbose: self._log("Loop 500 15: TEMP=> %s" % self.TEMP)
                        self.CA51 = -self.RHO2 * self.TEMP / self.WVNO2[J]
                        self.CA14 = (-self.WVNO2[J] * self.CQW + self.CPZ) / self.RHO[M]
                        self.CA21 = (-self.GAMM2[M][J] * self.CQW + self.GAM2[M][J] * self.CPZ / self.WVNO2[J]) * self.RHO[M]
                        self.CA23 = -(self.WVNO2[J] * self.GAMM1[M][J] * self.CQW - self.GAM[M][J] * self.CPZ) / self.WVNO2[J]
                        self.CA12 = (-self.CQX + self.WVNO2[J] * self.CPY) / self.RHO[M]
                        self.CA32 = self.WVNO2[J] * (self.GAM[M][J] * self.CQX / self.WVNO2[J] - self.GAMM1[M][J] * self.CPY) * 2.0
                        self.CA41 = (-self.GAM2[M][J] * self.CQX / self.WVNO2[J] + self.GAMM2[M][J] * self.CPY) * self.RHO[M]
                        self.CA22 = self.CPCQ
                        self.CA24 = -self.WZ
                        self.CA42 = -self.XY
                        self.CA25 = self.CA14
                        self.CA55 = self.CA11
                        self.CA54 = self.CA21
                        self.CA53 = -self.CA31 / (2.0 * self.WVNO2[J])
                        self.CA52 = self.CA41
                        self.CA43 =  -self.CA32 / (2.0 * self.WVNO2[J])
                        self.CA45 = self.CA12
                        self.CA44 = self.CA22
                        self.CA34 = -2.0 * self.WVNO2[J] * self.CA23
                        self.CA35 = -2.0 * self.WVNO2[J] * self.CA13

                        if self.verbose: self._log("Loop 500 15: CA51=> %s" % self.CA51)
                        if self.verbose: self._log("Loop 500 15: CA14=> %s" % self.CA14)
                        if self.verbose: self._log("Loop 500 15: CA21=> %s" % self.CA21)
                        if self.verbose: self._log("Loop 500 15: CA23=> %s" % self.CA23)
                        if self.verbose: self._log("Loop 500 15: CA12=> %s" % self.CA12)
                        if self.verbose: self._log("Loop 500 15: CA32=> %s" % self.CA32)
                        if self.verbose: self._log("Loop 500 15: CA41=> %s" % self.CA41)
                        if self.verbose: self._log("Loop 500 15: CA22=> %s" % self.CA22)
                        if self.verbose: self._log("Loop 500 15: CA24=> %s" % self.CA24)
                        if self.verbose: self._log("Loop 500 15: CA42=> %s" % self.CA42)
                        if self.verbose: self._log("Loop 500 15: CA25=> %s" % self.CA25)
                        if self.verbose: self._log("Loop 500 15: CA55=> %s" % self.CA55)
                        if self.verbose: self._log("Loop 500 15: CA54=> %s" % self.CA54)
                        if self.verbose: self._log("Loop 500 15: CA53=> %s" % self.CA53)
                        if self.verbose: self._log("Loop 500 15: CA52=> %s" % self.CA52)
                        if self.verbose: self._log("Loop 500 15: CA43=> %s" % self.CA43)
                        if self.verbose: self._log("Loop 500 15: CA45=> %s" % self.CA45)
                        if self.verbose: self._log("Loop 500 15: CA44=> %s" % self.CA44)
                        if self.verbose: self._log("Loop 500 15: CA34=> %s" % self.CA34)
                        if self.verbose: self._log("Loop 500 15: CA35=> %s" % self.CA35)

                        if self.verbose: self._log("")

                        self.EXE[J] = self.EXE[J] + self.EXA
                        if self.verbose: self._log("Loop 500: EXE[J]=> %s" % self.EXE[J])
                        self.ER1 = self.R1[J] * self.CA11 + self.R2[J] * self.CA21 + self.R3[J] * self.CA31 + self.R4[J] * self.CA41 + self.R5[J] * self.CA51
                        self.ER2 = self.R1[J] * self.CA12 + self.R2[J] * self.CA22 + self.R3[J] * self.CA32 + self.R4[J] * self.CA42 + self.R5[J] * self.CA52
                        self.ER3 = self.R1[J] * self.CA13 + self.R2[J] * self.CA23 + self.R3[J] * self.CA33 + self.R4[J] * self.CA43 + self.R5[J] * self.CA53
                        self.ER4 = self.R1[J] * self.CA14 + self.R2[J] * self.CA24 + self.R3[J] * self.CA34 + self.R4[J] * self.CA44 + self.R5[J] * self.CA54
                        self.ER5 = self.R1[J] * self.CA15 + self.R2[J] * self.CA25 + self.R3[J] * self.CA35 + self.R4[J] * self.CA45 + self.R5[J] * self.CA55
                        if self.verbose: self._log("Loop 500: ER1=> %s" % self.ER1)
                        if self.verbose: self._log("Loop 500: ER2=> %s" % self.ER2)
                        if self.verbose: self._log("Loop 500: ER3=> %s" % self.ER3)
                        if self.verbose: self._log("Loop 500: ER4=> %s" % self.ER4)
                        if self.verbose: self._log("Loop 500: ER5=> %s" % self.ER5)

                        self.TESTT = 0.0
                        if abs(self.ER1.real) > self.TESTT: self.TESTT = abs(self.ER1.real)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER2.real) > self.TESTT: self.TESTT = abs(self.ER2.real)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER3.real) > self.TESTT: self.TESTT = abs(self.ER3.real)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER4.real) > self.TESTT: self.TESTT = abs(self.ER4.real)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER5.real) > self.TESTT: self.TESTT = abs(self.ER5.real)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)

                        if abs(self.ER1.imag) > self.TESTT: self.TESTT = abs(self.ER1.imag)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER2.imag) > self.TESTT: self.TESTT = abs(self.ER2.imag)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER3.imag) > self.TESTT: self.TESTT = abs(self.ER3.imag)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER4.imag) > self.TESTT: self.TESTT = abs(self.ER4.imag)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if abs(self.ER5.imag) > self.TESTT: self.TESTT = abs(self.ER5.imag)
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)
                        if self.TESTT < 1.0e-30: self.TESTT = 1.0
                        if self.verbose: self._log("Loop 500: TESTT=> %s" % self.TESTT)

                        if self.verbose: self._log("")


                        self.TEST = 0.0
                        self.FAC = abs(self.ER1) / self.TESTT
                        if self.verbose: self._log("Loop 500: FAC=> %s" % self.FAC)
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER2) / self.TESTT
                        if self.verbose: self._log("Loop 500: FAC=> %s" % self.FAC)
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER3) / self.TESTT
                        if self.verbose: self._log("Loop 500: FAC=> %s" % self.FAC)
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER4) / self.TESTT
                        if self.verbose: self._log("Loop 500: FAC=> %s" % self.FAC)
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.FAC = abs(self.ER5) / self.TESTT
                        if self.verbose: self._log("Loop 500: FAC=> %s" % self.FAC)
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.TEST = self.TEST * self.TESTT
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        if self.TEST < 1.0e-30: self.TEST = 1.0
                        if self.verbose: self._log("Loop 500: TEST=> %s" % self.TEST)
                        self.XNORM = 1.0 / self.TEST
                        if self.verbose: self._log("Loop 500: XNORM=> %s" % self.XNORM)
                        self.EXA = -log(self.XNORM)
                        if self.verbose: self._log("Loop 500: EXA=> %s" % self.EXA)

                        self.R1[J] = self.ER1 * self.XNORM
                        self.R2[J] = self.ER2 * self.XNORM
                        self.R3[J] = self.ER3 * self.XNORM
                        self.R4[J] = self.ER4 * self.XNORM
                        self.R5[J] = self.ER5 * self.XNORM
                        if self.verbose: self._log("Loop 500: R1[J]=> %s" % self.R1[J])
                        if self.verbose: self._log("Loop 500: R2[J]=> %s" % self.R2[J])
                        if self.verbose: self._log("Loop 500: R3[J]=> %s" % self.R3[J])
                        if self.verbose: self._log("Loop 500: R4[J]=> %s" % self.R4[J])
                        if self.verbose: self._log("Loop 500: R5[J]=> %s" % self.R5[J])


                        self.EXE[J] = self.EXE[J] + self.EXA
                        self.EXEL[J] = self.EXEL[J] + self.EXB

                        if self.verbose: self._log("Loop 500: EXE[J]=> %s" % self.EXE[J])
                        if self.verbose: self._log("Loop 500: EXEL[J]=> %s" % self.EXEL[J])

                        #
                        #   self.  PROPAGATION OF 4X4 HASKEL MATRIX UPWARD ++++++++++++>
                        #
                        self.EXL[J] = self.EXL[J] + self.EX
                        if self.verbose: self._log("Loop 500: EXL=> %s" % self.EXL[J])
                        self.AA11 = self.GAM[M][J] * self.COSP - self.GAMM1[M][J] * self.COSQ
                        self.AA12 = -self.GAMM1[M][J] * self.W + self.GAM[M][J] * self.Z / self.WVNO2[J]
                        self.AA13 = -(self.COSP - self.COSQ) / self.RHO[M]
                        self.AA14 = (self.WVNO2[J] * self.W - self.Z) / self.RHO[M]
                        self.AA21 =  self.GAM[M][J] * self.X - self.WVNO2[J] * self.GAMM1[M][J] * self.Y
                        self.AA22 =  -self.GAMM1[M][J] * self.COSP + self.GAM[M][J] * self.COSQ
                        self.AA23 = (-self.X + self.WVNO2[J] * self.Y) / self.RHO[M]
                        self.AA24 =  -self.WVNO2[J] * self.AA13
                        self.AA31 =  self.RHO[M] * self.GAM[M][J] * self.GAMM1[M][J] * (self.COSP - self.COSQ)
                        self.AA32 = self.RHO[M] * ( -self.GAMM2[M][J] * self.W + self.GAM2[M][J] * self.Z / self.WVNO2[J])
                        self.AA34 = -self.WVNO2[J] * self.AA12
                        self.AA33 = self.AA22
                        self.AA41 = self.RHO[M] * (self.GAM2[M][J] * self.X / self.WVNO2[J] - self.GAMM2[M][J] * self.Y)
                        self.AA42 = -self.AA31 / self.WVNO2[J]
                        self.AA43 = -self.AA21 / self.WVNO2[J]
                        self.AA44 = self.AA11

                        if self.verbose: self._log("")
                        if self.verbose: self._log("Loop 500: J=> %s" % J)
                        if self.verbose: self._log("Loop 500: EXL[J]=> %s" % self.EXL[J])
                        if self.verbose: self._log("Loop 500: AA11=> %s" % self.AA11)
                        if self.verbose: self._log("Loop 500: AA12=> %s" % self.AA12)
                        if self.verbose: self._log("Loop 500: AA13=> %s" % self.AA13)
                        if self.verbose: self._log("Loop 500: AA14=> %s" % self.AA14)
                        if self.verbose: self._log("Loop 500: AA21=> %s" % self.AA21)
                        if self.verbose: self._log("Loop 500: AA22=> %s" % self.AA22)
                        if self.verbose: self._log("Loop 500: AA23=> %s" % self.AA23)
                        if self.verbose: self._log("Loop 500: AA24=> %s" % self.AA24)
                        if self.verbose: self._log("Loop 500: AA31=> %s" % self.AA31)
                        if self.verbose: self._log("Loop 500: AA32=> %s" % self.AA32)
                        if self.verbose: self._log("Loop 500: AA34=> %s" % self.AA34)
                        if self.verbose: self._log("Loop 500: AA33=> %s" % self.AA33)
                        if self.verbose: self._log("Loop 500: AA41=> %s" % self.AA41)
                        if self.verbose: self._log("Loop 500: AA42=> %s" % self.AA42)
                        if self.verbose: self._log("Loop 500: AA43=> %s" % self.AA43)
                        if self.verbose: self._log("Loop 500: AA44=> %s" % self.AA44)
                        #self._log("%s %s %s %s %s %s %s %s" %(self.ZX[J][0][0],self.AA11,self.ZX[J][0][1],self.AA21,self.ZX[J][0][2],self.AA31,self.ZX[J][0][3],self.AA41))
                        self.D11 = self.ZX[J][0][0] * self.AA11 + self.ZX[J][0][1] * self.AA21 + self.ZX[J][0][2] * self.AA31 + self.ZX[J][0][3] * self.AA41
                        self.D21 = self.ZX[J][1][0] * self.AA11 + self.ZX[J][1][1] * self.AA21 + self.ZX[J][1][2] * self.AA31 + self.ZX[J][1][3] * self.AA41
                        self.D31 = self.ZX[J][2][0] * self.AA11 + self.ZX[J][2][1] * self.AA21 + self.ZX[J][2][2] * self.AA31 + self.ZX[J][2][3] * self.AA41
                        self.D41 = self.ZX[J][3][0] * self.AA11 + self.ZX[J][3][1] * self.AA21 + self.ZX[J][3][2] * self.AA31 + self.ZX[J][3][3] * self.AA41
                        self.D12 = self.ZX[J][0][0] * self.AA12 + self.ZX[J][0][1] * self.AA22 + self.ZX[J][0][2] * self.AA32 + self.ZX[J][0][3] * self.AA42
                        self.D22 = self.ZX[J][1][0] * self.AA12 + self.ZX[J][1][1] * self.AA22 + self.ZX[J][1][2] * self.AA32 + self.ZX[J][1][3] * self.AA42
                        self.D32 = self.ZX[J][2][0] * self.AA12 + self.ZX[J][2][1] * self.AA22 + self.ZX[J][2][2] * self.AA32 + self.ZX[J][2][3] * self.AA42
                        self.D42 = self.ZX[J][3][0] * self.AA12 + self.ZX[J][3][1] * self.AA22 + self.ZX[J][3][2] * self.AA32 + self.ZX[J][3][3] * self.AA42
                        self.D13 = self.ZX[J][0][0] * self.AA13 + self.ZX[J][0][1] * self.AA23 + self.ZX[J][0][2] * self.AA33 + self.ZX[J][0][3] * self.AA43
                        self.D23 = self.ZX[J][1][0] * self.AA13 + self.ZX[J][1][1] * self.AA23 + self.ZX[J][1][2] * self.AA33 + self.ZX[J][1][3] * self.AA43
                        self.D33 = self.ZX[J][2][0] * self.AA13 + self.ZX[J][2][1] * self.AA23 + self.ZX[J][2][2] * self.AA33 + self.ZX[J][2][3] * self.AA43
                        self.D43 = self.ZX[J][3][0] * self.AA13 + self.ZX[J][3][1] * self.AA23 + self.ZX[J][3][2] * self.AA33 + self.ZX[J][3][3] * self.AA43
                        self.D14 = self.ZX[J][0][0] * self.AA14 + self.ZX[J][0][1] * self.AA24 + self.ZX[J][0][2] * self.AA34 + self.ZX[J][0][3] * self.AA44
                        self.D24 = self.ZX[J][1][0] * self.AA14 + self.ZX[J][1][1] * self.AA24 + self.ZX[J][1][2] * self.AA34 + self.ZX[J][1][3] * self.AA44
                        self.D34 = self.ZX[J][2][0] * self.AA14 + self.ZX[J][2][1] * self.AA24 + self.ZX[J][2][2] * self.AA34 + self.ZX[J][2][3] * self.AA44
                        self.D44 = self.ZX[J][3][0] * self.AA14 + self.ZX[J][3][1] * self.AA24 + self.ZX[J][3][2] * self.AA34 + self.ZX[J][3][3] * self.AA44

                        if self.verbose: self._log("")
                        if self.verbose: self._log("Loop 500: D11=> %s" % self.D11)
                        if self.verbose: self._log("Loop 500: D12=> %s" % self.D12)
                        if self.verbose: self._log("Loop 500: D13=> %s" % self.D13)
                        if self.verbose: self._log("Loop 500: D14=> %s" % self.D14)
                        if self.verbose: self._log("Loop 500: D21=> %s" % self.D21)
                        if self.verbose: self._log("Loop 500: D22=> %s" % self.D22)
                        if self.verbose: self._log("Loop 500: D23=> %s" % self.D23)
                        if self.verbose: self._log("Loop 500: D24=> %s" % self.D24)
                        if self.verbose: self._log("Loop 500: D31=> %s" % self.D31)
                        if self.verbose: self._log("Loop 500: D32=> %s" % self.D32)
                        if self.verbose: self._log("Loop 500: D33=> %s" % self.D33)
                        if self.verbose: self._log("Loop 500: D34=> %s" % self.D34)
                        if self.verbose: self._log("Loop 500: D41=> %s" % self.D41)
                        if self.verbose: self._log("Loop 500: D42=> %s" % self.D42)
                        if self.verbose: self._log("Loop 500: D43=> %s" % self.D43)
                        if self.verbose: self._log("Loop 500: D44=> %s" % self.D44)
                        if self.verbose: self._log("")

                        self.ZX[J][0][0] = self.D11
                        self.ZX[J][0][1] = self.D12
                        self.ZX[J][0][2] = self.D13
                        self.ZX[J][0][3] = self.D14
                        self.ZX[J][1][0] = self.D21
                        self.ZX[J][1][1] = self.D22
                        self.ZX[J][1][2] = self.D23
                        self.ZX[J][1][3] = self.D24
                        self.ZX[J][2][0] = self.D31
                        self.ZX[J][2][1] = self.D32
                        self.ZX[J][2][2] = self.D33
                        self.ZX[J][2][3] = self.D34
                        self.ZX[J][3][0] = self.D41
                        self.ZX[J][3][1] = self.D42
                        self.ZX[J][3][2] = self.D43
                        self.ZX[J][3][3] = self.D44

                    #}}} P
                #}}} O

                #for J in self.ZX:
                #    for L in self.ZX[J]:
                #        self._log(' [%s][%s] %s ' % (J,L,[self.ZX[J][L][W] for W in self.ZX[J][L]]))


                #{{{ Q
                if self.ISRC[2] == 1 or self.ISRC[3] == 1:
                    for K in range(self.IBLOCK):
                        self.Y11[K] =  self.X1[K] * self.ZX[K][1][0] + self.X2[K] * self.ZX[K][2][0] - self.WVNO2[K] * self.X3[K] *self.ZX[K][3][0]
                        self.Y12[K] =  self.X1[K] * self.ZX[K][1][1] + self.X2[K] * self.ZX[K][2][1] - self.WVNO2[K] * self.X3[K] * self.ZX[K][3][1]
                        self.Y31[K] = -self.X2[K] * self.ZX[K][0][0] - self.X3[K] * self.ZX[K][1][0] + self.X5[K] * self.ZX[K][3][0]
                        self.Y32[K] = -self.X2[K] * self.ZX[K][0][1] - self.X3[K] * self.ZX[K][1][1] + self.X5[K] * self.ZX[K][3][1]

                if self.ISRC[0] == 1 or self.ISRC[1] == 1 or self.ISRC[4] == 1 or self.ISRC[5] == 1 or self.ISRC[6] == 1 or self.ISRC[7] == 1:
                    for K in range(self.IBLOCK):
                        self.Y21[K] = -self.X1[K] * self.ZX[K][0][0] + self.X3[K] * self.ZX[K][2][0] + self.X4[K] * self.ZX[K][3][0]
                        self.Y22[K] = -self.X1[K] * self.ZX[K][0][1] + self.X3[K] * self.ZX[K][2][1] + self.X4[K] * self.ZX[K][3][1]
                        self.Y41[K] = self.WVNO2[K] * self.X3[K] * self.ZX[K][0][0] - self.X4[K] * self.ZX[K][1][0] - self.X5[K] * self.ZX[K][2][0]
                        self.Y42[K] = self.WVNO2[K] * self.X3[K] * self.ZX[K][0][1] - self.X4[K] * self.ZX[K][1][1] - self.X5[K] * self.ZX[K][2][1]
                        #self._log(" TEST  X1[K]=> %s" % self.X1[K])
                        #self._log(" TEST  ZX[K,0,0]=> %s" % self.ZX[K][0][0])
                        #self._log(" TEST  X3[K]=> %s" % self.X3[K])
                        #self._log(" TEST  ZX[K,2,0]=> %s" % self.ZX[K][2][0])
                        #self._log(" TEST  X4[K]=> %s" % self.X4[K])
                        #self._log(" TEST  ZX[K,3,0]=> %s" % self.ZX[K][3][0])
                        #self._log(" TEST  Y21[K]=> %s" % self.Y21[K])

                #exit()
                #self._log("")
                #self._log("Loop 500 16 17 : Y11=> %s" % self.Y11)
                #self._log("Loop 500 16 17 : Y12=> %s" % self.Y12)
                #self._log("Loop 500 16 17 : Y31=> %s" % self.Y31)
                #self._log("Loop 500 16 17 : Y32=> %s" % self.Y32)
                #self._log("Loop 500 16 17 : Y21=> %s" % self.Y21)
                #self._log("Loop 500 16 17 : Y22=> %s" % self.Y22)
                #self._log("Loop 500 16 17 : Y41=> %s" % self.Y41)
                #self._log("Loop 500 16 17 : Y42=> %s" % self.Y42)
                #sys.exit()
                #}}} Q

            #}}} BR

            if self.verbose: self._log("Loop 500 section CR => self.ISYST:%s" % self.ISYST)
            if self.ISYST != 1:
            #{{{ CR
                if self.verbose: self._log("\t ##### 300 #####")

                #
                # SH-MOTION STARTS HERE +++++>
                #
                for J in range(self.IBLOCK):
                    self.S32E[J] = 2.0 * self.XKA2 * self.WVNO[J] / self.FRHO
                    self.S34E[J] = 4.0 * self.WVNO[J] * self.XKA2 / (self.FOURPI * self.XKB2)
                    if self.verbose: self._log("")
                    if self.verbose: self._log("Loop 500 300 18:   S32E[J]=> %s" % self.S32E[J])
                    if self.verbose: self._log("Loop 500 300 18:   S34E[J]=> %s" % self.S34E[J])

                if self.JBDRY == 0:
                    for J in range(self.IBLOCK):
                        #if self.verbose: self._log("\nJ:%s \nMMAX:%s \nRB:%s \nRHO:%s" % (J,self.MMAX,self.RB,self.RHO))
                        self.R6[J] = (self.RHO[self.MMAX]) * self.RB[self.MMAX][J]
                        self.R7[J] = complex(1) / (self.B[self.MMAX] * self.ATNB[self.MMAX]) ** 2
                        if self.verbose: self._log("")
                        if self.verbose: self._log("Loop 500 300 19: R6[J]=> %s" % self.R6[J])
                        if self.verbose: self._log("Loop 500 300 19: R7[J]=> %s" % self.R7[J])

                elif self.JBDRY < 0:
                    for J in range(self.IBLOCK):
                        self.R6[J] = complex(1)
                        self.R7[J] = complex()
                        if self.verbose: self._log("")
                        if self.verbose: self._log("Loop 500 300 20: R6[J]=> %s" % self.R6[J])
                        if self.verbose: self._log("Loop 500 300 20: R7[J]=> %s" % self.R7[J])

                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R6[J] = complex()
                        self.R7[J] = complex(1)
                        if self.verbose: self._log("")
                        if self.verbose: self._log("Loop 500 300 21: R6[J]=> %s" % self.R6[J])
                        if self.verbose: self._log("Loop 500 300 21: R7[J]=> %s" % self.R7[J])
            #}}} CR


            if self.verbose: self._log("Loop 500 section DR => self.ISYST:%s" % self.ISYST)
            if self.ISYST != 3 and self.ISYST != 1:
            #{{{ DR
                for M in range(self.MMAX+1):
                    for J in range(self.IBLOCK):
                        #self._log("Loop 500:   self.YL=> %s" % self.YL)
                        #self._log("Loop 500:   [%s][%s]=> %s" % (M,J))
                        self.P = self.RA[M][J] * self.D[M]
                        self.Q = self.RB[M][J] * self.D[M]
                        self.A0 = 0.0
                        self.PR = self.P.real
                        self.PI = self.P.imag
                        self.QR = self.Q.real
                        self.QI = self.Q.imag
                        self.EPP = complex(cos(self.PI),sin(self.PI)) / 2.0
                        self.EPM = self.EPP.conjugate()
                        self.EQP = complex(cos(self.QI),sin(self.QI)) / 2.0
                        self.EQM = self.EQP.conjugate()
                        self.EXEL[J] = self.EXEL[J] + self.QR
                        if M > self.LMAX: self.EXLL[J] = self.EXLL[J] + self.QR

                        self.FAC = 0.0
                        if self.PR > 15.0: self.FAC = exp(-2.0 * self.PR)
                        self.COSP = self.EPP + self.FAC * self.EPM
                        self.SINP = self.EPP - self.FAC * self.EPM

                        self.FAC=0.0
                        if self.SNGL[self.QR] < 15.0: self.FAC = exp(-2.0 * self.QR)
                        self.COSQL[M][J] = self.EQP + self.FAC * self.EQM
                        self.SINQ = self.EQP - self.FAC * self.EQM
                        self.YL[M][J] = self.SINQ / self.RB[M][J]
                        self.ZL[M][J] = self.RB[M][J] * self.SINQ
                        #self._log("Loop 500 302 22:   P => %s" % self.P)
                        #self._log("Loop 500 302 22:   Q => %s" % self.Q)
                        #self._log("Loop 500 302 22:   A0 => %s" % self.A0)
                        #self._log("Loop 500 302 22:   PR => %s" % self.PR)
                        #self._log("Loop 500 302 22:   PI => %s" % self.PI)
                        #self._log("Loop 500 302 22:   QR => %s" % self.QR)
                        #self._log("Loop 500 302 22:   QI => %s" % self.QI)
                        #self._log("Loop 500 302 22:   EPP => %s" % self.EPP)
                        #self._log("Loop 500 302 22:   EPM => %s" % self.EPM)
                        #self._log("Loop 500 302 22:   EQP => %s" % self.EQP)
                        #self._log("Loop 500 302 22:   EQM => %s" % self.EQM)
                        #self._log("Loop 500 302 22:   EXEL[J] => %s" % self.EXEL[J])
                        #self._log("Loop 500 302 22:   EXLL[J] => %s" % self.EXLL[J])
                        #self._log("Loop 500 302 22:   COSP => %s" % self.COSP)
                        #self._log("Loop 500 302 22:   SINP => %s" % self.SINP)
                        #self._log("Loop 500 302 22:   COSQL[M][J] => %s" % self.COSQL[M][J])
                        #self._log("Loop 500 302 22:   SINQ => %s" % self.SINQ)
                        #self._log("Loop 500 302 22:   YL[M][J] => %s" % self.YL[M][J])
                        #self._log("Loop 500 302 22:   ZL[M][J] => %s" % self.ZL[M][J])

                        sys.exit()

            #}}} DR

            if self.verbose: self._log("Loop 500 section ER => self.ISYST:%s" % self.ISYST)
            if self.ISYST != 1:
            #{{{ ER
                if self.verbose: self._log(" ##### AFTER 302 ####")
                if self.verbose: self._log("Loop 500:   LMAX=> %s" % self.LMAX)
                if self.verbose: self._log("Loop 500:   MMAX=> %s" % self.MMAX)

                for I in range(self.LMAX, self.MMAX):

                    if self.verbose: self._log("Loop 500:   I => %s" % I)
                    M = self.MMAX + self.LMAX - I -1
                    self.HH = self.RHO[M] * (self.B[M] * self.ATNB[M])**2
                    if self.verbose: self._log("Loop 500:   M => %s" % M)
                    if self.verbose: self._log("Loop 500:   self.HH => %s" % self.HH)
                    for J in range(self.IBLOCK):
                        self.YLI = self.YL[M][J] / self.HH
                        self.ZLI = self.ZL[M][J] * self.HH
                        self.D11 = self.R6[J]
                        self.D12 = self.R7[J]
                        self.R6[J] = self.D11 * self.COSQL[M][J] + self.D12 * self.ZLI
                        self.R7[J] = self.D11 * self.YLI + self.D12 * self.COSQL[M][J]
                        if self.verbose: self._log("Loop 500:     J => %s" % J)
                        if self.verbose: self._log("Loop 500:     YLI=> %s" % self.YLI)
                        if self.verbose: self._log("Loop 500:     ZLI=> %s" % self.ZLI)
                        if self.verbose: self._log("Loop 500:     D11=> %s" % self.D11)
                        if self.verbose: self._log("Loop 500:     D12=> %s" % self.D12)
                        if self.verbose: self._log("Loop 500:     R6[J]=> %s" % self.R6[J])
                        if self.verbose: self._log("Loop 500:     R7[J]=> %s" % self.R7[J])


                if self.verbose: self._log("")
                for J in range(self.IBLOCK):
                    self.X6[J] = self.R6[J]
                    self.X7[J] = self.R7[J]
                    if self.verbose: self._log("Loop 500:     J => %s" % J)
                    if self.verbose: self._log("Loop 500:     X6[J]=> %s" % self.X6[J])
                    if self.verbose: self._log("Loop 500:     X7[J]=> %s" % self.X7[J])

                if self.verbose: self._log("")
                for I in range(self.LMAX):
                #for M in range(self.LMAX,-1,-1):
                    M = self.LMAX - I - 1
                    if self.verbose: self._log("Loop 500:   I => %s" % I)
                    if self.verbose: self._log("Loop 500:   M => %s" % M)
                    self.HH = self.RHO[M] * (self.B[M] * self.ATNB[M])**2
                    if self.verbose: self._log("Loop 500:   self.HH => %s" % self.HH)
                    for J in range(self.IBLOCK):
                        self.YLI = self.YL[M][J] / self.HH
                        self.ZLI = self.ZL[M][J] * self.HH
                        self.D11 = self.R6[J]
                        self.D12 = self.R7[J]
                        self.R6[J] = self.D11 * self.COSQL[M][J] + self.D12 * self.ZLI
                        self.R7[J] = self.D11 * self.YLI + self.D12 * self.COSQL[M][J]
                        #self._log("Loop 500:        J=> %s" % J
                        #self._log("Loop 500:        YLI=> %s" % self.YLI)
                        #self._log("Loop 500:        ZLI=> %s" % self.ZLI)
                        #self._log("Loop 500:        D11=> %s" % self.D11)
                        #self._log("Loop 500:        D12=> %s" % self.D12)
                        #self._log("Loop 500:        R6[J]=> %s" % self.R6[J])
                        #self._log("Loop 500:        R7[J]=> %s" % self.R7[J])
                #self._log("Loop 500:   R6=> %s" % self.R6)
                #if self.verbose: self._log("Loop 500:   M=> %s" % M)
                #if self.verbose: self._log("Loop 500:   HH=> %s" % self.HH)
                #if self.verbose: self._log("Loop 500:   YLI=> %s" % self.YLI)
                #if self.verbose: self._log("Loop 500:   ZLI=> %s" % self.ZLI)
                #if self.verbose: self._log("Loop 500:   D11=> %s" % self.D11)
                #if self.verbose: self._log("Loop 500:   D12=> %s" % self.D12)
                #if self.verbose: self._log("Loop 500:   R6=> %s" % self.R6)
                #if self.verbose: self._log("Loop 500:   R7=> %s" % self.R7)
            #}}} ER

            if self.verbose: self._log("Loop 500 section FR => self.ISYST:%s" % self.ISYST)

            if self.verbose: self._log("")
            if self.verbose: self._log("######## AFTER 301 #######")
            #{{{ FR
            if self.verbose: self._log("Loop 500 28: IBLOCK => %s" % self.IBLOCK)
            for J in range(self.IBLOCK):
                if self.verbose: self._log("Loop 500 28: J => %s" % J)
                self.FACT[J] = 0.0
                self.ELJ = self.EXE[J] - self.EXL[J]
                if self.ELJ < 55.0: self.FACT[J] = exp(-self.ELJ).real
                self.FACT0 = 0.0
                self.ELJ = self.EXLL[J] - self.EXEL[J]
                if self.SNGL[self.ELJ] > -40.0: self.FACT0 = exp(self.ELJ).real
                self.FACX[J] = self.FACT0 * self.FACXX
                if self.verbose: self._log("Loop 500:   FACT[J]=> %s" % self.FACT[J])
                if self.verbose: self._log("Loop 500:   ELJ=> %s" % self.ELJ)
                if self.verbose: self._log("Loop 500:   FACT0=> %s" % self.FACT0)
                if self.verbose: self._log("Loop 500:   FACX[J]=> %s" % self.FACX[J])

            if self.verbose: self._log("")

            if self.ISRC[0] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t29")
                for J in range(self.IBLOCK):
                    self.G = self.S32[J] * self.Y21[J] + self.S34[J] * self.Y41[J]
                    self.G1[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t29:     G = > %s" % self.G)
                    if self.verbose: self._log("\t29:     G1[J] = > %s" % self.G1[J])

            if self.ISRC[1] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t30")
                for J in range(self.IBLOCK):
                    self.G = self.S32[J] * self.Y22[J] + self.S34[J] * self.Y42[J]
                    self.G2[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t30:     G = > %s" % self.G)
                    if self.verbose: self._log("\t30:     G2[J] = > %s" % self.G2[J])

            if self.ISRC[2] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t31")
                for J in range(self.IBLOCK):
                    self.G = self.S21[J] * self.Y11[J]
                    self.G3[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t31:     G = > %s" % self.G)
                    if self.verbose: self._log("\t31:     G3[J] = > %s" % self.G3[J])

            if self.ISRC[3] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t32")
                for J in range(self.IBLOCK):
                    self.G = self.S21[J] * self.Y12[J]
                    self.G4[J]= -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t32:     G = > %s" % self.G)
                    if self.verbose: self._log("\t32:     G4[J] = > %s" % self.G4[J])

            if self.ISRC[4] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t33")
                for J in range(self.IBLOCK):
                    self.G = self.S14[J] * self.Y41[J]
                    self.G5[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t33:     G = > %s" % self.G)
                    if self.verbose: self._log("\t33:     G5[J] = > %s" % self.G5[J])

            if self.ISRC[5] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t34")
                for J in range(self.IBLOCK):
                    self.G = self.S14[J] * self.Y42[J]
                    self.G6[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t34:     G = > %s" % self.G)
                    if self.verbose: self._log("\t34:     G6[J] = > %s" % self.G6[J])

            if self.ISRC[6] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t35")
                for J in range(self.IBLOCK):
                    self.G = self.S32E[J] * self.Y21[J] + self.S34E[J] * self.Y41[J]
                    self.G7[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t35:     Y41[J] = > %s" % self.Y41[J])
                    if self.verbose: self._log("\t35:     S34E[J] = > %s" % self.S34E[J])
                    if self.verbose: self._log("\t35:     Y21[J] = > %s" % self.Y21[J])
                    if self.verbose: self._log("\t35:     S32E[J] = > %s" % self.S32E[J])
                    if self.verbose: self._log("\t35:     G = > %s" % self.G)
                    if self.verbose: self._log("\t35:     FACT[J] = > %s" % self.FACT[J])
                    if self.verbose: self._log("\t35:     R1[J] = > %s" % self.R1[J])
                    if self.verbose: self._log("\t35:     G7[J] = > %s" % self.G7[J])
            #exit()

            if self.ISRC[7] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t36")
                for J in range(self.IBLOCK):
                    self.G = self.S32E[J] * self.Y22[J] + self.S34E[J] * self.Y42[J]
                    self.G8[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.verbose: self._log("\t36:     G = > %s" % self.G)
                    if self.verbose: self._log("\t36:     FACT[J] = > %s" % self.FACT[J])
                    if self.verbose: self._log("\t36:     R1[J] = > %s" % self.R1[J])
                    if self.verbose: self._log("\t36:     G8[J] = > %s" % self.G8[J])

            if self.ISRC[8] == 1:
                if self.verbose: self._log("")
                if self.verbose: self._log("\t37")
                for J in range(self.IBLOCK):
                    self.G = 2.0 * self.X6[J] / self.RHO[self.LMAX]
                    self.G9[J] = self.G * self.FACX[J] / ( self.R6[J] * self.XKK )
                    if self.verbose: self._log("\t37:     G = > %s" % self.G)
                    if self.verbose: self._log("\t37:     G9[J] = > %s" % self.G9[J])


            if self.ISRC[9] == 1:
                for J in range(self.IBLOCK):
                    self.G = -2.0 * self.WVNO[J] * self.X7[J] * (self.B[self.LMAX] * self.ATNB[self.LMAX])**2
                    self.G10[J] = self.G * self.FACX[J] / (self.R6[J] * self.XKK)

            self.N11 = (self.NBLK - 1) * self.IB + 1
            self.N22 = self.N11 + self.IBLOCK - 1


            self.ICN = self.ICNT[1]
            if self.verbose: self._log("Loop 600 [ sub WVINT ]:   ICN=> %s" % self.ICN)

            if self.verbose: self._log("Loop 600:   WVNO=> %s" % self.WVNO)
            if self.verbose: self._log("Loop 600:   R=> %s" % self.DISTANCE)
            if self.verbose: self._log("Loop 600:   NBLK=> %s" % self.NBLOCK)
            if self.verbose: self._log("Loop 600:   IS=> %s" % self.NBLOCK)
            if self.verbose: self._log("Loop 600:   IBLOCK=> %s" % self.IBLOCK)
            if self.verbose: self._log("Loop 600:   ICN=> %s" % self.ICNT[2])
            if self.verbose: self._log("Loop 600:   ISRC=> %s" % self.ISRC)
            if self.verbose: self._log("Loop 600:   N11=> %s" % self.N11)
            if self.verbose: self._log("Loop 600:   N22=> %s" % self.N22)
            if self.verbose: self._log("Loop 600:   I=> %s" % I)

            if self.verbose: self._log("")

            #{{{ sub WVINT
            #
            #SUBROUTINE WVINT(WVNO,R,NBL,IS,IBL,ICNT,ISRC,N11,N22,II)
            #     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #
            #     R     = DISTANCES IN KM FOR WAVENUMBER INTEGRATION
            #     NBL   = INDEX AT WHICH [KR] BECOMES LESS TAN AND EQUAL TO 3.0
            #     IBL   = NO OF WAVENUMBERS TO BE INTEGRETAED. LESS THN IB
            #     ISRC  = CONTROL FOR REQUIRED COMPONENT
            #     N11   = INDEX FOR FIRST WAVENUMBER IN ARRAY
            #     N22   = INDEX FOR LAST  WAVENUMBER IN ARRAY
            #     II    = INDEX FOR DISTANCES
            #
            #     PROGRAMMED BY CHANDAN K. SAIKIA
            #     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            self.IB = 100
            self.NR = 100
            self.PI = acos(-1.0)
            self.NN1 = 1

            if self.verbose: self._log("")
            if self.verbose: self._log("++++++++++++++++++++ WVINT ++++++++++++++++++")
            for J in range(self.IB):
                self.SUME[J] = complex()
                self.SUMC[J] = complex()

            #{{{ INIT
            #self._log('NBL = %s' % self.NBLK)
            if self.verbose: self._log("\t WVINT      NBLK => %s " % self.NBLK)
            if self.NBLK == 1:
                self.NNN1 = 3
                if self.verbose: self._log("\t WVINT      NNN1 => %s " % self.NNN1)
                if self.verbose: self._log("\t WVINT      IASYMP => %s " % self.IASYMP)
                #     FIND THE ASYMPTOTIC COEFFICIENTS ++++ NECESSARY ONLY WHEN IASYMP
                #     IS .TRUE.                        +++++++++++++>
                if self.IASYMP:

                    #{{{ setup(R)
                    #----------------------------------------------------------
                    #       JNKM = R INTEGRAL EXP(-KH) KSUP M J SUB N (KR) DK
                    #       THE R IN FRONT TAKES INTO ACCOUNT THE 1/R IN THE
                    #       DO 300 OF SUBROUTINE WVINT
                    #----------------------------------------------------------

                    if self.verbose: self._log("")
                    if self.verbose: self._log("%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%")
                    if self.verbose: self._log("\t SETUP:   R       => %s" % self.DISTANCE    )
                    self.DIST  = csqrt(self.DISTANCE**2 + self.DEPTH**2).real
                    self.DIST3 = self.DIST**3
                    self.DIST5 = self.DIST**5
                    self.DIST7 = self.DIST**7
                    self.RZ    = self.DISTANCE*self.DEPTH
                    self.Z2    = self.DEPTH*self.DEPTH
                    self.R2    = self.DISTANCE*self.DISTANCE
                    self.R3    = self.DISTANCE*self.R2
                    self.Z3    = self.DEPTH*self.Z2
                    self.RZ2   = self.DISTANCE*self.Z2
                    self.RZ3   = self.DISTANCE*self.Z3
                    self.ZOR   = self.DEPTH/self.DIST
                    self.ZOR2  = self.ZOR*self.ZOR
                    self.ZOR3  = self.ZOR*self.ZOR2
                    self.J0K0  = self.DISTANCE/self.DIST
                    self.J0K1  = self.RZ/self.DIST3
                    self.J0K2  = (2.0 * self.RZ2 - self.R3) / self.DIST5
                    self.J0K3  = (6.0 * self.RZ3 - 9.0 * self.DEPTH * self.R3) / self.DIST7
                    self.J1K0  = 1.0 - self.DEPTH / self.DIST
                    self.J1K1  = self.R2 / self.DIST3
                    self.J1K2  = 3.0 * self.DEPTH * self.R2 / self.DIST5
                    self.J1K3  = 3.0 * self.R2 * (4.0 * self.Z2 - self.R2) / self.DIST7
                    self.J2K0  = ( 1.0 - 2.0 * self.ZOR + self.ZOR2 ) * (self.DIST / self.DISTANCE)
                    self.J2K1  = (2.0 - 3.0 * self.ZOR + self.ZOR3) / self.DISTANCE
                    self.J2K2  = 3.0 * self.R3 / self.DIST5
                    self.J2K3  = 15.0 * self.DEPTH * self.R3 / self.DIST7

                    if self.verbose: self._log("\t SETUP:   DIST     => %s" % self.DIST  )
                    if self.verbose: self._log("\t SETUP:   DIST3    => %s" % self.DIST3 )
                    if self.verbose: self._log("\t SETUP:   DIST5    => %s" % self.DIST5 )
                    if self.verbose: self._log("\t SETUP:   DIST7    => %s" % self.DIST7 )
                    if self.verbose: self._log("\t SETUP:   RZ       => %s" % self.RZ    )
                    if self.verbose: self._log("\t SETUP:   Z2       => %s" % self.Z2    )
                    if self.verbose: self._log("\t SETUP:   R2       => %s" % self.R2    )
                    if self.verbose: self._log("\t SETUP:   R3       => %s" % self.R3    )
                    if self.verbose: self._log("\t SETUP:   Z3       => %s" % self.Z3    )
                    if self.verbose: self._log("\t SETUP:   RZ2      => %s" % self.RZ2   )
                    if self.verbose: self._log("\t SETUP:   RZ3      => %s" % self.RZ3   )
                    if self.verbose: self._log("\t SETUP:   ZOR      => %s" % self.ZOR   )
                    if self.verbose: self._log("\t SETUP:   ZOR2     => %s" % self.ZOR2  )
                    if self.verbose: self._log("\t SETUP:   ZOR3     => %s" % self.ZOR3  )
                    if self.verbose: self._log("\t SETUP:   J0K0     => %s" % self.J0K0  )
                    if self.verbose: self._log("\t SETUP:   J0K1     => %s" % self.J0K1  )
                    if self.verbose: self._log("\t SETUP:   J0K2     => %s" % self.J0K2  )
                    if self.verbose: self._log("\t SETUP:   J0K3     => %s" % self.J0K3  )
                    if self.verbose: self._log("\t SETUP:   J1K0     => %s" % self.J1K0  )
                    if self.verbose: self._log("\t SETUP:   J1K1     => %s" % self.J1K1  )
                    if self.verbose: self._log("\t SETUP:   J1K2     => %s" % self.J1K2  )
                    if self.verbose: self._log("\t SETUP:   J1K3     => %s" % self.J1K3  )
                    if self.verbose: self._log("\t SETUP:   J2K0     => %s" % self.J2K0  )
                    if self.verbose: self._log("\t SETUP:   J2K1     => %s" % self.J2K1  )
                    if self.verbose: self._log("\t SETUP:   J2K2     => %s" % self.J2K2  )
                    if self.verbose: self._log("\t SETUP:   J2K3     => %s" % self.J2K3  )

                    if self.verbose: self._log("%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%")
                    if self.verbose: self._log("")


                    #}}} setup(R)

                    if I == 0:
                    #{{{
                        if self.verbose: self._log("")
                        self.W1 = self.WVNO[0].real
                        self.W2 = self.WVNO[1].real
                        if self.verbose: self._log("\t WVINT:   W1       => %s" % self.W1)
                        if self.verbose: self._log("\t WVINT:   W2       => %s" % self.W2)
                        if self.ISRC[0] == 1: self._solu300(self.G1[0],self.G1[1],self.W1,self.W2,0)
                        if self.ISRC[1] == 1: self._solu200(self.G2[0],self.G2[1],self.W1,self.W2,1)
                        if self.ISRC[2] == 1: self._solu300(self.G3[0],self.G3[1],self.W1,self.W2,2)
                        if self.ISRC[3] == 1: self._solu200(self.G4[0],self.G4[1],self.W1,self.W2,3)
                        if self.ISRC[4] == 1: self._solu300(self.G5[0],self.G5[1],self.W1,self.W2,4)
                        if self.ISRC[5] == 1: self._solu200(self.G6[0],self.G6[1],self.W1,self.W2,5)
                        if self.ISRC[6] == 1: self._solu300(self.G7[0],self.G7[1],self.W1,self.W2,6)
                        if self.ISRC[7] == 1: self._solu200(self.G8[0],self.G8[1],self.W1,self.W2,7)
                        if self.ISRC[8] == 1: self._solu200(self.G9[0],self.G9[1],self.W1,self.W2,8)
                        if self.ISRC[9] == 1: self._solu200(self.G10[0],self.G10[1],self.W1,self.W2,9)

                        if self.verbose: self._log("Done with SOLU")

                        for J in range(10):
                            if self.verbose: self._log("\t WVINT 101: J => %s" % J)
                            if self.verbose: self._log("\t WVINT 101: ISRC[J] => %s" % self.ISRC[J])
                            if self.ISRC[J] != 1: 
                                self.AA[J] = complex()
                                self.BB[J] = complex()
                                self.CC[J] = complex()
                            if self.verbose: self._log("WVINT  101 :   AA[J]=> %s" % self.AA[J])
                            if self.verbose: self._log("WVINT  101 :   BB[J]=> %s" % self.BB[J])
                            if self.verbose: self._log("WVINT  101 :   CC[J]=> %s" % self.CC[J])
                        #exit()
                    #}}}

                    # INITIALIZE THE INTEGRALS +++  NECESSARY FOR IASYMP.EQ.TRUE
                    #{{{
                    self.SMM[1][0] = self.AA[0] * self.J0K0 + self.BB[0] * self.J0K1 + self.CC[0] * self.J0K2
                    self.SMM[1][1] = self.AA[1] * self.J1K1 + self.BB[1] * self.J1K2 + self.CC[1] * self.J1K3
                    self.SMM[1][2] = self.AA[2] * self.J1K0 + self.BB[2] * self.J1K1 + self.CC[2] * self.J1K2 
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][0])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][1])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][2])

                    if self.ISRC[8] == 1 and self.ISRC[3] == 1: 
                        self.SUMD = (self.AA[3] + self.AA[8]) * self.J1K0 + (self.BB[3] + self.BB[8]) * self.J1K1 + (self.CC[3] + self.CC[8]) * self.J1K2
                    else:
                        self.SUMD=complex(0.0,0.0)

                    if self.verbose: self._log("\t WVINT after 101 SUMD() => %s" % self.SUMD)
                    self.SUMD= -self.SUMD / self.DISTANCE
                    if self.verbose: self._log("\t WVINT after 101 SUMD() => %s" % self.SUMD)
                    self.SMM[1][3] = self.SUMD + self.AA[3] * self.J0K1 + self.BB[3] * self.J0K2 + self.CC[3] * self.J0K3
                    self.SMM[1][4] = self.SUMD + self.AA[8] * self.J0K1 + self.BB[8] * self.J0K2 + self.CC[8] * self.J0K3
                    self.SMM[1][5] = self.AA[4] * self.J2K0 + self.BB[4] * self.J2K1 + self.CC[4] * self.J2K2
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][3])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][4])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][5])

                    if self.ISRC[5] == 1 and self.ISRC[9] == 1:
                        self.SUMD = (self.AA[5] + self.AA[9]) * self.J2K0 + (self.BB[5] + self.BB[9]) * self.J2K1 + (self.CC[5] + self.CC[9]) * self.J2K2
                    else:
                        self.SUMD = complex(0.0,0.0)

                    if self.verbose: self._log("\t WVINT after 101 SUMD() => %s" % self.SUMD)
                    self.SUMD = -2.0 * self.SUMD / self.DISTANCE
                    if self.verbose: self._log("\t WVINT after 101 SUMD() => %s" % self.SUMD)
                    self.SMM[1][6] = self.SUMD + self.AA[5] * self.J1K1 + self.BB[5] * self.J1K2 + self.CC[5] * self.J1K3
                    self.SMM[1][7] = self.SUMD + self.AA[9] * self.J1K1 + self.BB[9] * self.J1K2 + self.CC[9] * self.J1K3
                    self.SMM[1][8] = self.AA[6] * self.J0K0 + self.BB[6] * self.J0K1 + self.CC[6] * self.J0K2
                    self.SMM[1][9] = self.AA[7] * self.J1K1 + self.BB[7] * self.J1K2 + self.CC[7] * self.J1K3
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][6])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][7])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][8])
                    if self.verbose: self._log("\t WVINT after 101 SMM() => %s" % self.SMM[1][9])

                    if self.verbose: self._log("")
                    #}}}

                else:
                    for II in range(10):
                        self.SMM[1][II] = complex(0.0,0.0)


            #}}} INIT

            #     GET THE BESSEL'S FUNCTIONS +++++++++++++++++++++++++++++++>
            if self.verbose: self._log("")
            if self.verbose: self._log("***************** BESSEL FUNCTIONS *****************")
            if self.verbose: self._log("")

            if self.verbose: self._log("BESSEL FUNCTIONS: N22: %s" % self.N22)
            if self.verbose: self._log("BESSEL FUNCTIONS: ICN: %s" % self.ICN)
            if self.verbose: self._log("BESSEL FUNCTIONS: NNN1: %s" % self.NNN1)
            if self.verbose: self._log("BESSEL FUNCTIONS: IBLOCK: %s" % self.IBLOCK)
            #{{{
            if self.N22 < self.ICN:
                for J in range(self.NNN1-1,self.IBLOCK):
                    #if self.verbose: self._log("BESSEL FUNCTIONS: 1st J: %s" % J)
                    Z = self.WVNO[J].real * self.DISTANCE
                    X = (Z/3.0) * (Z/3.0)
                    AJ0 = 1.0-X*(2.2499997-X*(1.2656208-X*(.3163866-X*( .0444479-X*(.0039444-X*(.0002100))))))
                    A1Z = 0.5-X*(.56249985-X*(.21093573-X*(.03954289-X*( .00443319-X*(.00031761-X*(.00001109))))))
                    AJ1 = Z * A1Z
                    self.J2[J] = 2.0 * AJ1 / Z - AJ0
                    self.J1[J] = AJ1
                    self.J0[J] = AJ0
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: NNN1: %s" % self.NNN1)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: IBLOCK: %s" % self.IBLOCK)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: J: %s" % J)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: WVNO[J]: %s" % self.WVNO[J])
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: Z: %s" % Z)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: X: %s" % X)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: AJ0: %s" % AJ0)
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: AJ1: %s" % AJ1)
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J: %s" % J)
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J])
                    #exit()


            if self.verbose: self._log("BESSEL FUNCTIONS: N11: %s" % self.N11)
            if self.N11 > self.ICN:
                for J in range(self.NNN1-1,self.IBLOCK):
                    Z = self.WVNO[J].real * self.DISTANCE
                    X = 3.0 / Z
                    FAC = 1.0 / sqrt(Z)
                    F0 = 0.79788456+X*(-0.00000077 + X*(-0.00552740 + X*(-0.00009512+X*(0.00137237+X*(-0.00072805+X*(0.00014476))))))
                    T0 = Z - 0.78539816+X*(-0.04166397+X*(-0.00003954+X*(0.00262573+X*(-0.00054125+X*(-0.00029333+X*(0.00013558))))))
                    F1 = 0.79788456+X*(0.00000156+X*(0.01659667+X*(0.00017105+X*(-0.00249511+X*(0.00113653+X*(-0.00020033))))))
                    T1 = Z-2.35619449+X*(0.12499612+X*(0.00005650+X*(-0.00637879+X*(0.00074348+X*(0.00079824+X*(-0.00029166))))))
                    AJ0 = FAC * F0 * cos(T0)
                    AJ1 = FAC * F1 * cos(T1)
                    self.J2[J]= 2.0 * AJ1 / Z - AJ0
                    self.J1[J]= AJ1
                    self.J0[J]= AJ0
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: J: %s" % J)
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J])

            if self.N22 > self.ICN and self.N11 < self.ICN: 
                self.IM = self.ICN - self.N11 + 1
                #if self.verbose: self._log("BESSEL FUNCTIONS: IM: %s" % self.IM)
                for J in range(self.NNN1-1,self.IM):
                    #if self.verbose: self._log("BESSEL FUNCTIONS: 3rd J: %s" % J)
                    Z = self.WVNO[J].real * self.DISTANCE
                    X = (Z/3.0)*(Z/3.0)
                    AJ0 = 1.0-X*(2.2499997-X*(1.2656208-X*(0.3163866-X*(0.0444479-X*(0.0039444-X*(0.0002100))))))
                    A1Z = 0.5-X*(0.56249985-X*(0.21093573-X*(0.03954289-X*(0.00443319-X*(0.00031761-X*(0.00001109))))))
                    AJ1 = Z * A1Z
                    self.J2[J] = 2.0 * AJ1 / Z - AJ0
                    self.J1[J] = AJ1
                    self.J0[J] = AJ0
                    #if self.verbose: self._log("\t BESSEL FUNCTIONS: J: %s" % J)
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J])
                    #exit()

                for J in range(self.IM+1,self.IBLOCK):
                    #if self.verbose: self._log("BESSEL FUNCTIONS: 4th J: %s" % J)
                    Z = self.WVNO[J].real * self.DISTANCE
                    X = 3.0 / Z
                    FAC = 1.0 / sqrt(Z)
                    F0 = 0.79788456+X*(-0.00000077 + X*(-0.00552740 + X*(-0.00009512+X*(0.00137237+X*(-0.00072805+X*(0.00014476))))))
                    T0 = Z - 0.78539816+X*(-0.04166397+X*(-0.00003954+X*(0.00262573+X*(-0.00054125+X*(-0.00029333+X*(0.00013558))))))
                    F1 = 0.79788456+X*(0.00000156+X*(0.01659667+X*(0.00017105+X*(-0.00249511+X*(0.00113653+X*(-0.00020033))))))
                    T1 = Z-2.35619449+X*(0.12499612+X*(0.00005650+X*(-0.00637879+X*(0.00074348+X*(0.00079824+X*(-0.00029166))))))
                    AJ0 = FAC * F0 * cos(T0)
                    AJ1 = FAC * F1 * cos(T1)
                    self.J2[J] = 2.0 * AJ1 / Z - AJ0
                    self.J1[J] = AJ1
                    self.J0[J] = AJ0
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J])
                    if self.verbose: self._log("\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J])
                    #exit()

            if self.verbose: self._log("BESSEL FUNCTIONS:   J2 => %s" % self.J2)
            if self.verbose: self._log("BESSEL FUNCTIONS:   J1 => %s" % self.J1)
            if self.verbose: self._log("BESSEL FUNCTIONS:   J0 => %s" % self.J0)
            #}}}

            #     GET THE PHASE VELOCITY FILTER  +++++++++++++++++++++++++>
            #{{{
            if self.verbose: self._log('')
            if self.verbose: self._log('############### PHASE VELOCITY FILTER ##############')
            if self.verbose: self._log('')

            if self.CMAX >= 0.0:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.FACT[J] = 0.0
                    self.WVN = self.WVNO[J].real

                    if self.WVN > self.WVC1 and self.WVN < self.WVC2: self.FACT[J] = 1.0
                    if self.WVN > self.WVCM and self.WVN < self.WVC1: 
                        self.FACT[J] = (1.0 - cos(self.PI * (self.WVN-self.WVCM) / (self.WVC1-self.WVCM))).real / 2.0
                    if self.WVN > self.WVC2 and self.WVN < self.WVCN:
                        self.FACT[J] = (1.0 - cos(self.PI * (self.WVN - self.WVCN) / (self.WVC2 - self.WVCN))).real / 2.0

                    if self.verbose: self._log(' PHASE VELOCITY FILTER:    FACT[J] => %s' % self.FACT[J])
            #}}}

            #     REMOVAL OF ASYMPTOTIC TREND FROM THE INTEGRAND
            #{{{
            if self.verbose: self._log('')
            if self.verbose: self._log('############### ASYMP TREND REMOVAL ##############')
            if self.verbose: self._log('')

            if self.IASYMP:
                for K in range(10):
                    if self.ISRC[K] == 1:
                        for J in range(self.NNN1-1,self.IBLOCK):
                            self.AK = self.WVNO[J].real
                            self.EX = exp(-self.AK * self.DEPTH)
                            self.ASM = self.EX * (self.AA[K] + self.AK * (self.BB[K] + self.AK * self.CC[K]))
                            if self.verbose: self._log("    RM ASYMP TREND  10:    AK => %s  " % self.AK)
                            if self.verbose: self._log("    RM ASYMP TREND  10:    H => %s  " % self.DEPTH)
                            if self.verbose: self._log("    RM ASYMP TREND  10:    EX => %s  " % self.EX)
                            if self.verbose: self._log("    RM ASYMP TREND  10:    AA[K] => %s  " % self.AA[K])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    BB[K] => %s  " % self.BB[K])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    CC[K] => %s  " % self.CC[K])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    ASM => %s  " % self.ASM)
                            if self.verbose: self._log("    RM ASYMP TREND  10:    S1[J] => %s  " % self.G1[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    K => %s  " % K)
                            if self.verbose: self._log("    RM ASYMP TREND  10:    J => %s  " % J)
                            if K == 0: self.S1[J] = self.G1[J] - self.ASM
                            if K == 1: self.S2[J] = self.G2[J] - self.ASM
                            if K == 2: self.S3[J] = self.G3[J] - self.ASM
                            if K == 3: self.S4[J] = self.G4[J] - self.ASM
                            if K == 4: self.S5[J] = self.G5[J] - self.ASM
                            if K == 5: self.S6[J] = self.G6[J] - self.ASM
                            if K == 6: self.S7[J] = self.G7[J] - self.ASM
                            if K == 7: self.S8[J] = self.G8[J] - self.ASM
                            if K == 8: self.S9[J] = self.G9[J] - self.ASM
                            if K == 9: self.S10[J] = self.G10[J] - self.ASM
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G1[J] =>  %s " % self.S1[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G2[J] =>  %s " % self.S2[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G3[J] =>  %s " % self.S3[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G4[J] =>  %s " % self.S4[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G5[J] =>  %s " % self.S5[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G6[J] =>  %s " % self.S6[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G7[J] =>  %s " % self.S7[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G8[J] =>  %s " % self.S8[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G9[J] =>  %s " % self.S9[J])
                            if self.verbose: self._log("    RM ASYMP TREND  10:    G10[J] => %s  " % self.S10[J])

            if self.NBLK == self.NBLOCK:
                if self.IASYMP:
                    for K in range(10):
                        if self.ISRC[K] == 1:
                            #self._log("    ############ WVNO => %s  " % self.WVNO)
                            #self._log("    ############ len(WVNO) => %s  " % len(self.WVNO))
                            #self._log("    ############ IBLOCK => %s  " % self.IBLOCK)
                            #self._log("    ############ NK => %s  " % self.NK)
                            #self._log("    ############ IB => %s  " % self.IB)
                            self.AK = self.WVNO[self.IBLOCK].real
                            self.EX = exp(-self.AK )* self.DEPTH
                            self.ASM = self.EX * ( self.AA[K] + self.AK * (self.BB[K] + self.AK * self.CC[K]) )
                            if K == 0: self.S1[self.IBLOCK] = self.S1[self.IBLOCK] - self.ASM
                            if K == 1: self.S2[self.IBLOCK] = self.S2[self.IBLOCK] - self.ASM
                            if K == 2: self.S3[self.IBLOCK] = self.S3[self.IBLOCK] - self.ASM
                            if K == 3: self.S4[self.IBLOCK] = self.S4[self.IBLOCK] - self.ASM
                            if K == 4: self.S5[self.IBLOCK] = self.S5[self.IBLOCK] - self.ASM
                            if K == 5: self.S6[self.IBLOCK] = self.S6[self.IBLOCK] - self.ASM
                            if K == 6: self.S7[self.IBLOCK] = self.S7[self.IBLOCK] - self.ASM
                            if K == 7: self.S8[self.IBLOCK] = self.S8[self.IBLOCK] - self.ASM
                            if K == 8: self.S9[self.IBLOCK] = self.S9[self.IBLOCK] - self.ASM
                            if K == 9: self.S10[self.IBLOCK] = self.S10[self.IBLOCK] - self.ASM
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G1[J] =>  %s " % self.S1[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G2[J] =>  %s " % self.S2[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G3[J] =>  %s " % self.S3[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G4[J] =>  %s " % self.S4[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G5[J] =>  %s " % self.S5[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G6[J] =>  %s " % self.S6[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G7[J] =>  %s " % self.S7[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G8[J] =>  %s " % self.S8[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G9[J] =>  %s " % self.S9[1])
                            if self.verbose: self._log("    RM ASYMP TREND  666:    G10[J] => %s  " % self.S10[1])

            if not self.IASYMP:
                for K in range(10):
                    if self.ISRC[K] == 1:
                        for J in range(self.NNN1,self.IBLOCK):
                            if K == 0: self.S1[J] = self.G1[J]
                            if K == 1: self.S2[J] = self.G2[J]
                            if K == 2: self.S3[J] = self.G3[J]
                            if K == 3: self.S4[J] = self.G4[J]
                            if K == 4: self.S5[J] = self.G5[J]
                            if K == 5: self.S6[J] = self.G6[J]
                            if K == 6: self.S7[J] = self.G7[J]
                            if K == 7: self.S8[J] = self.G8[J]
                            if K == 8: self.S9[J] = self.G9[J]
                            if K == 9: self.S10[J] = self.G10[J]
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G1[J] =>  %s " % self.S1[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G2[J] =>  %s " % self.S2[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G3[J] =>  %s " % self.S3[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G4[J] =>  %s " % self.S4[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G5[J] =>  %s " % self.S5[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G6[J] =>  %s " % self.S6[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G7[J] =>  %s " % self.S7[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G8[J] =>  %s " % self.S8[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G9[J] =>  %s " % self.S9[1])
                            if self.verbose: self._log("    RM ASYMP TREND  777:    G10[J] => %s  " % self.S10[1])


            if self.verbose: self._log("    RM ASYMP TREND  13:    ISRC => %s  " %  self.ISRC)

            if self.ISRC[0] == 1:
                #if self.verbose: self._log("    RM ASYMP TREND  13:    IBL => %s  " %  self.IBLOCK)
                #if self.verbose: self._log("    RM ASYMP TREND  13:    N1 => %s  " %  self.NNN1)
                for J in range(self.IBLOCK-1,self.NNN1-2,-1):
                    self.SUMD = self.FACT[J] * self.S1[J] * self.J0[J]
                    self.SMM[1][0] = self.SMM[1][0] + self.SUMD
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    J => %s  " %  J)
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    FACT[J] => %s  " % self.FACT[J])
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    G1[J] => %s  " % self.S1[J])
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    J0[J] => %s  " % self.J0[J])
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    SUMD => %s  " % self.SUMD)
                    if self.verbose: self._log("    RM ASYMP TREND  13:    SMM[1][0] => %s  " % self.SMM[1][0])


            if self.ISRC[1] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.AK = self.WVNO[J].real
                    self.SUMD = self.FACT[J] * self.S2[J] * self.AK * self.J1[J]
                    self.SMM[1][1] = self.SMM[1][1] + self.SUMD
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    J => %s  " %  J)
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    AK => %s  " % self.AK)
                    #if self.verbose: self._log("    RM ASYMP TREND  13:    SUMD => %s  " % self.SUMD)
                    if self.verbose: self._log("    RM ASYMP TREND  14:    SMM[1][1] => %s  " % self.SMM[1][1])

            if self.ISRC[2] == 1: 
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.SUMD = self.FACT[J] * self.S3[J] * self.J1[J]
                    self.SMM[1][2] = self.SMM[1][2] + self.SUMD
                    if self.verbose: self._log("    RM ASYMP TREND  15:    SMM[1][2] => %s  " % self.SMM[1][2])
            #}}}

            #     INCLUDE NEAR FIELD TERM IF BOTH SH AND P-SV ARE COMPUTED
            #{{{
            if self.ISRC[3] == 1 and self.ISRC[8] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.SUMC[J] = (self.S4[J] + self.S9[J]) * self.J1[J]
                #    if self.verbose: self._log("    RM ASYMP TREND  16:    J => %s  " % J)
                #    if self.verbose: self._log("    RM ASYMP TREND  16:    G4[J] => %s  " % self.S4[J])
                #    if self.verbose: self._log("    RM ASYMP TREND  16:    G9[J] => %s  " % self.S9[J])
                #    if self.verbose: self._log("    RM ASYMP TREND  16:    J1[J] => %s  " % self.J1[J])
                    if self.verbose: self._log("    RM ASYMP TREND  16:    SUMC[J] => %s  " % self.SUMC[J])
                #self._log("")
                #self._log("\t\t\t G4 => %s " % self.S4)
                #self._log("\t\t\t G9 => %s " % self.S9)
                #self._log("")

            if self.ISRC[3] == 1:
                for J in range(self.NNN1-1, self.IBLOCK):
                    self.AK= self.WVNO[J].real
                    self.SUMD = self.S4[J] * self.AK * self.J0[J]
                    self.SMM[1][3] = self.SMM[1][3] + (self.SUMD - self.SUMC[J] / self.DISTANCE) * self.FACT[J]
                    if self.verbose: self._log("    RM ASYMP TREND  17:    SMM[1][3] => %s  " % self.SMM[1][3])

            if self.ISRC[8] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.AK = self.WVNO[J].real
                    self.SUMD = self.S9[J] * self.AK * self.J0[J]
                    self.SMM[1][4] = self.SMM[1][4] + (self.SUMD - self.SUMC[J] / self.DISTANCE) * self.FACT[J]
                    if self.verbose: self._log("   RM ASYMP TREND  18 :    SMM[1][4] => %s  " % self.SMM[1][4])

            if self.ISRC[4] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.SUMD = self.FACT[J] * self.S5[J] * self.J2[J]
                    self.SMM[1][5] = self.SMM[1][5] + self.SUMD
                    if self.verbose: self._log("    RM ASYMP TREND  19:    SMM[1][5] => %s  " % self.SMM[1][5])
            #}}}

            #     NEAR FIELD TERM IF BOTH SH AND P-SV TERMS ARE COMPUTED
            #{{{
            if self.ISRC[5] == 1 and self.ISRC[9] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.SUME[J] = (self.S6[J] + self.S10[J]) * self.J2[J]
                    if self.verbose: self._log("    RM ASYMP TREND  20:    SUME[J] => %s  " % self.SUME[J])

            if self.ISRC[5] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.AK = self.WVNO[J].real
                    self.SUMD = self.S6[J] * self.AK * self.J1[J]
                    self.SMM[1][6] = self.SMM[1][6] + (self.SUMD - 2.0 * self.SUME[J] /self.DISTANCE) * self.FACT[J]
                    if self.verbose: self._log("     RM ASYMP TREND  21:    SMM[1][6] => %s  " % self.SMM[1][6])

            if self.ISRC[9] == 1:
                for J in range(self.IBLOCK-1,self.NNN1-2,-1):
                    self.AK = self.WVNO[J].real
                    self.SUMD = self.S10[J] * self.AK * self.J1[J]
                    self.SMM[1][7] = self.SMM[1][7] + (self.SUMD - 2.0 * self.SUME[J] / self.DISTANCE) * self.FACT[J]
                    if self.verbose: self._log("     RM ASYMP TREND  22:    SMM[1][7] => %s  " % self.SMM[1][7])

            if self.ISRC[6] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.SUMD = self.FACT[J] * self.S7[J] * self.J0[J]
                    self.SMM[1][8] = self.SMM[1][8] + self.SUMD
                    if self.verbose: self._log("     RM ASYMP TREND  23:    SMM[1][8] => %s  " % self.SMM[1][8])

            if self.ISRC[7] == 1:
                for J in range(self.NNN1-1,self.IBLOCK):
                    self.AK = self.WVNO[J].real
                    self.SUMD = self.FACT[J] * self.S8[J] * self.AK * self.J1[J]
                    self.SMM[1][9] = self.SMM[1][9] + self.SUMD
                    if self.verbose: self._log("     RM ASYMP TREND  24:    SMM[1][9] => %s  " % self.SMM[1][9])

            if self.verbose: self._log("")
            if self.verbose: self._log(" END")
            if self.verbose: self._log("")

            #}}}

            #}}} sub WVINT

            #}}} FR

            self.NBLK += 1

        #}}} A

        #self._log(' verbose: T0 = %s' % self.T0)
        #{{{ B
        self.T0X = self.T0[1] + self.DISTANCE / 8
        self.FAT = self.OMEGA * self.T0X
        self.XR = cos(self.FAT)
        self.XI = sin(self.FAT)
        self.DKK = self.DK * self.DISTANCE
        self.SMM[1][1] *= -1
        self.SMM[1][9] *= -1
        for J in range(10):
            self.SMM[1][J] = self.SMM[1][J] * self.DKK
            self.SMM[1][J] = self.SMM[1][J] / self.DISTANCE

        for J in range(10):
            if self.ISRC[J] != 1: break
            self.GG[1][J] = self.SMM[1][J] * complex(self.XR,self.XI)
            if self.verbose: self._log("%s" % self.GG[1][J])
        #}}} B

        if self.verbose: self._log("%s" % self.GG)
        #return self.GG 

#}}}

    def _solu300(self,Y1,Y2,X1,X2,J):
        #{{{
        if self.verbose: self._log("")
        if self.verbose: self._log("&&&&&&&&&&&&&&&&&&& SOLU &&&&&&&&&&&&&&&")
        self.CC[J] = complex(0.0,0.0)
        U1 = X1 * self.DEPTH
        U2 = X2 * self.DEPTH
        DET = X2 - X1
        if self.verbose: self._log("SOLU 300: U1 =%s" % U1)
        if self.verbose: self._log("SOLU 300: U2 =%s" % U2)
        if self.verbose: self._log("SOLU 300: DET =%s" % DET)
        self.AA[J] = complex(0.0,0.0)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: X1 => %s" % X1)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: X2 => %s" % X2)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: Y1 => %s" % Y1)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: Y2 => %s" % Y2)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: U1 => %s" % U1)
        if self.verbose: self._log("^^^^^^^^ SOLU 300: U2 => %s" % U2)

        self.BB[J] = X2 * Y1 * exp(U1) / X1
        if self.verbose: self._log("^^^^^^^^ SOLU 300: BB[%s]:%s" % (J,self.BB[J]))
        self.BB[J] = X1 * Y2 * exp(U2) / X2
        if self.verbose: self._log("^^^^^^^^ SOLU 300: BB[%s]:%s" % (J,self.BB[J]))

        self.BB[J] = X2 * Y1 * exp(U1) / X1 - X1 * Y2 * exp(U2) / X2
        self.BB[J] = self.BB[J] / DET
        self.CC[J] = Y2 * exp(U2) / X2 - Y1 * exp(U1) / X1
        self.CC[J] = self.CC[J] / DET
        if self.verbose: self._log("SOLU 300: AA[%s]:%s" % (J,self.AA[J]))
        if self.verbose: self._log("SOLU 300: BB[%s]:%s" % (J,self.BB[J]))
        if self.verbose: self._log("SOLU 300: CC[%s]:%s" % (J,self.CC[J]))
        if self.verbose: self._log("")
        #exit()
        #}}}

    def _solu200(self,Y1,Y2,X1,X2,J):
        #{{{
        if self.verbose: self._log("")
        if self.verbose: self._log("&&&&&&&&&&&&&&&&&&& SOLU &&&&&&&&&&&&&&&")
        self.CC[J] = complex(0.0,0.0)
        U1 = X1 * self.DEPTH
        U2 = X2 * self.DEPTH
        DET = X2 - X1
        if self.verbose: self._log("SOLU 200: U1 =%s" % U1)
        if self.verbose: self._log("SOLU 200: U2 =%s" % U2)
        if self.verbose: self._log("SOLU 200: DET =%s" % DET)
        self.AA[J] = X2 * Y1 * exp(U1) - X1 * Y2 * exp(U2)
        self.AA[J] = self.AA[J] / DET
        self.BB[J] = Y2 * exp(U2) - Y1 * exp(U1)
        self.BB[J] = self.BB[J] / DET
        if self.verbose: self._log("SOLU 200: AA[%s]:%s" % (J,self.AA[J]))
        if self.verbose: self._log("SOLU 200: BB[%s]:%s" % (J,self.BB[J]))
        if self.verbose: self._log("")
        #}}}

    def _reformat(self):
#{{{
        for l in range(10):

            self._log("NOW: %s" % l)
            TEMPDATA = self.DATA[l]
            if l == 0:
                #self.ZDD = [x * -1 for x in self.TEMPDATA ]
                if not 'ZDD' in self.ELEMENTS: self.ELEMENTS['ZDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['ZDD']['data'] = TEMPDATA.values()
                #self.ELEMENTS['ZDD']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['ZDD']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZDD']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZDD']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZDD']['data'] ] # convert from cm/s to nm/s
                self._log("\tZDD: %s" % len(self.ELEMENTS['ZDD']['data']))
                #self.ELEMENTS['ZDD']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['ZDD']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 1:
                #self.XDD = [x for x in self.TEMPDATA]
                if not 'XDD' in self.ELEMENTS: self.ELEMENTS['XDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['XDD']['data'] = TEMPDATA.values()
                #self.ELEMENTS['XDD']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['XDD']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XDD']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XDD']['data'] = [ x * 10000000 for x in self.ELEMENTS['XDD']['data'] ] # convert from cm/s to nm/s
                self._log("\tXDD: %s" % len(self.ELEMENTS['XDD']['data']))
                #self.ELEMENTS['XDD']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['XDD']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 2:
                #self.ZDS = [x for x in self.TEMPDATA]
                if not 'ZDS' in self.ELEMENTS: self.ELEMENTS['ZDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['ZDS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['ZDS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['ZDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZDS']['data'] ] # convert from cm/s to nm/s
                self._log("\tZDS: %s" % len(self.ELEMENTS['ZDS']['data']))
                #self.ELEMENTS['ZDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['ZDS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 3:
                #self.XDS = [x * -1 for x in self.TEMPDATA ]
                if not 'XDS' in self.ELEMENTS: self.ELEMENTS['XDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['XDS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['XDS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['XDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['XDS']['data'] ] # convert from cm/s to nm/s
                self._log("\tXDS: %s" % len(self.ELEMENTS['XDS']['data']))
                #self.ELEMENTS['XDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['XDS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 4:
                #self.TDS = [x for x in self.TEMPDATA]
                if not 'TDS' in self.ELEMENTS: self.ELEMENTS['TDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['TDS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['TDS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['TDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['TDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['TDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['TDS']['data'] ] # convert from cm/s to nm/s
                self._log("\tTDS: %s" % len(self.ELEMENTS['TDS']['data']))
                #self.ELEMENTS['TDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['TDS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 5:
                #self.ZSS = [x * -1 for x in self.TEMPDATA ]
                if not 'ZSS' in self.ELEMENTS: self.ELEMENTS['ZSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['ZSS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['ZSS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['ZSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZSS']['data'] ] # convert from cm/s to nm/s
                self._log("\tZSS: %s" % len(self.ELEMENTS['ZSS']['data']))
                #self.ELEMENTS['ZSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['ZSS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 6:
                #self.XSS = [x for x in self.TEMPDATA]
                if not 'XSS' in self.ELEMENTS: self.ELEMENTS['XSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['XSS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['XSS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['XSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['XSS']['data'] ] # convert from cm/s to nm/s
                self._log("\tXSS: %s" % len(self.ELEMENTS['XSS']['data']))
                #self.ELEMENTS['XSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['XSS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 7:
                #self.TSS = [x * -1 for x in self.TEMPDATA ]
                if not 'TSS' in self.ELEMENTS: self.ELEMENTS['TSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['TSS']['data'] = TEMPDATA.values()
                #self.ELEMENTS['TSS']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['TSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['TSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['TSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['TSS']['data'] ] # convert from cm/s to nm/s
                self._log("\tTSS: %s" % len(self.ELEMENTS['TSS']['data']))
                #self.ELEMENTS['TSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['TSS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 8:
                #self.REX = [x for x in self.TEMPDATA]
                if not 'REX' in self.ELEMENTS: self.ELEMENTS['REX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['REX']['data'] = TEMPDATA.values()
                #self.ELEMENTS['REX']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['REX']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['REX']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['REX']['data'] = [ x * 10000000 for x in self.ELEMENTS['REX']['data'] ] # convert from cm/s to nm/s
                self._log("\tREX: %s" % len(self.ELEMENTS['REX']['data']))
                #self.ELEMENTS['REX']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['REX']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 9:
                #self.ZEX = [x for x in self.TEMPDATA]
                if not 'ZEX' in self.ELEMENTS: self.ELEMENTS['ZEX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                #self.ELEMENTS['ZEX']['data'] = TEMPDATA.values()
                #self.ELEMENTS['ZEX']['data'] = [TEMPDATA[0] for x in TEMPDATA ]
                self.ELEMENTS['ZEX']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZEX']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZEX']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZEX']['data'] ] # convert from cm/s to nm/s
                self._log("\tZEX: %s" % len(self.ELEMENTS['ZEX']['data']))
                #self.ELEMENTS['ZEX']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                #self.ELEMENTS['ZEX']['data'].extend( [x for x in self.TEMPDATA ])
            else:
                raise SystemExit('\n\nERROR: GreenFunctions(): generate(): Too many elements\n\n')


        self._save_to_db()
#}}}

    def _wvint9(self):
#{{{
        if self.verbose: self._log('N2:  %s' % self.N2)

        for i in range(self.N2):
            if self.verbose: self._log('OMEGA(%s):   %s' % ( i, self.ALLOMEGA[i]))
            for k in range(10):
                try:
                    if self.verbose: self._log('DATA(%s,%s,%s):   %s' % ( 1, i, k, self.DATA[i][k]))
                except Exception,e:
                    self.DATA[i][k] = complex(0.0,0.0)
                    if self.verbose: self._log('DATA(%s,%s,%s):   %s' % ( 1, i, k, Exception))

        self.REP = 'n'


        self.N = 2 * (self.NYQ-1)

        self.TAU = self.DT

        self.FMAX = self.FU
        self.inst = 0

        """
        Beginning of distance loop
        """

        self.PI = cos(-1.0)
        self.TWPI = 2.0 * self.PI

        #self.OMEGA



        T0X = self.DISTANCE/8.0
        yr = 0.0

        if self.verbose: self._log("range=%s yr=%s depth=%s" % (self.DISTANCE,yr,self.DEPTH))
        if self.verbose: self._log("npoint=%s t0x=%s dt=%s tau=%s" % (self.N,T0X,self.DT,self.TAU))

        if self.verbose: self._log("range=%s yr=%s depth=%s npoint=%s t0x=%s dt=%s tau=%s" % (self.DISTANCE,yr,self.DEPTH,self.N,T0X,self.DT,self.TAU))
        for l in range(10):

            self.TEMPDATA = defaultdict(complex)
            temp = defaultdict(complex)

            if self.ISRC[l] != 1:
                for j in range(self.N):
                    self.TEMPDATA[j] = complex(0.0,0.0)
                    #self._log("TEST 1: TEMPDATA[j] => %s" % self.TEMPDATA[j])
            else:

                for j in range(self.NYQ):

                    if j < self.N2:
                        #self._log("self.DATA= %s" % self.DATA)
                        #self._log("j = %s l = %s" % (j,l))
                        self.TEMPDATA[j]=self.DATA[j][l]
                        #self._log("TEMPDATA[%s] = %s" % (j,self.TEMPDATA[j]))
                        #self.TEMPDATA[self.N2-1-j]=self.DATA[nd][j][l].conjugate()

                    else:
                        self.TEMPDATA[j] = complex(0.0,0.0)
                        #self.TEMPDATA[self.N2-1-j] = complex(0.0,0.0).conjugate()

                for j in range(self.NYQ):
                    #self._log("TEMPTDATA:[%s]" % self.TEMPDATA)
                    #self._log("TEMPTDATA[%s]:[%s]" % (j,self.TEMPDATA[j]))
                    temp[self.NYQ-j] =  self.TEMPDATA[j].conjugate() 

                for j in range(self.NYQ):
                    self.TEMPDATA[self.NYQ+1+j] =  temp[j]



                """
                Correct for damping.
                """
                fac = exp(self.ALPHA  * T0X)
                dfac = exp(self.ALPHA*self.DT)

                for i in range(self.N2):
                    self.TEMPDATA[i] = self.TEMPDATA[i] * fac
                    fac = fac * dfac

                #convert to list
                self.TEMPDATA = [ x for x in self.TEMPDATA.values()]

                self.TEMPDATA = ifft(self.TEMPDATA)

                #Convert to real values
                self.TEMPDATA = self.TEMPDATA.real

            if l == 0:
                #self.ZDD = [x * -1 for x in self.TEMPDATA ]
                if not 'ZDD' in self.ELEMENTS: self.ELEMENTS['ZDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZDD']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['ZDD']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 1:
                #self.XDD = [x for x in self.TEMPDATA]
                if not 'XDD' in self.ELEMENTS: self.ELEMENTS['XDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XDD']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['XDD']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 2:
                #self.ZDS = [x for x in self.TEMPDATA]
                if not 'ZDS' in self.ELEMENTS: self.ELEMENTS['ZDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['ZDS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 3:
                #self.XDS = [x * -1 for x in self.TEMPDATA ]
                if not 'XDS' in self.ELEMENTS: self.ELEMENTS['XDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['XDS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 4:
                #self.TDS = [x for x in self.TEMPDATA]
                if not 'TDS' in self.ELEMENTS: self.ELEMENTS['TDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['TDS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['TDS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 5:
                #self.ZSS = [x * -1 for x in self.TEMPDATA ]
                if not 'ZSS' in self.ELEMENTS: self.ELEMENTS['ZSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['ZSS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 6:
                #self.XSS = [x for x in self.TEMPDATA]
                if not 'XSS' in self.ELEMENTS: self.ELEMENTS['XSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['XSS']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 7:
                #self.TSS = [x * -1 for x in self.TEMPDATA ]
                if not 'TSS' in self.ELEMENTS: self.ELEMENTS['TSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['TSS']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['TSS']['data'].extend( [x * -1 for x in self.TEMPDATA ])

            elif l == 8:
                #self.REX = [x for x in self.TEMPDATA]
                if not 'REX' in self.ELEMENTS: self.ELEMENTS['REX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['REX']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['REX']['data'].extend( [x for x in self.TEMPDATA ])

            elif l == 9:
                #self.ZEX = [x for x in self.TEMPDATA]
                if not 'ZEX' in self.ELEMENTS: self.ELEMENTS['ZEX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZEX']['data'] = [self.TEMPDATA[0] for x in self.TEMPDATA ]
                self.ELEMENTS['ZEX']['data'].extend( [x for x in self.TEMPDATA ])
            else:
                raise SystemExit('\n\nERROR: GreenFunctions(): generate(): Too many elements\n\n')


        self._save_to_db()
#}}}

    def _integrate(self):
#{{{
        """
        Time Domain Integration
        """
        for trace in self.ELEMENTS:
            if not self.ELEMENTS[trace]['data']:
                raise SystemExit('\n\nERROR: GreenFunctions(): Problems while integration(): no data for %s\n\n' % trace)

            if self.verbose: self._log('GreenFunctions(): integrat(%s)' % trace)
            self.ELEMENTS[trace]['data'] = cumtrapz(self.ELEMENTS[trace]['data'])

#}}}

    def plot(self,start=0,end=-1):
#{{{
        """ 
        Plot all traces in memory. They are containe in 
        the dictionary ELEMENTS in self.
        """
        total = len(self.ELEMENTS)

        half = int(total/2)

        if half != total/2: half += 1

        now = 0
        for trace in self.ELEMENTS:

            now += 1
            plt.subplot(5,2,now)
            data = self.ELEMENTS[trace]['data']
            plt.plot(data[start:end])
            plt.legend([trace])

        plt.suptitle("Green Functions: depth:%s distance:%s" % (self.DEPTH,self.DISTANCE))
        plt.show()
#}}}

    def __getitem__(self,i):
#{{{
        i.upper()
        if i == 'ALL':
            temp = {}
            for x in self.ELEMENTS:
                temp[x] = self.ELEMENTS[x]['data']

            return temp

        else:

            try:
                return self.ELEMENTS[i]['data']
            except:
                raise SystemExit('Wrong name of element. (%s)\n'% i)
#}}}

    def __call__(self,i):
#{{{
        i.upper()
        if i == 'ALL':
            temp = {}
            for x in self.ELEMENTS:
                temp[x] = self.ELEMENTS[x]['data']

            return temp

        else:

            try:
                return self.ELEMENTS[i]['data']
            except:
                raise SystemExit('Wrong name of element. (%s)\n'% i)
#}}}

    def __getattr__(self,i):
#{{{
        i.upper()
        if i == 'ALL':
            temp = {}
            for x in self.ELEMENTS:
                temp[x] = self.ELEMENTS[x]['data']

            return temp

        else:

            try:
                return self.ELEMENTS[i]['data']
            except:
                raise SystemExit('Wrong name of element. (%s)\n'% i)
#}}}

#}}}

"""
If we call this script directly, then output help.
"""
if __name__ == "__main__":
#{{{
    self._log("Green's Functions Library:")
    self._log("\t1) Import the library.  ")
    self._log("\t\timport moment_tensor.fkrprog as fkr")
    self._log("\t2) Build object for GF's.  ")
    self._log("\t\tgf = GreenFunctions(modle.pf)")
    self._log("\t3) Build for [distance,depth]")
    self._log("\t\tgf.build(depth=8,distance=1)")
    self._log("\t3) Get GF's matrix with elements")
    self._log("\t\tdata=gf.ALL")

    self._log("\n\n")
    self._log("No BRTT support.")
    self._log("Juan Reyes <reyes@ucsd.edu>")
    self._log("Rob Newman <rlnewman@ucsd.edu>")

#}}}

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
    def __init__(self,model):
#{{{
        log( "---------------------------------------" )
        log( "FKRPROG.PY: __init__(%s)" % model )

        self.pf = model

        self.IVEL = 'v'
        self.IB  = 100
        self.IFREQ = 0
        self.DATA = defaultdict(lambda: defaultdict(dict))
        self.ALLOMEGA = defaultdict(lambda: defaultdict(dict))
        self.ELEMENTS = defaultdict(dict)

        """ Read configuration from parameter file self.pf """
        self._read_pf()

        """ Return a normalized absolutized version of the path. """
        self.gf_db_path = os.path.abspath(self.gf_db_path)

        """ Verify database """
        db = self._open_db()

        db.close()

        log("FKRPROG.PY: Using database [%s]" % os.path.normpath(self.gf_db_path + '/' + self.gf_db + '.wfdisc'))

#}}}

    def _open_db(self):
        """  Open the database that will save our functions.  """
#{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: open_db()")

        """ Recursive directory creation function. """
        try:
            if not os.path.exists(self.gf_db_path): os.makedirs(self.gf_db_path)
        except Exception,e:
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot create direcotry (%s) %s => %s \n\n' % (self.gf_db_path,Exception,e))

        """ Verify the directory. """
        if not os.path.exists(self.gf_db_path):
            raise SystemExit('\n\nERROR: GreenFunctions(): Cannot find direcotry (%s)\n\n' % self.gf_db_path)

        log("FKRPROG.PY: open_db()  -   Open database: %s/%s\n" % (self.gf_db_path,self.gf_db))

        """ Open database and keep pointer in db object"""
        try: 
            db = datascope.dbopen( os.path.normpath(self.gf_db_path + '/' + self.gf_db), "r+" )
            db = db.lookup( table='wfdisc' )

        except Exception, e:
            raise SystemExit('\n\nERROR: dbopen() %s => %s\n\n' % (Exception,e))

        return db
#}}}

    def build(self,depth,distance,sps=1,type='v',filter=None):
        """ Main function to retrieve ( and produce if missing )
        the requested elements.

        """
#{{{

        """
        sps: The requested samplerate for the data
        type: Output velocity of displacepment
        filter: Apply this filter to the data 

        """

        log( "---------------------------------------" )
        log("FKRPROG.PY: build(%s,%s,sps=%s,type=%s,filter=%s)" % (depth,distance,sps,type,filter))

        self.DEPTH = depth
        self.DISTANCE = distance
        self.IVEL = type

        log("PKRPROG.PY: build - Get from database.")
        if not self._get_from_db(depth,distance,sps,type,filter): 
            log("PKRPROG.PY: build - Not in db. Generate new.")
            self._generate(depth,distance)
            log("PKRPROG.PY: build - Get from db newly created funcs.")
            self._get_from_db(depth,distance,sps,type,filter)

#}}}

    def _generate(self,depth,distance,type='v'):
#{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: generate(%s,%s,type=%s)" % (depth,distance,type))

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
        log("FKRPROG.PY: generate()  -  Build new model.")

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

        log("FKRPROG.PY: generate()  -  Setup vars.")
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


        self.MMAX = len(self.D) - 1

        #log("FKRPROG.PY: generate()  -  LMAX = %s" % self.LMAX)
        #log("FKRPROG.PY: generate()  -  MMAX = %s" % self.MMAX)
        #log("FKRPROG.PY: generate()  -  DEPTH = %s" % self.DEPTH)
        #log("FKRPROG.PY: generate()  -  RANGE = %s" % self.DISTANCE)


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

        log("FKRPROG.PY: generate()  -  MODEL: \n%s" % model,1)

        log("FKRPROG.PY: generate()  -  Running: fortran_fkrprog < ./greens_func")
        p = os.popen('fortran_fkrprog < ./.greens_funcs',"r")
        while 1:
            line = p.readline()
            log(line)
            if not line: break

        p.close()

        log("FKRPROG.PY: generate()  -  Running: wvint9")
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
                debug("DATA[%s][%s] = %s" % (l,n,m.group(3)))
            else: 
                debug("ERROR ** no match ** : %s" % line)

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
        """ Open the database and extract the traces for 
        requested depth and distance. 

        """
#{{{

        """ 

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

        log( "---------------------------------------" )
        log("FKRPROG.PY: get_from_db()")


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
            db = db.subset('sta =~ /%s/ && chan =~ /%s_.*/' % (distance,depth))
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

                debug('FKRPROG.PY: get_from_db()  -   getv()=> (%s,%s,%s,%s,%s,%s)'%(db_sta,db_chan,db_nsamp,db_time,db_endtime,db_samprate))

                """ Get steps needed to get to new samplerate. """
                decimate_files = self._decimate_file(db_samprate,sps)
                log('FKRPROG.PY: get_from_db()  -  decimate factor [%s] file [%s]' % ( db_samprate/sps,decimate_files))

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
                    log('FKRPROG.PY: get_from_db()  -  trloadchan(%s,%s,%s,%s)'% (db_time,db_endtime,db_sta,db_chan))
                    tr = datascope.trloadchan(db,db_time,db_endtime,db_sta,db_chan)
                    log('FKRPROG.PY: get_from_db()  -  tr.slice()')
                    tr.splice()
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): trace object %s: %s\n\n' % (Exception,e))

                try:
                    # integrate data
                    if type == 'D' or type == 'd':
                        log('FKRPROG.PY: get_from_db()  -  integrate for velocity ')
                        tr.filter('INT')
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): integrate %s: %s\n\n' % (Exception,e))

                try:
                    # filter data
                    if filter:
                        log('FKRPROG.PY: get_from_db()  -  filter [%s]' % filter)
                        tr.filter(filter)
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): filter %s: %s\n\n' % (Exception,e))


                try:
                    # decimate data
                    for f in decimate_files:
                        full_path = os.environ['ANTELOPE'] + '/data/responses/' + f
                        log('FKRPROG.PY: get_from_db()  -  filter(DECIMATE %s)' % full_path)
                        tr.filter('DECIMATE %s' % full_path)
                except Exception,e:
                    raise SystemExit('\n\nERROR: GreenFunctions(): decimate %s: %s\n\n' % (Exception,e))


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
            log('FKRPROG.PY: ERROR: No records for subset [sta =~ /%s/ && chan =~ /%s/].' % (distance,depth))
            return False

        try:
            db.close()
        except:
            pass

        return True
#}}}

    def _save_to_db(self):
        """ Open the database and save the new traces.

        """
#{{{
        """ 
        The functions are archived in ascii files referenced by Datascope using a simple 
        wfdisc table. The value for the station is our DISTANCE. The value for the channel  
        is our DEPTH and the element is specified in the location code. 
        i.e.   
            depth: 8
            distance: 10
            element: TDS
            => 10_8_TDS ( format: sta_chan_loc )

        """

        log( "---------------------------------------" )
        log('FKRPROG.PY: save_to_db()')

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

            log('FKRPROG.PY: save_to_db()  -  %s' % element)
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

#}}}

    def _decimate_file(self,have=1,want=1):
        """ Function to return needed response file 
        for proper decimation. 

        """
#{{{

        """
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

        log( "---------------------------------------" )
        log("FKRPROG.PY: decimate_file()")

        factor = int(have/want)

        log("FKRPROG.PY: decimate_file()  -  factor [%s]" % factor)

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
        """ Template for the earth model.

        """
    #{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: file_format()")

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
        """ Read parameters from 
        configuration file.

        """
    #{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: read_pf(%s)" % self.pf)

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

        log("FKRPROG.PY: read_pf()  -  DECAY=%s N1=%s N2=%s N=%s DT=%s " % (self.DECAY,self.N1,self.N2,self.N,self.DT))


        # ISRC, JBDRY 
        self.ISRC = [ 1 for x in range(10) ]
        self.JBDRY = 0
        log("FKRPROG.PY: read_pf()  -  ISRC=%s JBDRY=%s " % (self.ISRC,self.JBDRY))


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
                self.D[x]   = float(t[0])
                self.A[x]   = float(t[1])
                self.B[x]   = float(t[2])
                self.RHO[x] = float(t[3])
                self.QA[x]  = 1/float(t[4])
                self.QB[x]  = 1/float(t[5])

        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n RAW: %s'% (self.pf,Exception,e,temp))


        if debug:
            for I in range(len(self.D)):
                log("%s: %s\t%s\t%s\t%s\t%s\t%s" % (self.pf,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I]))


        #  NFC...
        self.XLENG  = 400 
        self.XFAC  = 1.5


        self.CMAX = stock.pfget_double(self.pf,'cmax')
        self.C1   = stock.pfget_double(self.pf,'c1')
        self.C2   = stock.pfget_double(self.pf,'c2')
        self.CMIN = stock.pfget_double(self.pf,'cmin')
        log("FKRPROG.PY: read_pf()  -  CMAX=%s C1=%s C2=%s CMIN=%s " % (self.ISRC,self.C1,self.C2,self.CMIN))

        # need to fix this while running
        self.T0    = defaultdict(int) 


        ##
        ## FIX VALUES TO MATCH PYTHON INDEXING
        ##
        self.LMAX = len(self.D)-1 # layer number below the source

        self.MMAX = len(self.D)-1 # Total number of layers

        return

    #}}}

    def _reformat(self):
        """ Take the GF elements and fix the 
        format to forward them to the dbmoment.py 
        tool.

        """
#{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: reformat()")

        for l in range(10):

            log("FKRPROG.PY: reformat()  -   %s" % l)
            TEMPDATA = self.DATA[l]
            if l == 0:
                if not 'ZDD' in self.ELEMENTS: self.ELEMENTS['ZDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZDD']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZDD']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZDD']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZDD']['data'] ] # convert from cm/s to nm/s
                debug("\tZDD: %s" % len(self.ELEMENTS['ZDD']['data']))

            elif l == 1:
                if not 'XDD' in self.ELEMENTS: self.ELEMENTS['XDD'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XDD']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XDD']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XDD']['data'] = [ x * 10000000 for x in self.ELEMENTS['XDD']['data'] ] # convert from cm/s to nm/s
                debug("\tXDD: %s" % len(self.ELEMENTS['XDD']['data']))

            elif l == 2:
                if not 'ZDS' in self.ELEMENTS: self.ELEMENTS['ZDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZDS']['data'] ] # convert from cm/s to nm/s
                debug("\tZDS: %s" % len(self.ELEMENTS['ZDS']['data']))

            elif l == 3:
                if not 'XDS' in self.ELEMENTS: self.ELEMENTS['XDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['XDS']['data'] ] # convert from cm/s to nm/s
                debug("\tXDS: %s" % len(self.ELEMENTS['XDS']['data']))

            elif l == 4:
                if not 'TDS' in self.ELEMENTS: self.ELEMENTS['TDS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['TDS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['TDS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['TDS']['data'] = [ x * 10000000 for x in self.ELEMENTS['TDS']['data'] ] # convert from cm/s to nm/s
                debug("\tTDS: %s" % len(self.ELEMENTS['TDS']['data']))

            elif l == 5:
                if not 'ZSS' in self.ELEMENTS: self.ELEMENTS['ZSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZSS']['data'] ] # convert from cm/s to nm/s
                debug("\tZSS: %s" % len(self.ELEMENTS['ZSS']['data']))

            elif l == 6:
                if not 'XSS' in self.ELEMENTS: self.ELEMENTS['XSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['XSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['XSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['XSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['XSS']['data'] ] # convert from cm/s to nm/s
                debug("\tXSS: %s" % len(self.ELEMENTS['XSS']['data']))

            elif l == 7:
                if not 'TSS' in self.ELEMENTS: self.ELEMENTS['TSS'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['TSS']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['TSS']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['TSS']['data'] = [ x * 10000000 for x in self.ELEMENTS['TSS']['data'] ] # convert from cm/s to nm/s
                debug("\tTSS: %s" % len(self.ELEMENTS['TSS']['data']))

            elif l == 8:
                if not 'REX' in self.ELEMENTS: self.ELEMENTS['REX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['REX']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['REX']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['REX']['data'] = [ x * 10000000 for x in self.ELEMENTS['REX']['data'] ] # convert from cm/s to nm/s
                debug("\tREX: %s" % len(self.ELEMENTS['REX']['data']))

            elif l == 9:
                if not 'ZEX' in self.ELEMENTS: self.ELEMENTS['ZEX'] = {'data':[], 'samplerate':self.DT, 'filter':filter} 
                self.ELEMENTS['ZEX']['data'] = [0 for x in TEMPDATA ]
                self.ELEMENTS['ZEX']['data'].extend( TEMPDATA.values() )
                self.ELEMENTS['ZEX']['data'] = [ x * 10000000 for x in self.ELEMENTS['ZEX']['data'] ] # convert from cm/s to nm/s
                debug("\tZEX: %s" % len(self.ELEMENTS['ZEX']['data']))
            else:
                raise SystemExit('\n\nERROR: GreenFunctions(): generate(): Too many elements\n\n')


        self._save_to_db()
#}}}

    def plot(self,start=0,end=-1):
        """ Plot all traces in memory. They are containe in 
        the dictionary ELEMENTS in self.

        """
#{{{

        log( "---------------------------------------" )
        log("FKRPROG.PY: plot()")

        total = len(self.ELEMENTS)

        half = int(total/2)

        if half != total/2: half += 1

        now = 0
        for trace in self.ELEMENTS:

            try:
                now += 1
                pyplot.subplot(5,2,now)
                data = self.ELEMENTS[trace]['data']
                pyplot.plot(data[start:end])
                pyplot.legend([trace])
            except Exception,e:
                sys.exit('ERROR: problem plotting green functions.[%s => %s]' % (Exception,e) )

        pyplot.suptitle("Green Functions: depth:%s distance:%s" % (self.DEPTH,self.DISTANCE))
        pyplot.show()
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

if __name__ == "__main__":
    """ If we call this script directly, then output help.  """
#{{{
    print "Green's Functions Library:"
    print "\n\t** Not to run directly!!!! **\n"
    print "\t1) Import the library.  "
    print "\t\timport moment_tensor.fkrprog as fkr"
    print "\t2) Build object for GF's.  "
    print "\t\tgf = GreenFunctions(modle.pf)"
    print "\t3) Build for [distance,depth]"
    print "\t\tgf.build(depth=8,distance=1)"
    print "\t3) Get GF's matrix with elements"
    print "\t\tdata=gf.ALL"

    print "\n\n"
    print "No BRTT support."
    print "Juan Reyes <reyes@ucsd.edu>"
    print "Rob Newman <rlnewman@ucsd.edu>"

#}}}

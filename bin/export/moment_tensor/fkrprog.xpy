import os
import sys
import re
import getopt
import numpy
import pprint 
from numpy.linalg import inv
from collections import defaultdict
from math import exp, log, sqrt, acos, asin, cos, sin
import  math 
from cmath import log as clog
from cmath import sqrt as csqrt
from cmath import exp as cexp


"""
      FREQUENCY WAVENUMBER-INTEGRATION PROGRAM.

    ORIGINAL FORTRAN CODE WRITTEN BY CHANDAN K. SAIKIA
    TRANSLATED TO PYTHON BY JUAN C. REYES 1/2011 <REYES@UCSD.EDU>

                 REFRENCES :

     1.  DUNKIN J. W.      [1965] BSSA 335-358
     2.  HASKELL N. A.     [1964] BSSA 337-393
     3.  WANG AND HERRMANN [1980] BSSA 1015-1036
     4.  WATSON T. H.      [1970] BSSA 161-166

              WRITTEN BY CHANDAN K. SAIKIA
"""



class new_dd():
#{{{ Test autovivification of elements in variable self.
    """
    Test autovivification of elements in variable self.
    """
    def __getitem__(self, key):
        print 'Getitem: key:%s ' % key
        try:
            return self.__dict__.get(attr)
        except:
            self.__dict__[key] = defaultdict(dict)

    def __setattr__(self,key,value):
        print 'Set: key:%s value:%s' % (key,value)
        self.__dict__[key] = defaultdict(dict)
        self.__dict__[key][value] = defaultdict(dict)
#   
#   # TEST 1
#   print '\n'
#   print 'defaultdict'
#   d = defaultdict(dict)
#   d[3] = 2
#   print '%s' % d
#   d[0][1] = [1,2,3]
#   print '%s' % d
#   d[1] = ('test','value')
#   print '%s' % d
#   
#   # TEST 2
#   print '\n'
#   print 'with self == test'
#   self = new_dd()
#   self.d = 9
#   self.d[3] = 2
#   print '%s' % self.d
#   self.d[0][1] = [1,2,3]
#   print '%s' % self.d
#   self.d[1] = ('test','value')
#   print '%s' % self.d
#   
#   # TEST 3
#   class dd():
#       def __init__(self):
#           print '\n'
#           print 'with self'
#           self.d = defaultdict(dict)
#           self.d[3] = 2
#           print '%s' % self.d
#           self.d[0][1] = [1,2,3]
#           print '%s' % self.d
#           self.d[1] = ('test','value')
#           print '%s' % self.d
#   
#   test = dd()
#   
#   exit()
#   """
#   End of test
#   """
#}}}

def configure():
#{{{
    """
    Get all information from command line 
    and verify the source file is accessible.
    """
    try:

        opts, pargs = getopt.getopt(sys.argv[1:], 'vd')

    except getopt.GetoptError:

        print "Usage: fkrprog [-v] [-d] file\n"
        sys.exit(-1)

    if( len(pargs) != 1):

        print "Usage: fkrprog [-v] [-d] file\n"
        sys.exit(-1)

    else:

        SOURCE_FILE = pargs[0]

    globals()['debug'] = False
    globals()['verbose'] = False

    for option, value in opts:

        if '-v' in option:
            globals()['verbose'] = True
            globals()['debug'] = False

        if '-d' in option:
            globals()['verbose'] = True
            globals()['debug'] = True

    #   
    # Sanity check
    #   
    if not os.path.isfile(SOURCE_FILE):
        raise SystemExit('\n\nCannot find specified file! (%s)\n'% SOURCE_FILE)

    return SOURCE_FILE

#}}}

class green_functions():
    " Build class to produce Green Functions "
    def __init__(self,source=None,verbose=False,debug=False):
#{{{
        if debug: print "green_functions(): init..."
        #self = defaultdict(lambda:defaultdict(list))
        test = new_dd()

        #
        # Constants
        #
        self.IB  = 100
        self.IFREQ = 0

        self.verbose = verbose
        self.debug = debug

        self.source = None

        if source: 
            if not os.path.isfile(source):
                raise SystemExit('\n\nCannot find specified file! (%s)\n'% source)

            self.source = self._read_file(source)

        if self.source:
            self._setup_vars()
            self._excit()


#}}}

    def file_format(self):
    #{{{
        print "File expected should have this format: \n \
                .F. \n \
                    0   64 \n \
                GREEN.X \n \
                    6.0      8.00       1  512 1024    0.500     5    1 \n \
                    1    1    1    1    1    1    1    1    1    1    0 \n \
                0.5500E+01 0.5500E+01 0.3180E+01 0.2400E+01   600.00    300.00 \n \
                0.2500E+01 0.6300E+01 0.3640E+01 0.2670E+01   600.00    300.00 \n \
                0.8000E+01 0.6300E+01 0.3640E+01 0.2670E+01   600.00    300.00 \n \
                0.1900E+02 0.6700E+01 0.3870E+01 0.2800E+01   600.00    300.00 \n \
                0.4000E+03 0.7800E+01 0.4500E+01 0.3300E+01   600.00    300.00 \n \
                    3 \n \
                0.4000000E+03  1.500000E+00         0   \n \
                    40  10000.0     30.0      2.9       2.5 \n \
                10.00      0.0      10.0 \n \
                20.00      0.0      10.0 \n \
                30.00      0.0      10.0 \n \
                40.00      0.0      10.0 \n \
                50.00      0.0      10.0 \n \
                60.00      0.0      10.0 \n \
                70.00      0.0      10.0 \n \
                80.00      0.0      10.0 \n \
                90.00      0.0      10.0 \n \
                100.00      0.0      10.0 \n \
                110.00      0.0      10.0 \n \
                120.00      0.0      10.0 \n \
                130.00      0.0      10.0 \n \
                140.00      0.0      10.0 \n \
                150.00      0.0      10.0 \n \
                160.00      0.0      10.0 \n \
                170.00      0.0      10.0 \n \
                180.00      0.0      10.0 \n \
                190.00      0.0      10.0 \n \
                200.00      0.0      10.0 \n \
                210.00      0.0      10.0 \n \
                220.00      0.0      10.0 \n \
                230.00      0.0      10.0 \n \
                240.00      0.0      10.0 \n \
                250.00      0.0      10.0 \n \
                260.00      0.0      10.0 \n \
                270.00      0.0      10.0 \n \
                280.00      0.0      10.0 \n \
                290.00      0.0      10.0 \n \
                310.00      0.0      10.0 \n \
                320.00      0.0      10.0 \n \
                330.00      0.0      10.0 \n \
                340.00      0.0      10.0 \n \
                350.00      0.0      10.0 \n \
                360.00      0.0      10.0 \n \
                370.00      0.0      10.0 \n \
                380.00      0.0      10.0 \n \
                390.00      0.0      10.0 \n \
                400.00      0.0      10.0 \n "
    #}}}

    def _read_file(self, file_name):
    #{{{
        #
        # Read data from file
        #
        if self.debug: print "%s: Open file." % file_name

        try: file = open(file_name)
        except Exception,e: raise SystemExit('\n\nCannot open input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        # First line. Verbosity
        if self.debug: print "%s: Read verbose flag."
        try: (PRNTA) = file.readline().strip()
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if PRNTA == '.F.' or PRNTA == 'False' or PRNTA == '.FALSE.'  or PRNTA == '0':
            self.PRNTA = False
        elif PRNTA == '.T.' or PRNTA == 'True' or PRNTA == '.TRUE.'  or PRNTA == '1':
            self.PRNTA = True
        else:
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        # M line
        if self.debug: print "%s: Read frequency range M1 and M2."
        try: (self.M1, self.M2) = [ int(s) for s in file.readline().strip().split() ]
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: (self.M1,self.M2)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: M1=%s, M2=%s" % (file_name,self.M1,self.M2)

        # Name2
        if self.debug: print "%s: Read name of output file NAME2"

        try: (self.NAME2) = file.readline().strip()
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: self.NAME2
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: NAME2=%s" % (file_name,self.NAME2)

        #DECAY,DEPTH,N1,N2,N,DT,MMAX,NSKIP 
        try: (self.DECAY,self.DEPTH,self.N1,self.N2,self.N,self.DT,self.MMAX,self.NSKIP) = file.readline().strip().split()
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try:
            self.DECAY  = float(self.DECAY)
            self.DEPTH  = float(self.DEPTH)
            self.N      = int(self.N)
            self.N1     = int(self.N1)
            self.N2     = int(self.N2)
            self.DT     = float(self.DT)
            self.MMAX   = int(self.MMAX)
            self.NSKIP   = int(self.NSKIP)
        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: (self.DECAY,self.DEPTH,self.N1,self.N2,self.N,self.DT,self.MMAX,self.NSKIP)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: DECAY=%s DEPTH=%s N1=%s N2=%s N=%s DT=%s MMAX=%s NSKIP=%s" % (file_name,self.DECAY,self.DEPTH,self.N1,self.N2,self.N,self.DT,self.MMAX,self.NSKIP)


        # ISRC, JBDRY 
        try:
            self.ISRC = [ int(x) for x in file.readline().strip().split() ]
            self.JBDRY = self.ISRC.pop()
            #self.JBDRY = 0
        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))
        
        try: (self.ISRC,self.JBDRY)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: ISRC=%s JBDRY=%s" % (file_name,self.ISRC,self.JBDRY)


        # D(I),A(I),B(I),RHO(I),QA(I),QB(I)
        self.D   = [0 for x in range(self.MMAX)]
        self.A   = [0 for x in range(self.MMAX)]
        self.B   = [0 for x in range(self.MMAX)]
        self.RHO = [0 for x in range(self.MMAX)]
        self.QA  = [0 for x in range(self.MMAX)]
        self.QB  = [0 for x in range(self.MMAX)]

        if self.debug: print "ORIGINAL:"
        try:
            for I in range(self.MMAX):
                (self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I])= [float(x) for x in file.readline().strip().split()]
                if self.debug: print "%s: %E %E %E %E %s %s" % (file_name,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I])

            self.QA = [ 1/x for x in self.QA ]
            self.QB = [ 1/x for x in self.QB ]
        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))


        try: (self.D,self.A,self.B,self.RHO,self.QA,self.QB)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        #if self.debug: 
            #print "CONVERTED:"
            #for I in range(self.MMAX):
                #print "%s: %E %E %E %E %s %s" % (file_name,self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I])

        # LMAX
        try: self.LMAX = int(file.readline().strip())
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: self.LMAX
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: LMAX=%s" % (file_name,self.LMAX)

        # XLENG XFAC XXX
        try: temp =  file.readline().strip().split()
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if len(temp) == 2:
            (self.XLENG,self.XFAC) = [ float(x) for x in temp ]
            self.XXX = 0
        elif len(temp) == 3: 
            (self.XLENG,self.XFAC,self.XXX) = [ float(x) for x in temp ]

        try: (self.XLENG,self.XFAC,self.XXX)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: XLENG=%s XFAC=%s XXX=%s" % (file_name,self.XLENG,self.XFAC,self.XXX)

        # NRANGE,CMAX,C1,C2,CMIN
        try: (self.NRANGE,self.CMAX,self.C1,self.C2,self.CMIN)= [ float(x) for x in file.readline().strip().split() ]
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: self.NRANGE = int(self.NRANGE)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        try: (self.NRANGE,self.CMAX,self.C1,self.C2,self.CMIN)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        if self.debug: print "%s: NRANGE=%s CMAX=%s C1=%s C2=%s CMIN=%s" % (file_name,self.NRANGE,self.CMAX,self.C1,self.C2,self.CMIN)

        # RANGE(I), T0(I), VRED(I)
        self.RANGE = [0 for x in range(self.NRANGE)]
        self.T0    = [0 for x in range(self.NRANGE)]
        self.VRED  = [0 for x in range(self.NRANGE)]

        try:
            for I in range(self.NRANGE):
                (self.RANGE[I],self.T0[I],self.VRED[I])= [ float(x) for x in file.readline().strip().split() ]
                if self.debug: print "RANGE(%s)=%s" % (I,self.RANGE[I])
                if self.debug: print "T0(%s)=%s" % (I,self.T0[I])
                if self.debug: print "VRED(%s)=%s" % (I,self.VRED[I])
        except Exception,e: 
            raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))


        try: (self.RANGE,self.T0,self.VRED)
        except Exception,e: raise SystemExit('\n\nWrong Format of input file[%s]. %s(%s) \n'% (file_name,Exception,e))

        #if self.debug: 
            #for I in range(self.NRANGE):
                #print "%s: %f %f %f" % (file_name,self.RANGE[I],self.T0[I],self.VRED[I])

        if self.debug: print '\n\n########## END OF FILE [%s] ############\n\n' % file_name


        ##
        ## FIX VALUES TO MATCH PYTHON INDEXING
        ##
        self.LMAX -= 1 # layer number below the source
        self.MMAX -= 1 # Total number of layers
        #self.N1   -= 1 # the starting frequency
        #self.N2   -= 1 # the total niumber of frequencies
        #self.N    -= 1 # the total number of time points
        return True

    #}}}

    def _setup_vars(self):
    #{{{

        # Fix values for C1, C2, CMAX and CMIN
        if self.verbose: print "CMAX=%s" % self.CMAX
        try: self.CMAX = 1/self.CMAX
        except Exception,e: raise SystemExit('\n\nCMAX. %s(%s) \n'% (Exception,e))
        try: self.C1 = 1/self.C1
        except Exception,e: raise SystemExit('\n\nC1. %s(%s) \n'% (Exception,e))
        try: self.C2 = 1/self.C2
        except Exception,e: raise SystemExit('\n\nC2. %s(%s) \n'% (Exception,e))
        try: self.CMIN = 1/self.CMIN
        except Exception,e: raise SystemExit('\n\nCMIN. %s(%s) \n'% (Exception,e))

        if self.debug: print "CMAX=%s CMIN=%s C1=%s C2=%s" % (self.CMAX,self.CMIN,self.C1,self.C2)

        # FACT(I)
        self.FACT = [0.0 for x in range(self.IB)]
        if self.CMAX < 0.0:
            for I in range(self.IB):
                self.FACT[I] = 1.0

        if self.debug: print "FACT=%s" % self.FACT

        #
        # CHECK IF DEPTH DOES CORRESPOND TO DEPTH INTERFACE
        #
        self.DPH = 0.0
        if self.debug: print "Adding layers on top of source[%s]" % self.LMAX

        if self.debug: print "\tLMAX= %s" % self.LMAX
        for I in range(self.LMAX):
            self.DPH += self.D[I]
            if self.debug: print "\tDPH= %s" % self.DPH

        if abs(self.DPH-self.DEPTH) <= 0.1: self.DEPTH = self.DPH

        if abs(self.DPH-self.DEPTH) > 0.1:
            raise SystemExit('\n\nCHECK TWO PARAMETERS <+> LMAX [%s] AND DEPTH [%s]\n'% (self.LMAX,self.DEPTH))

        #
        # CALCULATE ALPHA FROM DECAY
        #
        self.TLEN = self.N * self.DT
        if self.debug: print 'TLEN = %s' % self.TLEN

        try: self.DECAY = 1 / self.DECAY
        except Exception,e: raise SystemExit('\n\nDECAY. %s(%s) \n'% (Exception,e))

        if self.debug: print 'DECAY = %s' % self.DECAY

        try:
            self.ALPHA = -log(self.DECAY)
            self.ALPHA = self.ALPHA / self.TLEN
        except Exception,e: raise SystemExit('\n\nALPHA. %s(%s) \n'% (Exception,e))

        if self.debug: print 'ALPHA = %s' % self.ALPHA

        #
        # CHECK BOUCHON'S CRITERIA
        #
        try:
            self.RMAX  = max(self.RANGE)
            self.VAMIN = min(self.A)
            self.VAMAX = max(self.A)
            self.VBMIN = min(self.B)
            self.VBMAX = max(self.B)
        except Exception,e: raise SystemExit('\n\nRMAX VMIN VMAX VBMIN VBMAX. %s(%s) \n'% (Exception,e))

        #if self.debug: print 'A= %s' % self.A
        #if self.debug: print 'B= %s' % self.B
        #if self.debug: print 'RANGE= %s' % self.RANGE
        if self.debug: print 'RMAX = %s' % self.RMAX
        if self.debug: print 'VAMAX = %s' % self.VAMAX
        if self.debug: print 'VBMAX = %s' % self.VBMAX
        if self.debug: print 'VAMIN = %s' % self.VAMIN
        if self.debug: print 'VBMIN = %s' % self.VBMIN

        #for J in range(self.NRANGE): self.RMAX = max(self.RANGE[I])
        self.RMAX = max(self.RANGE)

        if self.debug: print 'RANGE= %s' % self.RANGE
        self.YMAX = 2.0 * self.RMAX
        if self.debug: print 'RMAX= %s' % self.RMAX
        if self.debug: print 'YMAX = %s' % self.YMAX

        if self.XLENG < self.YMAX: self.XLENG = self.YMAX
        if self.debug: print 'XLENG = %s' % self.XLENG

        self.RMAX = self.RMAX + sqrt( self.VAMAX * self.VAMAX * self.TLEN * self.TLEN - self.DEPTH * self.DEPTH )
        if self.debug: print 'RMAX = %s' % self.RMAX

        if self.XLENG < self.RMAX:
                self.XLENG=self.RMAX
                self.MM = int(self.XLENG / 100.0)
                self.XLL = int(self.MM) * 100.0
                if self.XLL != self.XLENG: 
                    self.XLENG = self.XLL + 100.0
        if self.debug: 
            print 'MM = %s' % self.MM
            print 'XLL = %s' % self.XLL
            print 'XLENG = %s' % self.XLENG

        try: self.DF = 1 / (self.N * self.DT)
        except Exception,e: raise SystemExit('\n\nDF. %s(%s) \n'% (Exception,e))

        if self.debug: print 'DF = %s' % self.DF

        try: self.NYQ = self.N / 2 + 1
        except Exception,e: raise SystemExit('\n\nNYQ. %s(%s) \n'% (Exception,e))

        if self.debug: print 'NYQ = %s' % self.NYQ

        self.FL =0.0

        if self.debug: print 'FL = %s' % self.FL

        try: self.FU = (self.NYQ - 1) * self.DF
        except Exception,e: raise SystemExit('\n\nFU. %s(%s) \n'% (Exception,e))

        if self.debug: print 'FU = %s' % self.FU


        #
        # Print values
        #
        #print ' FL =%-10.5f\tFU =%-10.5f\tDF =%-10.5f\n\tN1 =%-4s\tN2 =%-4s\tDEPTH =%-10.2f\tN =%-5s\tNSKIP =%-5s' % (self.FL,self.FU,self.DF,self.N1,self.N2,self.DEPTH,self.N,self.NSKIP)
        #print ' ALFPHA =%-10.5f\tDT =%-10.3f  XLENG= %-12.4E  XFAC = %-6.2f ' % (self.ALPHA,self.DT,self.XLENG,self.XFAC)
        #for I in range(self.MMAX+1):
            #print "%-11.4E %-11.4E %-11.4E %-11.4E %-11.4E %-11.4E" % (self.D[I],self.A[I],self.B[I],self.RHO[I],self.QA[I],self.QB[I])
        #print ' FREQUENCIES FOR WHICH RESPONSE COMPUTED'

        #
        # Test for file 
        #
        if os.path.isfile(self.NAME2):
            self.ISTAT3 = 'OLD'
            self.IXST3 = True
        else:
            self.ISTAT3 = 'NEW'
            self.IXST3 = False

        if self.IFREQ > self.N1: self.N1 = self.IFREQ + 1

        if not self.NSKIP: self.NSKIP = 1

        if self.debug: print 'N1 = %s' % self.N1
        if self.debug: print 'N2 = %s' % self.N2
        if self.debug: print 'NSKIP = %s' % self.NSKIP
        if self.debug: print 'DF = %s' % self.DF

        for x in range(self.N1,self.N2,self.NSKIP):
            if self.debug: print 'I=%s' % x
            self.FREQ = (x-1) * self.DF
            if self.FREQ < self.DF: self.FREQ = 0.01 * self.DF
            if self.debug: print 'FREQ => %s' % self.FREQ
            self.OMEGA = 2.0 * acos(-1.0) * self.FREQ
            if self.debug: print 'OMEGA => %s' % self.OMEGA

            #print 'X = %s' % x
            #print 'DF = %s' % self.DF
            #print 'FREQ => %s' % self.FREQ
            #print 'OMEGA => %s' % self.OMEGA

            self.WVCM = self.OMEGA * self.CMAX
            self.WVC1 = self.OMEGA * self.C1
            self.WVC2 = self.OMEGA * self.C2
            self.WVCN = self.OMEGA * self.CMIN

            if self.debug: print 'WVCM => %s' % self.WVCM
            if self.debug: print 'WVC1 => %s' % self.WVC1
            if self.debug: print 'WVC2 => %s' % self.WVC2
            if self.debug: print 'WVCN => %s' % self.WVCN

            self._excit()
            print "\t\t\t\t%s  F =  %12.4E NBLOCK =  %5s IBLOCK= %3s" % (x,self.FREQ,self.NBLOCK,self.IBLOCK)
            if x >= self.M1 and x <= self.M2:
                if self.PRNTA or self.verbose:
                    for y in range(self.NRANGE):
                        print "\tDIST =>    %s" % y
                        for z in range(10):
                            #print "\t\t%-10s%-16.9E%-4s%-16.9E" % self.GG[y][z]
                            print "\t\t%16.9E\t%16.9E" % (self.GG[y][z].real,self.GG[y][z].imag)


        self.XX = -999.0
        self.T0X = 0.0

        #print "%s %s" % (self.XX,self.T0X)


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
        #test = defaultdict(lambda:defaultdict(list))
        self.XKA = defaultdict(lambda: defaultdict(dict))
        self.XKB = defaultdict(lambda: defaultdict(dict))
        self.ATNA = defaultdict(lambda: defaultdict(dict))
        self.ATNB = defaultdict(lambda: defaultdict(dict))
        self.ICNT = defaultdict(lambda: defaultdict(dict))
        self.WVNO  = defaultdict(lambda: defaultdict(dict))
        self.WVNO2 = defaultdict(lambda: defaultdict(dict))
        self.WVNO3 = defaultdict(lambda: defaultdict(dict))
        self.SMM = defaultdict(lambda: defaultdict(dict))
        self.R1 = defaultdict(lambda: defaultdict(dict))
        self.R2 = defaultdict(lambda: defaultdict(dict))
        self.R3 = defaultdict(lambda: defaultdict(dict))
        self.R4 = defaultdict(lambda: defaultdict(dict))
        self.R5 = defaultdict(lambda: defaultdict(dict))
        self.R6 = defaultdict(lambda: defaultdict(dict))
        self.R7 = defaultdict(lambda: defaultdict(dict))
        self.X1 = defaultdict(lambda: defaultdict(dict))
        self.X2 = defaultdict(lambda: defaultdict(dict))
        self.X3 = defaultdict(lambda: defaultdict(dict))
        self.X4 = defaultdict(lambda: defaultdict(dict))
        self.X5 = defaultdict(lambda: defaultdict(dict))
        self.X6 = defaultdict(lambda: defaultdict(dict))
        self.X7 = defaultdict(lambda: defaultdict(dict))
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
        self.S14 = defaultdict(lambda: defaultdict(dict))
        self.S21 = defaultdict(lambda: defaultdict(dict))
        self.S32 = defaultdict(lambda: defaultdict(dict))
        self.S34 = defaultdict(lambda: defaultdict(dict))
        self.COSQL = defaultdict(lambda: defaultdict(dict))
        self.YL = defaultdict(lambda: defaultdict(dict))
        self.YLI = defaultdict(lambda: defaultdict(dict))
        self.ZL = defaultdict(lambda: defaultdict(dict))
        self.ZX = defaultdict(lambda: defaultdict(dict))
        self.K = defaultdict(lambda: defaultdict(dict))
        self.Y11 = defaultdict(lambda: defaultdict(dict))
        self.Y12 = defaultdict(lambda: defaultdict(dict))
        self.Y31 = defaultdict(lambda: defaultdict(dict))
        self.Y32 = defaultdict(lambda: defaultdict(dict))
        self.Y21 = defaultdict(lambda: defaultdict(dict))
        self.Y22 = defaultdict(lambda: defaultdict(dict))
        self.Y41 = defaultdict(lambda: defaultdict(dict))
        self.Y42 = defaultdict(lambda: defaultdict(dict))
        self.S32E = defaultdict(lambda: defaultdict(dict))
        self.S34E = defaultdict(lambda: defaultdict(dict))
        self.SNGL = defaultdict(lambda: defaultdict(dict))
        self.FACX = defaultdict(lambda: defaultdict(dict))
        self.SUME = defaultdict(lambda: defaultdict(dict))
        self.SUMC = defaultdict(lambda: defaultdict(dict))
        self.C = defaultdict(lambda: defaultdict(dict))
        self.S1 = defaultdict(lambda: defaultdict(dict))
        self.S2 = defaultdict(lambda: defaultdict(dict))
        self.S3 = defaultdict(lambda: defaultdict(dict))
        self.S4 = defaultdict(lambda: defaultdict(dict))
        self.S5 = defaultdict(lambda: defaultdict(dict))
        self.S6 = defaultdict(lambda: defaultdict(dict))
        self.S7 = defaultdict(lambda: defaultdict(dict))
        self.S8 = defaultdict(lambda: defaultdict(dict))
        self.S9 = defaultdict(lambda: defaultdict(dict))
        self.S10 = defaultdict(lambda: defaultdict(dict))
        self.G1 = defaultdict(lambda: defaultdict(dict))
        self.G2 = defaultdict(lambda: defaultdict(dict))
        self.G3 = defaultdict(lambda: defaultdict(dict))
        self.G4 = defaultdict(lambda: defaultdict(dict))
        self.G5 = defaultdict(lambda: defaultdict(dict))
        self.G6 = defaultdict(lambda: defaultdict(dict))
        self.G7 = defaultdict(lambda: defaultdict(dict))
        self.G8 = defaultdict(lambda: defaultdict(dict))
        self.G9 = defaultdict(lambda: defaultdict(dict))
        self.G10 = defaultdict(lambda: defaultdict(dict))

        self.GG = [complex(0) for J in range(100)]
        for J in range(100): self.GG[J] = [complex(0) for I in range(10)]

        self.S1 = [complex(0) for J in range(100)]
        self.S2 = [complex(0) for J in range(100)]
        self.S3 = [complex(0) for J in range(100)]
        self.S4 = [complex(0) for J in range(100)]
        self.S5 = [complex(0) for J in range(100)]
        self.S6 = [complex(0) for J in range(100)]
        self.S7 = [complex(0) for J in range(100)]
        self.S8 = [complex(0) for J in range(100)]
        self.S9 = [complex(0) for J in range(100)]
        self.S10 = [complex(0) for J in range(100)]

        self.G1 = [complex(0) for J in range(100)]
        self.G2 = [complex(0) for J in range(100)]
        self.G3 = [complex(0) for J in range(100)]
        self.G4 = [complex(0) for J in range(100)]
        self.G5 = [complex(0) for J in range(100)]
        self.G6 = [complex(0) for J in range(100)]
        self.G7 = [complex(0) for J in range(100)]
        self.G8 = [complex(0) for J in range(100)]
        self.G9 = [complex(0) for J in range(100)]
        self.G10 = [complex(0) for J in range(100)]

        self.J0 = [complex(0) for J in range(100)]
        self.J1 = [complex(0) for J in range(100)]
        self.J2 = [complex(0) for J in range(100)]
        self.SUME = [complex(0) for J in range(100)]
        self.SUMC = [complex(0) for J in range(100)]
        self.WVNO = [complex(0) for J in range(100)]


        self.FACT = [complex(0) for J in range(100)]
        self.FACX = [complex(0) for J in range(100)]
        self.S12 = [complex(0) for J in range(100)]
        self.S21 = [complex(0) for J in range(100)]
        self.S32 = [complex(0) for J in range(100)]
        self.S14 = [complex(0) for J in range(100)]
        self.S23 = [complex(0) for J in range(100)]
        self.S34 = [complex(0) for J in range(100)]
        self.S32E = [complex(0) for J in range(100)]
        self.S34E = [complex(0) for J in range(100)]

        self.R1 = [complex(0) for J in range(100)]
        self.R2 = [complex(0) for J in range(100)]
        self.R3 = [complex(0) for J in range(100)]
        self.R4 = [complex(0) for J in range(100)]
        self.R5 = [complex(0) for J in range(100)]
        self.R6 = [complex(0) for J in range(100)]
        self.R7 = [complex(0) for J in range(100)]

        self.X1 = [complex(0) for J in range(100)]
        self.X2 = [complex(0) for J in range(100)]
        self.X3 = [complex(0) for J in range(100)]
        self.X4 = [complex(0) for J in range(100)]
        self.X5 = [complex(0) for J in range(100)]
        self.X6 = [complex(0) for J in range(100)]
        self.X7 = [complex(0) for J in range(100)]

        self.Y11 = [complex(0) for J in range(100)]
        self.Y12 = [complex(0) for J in range(100)]
        self.Y21 = [complex(0) for J in range(100)]
        self.Y31 = [complex(0) for J in range(100)]
        self.Y41 = [complex(0) for J in range(100)]
        self.Y42 = [complex(0) for J in range(100)]
        self.Y32 = [complex(0) for J in range(100)]
        self.Y22 = [complex(0) for J in range(100)]


        self.WVNO = [complex(0) for J in range(100)]
        self.WVNO2 = [complex(0) for J in range(100)]
        self.WVNO3 = [complex(0) for J in range(100)]
        self.ATNA = [complex(0) for J in range(70)]
        self.ATNB = [complex(0) for J in range(70)]
        self.XKA = [complex(0) for J in range(70)]
        self.XKB = [complex(0) for J in range(70)]

        self.AA = [0 for J in range(10)]
        self.BB = [0 for J in range(10)]
        self.CC = [0 for J in range(10)]
#}}}

        self.IASYMP = True
        self.ID = 0
        self.ISYST = 0
        for I in range(8):
            if self.ISRC[I] == 1: self.ID = 1
        if self.debug: print 'ISRC => %s' % self.ISRC

        if self.ID == 1: self.ISYST = 1
        if self.debug: print 'ISYST => %s' % self.ISYST
        self.ID = 0
        if self.ISRC[8] == 1: self.ID=1
        if self.ISRC[9] == 1: self.ID=1
        if self.ID == 0 and self.ISYST == 0: return
        if self.ID == 1 and self.ISYST == 0: self.ISYST = 2
        if self.ID == 1 and self.ISYST == 1: self.ISYST = 3

        if self.debug: print 'ISYST => %s' % self.ISYST

        self.DK = 2.0 * acos(-1.0) / self.XLENG
        if self.debug: print 'DK => %s' % self.DK

        self.WVBM = self.OMEGA/self.VBMIN
        self.WVMM = (5.0/self.DEPTH)+self.XFAC*self.WVBM
        self.NK   = self.WVMM/self.DK
        self.NK   = int(self.NK+2)
        if self.debug: print 'WVBM => %s' % self.WVBM
        if self.debug: print 'WVMM => %s' % self.WVMM
        if self.debug: print 'NK => %s' % self.NK

        self.V = 6.0 / self.DEPTH

        self.WVNO[1]  = complex( self.V, 0 )
        self.WVNO2[1] = self.WVNO[1] * self.WVNO[1]
        self.WVNO3[1] = self.WVNO2[1] * self.WVNO[1]
        self.AKL      = ( self.NK - 3 ) * self.DK + 0.218 * self.DK

        if self.debug: print 'WVNO[1] => %s' % self.WVNO[1]
        if self.debug: print 'WVNO2[1] => %s' % self.WVNO2[1]
        if self.debug: print 'WVNO3[1] => %s' % self.WVNO3[1]
        if self.debug: print 'AKL => %s' % self.AKL


        if self.AKL > self.V: self.IASYMP = False
        if self.debug: print 'IASYMP => %s' % self.IASYMP

        self.V = 2.5 / self.DEPTH
        if self.debug: print 'V => %s' % self.V
        self.WVNO[0]  = complex( self.V , 0 )
        self.WVNO2[0] = self.WVNO[0] * self.WVNO[0]
        if self.debug: print 'WVNO[1] => %s' % self.WVNO[0]
        if self.debug: print 'WVNO2[1] => %s' % self.WVNO2[0]

        if self.CMAX > 0.0: self.MM = int(self.WVCN / self.DK) + 2
        if self.NK > self.MM: self.NK = self.MM
        if self.debug: print 'MM => %s' % self.MM
        if self.debug: print 'NK => %s' % self.NK
        if self.debug: print 'OMEGA => %s' % self.OMEGA
        if self.debug: print 'OMEGA:%s NK:%s' % (self.OMEGA,self.NK)

        #
        # DETERMINE INDEX FOR INTEGRATION  I.E., WHEN KR IS GT. 3.0
        #
#{{{
        
        if self.debug: print 'NRANGE => %s' % self.NRANGE
        for I in range(self.NRANGE):
            for J in range(1,50000):
                self.WVV = ( J - 1 ) * self.DK + 0.218 * self.DK
                self.AKR = self.WVV * self.RANGE[I]
                if self.debug: print ''
                if self.debug: print 'I => %s' % I
                if self.debug: print 'J => %s' % J
                if self.debug: print 'RANGE(I) => %s' % self.RANGE[I]
                if self.debug: print 'WVV => %s' % self.WVV
                if self.debug: print 'AKR => %s' % self.AKR
                if self.AKR <= 3.0: self.ICNT[I] = J + 2
                #if self.AKR <= 3.0: print '\tICNT[%s] = %s' %  (I,self.ICNT[I])
                if self.AKR > 3.0: break
        if self.debug: print 'ICNT => %s' % self.ICNT

        self.OM     = complex(self.OMEGA, -1 * self.ALPHA)
        self.FOURPI = 12.5663706 * self.OM * self.OM
        self.ZERO   = complex( 0.0 , 0 )
        self.ONE    = complex( 1.0 , 0 )
        self.TWO    = complex( 2.0 , 0 )
        self.PI     = 3.141592654
        self.OM1    = 6.283185307
        self.OML    = 0.06283185307

        if self.debug: print 'OM => %s' % self.OM
        if self.debug: print 'FOURPI => %s' % self.FOURPI
        if self.debug: print 'ZERO => %s' % self.ZERO
        if self.debug: print 'ONE => %s' % self.ONE
        if self.debug: print 'TWO => %s' % self.TWO
        if self.debug: print 'PI => %s' % self.PI
        if self.debug: print 'OM1 => %s' % self.OM1
        if self.debug: print 'OML => %s' % self.OML
#}}}

        #
        # COMPUTING THE ATTENUATION OF THE MDEIUM USING QA, QB
        #

        self.AT = self.ZERO
        #print 'OM => %s' % self.OM
        #print 'abs(OM) => %s' % abs(self.OM)
        #print 'OML => %s' % self.OML
        #print 'OM => %s' % self.OM
        #print 'OM1 => %s' % self.OM1
        if (abs(self.OM) > self.OML ): self.AT=log(self.OM.real/self.OM1.real)/self.PI.real
        if self.debug: print 'AT => %s' % self.AT
        if (abs(self.OM) < self.OML ):
            self.FA = sqrt(self.OML**2 + self.ALPHA**2)/self.OML
            self.FAC = float( self.FA )
            self.AT = clog(complex( self.OML, -self.ALPHA)/(self.OM1*self.FAC)) / self.PI

        if self.debug: print 'AT => %s' % self.AT
        if self.debug: print 'FA => %s' % self.FA
        if self.debug: print 'FAC => %s' % self.FAC

        for I in range(self.MMAX+1):
            self.ATNA[I] = self.ONE
            self.ATNB[I] = self.ONE
            if self.debug: print 'ATNA[I] => %s' % self.ATNA[I]
            if self.debug: print 'ATNB[I] => %s' % self.ATNB[I]
            if self.QA[I] > 0.0:
                self.PAL = complex(self.QA[I],0) * self.AT
                self.BA = float(0.5 * self.QA[I])
                self.QAL = complex(0, self.BA)
                self.ATNA[I] = self.ONE + self.PAL + self.QAL
                if self.debug: print 'PAL => %s' % self.PAL 
                if self.debug: print 'BA => %s' % self.BA 
                if self.debug: print 'QAL => %s' % self.QAL 
                if self.debug: print 'ATNA[I] => %s' % self.ATNA[I]

            if self.QB[I] > 0.0:
                self.PAL = complex(self.QB[I],0) * self.AT
                self.BA = float(0.5 * self.QB[I])
                self.QAL = complex(0, self.BA)
                self.ATNB[I] = self.ONE + self.PAL + self.QAL
                if self.debug: print 'PAL => %s' % self.PAL 
                if self.debug: print 'BA => %s' % self.BA 
                if self.debug: print 'QAL => %s' % self.QAL 
                if self.debug: print 'ATNB[I] => %s' % self.ATNB[I]

            self.XKA[I] = self.OM / ( float(self.A[I]) * self.ATNA[I] )
            self.XKB[I] = self.OM / ( float(self.B[I]) * self.ATNB[I] )
            if self.debug: print 'XKA[I] => %s' % self.XKA[I]
            if self.debug: print 'XKB[I] => %s' % self.XKB[I]

        if self.debug: print 'PAL => %s' % self.PAL
        if self.debug: print 'BA => %s' % self.BA
        if self.debug: print 'QAL => %s' % self.QAL
        if self.debug: print 'ATNA => %s' % self.ATNA
        if self.debug: print 'ATNB => %s' % self.ATNB
        if self.debug: print 'XKA => %s' % self.XKA
        if self.debug: print 'XKB => %s' % self.XKB


        if self.debug: print 'LMAX => %s' % self.LMAX
        if self.debug: print 'RHO => %s' % self.RHO
        if self.debug: print 'RHO[0] => %s' % self.RHO[0]
        if self.debug: print 'RHO[3] => %s' % self.RHO[3]
        if self.debug: print 'RHO[LMAX] => %s' % self.RHO[self.LMAX]
        if self.debug: print 'complex(RHO[LMAX]) => %s' % complex(self.RHO[self.LMAX])
        self.FRHO = complex( self.RHO[self.LMAX] ) * self.FOURPI
        self.XKA2 = self.XKA[self.LMAX] * self.XKA[self.LMAX]
        self.XKB2 = self.XKB[self.LMAX] * self.XKB[self.LMAX]
        self.XKK = self.ATNB[self.LMAX] * self.ATNB[self.LMAX]
        self.FACXX = 1.0 / ( 12.5663706 * self.B[self.LMAX] ** 2 )

        if self.debug: print 'FRHO => %s' % self.FRHO
        if self.debug: print 'XKA2 => %s' % self.XKA2
        if self.debug: print 'XKB2 => %s' % self.XKB2
        if self.debug: print 'XKK => %s' % self.XKK
        if self.debug: print 'FACXX => %s' % self.FACXX

        if self.NK < self.IB:
            self.IBLOCK = self.NK
            self.NBLOCK = 1
            self.LEFT = 0
        else:
            self.IBLOCK = self.IB
            self.NBLOCK = self.NK / self.IB
            self.LEFT = self.NK - self.NBLOCK * self.IB
            if self.LEFT: self.NBLOCK += 1

        if self.debug: print 'IBLOCK => %s' % self.IBLOCK
        if self.debug: print 'NBLOCK => %s' % self.NBLOCK
        if self.debug: print 'LEFT => %s' % self.LEFT
        if self.debug: print '\n\n'

        self.NBLK = 1
        self.LR = 0
        #}}}

        # Loop includes half-space also
        self.WVNO3[0] = self.WVNO2[0] * self.WVNO[0]
        self.NN1 = 2

        # 500 loop.
        if self.debug: print "Loop 500: while: %s <= %s" % (self.NBLK,self.NBLOCK)
        while self.NBLK <= self.NBLOCK:
            if self.debug: print "\tLoop 500: while: self.NBLK[%s] ?< self.NBOCLK[%s]" % (self.NBLK,self.NBLOCK)
        #{{{ A

            if self.debug: print "Loop 500 section AR => self.ISYST:%s" % self.ISYST
            #{{{ AR

            if self.NBLK == self.NBLOCK and self.LEFT != 0: self.IBLOCK = self.LEFT
            if self.NBLK > 1: self.NN1 = 0
            if self.debug: print "\tLoop 500: NBKL = %s" % self.NBLK
            if self.debug: print "\tLoop 500: IBLOCK = %s" % self.IBLOCK
            if self.debug: print "\tLoop 500: N1 = %s" % self.NN1

            if self.debug: print "\tLoop 500: self.NN1:[%s] self.IBLOCK[%s]" % (self.NN1,self.IBLOCK)
            #print "\tLoop 500: self.NN1 = %s" % self.NN1
            #print "\tLoop 500: self.IBLOCK = %s" % self.IBLOCK
            for J in range(self.NN1,self.IBLOCK):
            #{{{ D
                #print "\tLoop 500: J = %s" % J
                self.LR += 1
                self.WV = (self.LR - 1) * self.DK + 0.218 * self.DK
                self.WVNO[J] = complex(self.WV)
                self.WVNO2[J] = self.WVNO[J]**2
                self.WVNO3[J] = self.WVNO2[J] * self.WVNO[J]
                if self.debug: print "\tLoop 500: J = %s" % J
                if self.debug: print "\tLoop 500: LR = %s" % self.LR
                if self.debug: print "\tLoop 500: WVNO[J] = %s" % self.WVNO[J]
                if self.debug: print "\tLoop 500: WVNO2[J] = %s" % self.WVNO2[J]
                if self.debug: print "\tLoop 500: WVNO3[J] = %s" % self.WVNO3[J]
                if self.debug: print "\tLoop 500: WVNO2 = %s" % self.WVNO2
                #if J == 4: exit()
                #print "\tLoop 500: WVNO2[J] = %s" % self.WVNO2[J]
            #}}} D

            #if self.debug: print "self.WVNO2:%s" % self.WVNO2
            #exit()

            for J in range(self.IBLOCK):
            #{{{ E
                self.R1[J] = self.ZERO
                self.R2[J] = self.ZERO
                self.R3[J] = self.ZERO
                self.R4[J] = self.ZERO
                self.R5[J] = self.ZERO
                self.R6[J] = self.ZERO
                self.R7[J] = self.ZERO
                self.X1[J] = self.ZERO
                self.X2[J] = self.ZERO
                self.X3[J] = self.ZERO
                self.X4[J] = self.ZERO
                self.X5[J] = self.ZERO
                self.X6[J] = self.ZERO
                self.X7[J] = self.ZERO
            #}}} E

            #print "self.WVNO2:%s" % self.WVNO2
            if self.debug: print "self.WVNO2:%s" % self.WVNO2
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
                    self.GAMM1[I][J] =  self.GAM[I][J] - self.ONE
                    self.GAMM2[I][J] =  self.GAMM1[I][J] * self.GAMM1[I][J]
                    self.GAMM3[I][J] =  self.GAMM1[I][J] * self.GAMM2[I][J]
                    self.GAMX1[I][J] =  self.GAM[I][J] * self.GAMM1[I][J]
                    self.GAMX2[I][J] =  self.GAM2[I][J] * self.GAMM2[I][J]
                    if self.debug: print ""
                    if self.debug: print "Loop 500: RA[%s][%s]=> %s" % (I,J,self.RA[I][J])
                    if self.debug: print "Loop 500: RB[%s][%s]=> %s" % (I,J,self.RB[I][J])
                    if self.debug: print "Loop 500: GAM[%s][%s]=> %s" % (I,J,self.GAM[I][J])
                    if self.debug: print "Loop 500: GAM2[%s][%s]=> %s" % (I,J,self.GAM2[I][J])
                    if self.debug: print "Loop 500: GAM3[%s][%s]=> %s" % (I,J,self.GAM3[I][J])
                    if self.debug: print "Loop 500: GAM4[%s][%s]=> %s" % (I,J,self.GAM4[I][J])
                    if self.debug: print "Loop 500: GAMM1[%s][%s]=> %s" % (I,J,self.GAMM1[I][J])
                    if self.debug: print "Loop 500: GAMM2[%s][%s]=> %s" % (I,J,self.GAMM2[I][J])
                    if self.debug: print "Loop 500: GAMM3[%s][%s]=> %s" % (I,J,self.GAMM3[I][J])
                    if self.debug: print "Loop 500: GAMX1[%s][%s]=> %s" % (I,J,self.GAMX1[I][J])
                    if self.debug: print "Loop 500: GAMX2[%s][%s]=> %s" % (I,J,self.GAMX2[I][J])
                    if self.debug: print ""
                #}}} F

            #if self.debug: print "self.WVNO2:%s" % self.WVNO2
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

            if self.debug: print "Loop 500 section BR => self.ISYST:%s" % self.ISYST
            if self.ISYST != 2:
            #{{{ BR
                # Computation for P-SV system starts
                #{{{ H
                for J in range(self.IBLOCK):
                    self.S14[J] = -self.TWO * self.WVNO[J] / self.FOURPI
                    self.S21[J] =  self.TWO * self.XKB2 / self.FRHO
                    self.S32[J] = self.WVNO[J] * 4.0 * self.XKA2 / self.FRHO
                    self.TT = (2.0 * self.B[self.LMAX] / self.A[self.LMAX]) ** 2 - 3.
                    #print "\t\t\t B(LMAX)=>%s" % self.B[self.LMAX]
                    #print "\t\t\t A(LMAX)=>%s" % self.A[self.LMAX]
                    #print "\t\t\t TT=>%s" % self.TT
                    self.S34[J] = self.TWO * self.WVNO[J] * complex(self.TT) / self.FOURPI 
                    #print "\t\t\t TWO=>%s" % self.TWO
                    #print "\t\t\t WVNO[J]=>%s" % self.WVNO[J]
                    #print "\t\t\t TT=>%s" % self.TT
                    #print "\t\t\t complex(self.TT)=>%s" % complex(self.TT)
                    #print "\t\t\t FOURPI=>%s" % self.FOURPI
                    if self.debug: print ""
                    if self.debug: print "Loop 500: S14=> %s" % self.S14
                    if self.debug: print "Loop 500: S21=> %s" % self.S21
                    if self.debug: print "Loop 500: S32=> %s" % self.S32
                    if self.debug: print "Loop 500: S34=> %s" % self.S34
                    if self.debug: print ""
                #}}} H

                #     SET UP HALFSPACE BOUNDARY CONDITIONS
                #     JBDRY=-1 RIGID;  0= ELASTIC;  +1= FREE SURFACE
                #{{{ I 
                M = self.MMAX
                if self.debug: print "M:[%s] " % M
                if self.debug: print "RHO:%s" % self.RHO
                if self.JBDRY == 0:
                    for J in range(self.IBLOCK):
                        if self.debug: print "M:%s J:%s" % (M,J)
                        #if self.debug: print " self.GAM2[M][J] :%s" % self.GAM2[M][J]
                        #if self.debug: print " self.RA[M][J] :%s" % self.RA[M][J]
                        #if self.debug: print " self.RB[M][J] :%s" % self.RB[M][J]
                        #if self.debug: print " self.WVNO2[J] :%s" % self.WVNO2[J]
                        #if self.debug: print " self.GAMM2[M][J] :%s" % self.GAMM2[M][J]
                        self.RH = self.RHO[M]**2
                        if self.debug: print " self.RH :%s" % self.RH
                        self.R1[J] = self.RH * (-self.GAM2[M][J] * self.RA[M][J] * self.RB[M][J] + self.WVNO2[J] * self.GAMM2[M][J])
                        self.R2[J] = -1 * complex(self.RHO[M]) * self.WVNO2[J] * self.RA[M][J]
                        self.R3[J] = -1 * self.GAM[M][J] * self.RA[M][J] * self.RB[M][J] + self.WVNO2[J] * self.GAMM1[M][J]
                        self.R3[J] = complex(self.RHO[M]) * self.R3[J]
                        self.R4[J] = complex(self.RHO[M]) * self.WVNO2[J] * self.RB[M][J]
                        self.R5[J] = self.WVNO2[J] * ( self.WVNO2[J] - self.RA[M][J] * self.RB[M][J])
                        if self.debug: print "Loop 500: R1[J]=> %s" % self.R1[J]
                        if self.debug: print "Loop 500: R2[J]=> %s" % self.R2[J]
                        if self.debug: print "Loop 500: R3[J]=> %s" % self.R3[J]
                        if self.debug: print "Loop 500: R4[J]=> %s" % self.R4[J]
                        if self.debug: print "Loop 500: R5[J]=> %s" % self.R5[J]
                        if self.debug: print "\n"
                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R1[J] = self.ONE
                        self.R2[J] = self.ZERO
                        self.R3[J] = self.ZERO
                        self.R4[J] = self.ZERO
                        self.R5[J] = self.ZERO
                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R1[J] = self.ZERO
                        self.R2[J] = self.ZERO
                        self.R3[J] = self.ZERO
                        self.R4[J] = self.ZERO
                        self.R5[J] = self.ONE

                if self.debug: print "Loop 500: RH=> %s" % self.RH
                if self.debug: print "Loop 500: R1=> %s" % self.R1
                if self.debug: print "Loop 500: R2=> %s" % self.R2
                if self.debug: print "Loop 500: R3=> %s" % self.R3
                if self.debug: print "Loop 500: R4=> %s" % self.R4
                if self.debug: print "Loop 500: R5=> %s" % self.R5
                if self.debug: print ""
                #}}} I 

                # MATRIX-MULTIPLICATION - LAYERWISE ++++> FROM BOTTOM UPWARD TO LMAX
                #{{{ J 
                if self.debug: print "Loop 500: LMAX = >%s" % self.LMAX
                if self.debug: print "Loop 500: MMAX = >%s" % self.MMAX
                #for I in range(self.LMAX, self.MMAX):
                for I in range(self.LMAX, self.MMAX):
                    # MATRIX MULTIPLICATION FROM BOTTOM LAYER UPWARD
                    M = self.MMAX -1 + self.LMAX - I
                    self.DPTH = self.D[M]
                    if self.debug: print "Loop 500: I=> %s" % I
                    if self.debug: print "Loop 500: M=> %s" % M
                    if self.debug: print "Loop 500: DPTH=> %s" % self.DPTH
                    if self.debug: print "Loop 500: D=> %s" % self.D
                    if self.debug: print "Loop 500: IBLOCK=> %s" % self.IBLOCK
                    if self.debug: print "" 

                    if self.debug: print "IBLOCK=%s" % self.IBLOCK
                    #print "self.WVNO2 = %s " % self.WVNO2
                    for J in range(self.IBLOCK):
                        if self.debug: print "J=%s" % J
                        #print "self.WVNO2=%s" % self.WVNO2[J]
                        #{{{ K
                        self.P = self.RA[M][J] * complex(self.DPTH)
                        self.Q = self.RB[M][J] * complex(self.DPTH)
                        #if self.debug: print "Loop 500:   RA[%s][%s]=> %s" % (M,J,self.RA[M][J])
                        #if self.debug: print "Loop 500:   RB[%s][%s]=> %s" % (M,J,self.RB[M][J])
                        #if self.debug: print "Loop 500:   complex(DEPTH)=> %s" % complex(self.DEPTH)
                        if self.debug: print "Loop 500:   P=> %s" % self.P
                        if self.debug: print "Loop 500:   Q=> %s" % self.Q

                        self.A0 = 0.0
                        self.PR = self.P.real
                        self.PI = self.P.imag
                        self.QR = self.Q.real
                        self.QI = self.Q.imag

                        if self.debug: print "Loop 500:   A0=> %s" % self.A0
                        if self.debug: print "Loop 500:   PR=> %s" % self.PR
                        if self.debug: print "Loop 500:   PI=> %s" % self.PI
                        if self.debug: print "Loop 500:   QR=> %s" % self.QR
                        if self.debug: print "Loop 500:   QI=> %s" % self.QI

                        self.EPP = complex(cos(self.PI),sin(self.PI)) / 2.0
                        self.EPM = self.EPP.conjugate()
                        self.EQP = complex( cos(self.QI),sin(self.QI) ) / 2.0
                        self.EQM = self.EQP.conjugate()
                        self.EXB = self.QR
                        self.FAC = 0.0

                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   EPP=> %s" % self.EPP
                        if self.debug: print "Loop 500:   EPM=> %s" % self.EPM
                        if self.debug: print "Loop 500:   EQP=> %s" % self.EQP
                        if self.debug: print "Loop 500:   EQM=> %s" % self.EQM
                        if self.debug: print "Loop 500:   EXB=> %s" % self.EXB

                        if self.PR < 15.0: self.FAC = exp(-(2.0) * self.PR).real
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC


                        self.COSP = self.EPP + self.FAC * self.EPM
                        self.SINP = self.EPP - self.FAC * self.EPM
                        self.W = self.SINP / self.RA[M][J]
                        self.X = self.RA[M][J] * self.SINP
                        self.FAC = 0.0

                        if self.debug: print "Loop 500:   COSP=> %s" % self.COSP
                        if self.debug: print "Loop 500:   SINP=> %s" % self.SINP
                        if self.debug: print "Loop 500:   W=> %s" % self.W
                        if self.debug: print "Loop 500:   X=> %s" % self.X

                        if self.QR < 15.0: self.FAC = exp( -(2.0) * self.QR).real
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC

                        if self.debug: print ""

                        self.COSQL[M][J] = self.EQP + self.FAC * self.EQM
                        self.SINQ = self.EQP - self.FAC * self.EQM
                        self.YL[M][J] = self.SINQ / self.RB[M][J]
                        if self.debug: print "Loop 500:   M=> %s" % M
                        if self.debug: print "Loop 500:   J=> %s" % J
                        if self.debug: print "Loop 500:   COSQL[M][J]=> %s" % self.COSQL[M][J]
                        if self.debug: print "Loop 500:   SINQ=> %s" % self.SINQ
                        if self.debug: print "Loop 500:   YL[M][J]=> %s" % self.YL[M][J]
                        if self.debug: print ""
                        #print "Loop 500:   self.YL[%s][%s]=> %s" % (M,J,self.YL[M][J])
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
                        if self.debug: print "Loop 500:   ZL[M][J]=> %s" % self.ZL[M][J]
                        if self.debug: print "Loop 500:   EXA=> %s" % self.EXA
                        if self.debug: print "Loop 500:   CPCQ=> %s" % self.CPCQ
                        if self.debug: print "Loop 500:   CPY=> %s" % self.CPY
                        if self.debug: print "Loop 500:   CPZ=> %s" % self.CPZ
                        if self.debug: print "Loop 500:   CQW=> %s" % self.CQW
                        if self.debug: print "Loop 500:   CQX=> %s" % self.CQX
                        if self.debug: print "Loop 500:   XY=> %s" % self.XY
                        if self.debug: print "Loop 500:   XZ=> %s" % self.XZ
                        if self.debug: print "Loop 500:   WY=> %s" % self.WY
                        if self.debug: print "Loop 500:   WZ=> %s" % self.WZ
                        if self.debug: print ""
                        self.FAC = 0.0
                        self.QMP = self.QR - self.PR
                        if self.debug: print "Loop 500:   QMP=> %s" % self.QMP

                        if self.QMP > -40.0: self.FAC = exp(self.QMP)
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC

                        self.COSQ = self.COSQL[M][J] * self.FAC
                        if self.debug: print "Loop 500:   COSQ=> %s" % self.COSQ
                        self.Y = self.FAC * self.YL[M][J]
                        self.Z = self.FAC * self.ZL[M][J]
                        if self.debug: print "Loop 500:   Y=> %s" % self.Y
                        if self.debug: print "Loop 500:   Z=> %s" % self.Z
                        self.FAC = 0.0

                        if self.EXA < 60.0: self.A0 = exp(-self.EXA).real
                        if self.debug: print "Loop 500:   A0=> %s" % self.A0

                        self.RHO2  = self.RHO[M] * self.RHO[M]
                        if self.debug: print "Loop 500:   RHO2=> %s" % self.RHO2
                        if self.debug: print "Loop 500:   WVNO2=%s" % self.WVNO2[J]
                        if self.debug: print "Loop 500:   WVNO2=%s" % self.WVNO2
                        if self.debug: print "Loop 500:   J=%s" % J
                        self.A0C  = 2.0 * (self.A0 - self.CPCQ)
                        if self.debug: print "Loop 500:   A0C=> %s" % self.A0C
                        self.XZ2  = self.XZ / self.WVNO2[J]
                        if self.debug: print "Loop 500:   XZ2=> %s" % self.XZ2
                        self.WY2  = self.WY * self.WVNO2[J]
                        if self.debug: print "Loop 500:   WY2=> %s" % self.WY2
                        self.TEMP = self.A0C * self.WVNO2[J] + self.XZ + self.WY2 * self.WVNO2[J]
                        if self.debug: print "Loop 500:   TEMP=> %s" % self.TEMP

                        self.CA15 = -self.TEMP / self.RHO2
                        if self.debug: print "Loop 500:   CA15=> %s" % self.CA15

                        self.TEMP = 0.5 * self.A0C * ( self.GAM[M][J] + self.GAMM1[M][J] ) + self.GAM[M][J] * self.XZ2 + self.GAMM1[M][J] * self.WY2
                        if self.debug: print "Loop 500:   TEMP=> %s" % self.TEMP

                        if self.debug: print ""


                        self.CA13= -self.TEMP / self.RHO[M]
                        if self.debug: print "Loop 500:   CA13=> %s" % self.CA13
                        if self.debug: print "Loop 500:   GAMX1[M][J]=> %s" % self.GAMX1[M][J]
                        if self.debug: print "Loop 500:   GAM2[M][J]=> %s" % self.GAM2[M][J]
                        if self.debug: print "Loop 500:   GAMM2[M][J]=> %s" % self.GAMM2[M][J]

                        self.TEMP = self.A0C * self.GAMX1[M][J] + self.GAM2[M][J] * self.XZ2 + self.GAMM2[M][J] * self.WY2
                        if self.debug: print "Loop 500:   TEMP=> %s" % self.TEMP

                        self.CA33 = self.A0 + self.TEMP*2
                        self.CA11 = self.CPCQ - self.TEMP
                        if self.debug: print "Loop 500:   CA33=> %s" % self.CA33
                        if self.debug: print "Loop 500:   CA11=> %s" % self.CA11
                        #if self.debug: print "Loop 500:   CPCQ=> %s" % self.CPCQ

                        self.TEMP = 0.5 * self.A0C * self.GAMX1[M][J] * ( self.GAM[M][J] + self.GAMM1[M][J]) + self.GAM3[M][J] * self.XZ2 + self.GAMM3[M][J] * self.WY2
                        if self.debug: print "Loop 500:   TEMP=> %s" % self.TEMP

                        if self.debug: print ""
                        self.CA31 = 2.0 * self.TEMP * self.RHO[M]
                        self.TEMP = self.A0C * self.GAMX2[M][J] + self.GAM4[M][J] * self.XZ2 +self.GAMM3[M][J] * self.GAMM1[M][J] * self.WY2
                        if self.debug: print "Loop 500:   CA31=> %s" % self.CA31
                        if self.debug: print "Loop 500:   TEMP=> %s" % self.TEMP

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
                        if self.debug: print "Loop 500:   CA51=> %s" % self.CA51
                        if self.debug: print "Loop 500:   CA14=> %s" % self.CA14
                        if self.debug: print "Loop 500:   CA21=> %s" % self.CA21
                        if self.debug: print "Loop 500:   CA23=> %s" % self.CA23
                        if self.debug: print "Loop 500:   CA12=> %s" % self.CA12
                        if self.debug: print "Loop 500:   CA32=> %s" % self.CA32
                        if self.debug: print "Loop 500:   CA41=> %s" % self.CA41
                        if self.debug: print "Loop 500:   CA22=> %s" % self.CA22
                        if self.debug: print "Loop 500:   CA24=> %s" % self.CA24
                        if self.debug: print "Loop 500:   CA42=> %s" % self.CA42
                        if self.debug: print "Loop 500:   CA25=> %s" % self.CA25
                        if self.debug: print "Loop 500:   CA55=> %s" % self.CA55
                        if self.debug: print "Loop 500:   CA54=> %s" % self.CA54
                        if self.debug: print "Loop 500:   CA53=> %s" % self.CA53
                        if self.debug: print "Loop 500:   CA52=> %s" % self.CA52
                        if self.debug: print "Loop 500:   CA43=> %s" % self.CA43
                        if self.debug: print "Loop 500:   CA45=> %s" % self.CA45
                        if self.debug: print "Loop 500:   CA44=> %s" % self.CA44
                        if self.debug: print "Loop 500:   CA34=> %s" % self.CA34
                        if self.debug: print "Loop 500:   CA35=> %s" % self.CA35



                        self.EXE[J] = self.EXE[J] + self.EXA

                        self.ER1 = self.R1[J] * self.CA11 + self.R2[J] * self.CA21 + self.R3[J] * self.CA31 + self.R4[J] * self.CA41 + self.R5[J] * self.CA51
                        self.ER2 = self.R1[J] * self.CA12 + self.R2[J] * self.CA22 + self.R3[J] * self.CA32 + self.R4[J] * self.CA42 + self.R5[J] * self.CA52
                        self.ER3 = self.R1[J] * self.CA13 + self.R2[J] * self.CA23 + self.R3[J] * self.CA33 + self.R4[J] * self.CA43 + self.R5[J] * self.CA53
                        self.ER4 = self.R1[J] * self.CA14 + self.R2[J] * self.CA24 + self.R3[J] * self.CA34 + self.R4[J] * self.CA44 + self.R5[J] * self.CA54
                        self.ER5 = self.R1[J] * self.CA15 + self.R2[J] * self.CA25 + self.R3[J] * self.CA35 + self.R4[J] * self.CA45 + self.R5[J] * self.CA55

                        if self.debug: print ""
                        #print " TEST ER4 :  self.ER4      => %s" % self.ER4 
                        #print " TEST ER4 :  self.R1[J]    => %s" % self.R1[J] 
                        #print " TEST ER4 :  self.CA14     => %s" % self.CA14 
                        #print " TEST ER4 :  self.R2[J]    => %s" % self.R2[J] 
                        #print " TEST ER4 :  self.CA24     => %s" % self.CA24 
                        #print " TEST ER4 :  self.R3[J]    => %s" % self.R3[J] 
                        #print " TEST ER4 :  self.CA34     => %s" % self.CA34 
                        #print " TEST ER4 :  self.R4[J]    => %s" % self.R4[J] 
                        #print " TEST ER4 :  self.CA44     => %s" % self.CA44 
                        #print " TEST ER4 :  self.R5[J]    => %s" % self.R5[J] 
                        #print " TEST ER4 :  self.CA54     => %s" % self.CA54
                        #print ""



                        if self.debug: print "Loop 500:   ER1=> %s" % self.ER1
                        if self.debug: print "Loop 500:   ER2=> %s" % self.ER2
                        if self.debug: print "Loop 500:   ER3=> %s" % self.ER3
                        if self.debug: print "Loop 500:   ER4=> %s" % self.ER4
                        if self.debug: print "Loop 500:   ER5=> %s" % self.ER5
                        if self.debug: print ""

                        self.TESTT = 0.0
                        if abs(self.ER1.real) > self.TESTT: self.TESTT = abs(self.ER1.real)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER2.real) > self.TESTT: self.TESTT = abs(self.ER2.real)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER3.real) > self.TESTT: self.TESTT = abs(self.ER3.real)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER4.real) > self.TESTT: self.TESTT = abs(self.ER4.real)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER5.real) > self.TESTT: self.TESTT = abs(self.ER5.real)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER1.imag) > self.TESTT: self.TESTT = abs(self.ER1.imag)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER2.imag) > self.TESTT: self.TESTT = abs(self.ER2.imag)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER3.imag) > self.TESTT: self.TESTT = abs(self.ER3.imag)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER4.imag) > self.TESTT: self.TESTT = abs(self.ER4.imag)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT
                        if abs(self.ER5.imag) > self.TESTT: self.TESTT = abs(self.ER5.imag)
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT

                        if self.TESTT < 1.0E-30: self.TESTT = 1.0
                        if self.debug: print "Loop 500:   TESTT=> %s" % self.TESTT


                        self.TEST = 0.0

                        self.FAC = abs(self.ER1) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER2) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER3) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER4) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER5) / self.TESTT
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST
                        self.TEST = self.TEST * self.TESTT
                        if self.TEST < 1.0E-30: self.TEST = 1.0
                        if self.debug: print "Loop 500:   FAC=> %s" % self.FAC
                        if self.debug: print "Loop 500:   TEST=> %s" % self.TEST

                        if self.debug: print "self.TEST:%s" % self.TEST

                        self.XNORM = 1.0 / self.TEST
                        if self.debug: print "Loop 500:   XNORM=> %s" % self.XNORM
                        self.EXA  = -log(self.XNORM)
                        if self.debug: print "Loop 500:   EXA=> %s" % self.EXA
                        self.R1[J] = self.ER1 * self.XNORM
                        self.R2[J] = self.ER2 * self.XNORM
                        self.R3[J] = self.ER3 * self.XNORM
                        self.R4[J] = self.ER4 * self.XNORM
                        self.R5[J] = self.ER5 * self.XNORM
                        #print "Loop 500:        R1(%s)=> %s" % (J,self.R1[J])
                        #print "Loop 500:        R2(%s)=> %s" % (J,self.R2[J])
                        #print "Loop 500:        R3(%s)=> %s" % (J,self.R3[J])
                        #print "Loop 500:        R4(%s)=> %s" % (J,self.R4[J])
                        #print "Loop 500:        R5(%s)=> %s" % (J,self.R5[J])
                        #print "Loop 500:        R3(%s)=> %s" % (J,self.R3[J])
                        #print "Loop 500:        ER3=> %s" % self.ER3
                        #print "Loop 500:        XNORM=> %s" % self.XNORM

                        self.EXE[J] = self.EXE[J] + self.EXA
                        self.EXEL[J] = self.EXEL[J] + self.EXB
                        if self.debug: print "Loop 500:   EXE[%s]=> %s" % (J,self.EXE[J])
                        if self.debug: print "Loop 500:   EXEL[%s]=> %s" % (J,self.EXEL[J])
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
                if self.debug: print "Loop 500: X1=> %s" % self.X1
                if self.debug: print "Loop 500: X2=> %s" % self.X2
                if self.debug: print "Loop 500: X3=> %s" % self.X3
                if self.debug: print "Loop 500: X4=> %s" % self.X4
                if self.debug: print "Loop 500: X5=> %s" % self.X5
                if self.debug: print "Loop 500: EXL=> %s" % self.EXL
                if self.debug: print "Loop 500: EXLL=> %s" % self.EXLL
                if self.debug: print "Loop 500: ZX=> %s" % self.ZX
                #}}} L

                # MATRIX MULTIPLICATION FROM THE SOURCE INTERFACE UPWARD TO THE SURFACE
                #{{{ O
                for I in range(self.LMAX):
                    M = self.LMAX - I -1 
                    self.DPTH = self.D[M]
                    if self.debug: print "Loop 500 14: I=> %s" % I
                    if self.debug: print "Loop 500 14: LMAX=> %s" % self.LMAX
                    if self.debug: print "Loop 500 14: M=> %s" % M
                    if self.debug: print "Loop 500 14: DPTH=> %s" % self.DPTH

                    #{{{ P
                    for J in range(self.IBLOCK):
                        self.P = self.RA[M][J] * self.DPTH
                        self.Q = self.RB[M][J] * self.DPTH

                        if self.debug: print "Loop 500 15:   J=> %s" % J
                        if self.debug: print "Loop 500 15:   P=> %s" % self.P
                        if self.debug: print "Loop 500 15:   P=> %s" % self.Q

                        A0 = 0.0
                        self.PR = self.P.real
                        if self.debug: print "Loop 500 15: PR=> %s" % self.PR
                        self.PI = self.P.imag
                        if self.debug: print "Loop 500 15: PI=> %s" % self.PI
                        self.QR = self.Q.real
                        if self.debug: print "Loop 500 15: QR=> %s" % self.QR
                        self.QI = self.Q.imag
                        if self.debug: print "Loop 500 15: QI=> %s" % self.QI
                        self.EPP = complex(cos(self.PI),sin(self.PI)) / 2.0
                        if self.debug: print "Loop 500 15: EPP=> %s" % self.EPP
                        self.EPM = self.EPP.conjugate()
                        if self.debug: print "Loop 500 15: EPM=> %s" % self.EPM
                        self.EQP = complex(cos(self.QI),sin(self.QI)) / 2.0
                        if self.debug: print "Loop 500 15: EQP=> %s" % self.EQP
                        self.EQM = self.EQP.conjugate()
                        if self.debug: print "Loop 500 15: EQM=> %s" % self.EQM
                        self.EX = self.PR
                        if self.debug: print "Loop 500 15: EX=> %s" % self.EX
                        self.EXB = self.QR
                        if self.debug: print "Loop 500 15: EXB=> %s" % self.EXB
                        self.FAC = 0.0
                        if self.debug: print ""

                        if self.PR < 15.0: self.FAC = exp(-2*self.PR).real 
                        if self.debug: print "Loop 500 15: FAC=> %s" % self.FAC

                        self.COSP = self.EPP + self.FAC * self.EPM
                        if self.debug: print "Loop 500 15: COSP=> %s" % self.COSP
                        self.SINP = self.EPP - self.FAC * self.EPM
                        if self.debug: print "Loop 500 15: SINP=> %s" % self.SINP
                        self.W = self.SINP / self.RA[M][J]
                        if self.debug: print "Loop 500 15: W=> %s" % self.W
                        self.X = self.RA[M][J] * self.SINP
                        if self.debug: print "Loop 500 15: X=> %s" % self.X
                        self.FAC = 0.0
                        if self.QR < 15.0: self.FAC = exp(-2.0 * self.QR)
                        if self.debug: print "Loop 500 15: FAC=> %s" % self.FAC
                        self.COSQL[M][J] = self.EQP + self.FAC * self.EQM
                        if self.debug: print "Loop 500 15: COSQL[M][J]=> %s" % self.COSQL[M][J]
                        self.SINQ = self.EQP - self.FAC * self.EQM
                        if self.debug: print "Loop 500 15: SINQ=> %s" % self.SINQ
                        self.YL[M][J] = self.SINQ / self.RB[M][J]
                        if self.debug: print "Loop 500 15: YL[M][J]=> %s" % self.YL[M][J]
                        self.ZL[M][J] = self.RB[M][J] * self.SINQ
                        if self.debug: print "Loop 500 15: ZL[M][J]=> %s" % self.ZL[M][J]
                        self.EXA = self.PR + self.QR
                        if self.debug: print "Loop 500 15: EXA=> %s" % self.EXA
                        self.CPCQ = self.COSP * self.COSQL[M][J]
                        if self.debug: print "Loop 500 15: CPCQ=> %s" % self.CPCQ
                        self.CPY = self.COSP * self.YL[M][J]
                        if self.debug: print "Loop 500 15: CPY=> %s" % self.CPY
                        self.CPZ = self.COSP * self.ZL[M][J]
                        if self.debug: print "Loop 500 15: CPZ=> %s" % self.CPZ
                        self.CQW = self.COSQL[M][J] * self.W
                        if self.debug: print "Loop 500 15: CQW=> %s" % self.CQW
                        self.CQX = self.COSQL[M][J] * self.X
                        if self.debug: print "Loop 500 15: CQX=> %s" % self.CQX
                        self.XY = self.X * self.YL[M][J]
                        if self.debug: print "Loop 500 15: XY=> %s" % self.XY
                        self.XZ = self.X * self.ZL[M][J]
                        if self.debug: print "Loop 500 15: XZ=> %s" % self.XZ
                        self.WY = self.W * self.YL[M][J]
                        if self.debug: print "Loop 500 15: WY=> %s" % self.WY
                        self.WZ = self.W * self.ZL[M][J]
                        if self.debug: print "Loop 500 15: WZ=> %s" % self.WZ
                        if self.debug: print ""

                        self.FAC = 0.0
                        self.QMP = self.QR - self.PR
                        if self.debug: print "Loop 500 15: QMP=> %s" % self.QMP
                        if self.QMP > -40.0: self.FAC = exp(self.QMP).real
                        if self.debug: print "Loop 500 15: FAC=> %s" % self.FAC
                        self.COSQ = self.COSQL[M][J] * self.FAC
                        if self.debug: print "Loop 500 15: COSQ=> %s" % self.COSQ
                        self.Y = self.FAC * self.YL[M][J]
                        if self.debug: print "Loop 500 15: Y=> %s" % self.Y
                        self.Z = self.FAC * self.ZL[M][J]
                        if self.debug: print "Loop 500 15: Z=> %s" % self.Z

                        if self.debug: print ""

                        self.FAC=0.0
                        if self.EXA < 60.0: self.A0 = exp(-self.EXA)
                        self.RHO2 = self.RHO[M] ** 2 
                        self.A0C  = 2.0 * (self.A0 - self.CPCQ)
                        self.XZ2  = self.XZ / self.WVNO2[J]
                        self.WY2  = self.WY * self.WVNO2[J]
                        if self.debug: print "Loop 500 15: A0=> %s" % self.A0
                        if self.debug: print "Loop 500 15: RHO2=> %s" % self.RHO2
                        if self.debug: print "Loop 500 15: A0C=> %s" % self.A0C
                        if self.debug: print "Loop 500 15: XZ2=> %s" % self.XZ2
                        if self.debug: print "Loop 500 15: WY2=> %s" % self.WY2
                        self.TEMP = self.A0C * self.WVNO2[J] + self.XZ + self.WY2 * self.WVNO2[J]
                        if self.debug: print "Loop 500 15: TEMP=> %s" % self.TEMP
                        self.CA15 = -self.TEMP / self.RHO2
                        if self.debug: print "Loop 500 15: CA15=> %s" % self.CA15
                        self.TEMP = 0.5 * self.A0C * (self.GAM[M][J] + self.GAMM1[M][J]) + self.GAM[M][J] * self.XZ2 + self.GAMM1[M][J] * self.WY2
                        if self.debug: print "Loop 500 15: TEMP=> %s" % self.TEMP
                        self.CA13 = -self.TEMP / self.RHO[M]
                        if self.debug: print "Loop 500 15: CA13=> %s" % self.CA13
                        self.TEMP = self.A0C * self.GAMX1[M][J] + self.GAM2[M][J] * self.XZ2 + self.GAMM2[M][J] * self.WY2
                        if self.debug: print "Loop 500 15: TEMP=> %s" % self.TEMP
                        self.CA33 = self.A0 + self.TEMP + self.TEMP
                        if self.debug: print "Loop 500: CA33=> %s" % self.CA33
                        self.CA11 = self.CPCQ - self.TEMP
                        if self.debug: print "Loop 500 15: CA11=> %s" % self.CA11
                        self.TEMP = 0.5 * self.A0C * self.GAMX1[M][J] * (self.GAM[M][J] + self.GAMM1[M][J])+ self.GAM3[M][J] * self.XZ2 + self.GAMM3[M][J] * self.WY2
                        if self.debug: print "Loop 500 15: TEMP=> %s" % self.TEMP
                        self.CA31 = 2.0 * self.TEMP * self.RHO[M]
                        if self.debug: print "Loop 500 15: CA31=> %s" % self.CA31
                        self.TEMP = self.A0C * self.GAMX2[M][J] + self.GAM4[M][J] * self.XZ2 + self.GAMM3[M][J] * self.GAMM1[M][J] * self.WY2
                        #if self.debug: print "Loop 500 15: GAMX2=> %s" % self.GAMX2[M][J]
                        #if self.debug: print "Loop 500 15: GAM4=> %s" % self.GAM4[M][J]
                        #if self.debug: print "Loop 500 15: XZ2=> %s" % self.XZ2
                        #if self.debug: print "Loop 500 15: GAMM3=> %s" % self.GAMM3[M][J]
                        #if self.debug: print "Loop 500 15: GAMM1=> %s" % self.GAMM1[M][J]
                        #if self.debug: print "Loop 500 15: WY2=> %s" % self.WY2
                        if self.debug: print "Loop 500 15: TEMP=> %s" % self.TEMP
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

                        if self.debug: print "Loop 500 15: CA51=> %s" % self.CA51
                        if self.debug: print "Loop 500 15: CA14=> %s" % self.CA14
                        if self.debug: print "Loop 500 15: CA21=> %s" % self.CA21
                        if self.debug: print "Loop 500 15: CA23=> %s" % self.CA23
                        if self.debug: print "Loop 500 15: CA12=> %s" % self.CA12
                        if self.debug: print "Loop 500 15: CA32=> %s" % self.CA32
                        if self.debug: print "Loop 500 15: CA41=> %s" % self.CA41
                        if self.debug: print "Loop 500 15: CA22=> %s" % self.CA22
                        if self.debug: print "Loop 500 15: CA24=> %s" % self.CA24
                        if self.debug: print "Loop 500 15: CA42=> %s" % self.CA42
                        if self.debug: print "Loop 500 15: CA25=> %s" % self.CA25
                        if self.debug: print "Loop 500 15: CA55=> %s" % self.CA55
                        if self.debug: print "Loop 500 15: CA54=> %s" % self.CA54
                        if self.debug: print "Loop 500 15: CA53=> %s" % self.CA53
                        if self.debug: print "Loop 500 15: CA52=> %s" % self.CA52
                        if self.debug: print "Loop 500 15: CA43=> %s" % self.CA43
                        if self.debug: print "Loop 500 15: CA45=> %s" % self.CA45
                        if self.debug: print "Loop 500 15: CA44=> %s" % self.CA44
                        if self.debug: print "Loop 500 15: CA34=> %s" % self.CA34
                        if self.debug: print "Loop 500 15: CA35=> %s" % self.CA35

                        if self.debug: print ""

                        self.EXE[J] = self.EXE[J] + self.EXA
                        if self.debug: print "Loop 500: EXE[J]=> %s" % self.EXE[J]
                        self.ER1 = self.R1[J] * self.CA11 + self.R2[J] * self.CA21 + self.R3[J] * self.CA31 + self.R4[J] * self.CA41 + self.R5[J] * self.CA51
                        self.ER2 = self.R1[J] * self.CA12 + self.R2[J] * self.CA22 + self.R3[J] * self.CA32 + self.R4[J] * self.CA42 + self.R5[J] * self.CA52
                        self.ER3 = self.R1[J] * self.CA13 + self.R2[J] * self.CA23 + self.R3[J] * self.CA33 + self.R4[J] * self.CA43 + self.R5[J] * self.CA53
                        self.ER4 = self.R1[J] * self.CA14 + self.R2[J] * self.CA24 + self.R3[J] * self.CA34 + self.R4[J] * self.CA44 + self.R5[J] * self.CA54
                        self.ER5 = self.R1[J] * self.CA15 + self.R2[J] * self.CA25 + self.R3[J] * self.CA35 + self.R4[J] * self.CA45 + self.R5[J] * self.CA55
                        if self.debug: print "Loop 500: ER1=> %s" % self.ER1
                        if self.debug: print "Loop 500: ER2=> %s" % self.ER2
                        if self.debug: print "Loop 500: ER3=> %s" % self.ER3
                        if self.debug: print "Loop 500: ER4=> %s" % self.ER4
                        if self.debug: print "Loop 500: ER5=> %s" % self.ER5

                        self.TESTT = 0.0
                        if abs(self.ER1.real) > self.TESTT: self.TESTT = abs(self.ER1.real)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER2.real) > self.TESTT: self.TESTT = abs(self.ER2.real)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER3.real) > self.TESTT: self.TESTT = abs(self.ER3.real)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER4.real) > self.TESTT: self.TESTT = abs(self.ER4.real)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER5.real) > self.TESTT: self.TESTT = abs(self.ER5.real)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT

                        if abs(self.ER1.imag) > self.TESTT: self.TESTT = abs(self.ER1.imag)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER2.imag) > self.TESTT: self.TESTT = abs(self.ER2.imag)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER3.imag) > self.TESTT: self.TESTT = abs(self.ER3.imag)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER4.imag) > self.TESTT: self.TESTT = abs(self.ER4.imag)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if abs(self.ER5.imag) > self.TESTT: self.TESTT = abs(self.ER5.imag)
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT
                        if self.TESTT < 1.0e-30: self.TESTT = 1.0
                        if self.debug: print "Loop 500: TESTT=> %s" % self.TESTT

                        if self.debug: print ""


                        self.TEST = 0.0
                        self.FAC = abs(self.ER1) / self.TESTT
                        if self.debug: print "Loop 500: FAC=> %s" % self.FAC
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER2) / self.TESTT
                        if self.debug: print "Loop 500: FAC=> %s" % self.FAC
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER3) / self.TESTT
                        if self.debug: print "Loop 500: FAC=> %s" % self.FAC
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER4) / self.TESTT
                        if self.debug: print "Loop 500: FAC=> %s" % self.FAC
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.FAC = abs(self.ER5) / self.TESTT
                        if self.debug: print "Loop 500: FAC=> %s" % self.FAC
                        if self.TEST < self.FAC: self.TEST = self.FAC
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.TEST = self.TEST * self.TESTT
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        if self.TEST < 1.0e-30: self.TEST = 1.0
                        if self.debug: print "Loop 500: TEST=> %s" % self.TEST
                        self.XNORM = 1.0 / self.TEST
                        if self.debug: print "Loop 500: XNORM=> %s" % self.XNORM
                        self.EXA = -log(self.XNORM)
                        if self.debug: print "Loop 500: EXA=> %s" % self.EXA

                        self.R1[J] = self.ER1 * self.XNORM
                        self.R2[J] = self.ER2 * self.XNORM
                        self.R3[J] = self.ER3 * self.XNORM
                        self.R4[J] = self.ER4 * self.XNORM
                        self.R5[J] = self.ER5 * self.XNORM
                        if self.debug: print "Loop 500: R1[J]=> %s" % self.R1[J]
                        if self.debug: print "Loop 500: R2[J]=> %s" % self.R2[J]
                        if self.debug: print "Loop 500: R3[J]=> %s" % self.R3[J]
                        if self.debug: print "Loop 500: R4[J]=> %s" % self.R4[J]
                        if self.debug: print "Loop 500: R5[J]=> %s" % self.R5[J]


                        self.EXE[J] = self.EXE[J] + self.EXA
                        self.EXEL[J] = self.EXEL[J] + self.EXB

                        if self.debug: print "Loop 500: EXE[J]=> %s" % self.EXE[J]
                        if self.debug: print "Loop 500: EXEL[J]=> %s" % self.EXEL[J]

                        #
                        #   self.  PROPAGATION OF 4X4 HASKEL MATRIX UPWARD ++++++++++++>
                        #
                        self.EXL[J] = self.EXL[J] + self.EX
                        if self.debug: print "Loop 500: EXL=> %s" % self.EXL[J]
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

                        if self.debug: print ""
                        if self.debug: print "Loop 500: J=> %s" % J
                        if self.debug: print "Loop 500: EXL[J]=> %s" % self.EXL[J]
                        if self.debug: print "Loop 500: AA11=> %s" % self.AA11
                        if self.debug: print "Loop 500: AA12=> %s" % self.AA12
                        if self.debug: print "Loop 500: AA13=> %s" % self.AA13
                        if self.debug: print "Loop 500: AA14=> %s" % self.AA14
                        if self.debug: print "Loop 500: AA21=> %s" % self.AA21
                        if self.debug: print "Loop 500: AA22=> %s" % self.AA22
                        if self.debug: print "Loop 500: AA23=> %s" % self.AA23
                        if self.debug: print "Loop 500: AA24=> %s" % self.AA24
                        if self.debug: print "Loop 500: AA31=> %s" % self.AA31
                        if self.debug: print "Loop 500: AA32=> %s" % self.AA32
                        if self.debug: print "Loop 500: AA34=> %s" % self.AA34
                        if self.debug: print "Loop 500: AA33=> %s" % self.AA33
                        if self.debug: print "Loop 500: AA41=> %s" % self.AA41
                        if self.debug: print "Loop 500: AA42=> %s" % self.AA42
                        if self.debug: print "Loop 500: AA43=> %s" % self.AA43
                        if self.debug: print "Loop 500: AA44=> %s" % self.AA44
                        #print "%s %s %s %s %s %s %s %s" %(self.ZX[J][0][0],self.AA11,self.ZX[J][0][1],self.AA21,self.ZX[J][0][2],self.AA31,self.ZX[J][0][3],self.AA41)
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

                        if self.debug: print ""
                        if self.debug: print "Loop 500: D11=> %s" % self.D11
                        if self.debug: print "Loop 500: D12=> %s" % self.D12
                        if self.debug: print "Loop 500: D13=> %s" % self.D13
                        if self.debug: print "Loop 500: D14=> %s" % self.D14
                        if self.debug: print "Loop 500: D21=> %s" % self.D21
                        if self.debug: print "Loop 500: D22=> %s" % self.D22
                        if self.debug: print "Loop 500: D23=> %s" % self.D23
                        if self.debug: print "Loop 500: D24=> %s" % self.D24
                        if self.debug: print "Loop 500: D31=> %s" % self.D31
                        if self.debug: print "Loop 500: D32=> %s" % self.D32
                        if self.debug: print "Loop 500: D33=> %s" % self.D33
                        if self.debug: print "Loop 500: D34=> %s" % self.D34
                        if self.debug: print "Loop 500: D41=> %s" % self.D41
                        if self.debug: print "Loop 500: D42=> %s" % self.D42
                        if self.debug: print "Loop 500: D43=> %s" % self.D43
                        if self.debug: print "Loop 500: D44=> %s" % self.D44
                        if self.debug: print ""

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
                #        print ' [%s][%s] %s ' % (J,L,[self.ZX[J][L][W] for W in self.ZX[J][L]])


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
                        #print " TEST  X1[K]=> %s" % self.X1[K]
                        #print " TEST  ZX[K,0,0]=> %s" % self.ZX[K][0][0]
                        #print " TEST  X3[K]=> %s" % self.X3[K]
                        #print " TEST  ZX[K,2,0]=> %s" % self.ZX[K][2][0]
                        #print " TEST  X4[K]=> %s" % self.X4[K]
                        #print " TEST  ZX[K,3,0]=> %s" % self.ZX[K][3][0]
                        #print " TEST  Y21[K]=> %s" % self.Y21[K]

                #exit()
                #print ""
                #print "Loop 500 16 17 : Y11=> %s" % self.Y11
                #print "Loop 500 16 17 : Y12=> %s" % self.Y12
                #print "Loop 500 16 17 : Y31=> %s" % self.Y31
                #print "Loop 500 16 17 : Y32=> %s" % self.Y32
                #print "Loop 500 16 17 : Y21=> %s" % self.Y21
                #print "Loop 500 16 17 : Y22=> %s" % self.Y22
                #print "Loop 500 16 17 : Y41=> %s" % self.Y41
                #print "Loop 500 16 17 : Y42=> %s" % self.Y42
                #sys.exit()
                #}}} Q

            #}}} BR

            if self.debug: print "Loop 500 section CR => self.ISYST:%s" % self.ISYST
            if self.ISYST != 1:
            #{{{ CR
                if self.debug: print "\t ##### 300 #####"

                #
                # SH-MOTION STARTS HERE +++++>
                #
                for J in range(self.IBLOCK):
                    self.S32E[J] = 2.0 * self.XKA2 * self.WVNO[J] / self.FRHO
                    self.S34E[J] = 4.0 * self.WVNO[J] * self.XKA2 / (self.FOURPI * self.XKB2)
                    if self.debug: print ""
                    if self.debug: print "Loop 500 300 18:   S32E[J]=> %s" % self.S32E[J]
                    if self.debug: print "Loop 500 300 18:   S34E[J]=> %s" % self.S34E[J]

                if self.JBDRY == 0:
                    for J in range(self.IBLOCK):
                        #if self.debug: print "\nJ:%s \nMMAX:%s \nRB:%s \nRHO:%s" % (J,self.MMAX,self.RB,self.RHO)
                        self.R6[J] = (self.RHO[self.MMAX]) * self.RB[self.MMAX][J]
                        self.R7[J] = self.ONE / (self.B[self.MMAX] * self.ATNB[self.MMAX]) ** 2
                        if self.debug: print ""
                        if self.debug: print "Loop 500 300 19: R6[J]=> %s" % self.R6[J]
                        if self.debug: print "Loop 500 300 19: R7[J]=> %s" % self.R7[J]

                elif self.JBDRY < 0:
                    for J in range(self.IBLOCK):
                        self.R6[J] = self.ONE
                        self.R7[J] = self.ZERO
                        if self.debug: print ""
                        if self.debug: print "Loop 500 300 20: R6[J]=> %s" % self.R6[J]
                        if self.debug: print "Loop 500 300 20: R7[J]=> %s" % self.R7[J]

                elif self.JBDRY > 0:
                    for J in range(self.IBLOCK):
                        self.R6[J] = self.ZERO
                        self.R7[J] = self.ONE
                        if self.debug: print ""
                        if self.debug: print "Loop 500 300 21: R6[J]=> %s" % self.R6[J]
                        if self.debug: print "Loop 500 300 21: R7[J]=> %s" % self.R7[J]
            #}}} CR


            if self.debug: print "Loop 500 section DR => self.ISYST:%s" % self.ISYST
            if self.ISYST != 3 and self.ISYST != 1:
            #{{{ DR
                for M in range(self.MMAX+1):
                    for J in range(self.IBLOCK):
                        #print "Loop 500:   self.YL=> %s" % self.YL
                        #print "Loop 500:   [%s][%s]=> %s" % (M,J)
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
                        #print "Loop 500 302 22:   P => %s" % self.P
                        #print "Loop 500 302 22:   Q => %s" % self.Q
                        #print "Loop 500 302 22:   A0 => %s" % self.A0
                        #print "Loop 500 302 22:   PR => %s" % self.PR
                        #print "Loop 500 302 22:   PI => %s" % self.PI
                        #print "Loop 500 302 22:   QR => %s" % self.QR
                        #print "Loop 500 302 22:   QI => %s" % self.QI
                        #print "Loop 500 302 22:   EPP => %s" % self.EPP
                        #print "Loop 500 302 22:   EPM => %s" % self.EPM
                        #print "Loop 500 302 22:   EQP => %s" % self.EQP
                        #print "Loop 500 302 22:   EQM => %s" % self.EQM
                        #print "Loop 500 302 22:   EXEL[J] => %s" % self.EXEL[J]
                        #print "Loop 500 302 22:   EXLL[J] => %s" % self.EXLL[J]
                        #print "Loop 500 302 22:   COSP => %s" % self.COSP
                        #print "Loop 500 302 22:   SINP => %s" % self.SINP
                        #print "Loop 500 302 22:   COSQL[M][J] => %s" % self.COSQL[M][J]
                        #print "Loop 500 302 22:   SINQ => %s" % self.SINQ
                        #print "Loop 500 302 22:   YL[M][J] => %s" % self.YL[M][J]
                        #print "Loop 500 302 22:   ZL[M][J] => %s" % self.ZL[M][J]

                        sys.exit()

            #}}} DR

            if self.debug: print "Loop 500 section ER => self.ISYST:%s" % self.ISYST
            if self.ISYST != 1:
            #{{{ ER
                if self.debug: print " ##### AFTER 302 ####"
                if self.debug: print "Loop 500:   LMAX=> %s" % self.LMAX
                if self.debug: print "Loop 500:   MMAX=> %s" % self.MMAX

                for I in range(self.LMAX, self.MMAX):

                    if self.debug: print "Loop 500:   I => %s" % I
                    M = self.MMAX + self.LMAX - I -1
                    self.HH = self.RHO[M] * (self.B[M] * self.ATNB[M])**2
                    if self.debug: print "Loop 500:   M => %s" % M
                    if self.debug: print "Loop 500:   self.HH => %s" % self.HH
                    for J in range(self.IBLOCK):
                        self.YLI = self.YL[M][J] / self.HH
                        self.ZLI = self.ZL[M][J] * self.HH
                        self.D11 = self.R6[J]
                        self.D12 = self.R7[J]
                        self.R6[J] = self.D11 * self.COSQL[M][J] + self.D12 * self.ZLI
                        self.R7[J] = self.D11 * self.YLI + self.D12 * self.COSQL[M][J]
                        if self.debug: print "Loop 500:     J => %s" % J
                        if self.debug: print "Loop 500:     YLI=> %s" % self.YLI
                        if self.debug: print "Loop 500:     ZLI=> %s" % self.ZLI
                        if self.debug: print "Loop 500:     D11=> %s" % self.D11
                        if self.debug: print "Loop 500:     D12=> %s" % self.D12
                        if self.debug: print "Loop 500:     R6[J]=> %s" % self.R6[J]
                        if self.debug: print "Loop 500:     R7[J]=> %s" % self.R7[J]


                if self.debug: print ""
                for J in range(self.IBLOCK):
                    self.X6[J] = self.R6[J]
                    self.X7[J] = self.R7[J]
                    if self.debug: print "Loop 500:     J => %s" % J
                    if self.debug: print "Loop 500:     X6[J]=> %s" % self.X6[J]
                    if self.debug: print "Loop 500:     X7[J]=> %s" % self.X7[J]

                if self.debug: print ""
                for I in range(self.LMAX):
                #for M in range(self.LMAX,-1,-1):
                    M = self.LMAX - I - 1
                    if self.debug: print "Loop 500:   I => %s" % I
                    if self.debug: print "Loop 500:   M => %s" % M
                    self.HH = self.RHO[M] * (self.B[M] * self.ATNB[M])**2
                    if self.debug: print "Loop 500:   self.HH => %s" % self.HH
                    for J in range(self.IBLOCK):
                        self.YLI = self.YL[M][J] / self.HH
                        self.ZLI = self.ZL[M][J] * self.HH
                        self.D11 = self.R6[J]
                        self.D12 = self.R7[J]
                        self.R6[J] = self.D11 * self.COSQL[M][J] + self.D12 * self.ZLI
                        self.R7[J] = self.D11 * self.YLI + self.D12 * self.COSQL[M][J]
                        #print "Loop 500:        J=> %s" % J
                        #print "Loop 500:        YLI=> %s" % self.YLI
                        #print "Loop 500:        ZLI=> %s" % self.ZLI
                        #print "Loop 500:        D11=> %s" % self.D11
                        #print "Loop 500:        D12=> %s" % self.D12
                        #print "Loop 500:        R6[J]=> %s" % self.R6[J]
                        #print "Loop 500:        R7[J]=> %s" % self.R7[J]
                #print "Loop 500:   R6=> %s" % self.R6
                #if self.debug: print "Loop 500:   M=> %s" % M
                #if self.debug: print "Loop 500:   HH=> %s" % self.HH
                #if self.debug: print "Loop 500:   YLI=> %s" % self.YLI
                #if self.debug: print "Loop 500:   ZLI=> %s" % self.ZLI
                #if self.debug: print "Loop 500:   D11=> %s" % self.D11
                #if self.debug: print "Loop 500:   D12=> %s" % self.D12
                #if self.debug: print "Loop 500:   R6=> %s" % self.R6
                #if self.debug: print "Loop 500:   R7=> %s" % self.R7
            #}}} ER

            if self.debug: print "Loop 500 section FR => self.ISYST:%s" % self.ISYST

            if self.debug: print ""
            if self.debug: print "######## AFTER 301 #######"
            #{{{ FR
            if self.debug: print "Loop 500 28: IBLOCK => %s" % self.IBLOCK
            for J in range(self.IBLOCK):
                if self.debug: print "Loop 500 28: J => %s" % J
                self.FACT[J] = 0.0
                self.ELJ = self.EXE[J] - self.EXL[J]
                if self.ELJ < 55.0: self.FACT[J] = exp(-self.ELJ).real
                self.FACT0 = 0.0
                self.ELJ = self.EXLL[J] - self.EXEL[J]
                if self.SNGL[self.ELJ] > -40.0: self.FACT0 = exp(self.ELJ).real
                self.FACX[J] = self.FACT0 * self.FACXX
                if self.debug: print "Loop 500:   FACT[J]=> %s" % self.FACT[J]
                if self.debug: print "Loop 500:   ELJ=> %s" % self.ELJ
                if self.debug: print "Loop 500:   FACT0=> %s" % self.FACT0
                if self.debug: print "Loop 500:   FACX[J]=> %s" % self.FACX[J]

            if self.debug: print ""

            if self.ISRC[0] == 1:
                if self.debug: print ""
                if self.debug: print "\t29"
                for J in range(self.IBLOCK):
                    self.G = self.S32[J] * self.Y21[J] + self.S34[J] * self.Y41[J]
                    self.G1[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t29:     G = > %s" % self.G
                    if self.debug: print "\t29:     G1[J] = > %s" % self.G1[J]

            if self.ISRC[1] == 1:
                if self.debug: print ""
                if self.debug: print "\t30"
                for J in range(self.IBLOCK):
                    self.G = self.S32[J] * self.Y22[J] + self.S34[J] * self.Y42[J]
                    self.G2[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t30:     G = > %s" % self.G
                    if self.debug: print "\t30:     G2[J] = > %s" % self.G2[J]

            if self.ISRC[2] == 1:
                if self.debug: print ""
                if self.debug: print "\t31"
                for J in range(self.IBLOCK):
                    self.G = self.S21[J] * self.Y11[J]
                    self.G3[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t31:     G = > %s" % self.G
                    if self.debug: print "\t31:     G3[J] = > %s" % self.G3[J]

            if self.ISRC[3] == 1:
                if self.debug: print ""
                if self.debug: print "\t32"
                for J in range(self.IBLOCK):
                    self.G = self.S21[J] * self.Y12[J]
                    self.G4[J]= -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t32:     G = > %s" % self.G
                    if self.debug: print "\t32:     G4[J] = > %s" % self.G4[J]

            if self.ISRC[4] == 1:
                if self.debug: print ""
                if self.debug: print "\t33"
                for J in range(self.IBLOCK):
                    self.G = self.S14[J] * self.Y41[J]
                    self.G5[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t33:     G = > %s" % self.G
                    if self.debug: print "\t33:     G5[J] = > %s" % self.G5[J]

            if self.ISRC[5] == 1:
                if self.debug: print ""
                if self.debug: print "\t34"
                for J in range(self.IBLOCK):
                    self.G = self.S14[J] * self.Y42[J]
                    self.G6[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t34:     G = > %s" % self.G
                    if self.debug: print "\t34:     G6[J] = > %s" % self.G6[J]

            if self.ISRC[6] == 1:
                if self.debug: print ""
                if self.debug: print "\t35"
                for J in range(self.IBLOCK):
                    self.G = self.S32E[J] * self.Y21[J] + self.S34E[J] * self.Y41[J]
                    self.G7[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t35:     Y41[J] = > %s" % self.Y41[J]
                    if self.debug: print "\t35:     S34E[J] = > %s" % self.S34E[J]
                    if self.debug: print "\t35:     Y21[J] = > %s" % self.Y21[J]
                    if self.debug: print "\t35:     S32E[J] = > %s" % self.S32E[J]
                    if self.debug: print "\t35:     G = > %s" % self.G
                    if self.debug: print "\t35:     FACT[J] = > %s" % self.FACT[J]
                    if self.debug: print "\t35:     R1[J] = > %s" % self.R1[J]
                    if self.debug: print "\t35:     G7[J] = > %s" % self.G7[J]
            #exit()

            if self.ISRC[7] == 1:
                if self.debug: print ""
                if self.debug: print "\t36"
                for J in range(self.IBLOCK):
                    self.G = self.S32E[J] * self.Y22[J] + self.S34E[J] * self.Y42[J]
                    self.G8[J] = -self.G * self.FACT[J] / self.R1[J]
                    if self.debug: print "\t36:     G = > %s" % self.G
                    if self.debug: print "\t36:     FACT[J] = > %s" % self.FACT[J]
                    if self.debug: print "\t36:     R1[J] = > %s" % self.R1[J]
                    if self.debug: print "\t36:     G8[J] = > %s" % self.G8[J]

            if self.ISRC[8] == 1:
                if self.debug: print ""
                if self.debug: print "\t37"
                for J in range(self.IBLOCK):
                    self.G = 2.0 * self.X6[J] / self.RHO[self.LMAX]
                    self.G9[J] = self.G * self.FACX[J] / ( self.R6[J] * self.XKK )
                    if self.debug: print "\t37:     G = > %s" % self.G
                    if self.debug: print "\t37:     G9[J] = > %s" % self.G9[J]


            if self.ISRC[9] == 1:
                for J in range(self.IBLOCK):
                    self.G = -2.0 * self.WVNO[J] * self.X7[J] * (self.B[self.LMAX] * self.ATNB[self.LMAX])**2
                    self.G10[J] = self.G * self.FACX[J] / (self.R6[J] * self.XKK)

            self.N11 = (self.NBLK - 1) * self.IB + 1
            self.N22 = self.N11 + self.IBLOCK - 1

            #print "Loop 500:   G=> %s" % self.G
            #print "Loop 500:   G1=> %s" % self.G1
            #print "Loop 500:   G2=> %s" % self.G2
            #print "Loop 500:   G3=> %s" % self.G3
            #print "Loop 500:   G4=> %s" % self.G4
            #print "Loop 500:   G5=> %s" % self.G5
            #print "Loop 500:   G6=> %s" % self.G6
            #print "Loop 500:   G7=> %s" % self.G7
            #print "Loop 500:   G8=> %s" % self.G8
            #print "Loop 500:   G9=> %s" % self.G9
            #print "Loop 500:   G10=> %s" % self.G10
            #print "Loop 500:   N11=> %s" % self.N11
            #print "Loop 500:   N22=> %s" % self.N22
            #exit()

            if self.debug: print "self.NRANGE: %s" % self.NRANGE
            for I in range(self.NRANGE):
                if self.debug: print "loop WVINT I: %s" % I
                self.R = self.RANGE[I]
                self.ICN = self.ICNT[I]
                if self.debug: print "Loop 600 [ sub WVINT ]:   R=> %s" % self.R
                if self.debug: print "Loop 600 [ sub WVINT ]:   ICN=> %s" % self.ICN

                if self.debug: print "Loop 600:   WVNO=> %s" % self.WVNO
                if self.debug: print "Loop 600:   R=> %s" % self.RANGE[I]
                if self.debug: print "Loop 600:   NBLK=> %s" % self.NBLOCK
                if self.debug: print "Loop 600:   IS=> %s" % self.NBLOCK
                if self.debug: print "Loop 600:   IBLOCK=> %s" % self.IBLOCK
                if self.debug: print "Loop 600:   ICN=> %s" % self.ICNT[I]
                if self.debug: print "Loop 600:   ISRC=> %s" % self.ISRC
                if self.debug: print "Loop 600:   N11=> %s" % self.N11
                if self.debug: print "Loop 600:   N22=> %s" % self.N22
                if self.debug: print "Loop 600:   I=> %s" % I

                if self.debug: print ""

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
                self.ZERO = complex(0.0,0.0)
                self.PI = acos(-1.0)
                self.NN1 = 1

                if self.debug: print ""
                if self.debug: print "++++++++++++++++++++ WVINT ++++++++++++++++++"
                for J in range(self.IB):
                    self.SUME[J] = self.ZERO
                    self.SUMC[J] = self.ZERO

                #{{{ INIT
                #print 'NBL = %s' % self.NBLK
                if self.debug: print "\t WVINT      NBLK => %s " % self.NBLK
                if self.NBLK == 1:
                    self.NNN1 = 3
                    if self.debug: print "\t WVINT      NNN1 => %s " % self.NNN1
                    if self.debug: print "\t WVINT      IASYMP => %s " % self.IASYMP
                    #     FIND THE ASYMPTOTIC COEFFICIENTS ++++ NECESSARY ONLY WHEN IASYMP
                    #     IS .TRUE.                        +++++++++++++>
                    if self.IASYMP:

                        #{{{ setup(R)
                        #----------------------------------------------------------
                        #       JNKM = R INTEGRAL EXP(-KH) KSUP M J SUB N (KR) DK
                        #       THE R IN FRONT TAKES INTO ACCOUNT THE 1/R IN THE
                        #       DO 300 OF SUBROUTINE WVINT
                        #----------------------------------------------------------

                        if self.debug: print ""
                        if self.debug: print "%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%"
                        if self.debug: print "\t SETUP:   R       => %s" % self.R    
                        self.DIST  = csqrt(self.R**2 + self.DEPTH**2).real
                        self.DIST3 = self.DIST**3
                        self.DIST5 = self.DIST**5
                        self.DIST7 = self.DIST**7
                        self.RZ    = self.R*self.DEPTH
                        self.Z2    = self.DEPTH*self.DEPTH
                        self.R2    = self.R*self.R
                        self.R3    = self.R*self.R2
                        self.Z3    = self.DEPTH*self.Z2
                        self.RZ2   = self.R*self.Z2
                        self.RZ3   = self.R*self.Z3
                        self.ZOR   = self.DEPTH/self.DIST
                        self.ZOR2  = self.ZOR*self.ZOR
                        self.ZOR3  = self.ZOR*self.ZOR2
                        self.J0K0  = self.R/self.DIST
                        self.J0K1  = self.RZ/self.DIST3
                        self.J0K2  = (2.0 * self.RZ2 - self.R3) / self.DIST5
                        self.J0K3  = (6.0 * self.RZ3 - 9.0 * self.DEPTH * self.R3) / self.DIST7
                        self.J1K0  = 1.0 - self.DEPTH / self.DIST
                        self.J1K1  = self.R2 / self.DIST3
                        self.J1K2  = 3.0 * self.DEPTH * self.R2 / self.DIST5
                        self.J1K3  = 3.0 * self.R2 * (4.0 * self.Z2 - self.R2) / self.DIST7
                        self.J2K0  = ( 1.0 - 2.0 * self.ZOR + self.ZOR2 ) * (self.DIST / self.R)
                        self.J2K1  = (2.0 - 3.0 * self.ZOR + self.ZOR3) / self.R
                        self.J2K2  = 3.0 * self.R3 / self.DIST5
                        self.J2K3  = 15.0 * self.DEPTH * self.R3 / self.DIST7

                        if self.debug: print "\t SETUP:   DIST     => %s" % self.DIST  
                        if self.debug: print "\t SETUP:   DIST3    => %s" % self.DIST3 
                        if self.debug: print "\t SETUP:   DIST5    => %s" % self.DIST5 
                        if self.debug: print "\t SETUP:   DIST7    => %s" % self.DIST7 
                        if self.debug: print "\t SETUP:   RZ       => %s" % self.RZ    
                        if self.debug: print "\t SETUP:   Z2       => %s" % self.Z2    
                        if self.debug: print "\t SETUP:   R2       => %s" % self.R2    
                        if self.debug: print "\t SETUP:   R3       => %s" % self.R3    
                        if self.debug: print "\t SETUP:   Z3       => %s" % self.Z3    
                        if self.debug: print "\t SETUP:   RZ2      => %s" % self.RZ2   
                        if self.debug: print "\t SETUP:   RZ3      => %s" % self.RZ3   
                        if self.debug: print "\t SETUP:   ZOR      => %s" % self.ZOR   
                        if self.debug: print "\t SETUP:   ZOR2     => %s" % self.ZOR2  
                        if self.debug: print "\t SETUP:   ZOR3     => %s" % self.ZOR3  
                        if self.debug: print "\t SETUP:   J0K0     => %s" % self.J0K0  
                        if self.debug: print "\t SETUP:   J0K1     => %s" % self.J0K1  
                        if self.debug: print "\t SETUP:   J0K2     => %s" % self.J0K2  
                        if self.debug: print "\t SETUP:   J0K3     => %s" % self.J0K3  
                        if self.debug: print "\t SETUP:   J1K0     => %s" % self.J1K0  
                        if self.debug: print "\t SETUP:   J1K1     => %s" % self.J1K1  
                        if self.debug: print "\t SETUP:   J1K2     => %s" % self.J1K2  
                        if self.debug: print "\t SETUP:   J1K3     => %s" % self.J1K3  
                        if self.debug: print "\t SETUP:   J2K0     => %s" % self.J2K0  
                        if self.debug: print "\t SETUP:   J2K1     => %s" % self.J2K1  
                        if self.debug: print "\t SETUP:   J2K2     => %s" % self.J2K2  
                        if self.debug: print "\t SETUP:   J2K3     => %s" % self.J2K3  

                        if self.debug: print "%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%"
                        if self.debug: print ""


                        #}}} setup(R)

                        if I == 0:
                        #{{{
                            if self.debug: print ""
                            self.W1 = self.WVNO[0].real
                            self.W2 = self.WVNO[1].real
                            if self.debug: print "\t WVINT:   W1       => %s" % self.W1
                            if self.debug: print "\t WVINT:   W2       => %s" % self.W2
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

                            if self.debug: print "Done with SOLU"

                            for J in range(10):
                                if self.debug: print "\t WVINT 101: J => %s" % J
                                if self.debug: print "\t WVINT 101: ISRC[J] => %s" % self.ISRC[J]
                                if self.ISRC[J] != 1: 
                                    self.AA[J] = self.ZERO
                                    self.BB[J] = self.ZERO
                                    self.CC[J] = self.ZERO
                                if self.debug: print "WVINT  101 :   AA[J]=> %s" % self.AA[J]
                                if self.debug: print "WVINT  101 :   BB[J]=> %s" % self.BB[J]
                                if self.debug: print "WVINT  101 :   CC[J]=> %s" % self.CC[J]
                            #exit()
                        #}}}

                        # INITIALIZE THE INTEGRALS +++  NECESSARY FOR IASYMP.EQ.TRUE
                        #{{{
                        self.SMM[I][0] = self.AA[0] * self.J0K0 + self.BB[0] * self.J0K1 + self.CC[0] * self.J0K2
                        self.SMM[I][1] = self.AA[1] * self.J1K1 + self.BB[1] * self.J1K2 + self.CC[1] * self.J1K3
                        self.SMM[I][2] = self.AA[2] * self.J1K0 + self.BB[2] * self.J1K1 + self.CC[2] * self.J1K2 
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][0]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][1]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][2]

                        if self.ISRC[8] == 1 and self.ISRC[3] == 1: 
                            self.SUMD = (self.AA[3] + self.AA[8]) * self.J1K0 + (self.BB[3] + self.BB[8]) * self.J1K1 + (self.CC[3] + self.CC[8]) * self.J1K2
                        else:
                            self.SUMD=complex(0.0,0.0)

                        if self.debug: print "\t WVINT after 101 SUMD() => %s" % self.SUMD
                        self.SUMD= -self.SUMD / self.R
                        if self.debug: print "\t WVINT after 101 SUMD() => %s" % self.SUMD
                        self.SMM[I][3] = self.SUMD + self.AA[3] * self.J0K1 + self.BB[3] * self.J0K2 + self.CC[3] * self.J0K3
                        self.SMM[I][4] = self.SUMD + self.AA[8] * self.J0K1 + self.BB[8] * self.J0K2 + self.CC[8] * self.J0K3
                        self.SMM[I][5] = self.AA[4] * self.J2K0 + self.BB[4] * self.J2K1 + self.CC[4] * self.J2K2
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][3]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][4]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][5]

                        if self.ISRC[5] == 1 and self.ISRC[9] == 1:
                            self.SUMD = (self.AA[5] + self.AA[9]) * self.J2K0 + (self.BB[5] + self.BB[9]) * self.J2K1 + (self.CC[5] + self.CC[9]) * self.J2K2
                        else:
                            self.SUMD = complex(0.0,0.0)

                        if self.debug: print "\t WVINT after 101 SUMD() => %s" % self.SUMD
                        self.SUMD = -2.0 * self.SUMD / self.R
                        if self.debug: print "\t WVINT after 101 SUMD() => %s" % self.SUMD
                        self.SMM[I][6] = self.SUMD + self.AA[5] * self.J1K1 + self.BB[5] * self.J1K2 + self.CC[5] * self.J1K3
                        self.SMM[I][7] = self.SUMD + self.AA[9] * self.J1K1 + self.BB[9] * self.J1K2 + self.CC[9] * self.J1K3
                        self.SMM[I][8] = self.AA[6] * self.J0K0 + self.BB[6] * self.J0K1 + self.CC[6] * self.J0K2
                        self.SMM[I][9] = self.AA[7] * self.J1K1 + self.BB[7] * self.J1K2 + self.CC[7] * self.J1K3
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][6]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][7]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][8]
                        if self.debug: print "\t WVINT after 101 SMM() => %s" % self.SMM[I][9]

                        if self.debug: print ""
                        #}}}

                    else:
                        for II in range(10):
                            self.SMM[I][II] = complex(0.0,0.0)


                #}}} INIT

                #     GET THE BESSEL'S FUNCTIONS +++++++++++++++++++++++++++++++>
                if self.debug: print ""
                if self.debug: print "***************** BESSEL FUNCTIONS *****************"
                if self.debug: print ""

                if self.debug: print "BESSEL FUNCTIONS: N22: %s" % self.N22
                if self.debug: print "BESSEL FUNCTIONS: ICN: %s" % self.ICN
                if self.debug: print "BESSEL FUNCTIONS: NNN1: %s" % self.NNN1
                if self.debug: print "BESSEL FUNCTIONS: IBLOCK: %s" % self.IBLOCK
                #{{{
                if self.N22 < self.ICN:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        #if self.debug: print "BESSEL FUNCTIONS: 1st J: %s" % J
                        Z = self.WVNO[J].real * self.R
                        X = (Z/3.0) * (Z/3.0)
                        AJ0 = 1.0-X*(2.2499997-X*(1.2656208-X*(.3163866-X*( .0444479-X*(.0039444-X*(.0002100))))))
                        A1Z = 0.5-X*(.56249985-X*(.21093573-X*(.03954289-X*( .00443319-X*(.00031761-X*(.00001109))))))
                        AJ1 = Z * A1Z
                        self.J2[J] = 2.0 * AJ1 / Z - AJ0
                        self.J1[J] = AJ1
                        self.J0[J] = AJ0
                        #if self.debug: print "\t BESSEL FUNCTIONS: NNN1: %s" % self.NNN1
                        #if self.debug: print "\t BESSEL FUNCTIONS: IBLOCK: %s" % self.IBLOCK
                        #if self.debug: print "\t BESSEL FUNCTIONS: J: %s" % J
                        #if self.debug: print "\t BESSEL FUNCTIONS: R: %s" % self.R
                        #if self.debug: print "\t BESSEL FUNCTIONS: WVNO[J]: %s" % self.WVNO[J]
                        #if self.debug: print "\t BESSEL FUNCTIONS: Z: %s" % Z
                        #if self.debug: print "\t BESSEL FUNCTIONS: X: %s" % X
                        #if self.debug: print "\t BESSEL FUNCTIONS: AJ0: %s" % AJ0
                        #if self.debug: print "\t BESSEL FUNCTIONS: AJ1: %s" % AJ1
                        if self.debug: print "\t BESSEL FUNCTIONS: J: %s" % J
                        if self.debug: print "\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J]
                        #exit()


                if self.debug: print "BESSEL FUNCTIONS: N11: %s" % self.N11
                if self.N11 > self.ICN:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        Z = self.WVNO[J].real * self.R
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
                        #if self.debug: print "\t BESSEL FUNCTIONS: J: %s" % J
                        if self.debug: print "\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J]

                if self.N22 > self.ICN and self.N11 < self.ICN: 
                    self.IM = self.ICN - self.N11 + 1
                    #if self.debug: print "BESSEL FUNCTIONS: IM: %s" % self.IM
                    for J in range(self.NNN1-1,self.IM):
                        #if self.debug: print "BESSEL FUNCTIONS: 3rd J: %s" % J
                        Z = self.WVNO[J].real * self.R
                        X = (Z/3.0)*(Z/3.0)
                        AJ0 = 1.0-X*(2.2499997-X*(1.2656208-X*(0.3163866-X*(0.0444479-X*(0.0039444-X*(0.0002100))))))
                        A1Z = 0.5-X*(0.56249985-X*(0.21093573-X*(0.03954289-X*(0.00443319-X*(0.00031761-X*(0.00001109))))))
                        AJ1 = Z * A1Z
                        self.J2[J] = 2.0 * AJ1 / Z - AJ0
                        self.J1[J] = AJ1
                        self.J0[J] = AJ0
                        #if self.debug: print "\t BESSEL FUNCTIONS: J: %s" % J
                        if self.debug: print "\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J]
                        #exit()

                    for J in range(self.IM+1,self.IBLOCK):
                        #if self.debug: print "BESSEL FUNCTIONS: 4th J: %s" % J
                        Z = self.WVNO[J].real * self.R
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
                        if self.debug: print "\t BESSEL FUNCTIONS: J2[J]: %s" % self.J2[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J1[J]: %s" % self.J1[J]
                        if self.debug: print "\t BESSEL FUNCTIONS: J0[J]: %s" % self.J0[J]
                        #exit()

                if self.debug: print "BESSEL FUNCTIONS:   J2 => %s" % self.J2
                if self.debug: print "BESSEL FUNCTIONS:   J1 => %s" % self.J1
                if self.debug: print "BESSEL FUNCTIONS:   J0 => %s" % self.J0
                #}}}

                #     GET THE PHASE VELOCITY FILTER  +++++++++++++++++++++++++>
                #{{{
                if self.debug: print ''
                if self.debug: print '############### PHASE VELOCITY FILTER ##############'
                if self.debug: print ''

                if self.CMAX >= 0.0:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.FACT[J] = 0.0
                        self.WVN = self.WVNO[J].real

                        if self.WVN > self.WVC1 and self.WVN < self.WVC2: self.FACT[J] = 1.0
                        if self.WVN > self.WVCM and self.WVN < self.WVC1: 
                            self.FACT[J] = (1.0 - cos(self.PI * (self.WVN-self.WVCM) / (self.WVC1-self.WVCM))).real / 2.0
                        if self.WVN > self.WVC2 and self.WVN < self.WVCN:
                            self.FACT[J] = (1.0 - cos(self.PI * (self.WVN - self.WVCN) / (self.WVC2 - self.WVCN))).real / 2.0

                        if self.debug: print ' PHASE VELOCITY FILTER:    FACT[J] => %s' % self.FACT[J]
                #}}}

                #     REMOVAL OF ASYMPTOTIC TREND FROM THE INTEGRAND
                #{{{
                if self.debug: print ''
                if self.debug: print '############### ASYMP TREND REMOVAL ##############'
                if self.debug: print ''

                if self.IASYMP:
                    for K in range(10):
                        if self.ISRC[K] == 1:
                            for J in range(self.NNN1-1,self.IBLOCK):
                                self.AK = self.WVNO[J].real
                                self.EX = exp(-self.AK * self.DEPTH)
                                self.ASM = self.EX * (self.AA[K] + self.AK * (self.BB[K] + self.AK * self.CC[K]))
                                if self.debug: print "    RM ASYMP TREND  10:    AK => %s  " % self.AK
                                if self.debug: print "    RM ASYMP TREND  10:    H => %s  " % self.DEPTH
                                if self.debug: print "    RM ASYMP TREND  10:    EX => %s  " % self.EX
                                if self.debug: print "    RM ASYMP TREND  10:    AA[K] => %s  " % self.AA[K]
                                if self.debug: print "    RM ASYMP TREND  10:    BB[K] => %s  " % self.BB[K]
                                if self.debug: print "    RM ASYMP TREND  10:    CC[K] => %s  " % self.CC[K]
                                if self.debug: print "    RM ASYMP TREND  10:    ASM => %s  " % self.ASM
                                if self.debug: print "    RM ASYMP TREND  10:    S1[J] => %s  " % self.G1[J]
                                if self.debug: print "    RM ASYMP TREND  10:    K => %s  " % K
                                if self.debug: print "    RM ASYMP TREND  10:    J => %s  " % J
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
                                if self.debug: print "    RM ASYMP TREND  10:    G1[J] =>  %s " % self.S1[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G2[J] =>  %s " % self.S2[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G3[J] =>  %s " % self.S3[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G4[J] =>  %s " % self.S4[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G5[J] =>  %s " % self.S5[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G6[J] =>  %s " % self.S6[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G7[J] =>  %s " % self.S7[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G8[J] =>  %s " % self.S8[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G9[J] =>  %s " % self.S9[J]
                                if self.debug: print "    RM ASYMP TREND  10:    G10[J] => %s  " % self.S10[J]

                if self.NBLK == self.NBLOCK:
                    if self.IASYMP:
                        for K in range(10):
                            if self.ISRC[K] == 1:
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
                                if self.debug: print "    RM ASYMP TREND  666:    G1[J] =>  %s " % self.S1[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G2[J] =>  %s " % self.S2[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G3[J] =>  %s " % self.S3[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G4[J] =>  %s " % self.S4[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G5[J] =>  %s " % self.S5[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G6[J] =>  %s " % self.S6[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G7[J] =>  %s " % self.S7[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G8[J] =>  %s " % self.S8[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G9[J] =>  %s " % self.S9[1]
                                if self.debug: print "    RM ASYMP TREND  666:    G10[J] => %s  " % self.S10[1]

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
                                if self.debug: print "    RM ASYMP TREND  777:    G1[J] =>  %s " % self.S1[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G2[J] =>  %s " % self.S2[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G3[J] =>  %s " % self.S3[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G4[J] =>  %s " % self.S4[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G5[J] =>  %s " % self.S5[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G6[J] =>  %s " % self.S6[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G7[J] =>  %s " % self.S7[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G8[J] =>  %s " % self.S8[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G9[J] =>  %s " % self.S9[1]
                                if self.debug: print "    RM ASYMP TREND  777:    G10[J] => %s  " % self.S10[1]


                if self.debug: print "    RM ASYMP TREND  13:    ISRC => %s  " %  self.ISRC

                if self.ISRC[0] == 1:
                    #if self.debug: print "    RM ASYMP TREND  13:    IBL => %s  " %  self.IBLOCK
                    #if self.debug: print "    RM ASYMP TREND  13:    N1 => %s  " %  self.NNN1
                    for J in range(self.IBLOCK-1,self.NNN1-2,-1):
                        self.SUMD = self.FACT[J] * self.S1[J] * self.J0[J]
                        self.SMM[I][0] = self.SMM[I][0] + self.SUMD
                        #if self.debug: print "    RM ASYMP TREND  13:    J => %s  " %  J
                        #if self.debug: print "    RM ASYMP TREND  13:    FACT[J] => %s  " % self.FACT[J]
                        #if self.debug: print "    RM ASYMP TREND  13:    G1[J] => %s  " % self.S1[J]
                        #if self.debug: print "    RM ASYMP TREND  13:    J0[J] => %s  " % self.J0[J]
                        #if self.debug: print "    RM ASYMP TREND  13:    SUMD => %s  " % self.SUMD
                        if self.debug: print "    RM ASYMP TREND  13:    SMM[I][0] => %s  " % self.SMM[I][0]


                if self.ISRC[1] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.AK = self.WVNO[J].real
                        self.SUMD = self.FACT[J] * self.S2[J] * self.AK * self.J1[J]
                        self.SMM[I][1] = self.SMM[I][1] + self.SUMD
                        #if self.debug: print "    RM ASYMP TREND  13:    J => %s  " %  J
                        #if self.debug: print "    RM ASYMP TREND  13:    AK => %s  " % self.AK
                        #if self.debug: print "    RM ASYMP TREND  13:    SUMD => %s  " % self.SUMD
                        if self.debug: print "    RM ASYMP TREND  14:    SMM[I][1] => %s  " % self.SMM[I][1]

                if self.ISRC[2] == 1: 
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.SUMD = self.FACT[J] * self.S3[J] * self.J1[J]
                        self.SMM[I][2] = self.SMM[I][2] + self.SUMD
                        if self.debug: print "    RM ASYMP TREND  15:    SMM[I][2] => %s  " % self.SMM[I][2]
                #}}}

                #     INCLUDE NEAR FIELD TERM IF BOTH SH AND P-SV ARE COMPUTED
                #{{{
                if self.ISRC[3] == 1 and self.ISRC[8] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.SUMC[J] = (self.S4[J] + self.S9[J]) * self.J1[J]
                    #    if self.debug: print "    RM ASYMP TREND  16:    J => %s  " % J
                    #    if self.debug: print "    RM ASYMP TREND  16:    G4[J] => %s  " % self.S4[J]
                    #    if self.debug: print "    RM ASYMP TREND  16:    G9[J] => %s  " % self.S9[J]
                    #    if self.debug: print "    RM ASYMP TREND  16:    J1[J] => %s  " % self.J1[J]
                        if self.debug: print "    RM ASYMP TREND  16:    SUMC[J] => %s  " % self.SUMC[J]
                    #print ""
                    #print "\t\t\t G4 => %s " % self.S4
                    #print "\t\t\t G9 => %s " % self.S9
                    #print ""

                if self.ISRC[3] == 1:
                    for J in range(self.NNN1-1, self.IBLOCK):
                        self.AK= self.WVNO[J].real
                        self.SUMD = self.S4[J] * self.AK * self.J0[J]
                        self.SMM[I][3] = self.SMM[I][3] + (self.SUMD - self.SUMC[J] / self.R) * self.FACT[J]
                        if self.debug: print "    RM ASYMP TREND  17:    SMM[I][3] => %s  " % self.SMM[I][3]

                if self.ISRC[8] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.AK = self.WVNO[J].real
                        self.SUMD = self.S9[J] * self.AK * self.J0[J]
                        self.SMM[I][4] = self.SMM[I][4] + (self.SUMD - self.SUMC[J] / self.R) * self.FACT[J]
                        if self.debug: print "   RM ASYMP TREND  18 :    SMM[I][4] => %s  " % self.SMM[I][4]

                if self.ISRC[4] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.SUMD = self.FACT[J] * self.S5[J] * self.J2[J]
                        self.SMM[I][5] = self.SMM[I][5] + self.SUMD
                        if self.debug: print "    RM ASYMP TREND  19:    SMM[I][5] => %s  " % self.SMM[I][5]
                #}}}

                #     NEAR FIELD TERM IF BOTH SH AND P-SV TERMS ARE COMPUTED
                #{{{
                if self.ISRC[5] == 1 and self.ISRC[9] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.SUME[J] = (self.S6[J] + self.S10[J]) * self.J2[J]
                        if self.debug: print "    RM ASYMP TREND  20:    SUME[J] => %s  " % self.SUME[J]

                if self.ISRC[5] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.AK = self.WVNO[J].real
                        self.SUMD = self.S6[J] * self.AK * self.J1[J]
                        self.SMM[I][6] = self.SMM[I][6] + (self.SUMD - 2.0 * self.SUME[J] /self.R) * self.FACT[J]
                        if self.debug: print "     RM ASYMP TREND  21:    SMM[I][6] => %s  " % self.SMM[I][6]

                if self.ISRC[9] == 1:
                    for J in range(self.IBLOCK-1,self.NNN1-2,-1):
                        self.AK = self.WVNO[J].real
                        self.SUMD = self.S10[J] * self.AK * self.J1[J]
                        self.SMM[I][7] = self.SMM[I][7] + (self.SUMD - 2.0 * self.SUME[J] / self.R) * self.FACT[J]
                        if self.debug: print "     RM ASYMP TREND  22:    SMM[I][7] => %s  " % self.SMM[I][7]

                if self.ISRC[6] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.SUMD = self.FACT[J] * self.S7[J] * self.J0[J]
                        self.SMM[I][8] = self.SMM[I][8] + self.SUMD
                        if self.debug: print "     RM ASYMP TREND  23:    SMM[I][8] => %s  " % self.SMM[I][8]

                if self.ISRC[7] == 1:
                    for J in range(self.NNN1-1,self.IBLOCK):
                        self.AK = self.WVNO[J].real
                        self.SUMD = self.FACT[J] * self.S8[J] * self.AK * self.J1[J]
                        self.SMM[I][9] = self.SMM[I][9] + self.SUMD
                        if self.debug: print "     RM ASYMP TREND  24:    SMM[I][9] => %s  " % self.SMM[I][9]

                if self.debug: print ""
                if self.debug: print " END"
                if self.debug: print ""

                #}}}

                #}}} sub WVINT

            #}}} FR

            self.NBLK += 1

        #}}} A

        #print ' DEBUG: T0 = %s' % self.T0
        for I in range(self.NRANGE):
        #{{{ B
            #print ' DEBUG: I = %s' % I
            #print ' DEBUG: T0[%s] = %s' % (I,self.T0[I])
            #print ' DEBUG: RANGE[%s] = %s' % (I,self.RANGE[I])
            #print ' DEBUG: VRED[%s] = %s' % (I,self.VRED[I])
            if self.VRED[I] >= 0.0: self.VRED[I] = 10.0
            self.T0X = self.T0[I] + self.RANGE[I] / self.VRED[I]
            self.FAT = self.OMEGA * self.T0X
            self.XR = cos(self.FAT)
            self.XI = sin(self.FAT)
            self.DKK = self.DK * self.RANGE[I]
            #self.SMM[I][1] *= -1 * self.SMM[I][1]
            #self.SMM[I][9] *= -1 * self.SMM[I][9]
            self.SMM[I][1] *= -1
            self.SMM[I][9] *= -1
            for J in range(10):
                self.SMM[I][J] = self.SMM[I][J] * self.DKK
                self.SMM[I][J] = self.SMM[I][J] / self.RANGE[I]

            for J in range(10):
                if self.ISRC[J] != 1: break
                self.GG[I][J] = self.SMM[I][J] * complex(self.XR,self.XI)
                if self.debug: print "%s" % self.GG[I][J]
        #}}} B

        if self.debug: print "%s" % self.GG
        return 

#}}}

    def _solu300(self,Y1,Y2,X1,X2,J):
        #{{{
        if self.debug: print ""
        if self.debug: print "&&&&&&&&&&&&&&&&&&& SOLU &&&&&&&&&&&&&&&"
        self.CC[J] = complex(0.0,0.0)
        U1 = X1 * self.DEPTH
        U2 = X2 * self.DEPTH
        DET = X2 - X1
        if self.debug: print "SOLU 300: U1 =%s" % U1
        if self.debug: print "SOLU 300: U2 =%s" % U2
        if self.debug: print "SOLU 300: DET =%s" % DET
        self.AA[J] = complex(0.0,0.0)
        if self.debug: print "^^^^^^^^ SOLU 300: X1 => %s" % X1
        if self.debug: print "^^^^^^^^ SOLU 300: X2 => %s" % X2
        if self.debug: print "^^^^^^^^ SOLU 300: Y1 => %s" % Y1
        if self.debug: print "^^^^^^^^ SOLU 300: Y2 => %s" % Y2
        if self.debug: print "^^^^^^^^ SOLU 300: U1 => %s" % U1
        if self.debug: print "^^^^^^^^ SOLU 300: U2 => %s" % U2

        self.BB[J] = X2 * Y1 * exp(U1) / X1
        if self.debug: print "^^^^^^^^ SOLU 300: BB[%s]:%s" % (J,self.BB[J])
        self.BB[J] = X1 * Y2 * exp(U2) / X2
        if self.debug: print "^^^^^^^^ SOLU 300: BB[%s]:%s" % (J,self.BB[J])

        self.BB[J] = X2 * Y1 * exp(U1) / X1 - X1 * Y2 * exp(U2) / X2
        self.BB[J] = self.BB[J] / DET
        self.CC[J] = Y2 * exp(U2) / X2 - Y1 * exp(U1) / X1
        self.CC[J] = self.CC[J] / DET
        if self.debug: print "SOLU 300: AA[%s]:%s" % (J,self.AA[J])
        if self.debug: print "SOLU 300: BB[%s]:%s" % (J,self.BB[J])
        if self.debug: print "SOLU 300: CC[%s]:%s" % (J,self.CC[J])
        if self.debug: print ""
        #exit()
        #}}}

    def _solu200(self,Y1,Y2,X1,X2,J):
        #{{{
        if self.debug: print ""
        if self.debug: print "&&&&&&&&&&&&&&&&&&& SOLU &&&&&&&&&&&&&&&"
        self.CC[J] = complex(0.0,0.0)
        U1 = X1 * self.DEPTH
        U2 = X2 * self.DEPTH
        DET = X2 - X1
        if self.debug: print "SOLU 200: U1 =%s" % U1
        if self.debug: print "SOLU 200: U2 =%s" % U2
        if self.debug: print "SOLU 200: DET =%s" % DET
        self.AA[J] = X2 * Y1 * exp(U1) - X1 * Y2 * exp(U2)
        self.AA[J] = self.AA[J] / DET
        self.BB[J] = Y2 * exp(U2) - Y1 * exp(U1)
        self.BB[J] = self.BB[J] / DET
        if self.debug: print "SOLU 200: AA[%s]:%s" % (J,self.AA[J])
        if self.debug: print "SOLU 200: BB[%s]:%s" % (J,self.BB[J])
        if self.debug: print ""
        #}}}

    def wvint9(self):
#{{{
        self.NM = self.N2
        self.REP = 'n'

        """
        Change ivel to "v" to output velocity.
        """
        self.IVEL = 'd'


        self.N = 2 * (self.NYQ-1)

        self.NPOINT = self.N
        self.TAU = self.DT

        self.FMAX = self.FU
        self.inst = 0

        """
        Beginning of distance loop
        """

        self.PI = cos(-1.0)
        self.TWPI = 2.0 * self.PI

        #self.OMEGA



        for nd in range(self.NRANGE):
            T0X = self.RANGE(nd)/self.VRED(nd)
            yr = 0.0

            if self.debug: print "range=%s yr=%s depth=%s npoint=%s t0x=%s dt=%s tau=%s" % (self.RANGE(nd),yr,self.DEPTH,self.NPOINT,T0X,self.DT,self.TAU))
            #write(11) range(nd),yr,depth,npoint,t0x,dt,tau

            for l in range(10):

                if self.ISRC[l] != 1:
                    for j in range(n):
                        self.DATA[j] = complex(0.0,0.0)
                else:
                    #if(isrc(l).ne.1) go to 444
                    for j in range(self.NM):
                        self.DATA[j]=self.GG(nd,j,l)

                    for j in range(self.NM+1,self.NYQ):
                        self.DATA[j]=complex(0.0,0.0)

                    for j in range(self.NM):
                        self.FREQ=(j-1)*self.DF
                        if self.FREQ > self.DF: 
                            freq=0.01*df

                        self.OM = self.TWPI * self.FREQ

                        if j != 1: 
                            seld.DATA[n+2-j] = self.DATA[j].conjugate()

                    for i in range(self.NM+1,self.NYQ):
                        self.DATA[i] = complex(0.,0.)
                        self.DATA[n+2-i] = self.DATA[i].conjugate()

                    #self.DATA[1] = complex(self.DATA[1].real,0.0)
                    #self.DATA[self.NYQ] = complex(seld.DATA[self.NYQ].real,0.0)
                    self.DATA[1] = complex(0.0,0.0)
                    self.DATA[self.NYQ] = complex(0.0,0.0)

                    #
                    #call four(data,n,+1,dt,df)
                    #
                    self._four()

                    """
                    Correct for damping.
                    """
                    fac = exp(self.ALPHA  * T0X)
                    dfac = exp(self.ALPHA*self.DT)

                    for i in range(npoints):
                        self.DATA[i] = self.DATA[i] * fac
                        fac = fac * dfac
                """
                Time Domain Integration
                """
                if self.IVEL == 'd':
                    totint = 0.0
                    for i in range(self.NPOINT):
                        temp[i] = self.DATA[i].real()
                    for i in range(self.NPOINT-1):
                        prtint = 0.5 * self.DT * (temp[i] + temp[i+1])
                        totint = totint+printt
                        temp[i] = totint
                    for i in range(self.NPOINT-1):
                        self.DATA[i] = temp[i] - temp[1]

                if l == 1:
                    for i in range(self.NPOINT):
                       self.DATA[i] = self.DATA[i] * -1.0 
                    print "\nNPTS %s    DT %s\n" % (self.NPOINT,self.DT)

                if l == 4:
                    for i in range(self.NPOINT):
                        self.DATA[i] = self.DATA[i] * -1.0

                if l == 6:
                    for i in range(self.NPOINT):
                        self.DATA[i] = self.DATA[i] * -1.0

                if l == 8:
                    for i in range(self.NPOINT):
                        self.DATA[i] = self.DATA[i] * -1.0

                for i in range(self.NPOINT):
                    dataout = self.DATA[i].real()

                print "DATAOUT: %s\n" dataout
                print "NPOINT: %s\n" self.NPOINT
                print "ND: %s\n" self.ND
                print "l: %s\n" l

#}}}

    def _four(self):
#{{{
        """ 
        subroutine four(data,nn,isign,dt,df)
        """ 

        n = 2 * self.NN

        if self.DT == 0.0: self.DT = 1/(self.NN*self.DF)
        if self.DF == 0.0: self.DF = 1/(self.NN*self.DT)
        if self.DT != (self.NN * self.DF): self.DF = 1/(self.NN*self.DT)

        j = 1
        for i in range(1,n,2):
            if (i - j) < 0:
                tempr = self.DATA[j]
                tempi = self.DATA[j+1]
                self.DATA[j] = self.DATA[i]
                self.DATA[j+1] = self.DATA[i+1]
                self.DATA[i] = tempr
                self.DATA[i] = tempi

            m = n / 2
            while (m - 2) >= 0:
                if (j - m) <= 0:
                    break
                j = j - m
                m = m / 2

            j = j + m
            mmax = 2

            while (mmax - n) < 0:
                istep = 2 * mmax
                theta = 6.283185307/mmax
                sinth = sin(theta/2)
                wstpr = -2 * sinth * sinth
                wstpi = sin(theta)
                wr = 1
                wi = 0
                for m in range(1,mmax,2):
                    for i in range(m,n,istep):
                        j = i + mmax
                        tempr = wr * self.DATA[j] - wi * self.DATA[j+1]
                        tempi = wr * self.DATA[j+1] + wi * self.DATA[j]
                        self.DATA[j] = self.DATA[i] - tempr
                        self.DATA[j+1] = self.DATA[i+1] - tempi
                        self.DATA[i] = self.DATA[i] + tempr
                        self.DATA[i+1] = self.DATA[i+1] + tempi 

                    tempr = wr
                    wr = wr * wstpr - wi * wstpi + wr
                    wi = wi * wstpr + tempr * wstpi + wi
                mmax = istep

            #
            # time to frequency domain
            #
            for h in range(n):
                self.DATA[h] = self.DATA[h] * self.DT
#}}}


if __name__ == "__main__":
    """
    If we call this script directly, then...
    """

    file = configure()
    if verbose: print 'fkrprog done configure()'
    if verbose: print "fkrprog start green_functions(%s,True,True)" % file
    GF = green_functions(file,False,False)
    if verbose: print 'Done!'


import os
import re
import sys
import getopt
from collections import defaultdict

# ANTELOPE
try:
    import antelope.stock as stock
    import antelope.datascope as datascope
except Exception,e:
    sys.exit("Antelope Import Error: [%s] => [%s]" % (Exception,e))

# Needed for fkrprog.py
import matplotlib as mpl
mpl.use('tkagg')
import pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from numpy.fft  import ifft
from numpy import trapz

from math import exp, log, sqrt, acos, asin, cos, sin
from cmath import log as clog
from cmath import sqrt as csqrt
from cmath import exp as cexp

import moment_tensor.fkrprog as fkr

"""
Get all information from command line 
and verify the source file is accessible.
"""
try:

    opts, pargs = getopt.getopt(sys.argv[1:], 'vtp:')

except getopt.GetoptError:

    raise SystemExit('\n\n\tUSAGE: fkrprog [-v] [-p model_pf_file]\n')

if len(pargs) != 0: raise SystemExit('\n\n\tUSAGE: fkrprog [-v] [-p model_pf_file]\n')


verbose = False
pf_file = 'SOCAL_MODEL'

for option, value in opts:

    if '-p' in option:
        pf_file = str(value)

    if '-v' in option:
        verbose = True

"""
Build object for GFs.
"""
if verbose: print 'Load lib: GreenFunctions(%s)' % pf_file
GF =  fkr.GreenFunctions(pf_file,verbose)

#"""
#Generate GFs for depth of 8km and distance of 10km.
#"""
#if verbose: print 'generate(depth=%s,distance=%s)' % (8,10)
#GF.build(depth=8,distance=10,sps=1,type='v',filter='BW 0.01 5 0.05 5')

#"""
#Plot the GFs
#"""
#if verbose: print 'plot()'
#GF.plot()

if verbose: print 'generate(depth=%s,distance=%s)' % (8,10)
#GF.build(depth=30,distance=200,sps=1,type='v',filter='BW 0.01 5 0.05 5')
GF.build(depth=5,distance=190,sps=1,type='v',filter='BW 0.01 5 0.05 5')

"""
Plot the GFs
"""
if verbose: print 'plot()'
GF.plot()

if verbose: print 'Done!'



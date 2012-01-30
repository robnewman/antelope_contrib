import getopt
import moment_tensor.fkrprog as fkr
try:
    import antelope.stock as stock
except Exception,e:
    print  "Import Error: Cannot import Antelope's stock lib."
    exit()


"""
Get all information from command line 
and verify the source file is accessible.
"""
try:

    opts, pargs = getopt.getopt(sys.argv[1:], 'vfp:')

except getopt.GetoptError:

    print "Usage: fkrprog [-v] [-f] [-p model_pf_file]\n"
    sys.exit(-1)

if( len(pargs) != 0):

    print "Usage: fkrprog [-v] [-f] [-p model_pf_file]\n"
    sys.exit(-1)


verbose = False
fortran = False
pf_file = 'SOCAL_MODEL'

for option, value in opts:

    if '-p' in option:
        pf_file = str(value)

    if '-v' in option:
        verbose = True

    if '-f' in option:
        fortran = True

"""
Build object for GFs.
"""
if verbose: print 'Load lib: GreenFunctions(%s)' % pf_file
GF =  fkr.GreenFunctions(pf_file)

"""
Generate GFs for depth of 8km and distance of 1km.
"""
if fortran:
    if verbose: print 'generate_fortran(depth=%s,distance=%s)' % (8,10)
    GF.generate_fortran(depth=8,distance=10)
else:
    if verbose: print 'generate(depth=%s,distance=%s)' % (8,10)
    GF.generate(depth=8,distance=10)

"""
Plot the GFs
"""
if verbose: print 'plot()'
GF.plot()

if verbose: print 'Done!'



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

    opts, pargs = getopt.getopt(sys.argv[1:], 'vp:')

except getopt.GetoptError:

    print "Usage: fkrprog [-v] [-p model_pf_file]\n"
    sys.exit(-1)

if( len(pargs) != 0):

    print "Usage: fkrprog [-v] [-p model_pf_file]\n"
    sys.exit(-1)


verbose = False
pf_file = 'SOCAL_MODEL'

for option, value in opts:

    if '-p' in option:
        pf_file = str(value)

    if '-v' in option:
        verbose = True

"""
Build object for GF's.
"""
#GF =  fkr.GreenFunctions(pf_file,verbose=verbose)
GF =  fkr.GreenFunctions(pf_file)
"""
Generate GF's for depth of 8km and distance of 1km.
"""
GF.generate_fortran(depth=8,distance=10)
#GF.generate(depth=8,distance=10)
GF.test_plot()
#GF.plot()
#if verbose: print 'Done!'



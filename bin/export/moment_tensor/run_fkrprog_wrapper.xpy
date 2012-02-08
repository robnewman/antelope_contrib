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

"""
Generate GFs for depth of 8km and distance of 10km.
"""
if verbose: print 'generate(depth=%s,distance=%s)' % (8,10)
GF.build(depth=8,distance=10)

"""
Plot the GFs
"""
if verbose: print 'plot()'
GF.plot()

"""
Generate GFs for depth of 25 and distance of 20km.
"""
if verbose: print 'generate(depth=%s,distance=%s)' % (25,20)
GF.build(depth=25,distance=20)

"""
Plot the GFs
"""
if verbose: print 'plot()'
GF.plot()

if verbose: print 'Done!'



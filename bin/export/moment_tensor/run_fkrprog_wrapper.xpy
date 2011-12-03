import getopt
import moment_tensor.fkrprog as fkr


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

debug = False
verbose = False

for option, value in opts:

    if '-v' in option:
        debug = False
        verbose = True

    if '-d' in option:
        debug = True
        verbose = True

if not os.path.isfile(SOURCE_FILE):
    raise SystemExit('\n\nCannot find specified file! (%s)\n'% SOURCE_FILE)

"""
Build object for GF's.
"""
GF =  fkr.GreenFunctions(SOURCE_FILE)
"""
Generate GF's for depth of 8km and distance of 1km.
"""
GF.generate(depth=8,distance=80)
GF.plot()
#print GF('TSS')
#print GF['XDS']
#print GF.ALL
if verbose: print 'Done!'



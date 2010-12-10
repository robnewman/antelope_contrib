#! /opt/antelope/5.0-64/local/bin/python 


# Import statements
import sys
import os
import getopt
import re
import numpy as np
import math as math
from datetime import datetime
from collections import defaultdict
from time import gmtime

sys.path.append( os.environ['ANTELOPE'] + '/local/data/python' )
import antelope.stock as stock
from antelope.datascope import *

# Set defaults

pfname  = 'dbmoment'
chan_to_use = 'LH*'
model_type = 'v'
isoflag = 5
dw = False
use_incidence = 0
trim_value = 0.6
statmax = 20
event_db = ''
green_db = ''
wave_db = ''
resp_db = ''
filters = {}
mag_filters = {}
fstring = None
always = True
verbose = False
debug = False
flag = 0

def usage():
#{{{
#
# Message printed whenever the command-line input is incorrect
#
    print 'Usage: dbmoment [-vV] [-p pfname] orid'
#}}}
def logmt(flag,message):
#{{{
#
# Function to handle log output, which depends on the verbosity setting.
# Prints messages with timestamp. Whenever an error occurs it also exits.
#
    from time import time
    curtime = stock.strtime(time())
    if not flag:
        return
    elif flag < 3:
        print curtime,message
    else:
        print curtime,'ERROR:',message,'--> exiting'
        sys.exit(-1)
#}}}
def configure():
#{{{
#
# Function to configure paramters from command-line and the pf-file
# Check command line options
    try:
        opts,pargs = getopt.getopt(sys.argv[1:],'vVp:')
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    if(len(pargs) != 1):
        usage()
        sys.exit(-1)
    else:
        orid = pargs[0]
# Get command line options
    for option, value in opts:
        if option in ('-v'):
            globals()['verbose'] = True
        if option in ('-V'):
            globals()['verbose'] = True
            globals()['debug'] = True
        if option in ('-p'):
            globals()['pfname'] = value
# Get parameters from pf-file
    globals()['chan_to_use'] = stock.pfget_string(pfname,'chan_to_use')
    globals()['trim_value']  = stock.pfget_string(pfname,'trim_value')
    if(stock.pfget_string(pfname,'isoflag') == 1): globals()['isoflag'] = 6
    globals()['model_type']  = stock.pfget_string(pfname,'model_type')
    globals()['statmax']     = stock.pfget_int(pfname,'statmax')
    globals()['use_inc']     = stock.pfget_string(pfname,'use_inc')
    globals()['event_db']    = stock.pfget_string(pfname,'event_db')
    try:globals()['wave_db'] = stock.pfget_string(pfname,'wave_db')
    except:globals()['wave_db'] = event_db
    globals()['green_db']    = stock.pfget_string(pfname,'green_db')
    try:globals()['resp_db']     = stock.pfget_string(pfname,'resp_db')
    except:globals()['resp_db'] = event_db        
    globals()['mag_filters'] = stock.pfget_arr(pfname,'mag_filters')
    if(stock.pfget_string(pfname,'distance_weighting') == 'on'): globals()['dw'] = True
    if not wave_db: globals()['wave_db'] = event_db
    if not resp_db: globals()['resp_db'] = event_db

    return orid
#}}}
def get_view_from_db(orid):
#{{{
#
# Open event database and get event parameters for given origin-id
# Returns event parameters or exits when an error occurs
#
    if not os.path.isfile(event_db):
        logmt(3,'Database (%s) does not exist' % event_db)
    evdb = dbopen(event_db,'r')
    evdb = dblookup(evdb,'','origin','','')
    logmt(verbose,'Processing origin %s' %orid)
    evdb = dbsubset(evdb,'orid == %s' % orid)
    if evdb.nrecs() == 0:
        logmt(3,'Origin id (%s) does not exist in origintable'%orid)
    evdb = dbjoin(evdb,'netmag')
    if evdb.nrecs() == 0:
        logmt(3,'Could not join netmag tabel for orid %s' % orid)
    evparams = {}
    evdb[3] = 0
    for field in dbquery(evdb,'dbTABLE_FIELDS'):
        try: 
            evparams[field] = evdb.getv(field)
        except: logmt(3,'Could not find field %s in join of origin and netmag table for orid %s' % (field,orid))

    evdb = dbjoin(evdb,'assoc')
    evdb = dbjoin(evdb,'arrival')
    evdb = dbsort(evdb,'delta')
    evdb = dbsubset(evdb,'iphase=~/.*P.*|.*p.*/')
    if evdb.nrecs() == 0:
        logmt(3,'No arrivals for selected origin')
    for line in mag_filters:
        splitline = line.split()
        if not len(splitline)==3:
            lm,flt=splitline
        else:
            lm,um,flt = splitline
        if (float(evparams['magnitude'][0]) > float(lm) and float(evparams['magnitude'][0]) < float(um)):
            flt = flt.replace('_',' ')
            globals()['fstring'] = flt
    return(evdb,evparams)
#}}}
def choose_chan(wvstadb):
#{{{
#
# Get all channels for a given station. Returns the channels with 
# samplerate closest to and bigger than 1
#
    view = dbsort(wvstadb,'samprate',unique=True)
    for i in range(view.nrecs()):
        view[3]=i
        chan = dbgetv(view,'chan')
        chanview = dbsubset(wvstadb,'chan =~ /%s.*/' % chan[0][:2])
        chanview = dbsort(chanview,'chan',unique=True)
        channels = []
        if chanview.nrecs() == 3:
            for j in range(chanview.nrecs()):
                chanview[3] = j
                chan = dbgetv(chanview,'chan')
                channels.append(chan[0])
        else:
            logmt(verbose,'Did not find 3 components for %s, trying higher sample rate' % chan)
            channels = []
            continue
        if len(channels) == 3:
            break
    return channels
#}}}
def get_chan_data(wave_db,evdb):
#{{{
#
# Open wave_form database and returns trace objects based on sta_chan
# Applies calibration, splice, filter, decimation and rotation to channels
#
    try:
        wvdb = dbopen(wave_db,'r')
    except:
        logmt(3,'Could not open waveform database %s' % wave_db)
    wvdb = dblookup(wvdb,'','wfdisc','','')
    wvdb = dbsubset(wvdb,'samprate>=0.9')
    try:
        wvdb = dbsubset(wvdb,'chan =~/%s/' % chan_to_use)
    except:
        logmt(3,'No records found in wfdisc table (%s.wfdisc)' %wave_db) 
    numrec = evdb.nrecs()
    if numrec > int(statmax): numrec = int(statmax)
    logmt(verbose,'Processing %s stations for orid %s' % (numrec,orid))
    stachan_traces = {}
    counter = 0
    for i in range(evdb.nrecs()):
        evdb[3] = i
        (sta,at,esaz) = dbgetv(evdb,'sta','arrival.time','esaz')
        st = at - 12 
        et = at + 12
        wvstadb = dbsubset(wvdb,'sta=~/^%s$/' % sta)
        logmt(debug,'Looking for channels with a sample rate >= 1 Hz for sta = %s' % sta)
        try:
            chans = choose_chan(wvstadb)
        except:
            logmt(verbose,'No channels found with a sample-rate >= 1 Hz for sta = %s' % sta)
        logmt(debug,'Channels found: %s %s %s' % (chans[0],chans[1],chans[2]))
        wvstadb.subset('chan =~ /%s|%s|%s/' % (chans[0],chans[1],chans[2]))
        resample = 0
        for j in range(wvstadb.nrecs()):
            wvstadb[3] = j
            if wvstadb.getv('samprate')[0] != 1:
                logmt (verbose,'Samplerate higher than 1.0 --> resample')
                resample = 1
            break
        if resample == 1:
            logmt(always,'Resampling to be implemented, skipping station %s for now' % sta)
            continue
        vang = 0.0
        if use_inc == 1:
            vang = evdb.getv('vang')[0]
        trace = wvstadb.load_css(st,et)
        trace.apply_calib()
        trsplice(trace)
        trfilter(trace,fstring)
        rotchan = ('R','T','Z')
        trace.rotate(esaz,vang,rotchan)
        trace.subset('chan =~ /R|T|Z/')
        foundchans = []
        for j in range(trace.nrecs()):
            trace[3] = j
            (sta,chan,ns,sr) = trace.getv('sta','chan','nsamp','samprate')
            stachan = '%s_%s' % (sta,chan)
            ns_tra = 0
            ns_req = 0
            samples = trace.data()
            for k in range(len(samples)):
                ns_tra += ns
                ns_req += int((et-st)*sr)
            if not ns_tra == ns_req:
                logmt(verbose,'Incorrect number of samples for %s, skipping %s' % (stachan,sta))
                for delchan in foundchans:
                    if stachan_traces[delchan]:
                        del stachan_traces[delchan]
                        logmt(debug,'Deleting channel %s' % delchan)
                        counter -= 1
            stachan_traces[stachan] = samples  
            logmt(debug,'Trace extracted for %s' % stachan)
            foundchans.append(stachan)
            counter += 1
        if counter == statmax*3:
            break
    return (stachan_traces)
#}}}
def construct_data_matrix(stachan_traces):
#{{{
# 
# Construct data matrix from dictionary containing trace_objects per sta_chan
# Returns matrix s, which is a 3D matrix, containing all data per channel
#
    numsta = 0
    numchan = 0
    for stachan,tr in sorted(stachan_traces.items()):
        sta,chan = stachan.split('_')
        if dw == True:
            logmt(debug,'Apply distance weighting')
        if chan == 'R':
            numchan += 1
            for j in range(len(tr)):
                s['R'][numsta][j] = tr[j]
        if chan == 'T':
            numchan += 1
            for j in range(len(tr)):
                s['T'][numsta][j] = tr[j]
        if chan == 'Z':
            numchan += 1
            for j in range(len(tr)):
                s['Z'][numsta][j] = tr[j]
        if numchan == 3:
            numchan = 0
            numsta += 1
    return(s)
#}}}
def cross_cor(a,b):
#{{{
#
# Calculate cross correlation between data and accompanying 
# Green's function. Returns the max_cor and the shift
#
    xcor  = np.correlate(a.values(),b.values(),'full')
    pr = len(xcor)/2 
    maxval  = 0
    maxcoef = 0
    for j in range(pr):
        if abs(xcor[j+pr]) > maxval:
            maxval = abs(xcor[j+pr])
            maxcoef = j+1
    return (maxval,maxcoef)
#}}}

orid = configure()
(evdb,evparams) = get_view_from_db(orid)
stachan_traces = get_chan_data(wave_db,evdb)

# Declare matrices (data and Green's functions
s = defaultdict(lambda: defaultdict(defaultdict))
g = defaultdict(lambda: defaultdict(defaultdict))

s = construct_data_matrix(stachan_traces)
if len(s) != 0:
    logmt(debug,'Data matrix S created --> %s stations used' % len(s))

def get_greens_functions(evdb):
#{{{
#
# Opens Green's database and returns a dictionary which describes
# the path to the Green's function per station. Based on azimuth
# and distance. Some settings required from parameter-file for 
# the Green's function creation.
#
    tmp = {}
    green_pf = 'pf/create_green'
    ddist = stock.pfget_string(green_pf,'ddist')
    dazim = stock.pfget_string(green_pf,'dazim')
    ddip = stock.pfget_string(green_pf,'ddip')
    ddepth = stock.pfget_string(green_pf,'ddepth')
    if not os.path.isfile(green_db):
        logmt(3,'Database (%s) does not exist' % green_db)

    gdb = dbopen(green_db,'r')
    gdb = dblookup(gdb,'','moment_tensor_greensfuncs','','')
    for i in range(evdb.nrecs()):

        evdb[3] = i
        (sta,delta,seaz) = evdb.getv('sta','delta','seaz')
        ddist = float(ddist)
        dazim = float(dazim)
        expr = 'delta > %s && delta <= %s ' % (delta-ddist/2,delta+ddist/2)
        expr += '&& azimuth > %s && azimuth <= %s ' % (seaz-dazim/2,seaz+dazim/2)
        subgdb = dbsubset(gdb,expr)
        try:
            subgdb = dbsubset(gdb,expr)
            subgdb[3] = 0
            (dir,file) = subgdb.getv('dir','dfile')       
            dir += '/%s' % file 
            tmp[sta] = dir
        except:
            logmt(verbose,'No Green\'s function found for %s' % sta)
    return(tmp)
#}}}

green_funcs  = {}
#greens_funcs = get_greens_functions(evdb)




# read the Green's function from database......based on distance and depth.....
# Only used in testing, should be removed before committing first working version
######################################################################################################
# Reading test data and Green;s functions....from file                                               #
#{{{
az = [10,40,50]
for j in range(3):
    j += 1
    data =  open('data/data%d' % j,'r')
    samples = []
    for line in data:
        tmp = line.split()
        for samp in tmp:
            samples.append(float(samp))
    data.close()
    for i in range(200):
        s['T'][j-1][i] = samples.pop(0)
    for i in range(200):
        s['R'][j-1][i] = samples.pop(0)
    for i in range(200):
        s['Z'][j-1][i] = samples.pop(0)


greens = []
green = open('data/green','r')
for line in green:
    for j in range(len(line)/12):
        greens.append(line[j*12:j*12+12]) 
green.close()

for i in range(3):
    for k in range(8):
        for j in range(200):
            g[k][i][j] = float(greens[j + k*200])
            if k in [5,6,7]:
                g[k][i][j] *= -1
    if isoflag == 5:
        for k in [8,9]:
            for j in range(200):
                g[k][i][j] = 0.0
    if isoflag == 6:
        for k in [8,9]:
            for j in range(200):
                g[k][i][j] = float(greens[j + k*200])
#}}}
######################################################################################################
num_g = 0
for i in range(9):
    num_g += len(g[i])
if int(num_g/9) != len(s):
    logmt(3,'Number of Green\'s function does not match the number of stations')
else:
    logmt(verbose,'Matrix S (data) and G (synthetics) created')

W = []
for i in range(len(s)):
    W.append(1.0)

def get_time_shift():
#{{{
    for i in range(len(s)):
        shift = 0
        xcor  = 0
        if cross_cor(s['T'][i],g[0][i])[0] > xcor:
            xcor  = cross_cor(s['T'][i],g[0][i])[0]
            shift = cross_cor(s['T'][i],g[0][i])[1]
        if cross_cor(s['T'][i],g[1][i])[0] > xcor:
            xcor = cross_cor(s['T'][i],g[1][i])[0]
            shift = cross_cor(s['T'][i],g[1][i])[1]
        if cross_cor(s['R'][i],g[2][i])[0] > xcor:
            xcor = cross_cor(s['R'][i],g[2][i])[0]
            shift = cross_cor(s['R'][i],g[2][i])[1]
        if cross_cor(s['R'][i],g[3][i])[0] > xcor:
            xcor = cross_cor(s['R'][i],g[3][i])[0]
            shift = cross_cor(s['R'][i],g[3][i])[1]
        if cross_cor(s['R'][i],g[4][i])[0] > xcor:
            xcor = cross_cor(s['R'][i],g[4][i])[0]
            shift = cross_cor(s['R'][i],g[4][i])[1]
        if cross_cor(s['Z'][i],g[5][i])[0] > xcor:
            xcor = cross_cor(s['Z'][i],g[5][i])[0]
            shift = cross_cor(s['Z'][i],g[5][i])[1]
        if cross_cor(s['Z'][i],g[6][i])[0] > xcor:
            xcor = cross_cor(s['Z'][i],g[6][i])[0]
            shift = cross_cor(s['Z'][i],g[6][i])[1]
        if cross_cor(s['Z'][i],g[7][i])[0] > xcor:
            xcor = cross_cor(s['Z'][i],g[7][i])[0]
            shift = cross_cor(s['Z'][i],g[7][i])[1]
        timeshift.append(shift)
    return(timeshift)
#}}}

timeshift = []
timeshift = get_time_shift()
if len(timeshift) != len(s):
    logmt(3,'Error in correlating synthetic and measured data')
else:
    logmt(debug,'Timeshift between data and synthetics computed for %s stations' % len(timeshift))

def matrix_AIV_B(g,s):
#{{{
#
# Construct matrices AIV and B using the data and Green's function matrices
# Return AIV and B for further processingi (inversion)
#
    AJ = defaultdict(dict) 
    trim = 0
    cnt1=cnt2=cnt3 = 0
    trim  = int(len(s['T'][0])*float(trim_value))
    for i in range(len(s)):
        cnt1=cnt2=cnt3
        cnt2 += len(s['T'][0])-trim
        cnt3 += 2*len(s['T'][0])-2*trim
        az[i] *= math.pi/180
        for j in range(len(s['T'][0])-trim):
            AJ[0][cnt1] =  math.sin(2*az[i])*g[0][i][j]/2
            AJ[1][cnt1] = -math.sin(2*az[i])*g[0][i][j]/2
            AJ[2][cnt1] = -math.cos(2*az[i])*g[0][i][j]
            AJ[2][cnt2] = -math.sin(2*az[i])*g[2][i][j]
            AJ[2][cnt3] = -math.sin(2*az[i])*g[5][i][j]
            AJ[3][cnt1] = -math.sin(az[i])*g[1][i][j]
            AJ[3][cnt2] =  math.cos(az[i])*g[3][i][j]
            AJ[3][cnt3] =  math.cos(az[i])*g[6][i][j]
            AJ[4][cnt1] =  math.cos(az[i])*g[1][i][j]
            AJ[4][cnt2] =  math.sin(az[i])*g[3][i][j]
            AJ[4][cnt3] =  math.sin(az[i])*g[6][i][j]
            if isoflag == 5:
                AJ[0][cnt2] = (g[4][i][j])/2 - (math.cos(2*az[i])*g[2][i][j])/2
                AJ[0][cnt3] = (g[7][i][j])/2 - (math.cos(2*az[i])*g[5][i][j])/2
                AJ[1][cnt2] = (g[4][i][j])/2 + (math.cos(2*az[i])*g[2][i][j])/2
                AJ[1][cnt3] = (g[7][i][j])/2 + (math.cos(2*az[i])*g[5][i][j])/2
            if isoflag == 6:
                AJ[0][cnt2] = (g[4][i][j])/6 - (math.cos(2*az[i])*g[2][i][j])/2 + (g[8][i][j])/3
                AJ[0][cnt3] = (g[7][i][j])/6 - (math.cos(2*az[i])*g[5][i][j])/2 + (g[9][i][j])/3
                AJ[1][cnt2] = (g[4][i][j])/6 + (math.cos(2*az[i])*g[2][i][j])/2 + (g[8][i][j])/3
                AJ[1][cnt3] = (g[7][i][j])/6 + (math.cos(2*az[i])*g[5][i][j])/2 + (g[9][i][j])/3
                AJ[5][cnt1] = 0.0
                AJ[5][cnt2] = (g[8][i][j])/3  - (g[4][i][j])/3
                AJ[5][cnt3] = (g[9][i][j])/3 - (g[7][i][j])/3
            cnt1 += 1
            cnt2 += 1
            cnt3 += 1
    AIV = defaultdict(dict) 
    for i in range(isoflag):
        for j in range(isoflag):
            AIV[i][j] = 0.0
    for i in range(5):
        for j in range(5):
            for k in range(cnt3):
                AIV[i][j] += AJ[i][k]*AJ[j][k]

    B = defaultdict(dict) 
    for i in range(isoflag):
        B[i][0] = 0.0
    cnt1 = cnt2 = cnt3 = 0
    tmp = defaultdict(dict) 
    for i in range(len(s)):
    
        cnt1 = cnt2 = cnt3
        cnt2 += len(s['T'][0])-trim
        cnt3 += 2*(len(s['T'][0])-trim)
        for j in range(len(s['T'][0])-trim):
            tmp[cnt1] = s['T'][i][j+timeshift[i]]
            tmp[cnt2] = s['R'][i][j+timeshift[i]]
            tmp[cnt3] = s['Z'][i][j+timeshift[i]]
            cnt1 += 1
            cnt2 += 1
            cnt3 += 1
    
    for i in range(isoflag):
        for j in range(cnt3):
            B[i][0] += AJ[i][j]*tmp[j]
    return(AIV,B)
#}}}

AIV = defaultdict(dict)
B = defaultdict(dict) 
(AIV,B) = matrix_AIV_B(g,s)
if len(AIV) != isoflag or len(AIV[0]) != isoflag:
    logmt(3,'Matrix AIV has dimension [%s,%s] should be [%s,%s]' % (len(AIV),len(AIV),isoflag,isoflag))
elif len(B) != isoflag:
    logmt(3,'Matrix AIV has dimension [%s,i%s] should be [%s,1]' % (len(B),len(B[0]),isoflag))
else:
    logmt(debug,'Matrices AIV and B created and have correct dimensions')

def swap(a,b):
#{{{
#
# Swap function....
#
    tmp = a
    a = b
    b = tmp
    return(a,b)

#}}}
def dyadic(v,n1,n2,c):
#{{{
#
# Compute the dyadic matrix of vector v
#
    tmp = np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
    for i in range(3):
        for j in range(3):
            tmp[i,j] = v[i,n1]*v[j,n2]*c
    return(tmp)
#}}}
def det_solution_vec(AIV,B):
#{{{
#
# Solve the inversion problem, returning the moment tensor
#
    ipiv = defaultdict(dict)
    for i in range(isoflag):
        ipiv[i] = 0
    for i in range(isoflag):
        big = 0.0
        for j in range(isoflag):
            if ipiv[j] != 1:
                for k in range(isoflag):
                    if ipiv[k] == 0:
                        if abs(AIV[j][k]) >= big:
                                big = abs(AIV[j][k])
                                irow = j
                                icol = k
                    elif ipiv[k] > 1:
                        print 'error......1'
        ipiv[icol] += 1
        if not irow == icol:
            for l in range(isoflag):
                (AIV[irow][l],AIV[icol][l]) = swap(AIV[irow][l],AIV[icol][l])
            for l in range(1):
                (B[irow][l],B[icol][l]) = swap(B[irow][l],B[icol][l])
        if AIV[icol][icol] == 0.0:
            print 'error.....2'
        pivinv = 1.0/AIV[icol][icol]
        AIV[icol][icol] = 1.0
        for l in range(isoflag):
            AIV[icol][l] *= pivinv
        for l in range(1):
            B[icol][l] *= pivinv
        for h in range(isoflag):
            if h != icol:
                dum = AIV[h][icol]
                AIV[h][icol] = 0.0
                for l in range(isoflag):
                    AIV[h][l] -= AIV[icol][l]*dum
                for l in range(1):
                    B[h][l] -= B[icol][l]*dum
    if isoflag == 6:
        M = np.matrix([[B[0][0],B[2][0],B[3][0]],[B[2][0],B[1][0],B[4][0]],[B[3][0],B[4][0],B[5][0]]])
    if isoflag == 5:
        M = np.matrix([[B[0][0],B[2][0],B[3][0]],[B[2][0],B[1][0],B[4][0]],[B[3][0],B[4][0],-(B[0][0]+B[1][0])]])
    
#    Coded in the original code by Doug Dreger, however, AIV is not used further in the program, so not needed???
#
#    indc = defaultdict(dict)
#    indr = defaultdict(dict)
#    for i in range(isoflag-1,0,-1):
#        if indr[i] != indc[i]:
#            for j in range(isoflag):
#                (AIV[j][indr[i]],AIV[j][indxc[i]]) = swap(AIV[j][indr[i]],AIV[j][indxc[i]])
    return(M)
#}}}

M = det_solution_vec(AIV,B)
if len(M) != 3:
    logmt(3,'Wrong dimension returned for matrix M after inversion')
else:
    logmt(debug,'Moment tensor matrix M[3,3] created')
gfscale = -1.0e+20
gfscale = -1.0
    
# Get eigenvalues and eigen vectors of moment tensor
def decompose_moment_tensor(M):
#{{{
#
# Decompose moment tensor into eigenvector/values.
# Calculate moment tensor parameters from eigenvalues and vectors. 
# Return M0, Mw, strike, slip, rake and precentages of present
# source characteristics
#
    M *= -1.0e+20
    trace = 0
    for i in range(3):
        trace += M[i,i]
    trace /= 3
    
    for i in range(3):
        M[i,i] -= trace
    miso = np.matrix([[trace,0,0],[0,trace,0],[0,0,trace]])
    (eval,evec) = np.linalg.eig(M)
    
    for i in (0,1):
        k = i
        p = eval[i]
        for j in (1,2):
            if abs(eval[j]) < abs(p):
                k = j
                p = eval[j]
        if k != i:
            eval[k] = eval[i]
            eval[i] = p
            for j in range(3):
                p=evec[j,i]
                evec[j,i] = evec[j,k]
                evec[j,k] = p
    
    f = -eval[0]/eval[2]
    c = eval[2]*(1-2*f)
    a2a2 = dyadic(evec,2,2,c)
    c *= -1
    a1a1 = dyadic(evec,1,1,c)
    mdc = a2a2+ a1a1
    c = 2*eval[2]*f
    a2a2 = dyadic(evec,2,2,c)
    c = -eval[2]*f
    a1a1 = dyadic(evec,1,1,c)
    a0a0 = dyadic(evec,0,0,c)
    mclvd = a2a2+a1a1+a0a0
    
    dd = []
    for i in range(3):
        dd.append(abs(eval[i]))
    for i in range(3):
        for j in range(i,3):
            if dd[j] < dd[i]:
                dd[i],dd[j] = swap(dd[i],dd[j])
    eps = dd[0]/dd[2]
    pcdc = 100*(1-2*eps)
    pcclvd = 200*eps
    
    for i in range(3):
        if evec[2,i] < 0:
            for j in range(3):
                evec[j,i] *= -1
    
    azimuth = []
    plunge = [] 
    for i in range(3):
        if evec[1,i] == 0 and evec[0,i] == 0:
            azimuth.append(0.0)
        else:
            tmp = math.degrees(math.atan2(evec[1,i],evec[0,i]))
            if tmp < 0:
                tmp += 360
            azimuth.append(tmp)
        r = math.sqrt(evec[0,i]*evec[0,i] + evec[1,i]*evec[1,i])
        if evec[2,i] == 0 and r == 0:
            plunge.append(0.0)
        else:
            plunge.append(math.atan2(evec[2,i],r))
    
    axis = []
    axis.append('N')
    if eval[1] > eval[2]:
        axis.append('T')
        axis.append('P')
    else:
        axis.append('P')
        axis.append('T')
    
    p = []
    t = []
    for i in range(3):
        if axis[i] == 'P':
            for j in range(3):
                p.append(evec[j,i])
        elif axis[i] == 'T':
            for j in range(3):
                t.append(evec[j,i])
    con = 1/math.sqrt(2)
    tmp1 = []
    tmp2 = []
    for i in range(3):
        tmp1.append(con*(t[i]+p[i]))
        tmp2.append(con*(t[i]-p[i]))
    u  = np.matrix([[tmp1[0],tmp1[1],tmp1[2]],[tmp2[0],tmp2[1],tmp2[2]]])
    nu = np.matrix([[tmp2[0],tmp2[1],tmp2[2]],[tmp1[0],tmp1[1],tmp1[2]]])
    
    dip = []
    slip = []
    strike = []
    for i in range(2):
        dip.append(math.acos(-nu[i,2]))
        if nu[i,0] == 0 and nu[i,1] == 0:
            strike.append(0)
        else:
            strike.append(math.atan2(-nu[i,0],nu[i,1]))
    for i in range(2):
        sstr = math.sin(strike[i])
        cstr = math.cos(strike[i])
        sdip = math.sin(dip[i])
        cdip = math.cos(dip[i])
        if abs(sdip) > 0:
            lamb = math.asin(-u[i,2])/math.sin(dip[i])
        else:
            if u[i,2] > 0:
                arg = 1
            else:
                arg = -1
            if arg < 0:
                lamb = math.pi
            else:
                lamb = 0
        slamb = math.sin(lamb)
        cdsl = cdip*slamb
        if abs(sstr) > abs(cstr):
            clamb = (u[i,1]+cdsl*cstr)/sstr
        else:
            clamb = (u[i,0]-cdsl*sstr)/cstr
        if slamb == 0 and clamb == 0:
            slip.append(0)
        else:
            slip.append(math.atan2(slamb,clamb))
        if dip[i] > math.pi/2:
            dip[i] = math.pi - dip[i]
            strike[i] += math.pi
            slip[i] = 2*math.pi - slip[i]
        if strike[i] < 0:
            strike[i] += 2*math.pi
        if slip[i] > math.pi:
            slip[i] -= 2*math.pi
    
    m0 = abs(miso[0,0]) + abs(eval[2])
    Mw = math.log10(m0)/1.5 - 10.7
    pciso = abs(miso[0,0])/m0
    return (m0,Mw,strike,dip,slip,pcdc,pcclvd,pciso)
#}}}

strike = []
dip = []
rake = []
(m0,Mw,strike,dip,rake,pcdc,pcclvd,pciso) = decompose_moment_tensor(M)
logmt(verbose,'Decomposition of moment tensor succesfull')

def fitcheck(g,s,W,M,m0,isoflag,timeshift,trim_value,az):
#{{{
# 
# Calculate the variance and variance reduction and flag bad stations.
#
    M /= -1.0e+20
    cnt = 0
    wsum = etot = var = dtot = dvar = 0
    svar = []
    sdpower = []
    for i in range(len(s)):
        dpower = 0
        e = 0
        trimrange = int(len(s['T'][0])*float(trim_value))
        for j in range(trimrange):
            etmp  = s['T'][i][j+timeshift[i]] 
            etmp -= M[0,0]*0.5*g[0][i][j]*math.sin(2*az[i]) 
            etmp += M[1,1]*0.5*g[0][i][j]*math.sin(2*az[i]) 
            etmp += M[0,1]*g[0][i][j]*math.cos(2*az[i]) 
            etmp += M[0,2]*g[1][i][j]*math.sin(az[i]) 
            etmp -= M[1,2]*g[1][i][j]*math.cos(az[i])
            
            e += etmp*etmp
            
            etmp  = s['R'][i][j+timeshift[i]] 
            etmp -= M[0,0]*(0.5*g[4][i][j] - 0.5*g[2][i][j]*math.cos(2*az[i]) + g[8][i][j]/3)
            etmp -= M[1,1]*(0.5*g[4][i][j] + 0.5*g[2][i][j]*math.cos(2*az[i]) + g[8][i][j]/3)
            etmp -= M[2,2]*g[8][i][j]/3 
            etmp += M[0,1]*g[2][i][j]*math.sin(2*az[i]) 
            etmp -= M[0,2]*g[3][i][j]*math.cos(az[i])
            etmp -= M[1,2]*g[3][i][j]*math.sin(az[i])
            
            e += etmp*etmp
            
            etmp  = s['Z'][i][j+timeshift[i]] 
            etmp -= M[0,0]*(0.5*g[7][i][j] - 0.5*g[5][i][j]*math.cos(2*az[i]) + g[9][i][j]/3)
            etmp -= M[1,1]*(0.5*g[7][i][j] + 0.5*g[5][i][j]*math.cos(2*az[i]) + g[9][i][j]/3)
            etmp -= M[2,2]*g[9][i][j]/3 
            etmp += M[0,1]*g[5][i][j]*math.sin(2*az[i]) 
            etmp -= M[0,2]*g[6][i][j]*math.cos(az[i])
            etmp -= M[1,2]*g[6][i][j]*math.sin(az[i])
            
            e += etmp*etmp
            dpower += s['T'][i][j+timeshift[i]]*s['T'][i][j+timeshift[i]]
            dpower += s['R'][i][j+timeshift[i]]*s['R'][i][j+timeshift[i]]
            dpower += s['Z'][i][j+timeshift[i]]*s['Z'][i][j+timeshift[i]]
            cnt += 1
        
        wsum += W[i]
        etot += e
        var += W[i]*e
        dtot += dpower
        dvar += W[i]*dpower
        e /= dpower
        svar.append((1.0 - e)*100.0)
        sdpower.append(dpower)
    pvar = etot/(3*cnt - isoflag - 1.0)
    etot /= dtot
    pvred = (1-etot)*100
    var /= wsum
    dvar /= wsum
    var /= dvar
    var = (1-var)*100
    return(pvar,pvred,var,svar,sdpower)
#}}}
            
svar = []
sdpower = []
(E,VR,VAR,svar,sdpower) = fitcheck(g,s,W,M,m0,isoflag,timeshift,trim_value,az)

if VR < 100:qlt=4
if VR < 80:qlt=3
if VR < 60:qlt=2
if VR < 40:qlt=1
if VR < 20:qlt=0

print 'M0      =',m0
print 'Mw      =',Mw 
print 'Strike  =',math.degrees(strike[0])
print 'Rake    =',math.degrees(rake[0])
print 'Dip     =',math.degrees(dip[0])
print 'Strike2 =',math.degrees(strike[1])
print 'Rake2   =',math.degrees(rake[1])
print 'Dip2    =',math.degrees(dip[1])
print 'Pdc     =',pcdc
print 'Pclvd   =',pcclvd
print 'Piso    =',pciso
for i in range(len(svar)):
    print 'Stat(%d) = %s  %s' % (i,svar[i],sdpower[i])
print 'Var     =',E
print 'VR      =',VR,'(UNWEIGHTED)'
print 'VR      =',VR,'(WEIGHTED)'
print 'Var/Pdc =',E/pcdc
print 'Quality =',qlt

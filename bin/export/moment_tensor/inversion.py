from __main__ import *      # Get all the libraries from parent

class MomentTensor():
    # {{{
    """
    Class for building moment tensors and doing the inversion

    """

    def __init__(self, distance_weighting, isoflag, trim_value, verbosity=0):
        """Initialize"""
        self.distance_weighting = distance_weighting
        self.isoflag = isoflag
        self.trim_value = trim_value
        self.verbosity = verbosity

    def construct_data_matrix(self, stachan_traces):
        """Construct data matrix from dictionary 
        containing trace_objects per sta_chan
        Returns rearranged 3D matrix, 
        containing all data per channel
        """
        this_dict = defaultdict(lambda: defaultdict(defaultdict))

        if self.verbosity > 0:
            logmt(1, 'Constructing data matrix from trace objects')
        numsta = 0
        numchan = 0
        for stachan, tr in sorted(stachan_traces.items()):
            sta, chan = stachan.split('_')

            if self.distance_weighting == True:
                logmt(1, 'Apply distance weighting')

            if chan == 'R':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['R'][numsta][j] = tr[j]
            if chan == 'T':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['T'][numsta][j] = tr[j]
            if chan == 'Z':
                numchan += 1
                for j in range(len(tr)):
                    this_dict['Z'][numsta][j] = tr[j]

            if numchan == 3:
                numchan = 0
                numsta += 1

        return this_dict

    def plot_cross_cor(self, a, b, xcor, components):
        """Plot the cross 
        correlation for 
        debugging purposes
        """
        stacode, component, shift = components
        title_string = "Station %s. Filtered data for %s. Shift = %s" % (stacode, component, shift)
        fig = plt.figure()
        fig.subplots_adjust(left=0.04, right=0.96, hspace=0.3)
        ax = fig.add_subplot(311)
        ax.set_title(title_string)
        ax.plot(self.normalize(a), 'b-')
        ax.plot(self.normalize(b), 'g-')
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.legend(('Station', 'Synthetic'), 'upper right', prop={'size':'10'})
        # Print shifted time series
        bx = fig.add_subplot(312)
        bx.set_title('Shifted time series')
        if shift < 0:
            try:
                [a.insert(0, 0) for x in range(abs(shift))]
            except Exception, e:
                logmt(1, "   - plot_cross_cor(): Could not insert values at start of real data. Error: %s" % e)
                print a
        elif shift > 0:
            try:
                [b.insert(0, 0) for x in range(abs(shift))]
            except Exception, e:
                logmt(1, "   - plot_cross_cor(): Could not insert values at start of synthetic data. Error: %s" % e)
                print b
        bx.plot(self.normalize(a), 'r-')
        bx.plot(self.normalize(b), 'c-')
        plt.setp(bx.get_yticklabels(), visible=False)
        bx.legend(('Station', 'Synthetic'), 'upper right', prop={'size':'10'})
        # Print xcor
        cx = fig.add_subplot(313)
        cx.set_title('Cross correlation results')
        cx.plot(xcor)
        plt.setp(cx.get_yticklabels(), visible=False)
        cx.set_ylabel('X-cor.')
        cx.set_xlabel('Point shift (in samples)')
        plt.show()
        return True
 
    def cross_cor(self, a, b, component=False, stacode=False):
        """Calculate cross correlation 
        between data and accompanying 
        Green's function. Returns the 
        max_cor and the shift
        """
        xcor = np.correlate(a, b, 'full')
        maxval = np.amax(xcor)
        maxshift = np.argmax(xcor)
        shift = maxshift - (len(xcor)-1)/2
        if self.verbosity > 1:
            logmt(1, '  - cross_cor(): Component: %s, max val: %s, timeshift (in samples): %s' % (component, maxval, shift))
            self.plot_cross_cor(a, b, xcor, [stacode, component, shift])
        return maxval, shift

    def normalize(self, data_as_list):
        """Determine the 
        normalization factor
        for all the data
        """
        try:
            normalizer = max([abs(x) for x in data_as_list])
        except Exception, e:
            logmt(1, "  - Exception encountered: %s" % e)
        else:
            if self.verbosity > 1:
                logmt(1, "  - Normalization factor: %s" % normalizer)
            return [x / normalizer for x in data_as_list]

    def get_time_shift(self, data_dict, greens_dict, stacode, size):
        """Get time shift for each station
        using cross correlation of station 
        data points and the Green's Function
        data points. Make sure both matrices
        are the same size for the comparison.
        Return a list of tuples (stacode, timeshift).

        !!! NOTE: They don't have to be the same
                  size for a successful cross 
                  correlation
        """
        if self.verbosity > 0:
            logmt(1, '  - Get the time shift for all components & return the max cross-cor and shift in samples')
        '''
        !!! NOTE: Loop over each stations data points.
                  Remember that the Green's Function 
                  matrix is just two dimensional - not 
                  replicating the number of stations
                  (3D matrix). RLN (2011-12-02)
        '''
        size = len(greens_dict['XDS'])
        max_shift = 0
        max_xcor = 0
        for c in greens_dict:
            if c == 'REX' or c == 'ZEX':
                continue
            else:
                if c[:1] == 'T':
                    data_key = 'T'
                elif c[:1] == 'X':
                    data_key = 'R'
                elif c[:1] == 'Z':
                    data_key = 'Z'
                this_xcor, this_shift = self.cross_cor(data_dict[data_key][0:size], greens_dict[c][0:size], c, stacode)
                if this_xcor > max_xcor:
                    max_xcor = this_xcor
                    max_shift = this_shift
        if self.verbosity > 0:
            logmt(1, '  - Cross-correlation value: %s, timeshift value: %s' % (max_xcor, max_shift))
        return max_xcor, max_shift

    # def matrix_AIV_B(self, dict_s, dict_g, list_az, timeshift, size):
    def matrix_AIV_B(self, dict_s, dict_g, ev2sta, size):
        """INVERSION ROUTINE

        Construct matrices AIV and B using 
        the data and Green's function matrices
        Return AIV and B for further processing 
        This is to normalize the two signals
        """
        if self.verbosity > 0:
            logmt(1, 'INVERSION ROUTINE: Construct matrices AIV and B using the data and Greens function matrices')

        AJ = defaultdict(dict)
        trim = 0
        cnt1 = cnt2 = cnt3 = 0

        if self.verbosity > 0:
            # logmt(1, ' - Timeshift: %s' % timeshift)
            # logmt(1, ' - Azimuthal distances tuple: %s' % list_az)
            logmt(1, ' - Number of stations used in calculation: %s' % len(dict_s))

        # Iterate over number of stations
        for i in range(len(ev2sta)):
            if self.verbosity > 0:
                logmt(1, '  - Working on inversion for station (%s)' % ev2sta[i][0])
            # Allocate the dictionary space for each section
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            # Convert azimuth to radians
            this_azimuth = ev2sta[i][1] * math.pi/180
            '''
            !!! NOTE: len(dict_s['T'][0]) is the number of 
                      the datapoints for the first station. 
                      Use var size instead because of timeshift.
                      RLN (2011-12-06)

            '''
            for j in range(size - trim):
                # Mxx term
                AJ[0][cnt1] =  math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Myy term
                AJ[1][cnt1] = -math.sin(2*this_azimuth)*dict_g[i]['TSS'][j]/2

                # Mxy term
                AJ[2][cnt1] = -math.cos(2*this_azimuth)*dict_g[i]['TSS'][j]
                AJ[2][cnt2] = -math.sin(2*this_azimuth)*dict_g[i]['XSS'][j]
                AJ[2][cnt3] = -math.sin(2*this_azimuth)*dict_g[i]['ZSS'][j]

                # Mxz term
                AJ[3][cnt1] = -math.sin(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[3][cnt2] =  math.cos(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[3][cnt3] =  math.cos(this_azimuth)*dict_g[i]['ZDS'][j]

                # Myz term
                AJ[4][cnt1] =  math.cos(this_azimuth)*dict_g[i]['TDS'][j]
                AJ[4][cnt2] =  math.sin(this_azimuth)*dict_g[i]['XDS'][j]
                AJ[4][cnt3] =  math.sin(this_azimuth)*dict_g[i]['ZDS'][j]

                # Vary the other values depending on isoflag value
                if self.isoflag == 5:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/2 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/2 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2
                if self.isoflag == 6:
                    # Mxx term
                    AJ[0][cnt2] = (dict_g[i]['XDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[0][cnt3] = (dict_g[i]['ZDD'][j])/6 - (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Myy term
                    AJ[1][cnt2] = (dict_g[i]['XDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['XSS'][j])/2 + (dict_g[i]['REX'][j])/3
                    AJ[1][cnt3] = (dict_g[i]['ZDD'][j])/6 + (math.cos(2*this_azimuth)*dict_g[i]['ZSS'][j])/2 + (dict_g[i]['ZEX'][j])/3
                    # Explosion-related values
                    AJ[5][cnt1] = 0.0
                    AJ[5][cnt2] = (dict_g[i]['REX'][j])/3  - (dict_g[i]['XDD'][j])/3
                    AJ[5][cnt3] = (dict_g[i]['ZEX'][j])/3 - (dict_g[i]['ZDD'][j])/3
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1
        if self.verbosity > 0:
            logmt(1, ' - Created matrix AJ with length: %s' % len(AJ))
            logmt(1, ' - Final counts: cnt1=%s, cnt2=%s, cnt3=%s' % (cnt1, cnt2, cnt3))
            logmt(1, ' - Now apply the timeshift if needed')

        '''
        !!! NOTE: Make placeholder AIV 
                  (or reset it like Dreger?)
                  RLN (2011-12-02)
        !!! NOTE: Attempt to replace this with scipy.linalg.lstsq
                  http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html#scipy.linalg.lstsq
                  RLN & MK (2012-02-07)
        '''

        AIV = defaultdict(dict) 
        for i in range(self.isoflag):
            for j in range(self.isoflag):
                AIV[i][j] = 0.0
        # Compute AtA
        for i in range(5):
            for j in range(5):
                for k in range(cnt3):
                    AIV[i][j] += AJ[i][k] * AJ[j][k]

        B = defaultdict(dict) 
        for i in range(self.isoflag):
            B[i][0] = 0.0

        cnt1 = cnt2 = cnt3 = 0
        tmp = defaultdict(dict) 

        # Iterate over the number of stations
        for i in range(len(ev2sta)):
            if self.verbosity > 0:
                logmt(1, '  - Applying timeshift to station (%s): Trim: %s, timeshift: %s' % (ev2sta[i][0], trim, ev2sta[i][3]))
            '''
            !!! FIX: Direct copy of Dregers code - don't need to do this
                     RLN (2011-12-08)
            Pre-populate dictionary of values
            |cnt1           |cnt2        |cnt3       |
            |*****----------|*****-------|*****------|
            '''
            cnt1 = cnt2 = cnt3
            cnt2 += size - trim
            cnt3 += 2*size - 2*trim
            for j in range(size - trim):
                correction_iter = j + ev2sta[i][3]
                try:
                    tmp[cnt1] = dict_s[i]['T'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s with index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['T'])
                try:
                    tmp[cnt2] = dict_s[i]['R'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s with index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['R'])
                try:
                    tmp[cnt3] = dict_s[i]['Z'][correction_iter]
                except KeyError as k:
                    logmt(1, 'KeyError (%s) for station: %s at index: %s' % (k, ev2sta[i][0], correction_iter))
                    pprint(dict_s[i]['Z'])
                cnt1 += 1
                cnt2 += 1
                cnt3 += 1

        # Calculate Righthand Side
        for i in range(self.isoflag):
            for j in range(cnt3):
                B[i][0] += AJ[i][j] * tmp[j]

        if len(AIV) != self.isoflag or len(AIV[0]) != self.isoflag:
            logmt(3, 'Matrix AIV has dimension [%s,%s] should be [%s,%s]' % (len(AIV), len(AIV), self.isoflag, self.isoflag))
        elif len(B) != self.isoflag:
            logmt(3, 'Matrix AIV has dimension [%s,i%s] should be [%s,1]' % (len(B), len(B[0]), self.isoflag))
        elif self.verbosity > 0:
            logmt(1, 'Matrices AIV and B created and have correct dimensions')

        return AIV, B

    def swap(self, a, b):
        """Remap list
        """
        tmp = a
        a = b
        b = tmp
        return(a, b)

    def dyadic(self, v, n1, n2, c):
        """Calculate the dyadic 
        matrix of eigenvectors v
        """
        if self.verbosity > 0:
            logmt(1, ' - Compute the dyadic matrix of vector v')
        tmp = np.matrix([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        for i in range(3):
            for j in range(3):
                tmp[i,j] = v[i, n1]*v[j, n2]*c
        return tmp

    def determine_solution_vector(self, dict_AIV, dict_B):
        """Determine the solution vector
        Solve the inversion problem, 
        returning the moment tensor
        From Dreger:
            Call Bobs MT decomposition routines
            The minus one is needed to map Helmbergers convention into Aki's
            Jost and Hermann (1989) state that AKI's convention is -1*LANGSTONS

        !!! NOTE: Matt thinks we can replace this with a library in SciPy
            RLN (2012-02-07)
        """
        if self.verbosity > 0:
            logmt(1, 'Solve the inversion problem and return the moment tensor')
        ipiv = defaultdict(dict)
        for i in range(self.isoflag):
            ipiv[i] = 0
        for i in range(self.isoflag):
            big = 0.0
            for j in range(self.isoflag):
                if ipiv[j] != 1:
                    for k in range(self.isoflag):
                        if ipiv[k] == 0:
                            if abs(dict_AIV[j][k]) >= big:
                                    big = abs(dict_AIV[j][k])
                                    irow = j
                                    icol = k
                        elif ipiv[k] > 1:
                            logmt(3, 'determine_solution_vector(): ERROR... 1')
            ipiv[icol] += 1
            if not irow == icol:
                for l in range(self.isoflag):
                    (dict_AIV[irow][l],dict_AIV[icol][l]) = swap(dict_AIV[irow][l],dict_AIV[icol][l])
                for l in range(1):
                    (dict_B[irow][l],dict_B[icol][l]) = swap(dict_B[irow][l],dict_B[icol][l])
            if dict_AIV[icol][icol] == 0.0:
                logmt(3, 'determine_solution_vector(): ERROR... 2')
            pivinv = 1.0/dict_AIV[icol][icol]
            dict_AIV[icol][icol] = 1.0
            for l in range(self.isoflag):
                dict_AIV[icol][l] *= pivinv
            for l in range(1):
                dict_B[icol][l] *= pivinv
            for h in range(self.isoflag):
                if h != icol:
                    dum = dict_AIV[h][icol]
                    dict_AIV[h][icol] = 0.0
                    for l in range(self.isoflag):
                        dict_AIV[h][l] -= dict_AIV[icol][l]*dum
                    for l in range(1):
                        dict_B[h][l] -= dict_B[icol][l]*dum
        if self.isoflag == 6:
            M = np.matrix([[dict_B[0][0], dict_B[2][0], dict_B[3][0]], [dict_B[2][0], dict_B[1][0], dict_B[4][0]], [dict_B[3][0], dict_B[4][0], dict_B[5][0]]])
        if self.isoflag == 5:
            M = np.matrix([[dict_B[0][0], dict_B[2][0], dict_B[3][0]], [dict_B[2][0], dict_B[1][0], dict_B[4][0]], [dict_B[3][0], dict_B[4][0], -(dict_B[0][0]+dict_B[1][0])]])
        
        """
        # Coded in the original code by Doug Dreger, however, AIV is not used further in the program, so not needed???
        indc = defaultdict(dict)
        indr = defaultdict(dict)
        for i in range(self.isoflag-1,0,-1):
            if indr[i] != indc[i]:
                for j in range(self.isoflag):
                    (dict_AIV[j][indr[i]],dict_AIV[j][indxc[i]]) = swap(dict_AIV[j][indr[i]],dict_AIV[j][indxc[i]])
        """
        if len(M) != 3:
            logmt(3, 'Wrong dimension returned for matrix M after inversion')
        elif self.verbosity > 0:
            logmt(1, ' - Moment tensor matrix M[3,3] created')
        return M

    def decompose_moment_tensor(self, matrix_M):
        """Decompose moment tensor into eigenvector/values.
        Calculate moment tensor parameters from eigenvalues and vectors. 
        Return M0, Mw, strike, slip, rake and precentages of present
        source characteristics
        From Dreger source: fmap_subs_linux.f
        Mxx = matrix_M[0, 0]
        Mxy = matrix_M[0, 1]
        Mxz = matrix_M[0, 2]
        Myy = matrix_M[1, 1]
        Myz = matrix_M[1, 2]
        Mzz = matrix_M[2, 2]
        """
        if self.verbosity > 0:
            logmt(1, 'Decompose moment tensor into eigenvector/values')
        matrix_M *= -1.0e+20
        trace = 0
        for i in range(3):
            trace += matrix_M[i,i]
        trace /= 3
        
        for i in range(3):
            matrix_M[i,i] -= trace
        miso = np.matrix([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
        eval, evec = np.linalg.eig(matrix_M)
        
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
        a2a2 = self.dyadic(evec, 2, 2, c)
        c *= -1
        a1a1 = self.dyadic(evec, 1, 1, c)
        mdc = a2a2+ a1a1
        c = 2*eval[2]*f
        a2a2 = self.dyadic(evec, 2, 2, c)
        c = -eval[2]*f
        a1a1 = self.dyadic(evec, 1, 1, c)
        a0a0 = self.dyadic(evec, 0, 0, c)
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

        if self.verbosity > 0:
            logmt(1, ' - Decomposition of moment tensor succesful!')

        return m0, Mw, strike, dip, slip, pcdc, pcclvd, pciso

    def fitcheck(self, dict_g, dict_s, matrix_M, m0, ev2sta, size):
        """Calculate the variance, 
        variance reduction
        and flag bad stations
        """
        if self.verbosity > 0:
            logmt(1, 'Calculate the variance, variance reduction & flag bad stations')
        matrix_M /= -1.0e+20
        cnt = 0
        wsum = etot = var = dtot = dvar = 0
        svar = []
        sdpower = []
        # Loop over the number of stations
        for i in range(len(ev2sta)):
            dpower = 0
            e = 0
            trimrange = int(int(size) * float(self.trim_value))
            for j in range(trimrange):
                correction_iter = j + ev2sta[i][3]
                etmp  = dict_s[i]['T'][correction_iter] 
                etmp -= matrix_M[0,0]*0.5*dict_g[i]['TSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp += matrix_M[1,1]*0.5*dict_g[i]['TSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp += matrix_M[0,1]*dict_g[i]['TSS'][j]*math.cos(2*ev2sta[i][1]) 
                etmp += matrix_M[0,2]*dict_g[i]['TDS'][j]*math.sin(ev2sta[i][1]) 
                etmp -= matrix_M[1,2]*dict_g[i]['TDS'][j]*math.cos(ev2sta[i][1])
                
                e += etmp*etmp
                
                etmp  = dict_s[i]['R'][correction_iter] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[i]['XDD'][j] - 0.5*dict_g[i]['XSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['REX'][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[i]['XDD'][j] + 0.5*dict_g[i]['XSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['REX'][j]/3)
                etmp -= matrix_M[2,2]*dict_g[i]['REX'][j]/3 
                etmp += matrix_M[0,1]*dict_g[i]['XSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp -= matrix_M[0,2]*dict_g[i]['XDS'][j]*math.cos(ev2sta[i][1])
                etmp -= matrix_M[1,2]*dict_g[i]['XDS'][j]*math.sin(ev2sta[i][1])
                
                e += etmp*etmp
                
                etmp  = dict_s[i]['Z'][correction_iter] 
                etmp -= matrix_M[0,0]*(0.5*dict_g[i]['ZDD'][j] - 0.5*dict_g[i]['ZSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['ZEX'][j]/3)
                etmp -= matrix_M[1,1]*(0.5*dict_g[i]['ZDD'][j] + 0.5*dict_g[i]['ZSS'][j]*math.cos(2*ev2sta[i][1]) + dict_g[i]['ZEX'][j]/3)
                etmp -= matrix_M[2,2]*dict_g[i]['ZEX'][j]/3 
                etmp += matrix_M[0,1]*dict_g[i]['ZSS'][j]*math.sin(2*ev2sta[i][1]) 
                etmp -= matrix_M[0,2]*dict_g[i]['ZDS'][j]*math.cos(ev2sta[i][1])
                etmp -= matrix_M[1,2]*dict_g[i]['ZDS'][j]*math.sin(ev2sta[i][1])
                
                e += etmp*etmp
                dpower += dict_s[i]['T'][correction_iter]*dict_s[i]['T'][correction_iter]
                dpower += dict_s[i]['R'][correction_iter]*dict_s[i]['R'][correction_iter]
                dpower += dict_s[i]['Z'][correction_iter]*dict_s[i]['Z'][correction_iter]
                cnt += 1
            
            # wsum += list_W[i]
            wsum += 1
            etot += e
            # var += list_W[i]*e
            var += 1*e
            dtot += dpower
            # dvar += list_W[i]*dpower
            dvar += 1*dpower
            e /= dpower
            svar.append((ev2sta[i][0], (1.0 - e)*100.0))
            sdpower.append(dpower)
        pvar = etot/(3*cnt - self.isoflag - 1.0)
        etot /= dtot
        pvred = (1-etot)*100
        var /= wsum
        dvar /= wsum
        var /= dvar
        var = (1-var)*100
        return pvar, pvred, var, svar, sdpower

    def quality_check(self, vr):
        """Check the quality of the result
        """
        if vr < 100:
            qlt = 4
        if vr < 80:
            qlt = 3
        if vr < 60:
            qlt = 2
        if vr < 40:
            qlt = 1
        if vr < 20:
            qlt = 0
        return qlt
    # }}}


"""
If we call this script directly, then output help.
"""
if __name__ == "__main__":
#{{{
    print "Moment Tensor Inversion Library:"
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

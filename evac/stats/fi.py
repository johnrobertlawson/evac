import os
import pdb

class FI:
    def __init__(self,xa,xfs,thresholds,decompose=True,
                    neighborhoods='auto',temporal_windows=[1,],
                    tidxs='all'):
        """Fractional Ignorance.

        The "continuous" part of this scheme is not over thresholds in the raw
            data (e.g., mm for rainfall), but for fractional values. There's
            a good reason for being able to choose thresholds: the user
            may only be interested in running FI for flash-flood cases, for
            example.

        Returns:
            dictionary of the form:
            
            dict[threshold][neighborhood][temporal_window][timestep][score] = x

            There is N.NaN returned for the first and last (W-1)/2 timesteps
            for temporal window W. The score keys are strings of FI, UNC, 
            REL, RES.

        Args:
            xa: 3-D observation array of size (ntimes,nlat, nlon).
            xfs: 4-D forecast array of size (N,ntimes,nlat,nlon) where N is 
                ensemble size and ntimes is the number of time steps.
            thresholds (tuple,list): thresholding to apply to input data.
            neighborhoods: if 'auto', use smallest to largest spatial
                kernels. These can accessed after processing with the
                self.neighborhoods attribute, or the keys of the 
                returned dictionary.
            temporal_windows (tuple,list): by default, no temporal
                windowing is performed. If a list is given, elements
                must be integers of odd numbers only.
            decompose (bool): if True, return a named tuple of 
                uncertainty (UNC), reliability (REL), and resolution (RES)
                dictionaries. The value of FI is the sum of all three.
            tidxs: if "all", then all timesteps are used. Otherwise, only
                these indices are used.


        Returns:
            Either FI or the decomposition as a dictionary.
        """

        assert xfs.shape[1:] == xa.shape 
        self.xa = xa
        self.xfs = xfs

        self.nens = self.xfs.shape[0]
        self.probthreshs = N.linspace(0,1,nens+1)

        check = _assert_list(thresholds)
        self.thresholds = thresholds
        
        if neighborhoods == 'auto':
            # The biggest box we can fit in the domain
            maxwindow = N.floor*(min(self.nlats,self.nlons)/2)
            neighborhoods = N.arange(1,maxwindow,2)
        check = _assert_list(neighborhoods)
        self.neighborhoods = neighborhoods

        check = _assert_list(temporal_windows)
        self.temporal_windows = temporal_windows

        check = _assert_list(tidxs)
        self.tidxs = tidxs
        
        self.decompose = decompose

        # Create results dictionary
        self.results = {th: {n: { tw: {t: {} for t in self.tidxs} for tw in 
                        self.temporal_windows} for n in self.neighborhoods}
                        for th in self.thresholds}

        namedtup = self.parallelise()
        return self.results

    def parallelise(self):
        # Loop over all thresholds, neighbourhood sizes, and temporal windows
        loop = self.generate_loop()

        if self.ncpus == 1:
            for i in loop:
                thresh, neigh, tempwindow, tidx = i
                result = self.compute_fi(i)
                self.update_results(result)
        else:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.compute_fi,loop)
            for result in results:
                self.update_results(result)
                # scores, thresh, neigh, tempwindow, tidx = tup
                # self.result[thresh][neigh][tempwindow][tidx] = fi

    def update_results(self,result):
        scores, thresh, neigh, tempwindow, tidx = result
            for score in ("REL","RES","UNC","FI"):
                self.result[thresh][neigh][tempwindow][tidx][score] = result[score]
                                # thresh=thresh,neigh=neigh,
                                # tempwindow=tempwindow,tidx=tidx,)
            return


    def compute_fi(self,i):
        """
        Consider coding a continuous version (loop over all possible 
        fractions possible in ensemble members) or simply ranked
        (where thresholding is done earlier).
        """
        thresh, neigh, tempwindow, tidx = i
    
        # First threshold the data into a binary array
        # Whether each gridpoint is over (1) or under (0) the threshold
        Im = N.where(self.xfs[tidx,:,:,:] > thresh, 1, 0)
        Io = N.where(self.xa[tidx,:,:] > thresh, 1, 0)

        # Next, apply the 3D kernel
        kernel = N.ones([neigh,neigh,tempwindow])
        M = signal.fftconvolve(Im,kernel)
        O = signal.fftconvolve(Io,kernel)

        # Next, loop over all possible fractional values
        # Calculate REL, RES, UNC for each 
        num = (neigh**2)+1
        REL = N.zeros([num-1])
        RES = N.copy(REL)
        UNC = N.copy(REL)

        # As we're using fractional arrays, the distance is constant,
        # see eq. 27/28 in Toedter and Ahrens 2010
        distance = 1/(neigh**2)

        for fidx, fracval in enumerate(N.linspace(0,1,num=num)):
            REL[fidx] = self.compute_rel(M,O,f)
            RES[fidx] = self.compute_res(M,O,f)
            UNC[fidx] = self.compute_unc(O,f)

        rel = distance*N.sum(REL)
        res = distance*N.sum(RES)
        unc = distance*N.sum(UNC)

        fi = rel - res + unc

        scores = dict(REL=rel,RES=res,UNC=unc,FI=fi)
        return scores, thresh, neigh, tempwindow, tidx
            
    def eq_14_16(self,M,O,f,func):   
        """ Eqs 14 and 16 are very similar in Toedter and Ahrens
            (2012, MWR), so this combines the logic.
        """
        Mf = utils.exceed_probs_2d(M,f,overunder='over',fmt='decimal'):
        results = N.zeros_like(self.probthreshs)
        for iidx,yi in enumerate(self.probthreshs):
            widx = N.where(Mf == yi)

            # pyi = N.count(Mf[widx] == yi)
            pyi = self.compute_pyi(Mf,yi)
            zi = self.compute_zi(O[widx],f,yi)
            z = self.compute_z(O,f)
            results[iidx] = func(pyi=pyi,yi=yi,zi=zi,z=z)
        return N.sum(results)

    def compute_pyi(probs,yi):
        """ Freq. of probs equal to yi.
        """
        pyi = (N.where(probs == yi).size)/probs.size
        return pyi

    def compute_z(self,O,f):
        """

        z is the prob of observed occurrence in the
        sample (0-1 frequency).
        """
        # Not sure which z or zi is correct.
        # zi = N.count(o[widx] > f)
        # zi = N.count(O[widx] == f)

        # z = N.count(O > f)
        # z = N.count(O > f)/O.size
        # z = N.count(O==f)/O.size
        # z = N.count(O[widx])

        z = (N.where(O == f).size)/O.size
        return z

    def compute_zi(self,o,f,yi):
        """ 
        Args:
            o: subset of fractional array O that were
                points with forecast prob yi.
            f: fraction value of interest (later, integrated)
            yi: forecast prob of interest (later, integrated)

        A well calibrated ensemble has zi == yi.
        It is the observed (0-1) frequency of a ob
        that occurred for prob forecasts for the ith
        probability (between 0 and 1).

        For a given yi (e.g. 70% or 0.7), what percentage
        of those points verified?

        Conditional freq of occurrence on all occasions where
        yi was forecasted.
        """
        zi = (N.where(o > f).size)/o.size
        return zi

    def compute_rel(self,M,O,f):
        """ 
        Loop over all bins of probability
        For each, if the M grid point

        Generate probs of exceedence for the given f.
        bar zi = yi ... bar zi is frequency of 
        obs for a given bin yi.
        """
        def rel_eq(pyi,yi,zi,**kwargs):
            """ Eq. 14 in Toedter and Ahrens 2012 MWR,
            without the summation.
            """
            rel = pyi * (zi*N.log(zi/yi) + 
                    (1-zi)*N.log((1-zi)/(1-yi)))
            return rel
        return self.eq_14_16(M,O,f,rel_eq)

    def compute_res(self,M,O,f):
        def res_eq(pyi,zi,z,**kwargs):
            """ Eq. 16 in Toedter and Ahrens 2012 MWR,
            without the summation.
            """
            res = pyi * (zi*N.log(zi/z) + 
                    (1-zi)*N.log((1-zi)/(1-z)))
            return res
        return self.eq_14_16(M,O,f,res_eq)

    def compute_unc(self,O,f):
        """ 
        z is the frequency of f in the sample O.
        """
        # z = N.count(O==f)/O.size
        z = self.compute_z(O,f)
        UNC = -z*N.log(z) - (1-z)*N.log(1-z)
        return UNC
    
    def generate_loop(self):
        for thresh, neigh, tempwindow, tidx in itertools.product(self.thresholds,
            self.neighborhoods, self.temporal_windows,self.tidxs):
            yield thresh, neigh, tempwindow, tidx

    @staticmethod
    def _assert_list(l):
        if isinstance(l,(N.array,list,tuple,set)):
            return True
        else:
            raise Exception("This is not a valid list, tuple etc.)



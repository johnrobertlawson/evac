import os
import pdb

class FIbST:
    def __init__(self,xa,xfs,thresholds,decompose=True,
                    neighborhoods='auto',temporal_windows=[1,],
                    tidxs='all'):
        """Fractional Ignorance by Scale/Time (FIbST).

        Returns:
            dictionary of the form:
            
            dict[threshold][neighborhood][temporal_window][timestep][score] = x

            There is N.NaN returned for the first and last (W-1)/2 timesteps
            for temporal window W. The score keys are strings of FIBST, UNC, 
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
                dictionaries. The value of FIbST is the sum of all three.
            tidxs: if "all", then all timesteps are used. Otherwise, only
                these indices are used.


        Returns:
            Either FIbST or the decomposition as a dictionary.
        """

        assert xfs.shape[1:] == xa.shape 
        self.xa = xa
        self.xfs = xfs

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
                result = self.compute_fibst(i)
                self.update_results(result)
        else:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.compute_fibst,loop)
            for result in results:
                self.update_results(result)
                # scores, thresh, neigh, tempwindow, tidx = tup
                # self.result[thresh][neigh][tempwindow][tidx] = fibst

    def update_results(self,result):
        scores, thresh, neigh, tempwindow, tidx = result
            for score in ("REL","RES","UNC","FIBST"):
                self.result[thresh][neigh][tempwindow][tidx][score] = result[score]
                                # thresh=thresh,neigh=neigh,
                                # tempwindow=tempwindow,tidx=tidx,)
            return

    def compute_fibst(self,i):
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

        fibst = rel - res + unc

        scores = dict(REL=rel,RES=res,UNC=unc,FIBST=fibst)
        return scores, thresh, neigh, tempwindow, tidx
            

    def compute_rel(self,M,O,f):
        """ 
        Loop over all bins of probability
        For each, if the M grid point

        Generate probs of exceedence for the given f.
        bar zi = yi ... bar zi is frequency of 
        obs for a given bin yi.
        """
        Mf = probexceed(f)
        for iidx,yi in self.probbins:
            widx = N.where(Mf == yi)
            pyi = N.count(Mf[widx] == yi)
            zi = N.count(O[widx] > f)
        pass

    def compute_res(self,M,O,f):
        pass

    def compute_unc(self,O,f):
        """ 
        z is the frequency of f in the sample O.
        """
        z = N.count(O==f)/O.size
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



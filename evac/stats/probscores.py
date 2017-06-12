import pdb
import itertools

import numpy as N

"""Probabilistic scores, using Buizza 2001, MWR.
"""

class ProbScores:
    def __init__(self,og=None,pfg=None,xa=None,xfs=None):
        """
        Pick from:
        self.og         :   True/False field (1D)
        self.pfg        :   Prob. forecast of event (0.0-1.0)
        --- or ---
        self.xa         :   observation field (2D)
        self.xfs        :   forecast fields (3D)
        """
        self.og = og
        self.pfg = pfg
        self.xa = xa
        self.xfs = xfs
        if (og is not None) and (pfg is not None):
            assert self.og.size == self.pfg.size
        elif (xa is not None) and (xfs is not None):
            assert self.xa.shape == self.xfs[0,:,:].shape
        else:
            raise Exception

    def compute_briar(self,):
        Ng = self.og.size
        self.og = self.og.flatten()
        self.pfg = self.pfg.flatten()
        BS = (1/Ng) * N.sum((self.pfg-self.og.astype(float))**2)
        return BS
    
    def compute_bss(self):
        """TODO
        """
        BS = self.compute_briar()
        BSS = 0
        return

    def compute_rps(self,thresholds):
        """
        Arguments:
            thresholds  : 1D list/array, and in same order
                            of first dimension of pfg.
            self.pfg    : 3D array: [probability index,lat,lon]
                    
        Note:
        RPSg is the grid-point RPS.
        RPS is the area-average (ranked probability score)
        """
        assert len(self.og.shape) == 3
        Ng = self.og[0,:,:].size
        # evs = N.array((sorted(thresholds) )
        Jev = len(thresholds)
        RPSgs = N.zeros_like(self.og)

        # RPS for every grid point, for each threshold
        for thidx,th in enumerate(thresholds):
            RPSgs[thidx,:,:] = (N.sum(pfg[:th+1,:,:],axis=0)-
                                N.sum(og[:th+1,:,:],axis=0))**2

        # Sum up over all thresholds for RPS at each point
        RPSg = N.sum(RPSgs,axis=0)

        # Average over all grid points for RPS.
        RPS1 = (1/Ng) * N.sum(RPSg)
        RPS2 = N.mean(RPSg)
        assert RPS1 == RPS2
        return RPS2


    def compute_poutl(self,E,obs):
        """Percentage of outliers.

        E       :   Ensemble class instance
        obs     :   observation array the same 2D size.
        """
        pass

    # def compute_crps(self,xfs,xa):
    def compute_crps(self,mean=True):
        """
        Some inspiration from:
        https://github.com/TheClimateCorporation/properscoring.git
        CRPS (Continuous Ranked Probability Score)
        Measures distance between prob forecast and truth.
        

        Needs observation array (m x n)
        Needs ensemble forecast array (e x m x n) where
                e = number of ensemble members

        xa     :   observation 
        xfs    :   forecasts for N ensemble members
        rho    :   PDF forecast
        """
        from scipy.stats import norm

        # from scipy.integrate import quad,trapz

        # To do: get rid of NaNs or stupid numbers?

        # x is x-axis (spacings for integration)
        # x = N.linspace(self.xfs.min(),xfs.max(),50)

        s1,s2 = self.xa.shape
        ss = self.xa.size
        count = 0
        crps = N.zeros_like(self.xa)
        print("Starting CRPS calculation over the whole grid.")
        for x,y in itertools.product(range(s1),range(s2)):
            count += 1
            if not count % 1000:
                print("{} of {}".format(count,ss))

            # Px is the % that self.xa will be smaller than x.
            # Sort ensemble forecasts 
            xs = N.sort(self.xfs[:,x,y])
            # Generate CDF
            # Px = scipy.stats.norm.pdf(x)
            Px = norm.cdf(xs)
           
            # CDF of observation (0 or 1)
            Pax = self.heaviside(xs-self.xa[x,y])


            # dx = N.zeros(len(xs)+1)
            dx = N.zeros(len(xs))

            # WHICH:   ???
            # dx[1:-1] = N.diff(xs)
            # dx[:-1] = N.diff(xs)
            dx[1:] = N.diff(xs)

            # Px[0] should be 0
            # Px[-1] should be 1.

            integrand = dx*(Px-Pax)**2

            crps[x,y] = N.sum(integrand)
            # if N.any(xs > 0.0):
                # pdb.set_trace()

        if mean:
            return crps.mean()
        else:
            return crps
        # def integrand(Px0,Pax0):
            # return (Px0 - Pax0)**2
        # def crps(Px0,Pax0):
            # return quad(integrand,-N.inf,N.inf,args=(Px0,Pax0))
        # crps_vec = N.vectorize(crps)
        # CRPS = (Px(x) - Pax)**2
        # CRPS = crps_vec(Px,Pax)
        pdb.set_trace()
         
    def heaviside(self,x):
        # if x < 0:
            # return 0
        # else:
            # return 1
        return x >= 0

    def cdf_(self,):
        pass

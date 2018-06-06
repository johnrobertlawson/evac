import itertools
import pdb

import numpy as N
import scipy as S

class ProbScores:
    """Probabilistic scores, using Buizza 2001, MWR.

    Note:
        User must specify (og and pfg) or (xa and xfs).

    Args:
        self.og         :   (numpy.array) - True/False field (1D)
        self.pfg        :   Prob. forecast of event (0.0-1.0)
        self.xa         :   observation field (2D)
        self.xfs        :   forecast fields (3D)
    """
    def __init__(self,og=None,pfg=None,xa=None,xfs=None,
                    allow_1d=False):
        self.og = og
        self.pfg = pfg
        self.xa = xa
        self.xfs = xfs
        if (og is not None) and (pfg is not None):
            assert self.og.size == self.pfg.size
        elif (xa is not None) and (xfs is not None):
            if allow_1d:
                assert self.xfs.ndim == 1
            else:
                assert self.xa.shape == self.xfs[0,:,:].shape
        else:
            raise Exception

    def compute_briar(self,decompose=False):
        """ Briar Score.

        Todo:
            * Consider moving nested functions to methods or classfunctions.
        """

        def compute_briar_reliability(yi,pyi,zi):
            """ Reliability is the correspondence of observed to forecast frequencies. 

            Args:
                yi: set of possible values
                pyi: frequency with which each forecast value from yi is forecasted
                zi: probability of occurrence given forecast yi
            """
            return N.sum(pyi*((yi-zi)**2))

        def compute_briar_resolution(pyi,zi,z):
            """ Resolution is ability of a forecast to distinguish between events.

            Args:
                pyi:  ...
                zi: ...
                z: climatological probability of observed event in same sample
            """
            return N.sum(pyi*((zi-z)**2))

        def compute_briar_uncertainty(z):
            """ Uncertainty is the variance of observations.

            This is the score achieved if the climatological probability is
            forecasted each time (persistence?)

            Args:
                z: ...
            """
            return z(1-z)

        Ng = self.og.size
        self.og = self.og.flatten()
        self.pfg = self.pfg.flatten()
        if decompose:
            z = 0
            pyi = 0
            yi = 0
            zi = 0

            UNC = compute_briar_uncertainty(z)
            REL = compute_briar_reliability(yi,pyi,zi)
            RES = compute_briar_resolution(pyi,zi,z)
            
            BS = REL - RES + UNC
        else:
            BS = (1/Ng) * N.sum((self.pfg-self.og.astype(float))**2)
        return BS

    def compute_ignorance(self,):
        """ Compute the Ignorance score.

        Using Toeter and Ahrens 2012 paper, MWR

        Todo:
            * Decompose

        Args:
        """
        pass 

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

    def rank_hist(obdata,ensdata):
        """
        Wrapper for rank histogram over 3D.
        """
        pass


    def compute_poutl(self,E,obs):
        """Percentage of outliers.

        E       :   Ensemble class instance
        obs     :   observation array the same 2D size.
        """
        pass

    def compute_crps(self,threshs,mean=True,QC=False):
        """
        Notes:
            This method needs the followed during class instantiation:
            * Needs `xa`, observation array (m x n)
            * Needs `xfs`, ensemble forecast array (e x m x n) where
                e = number of ensemble members
            * Px is the probability that ob is less or equal to fcst.
                Loop through all probs from 0 to 100 %

        Todo:
            * Decompose CRPS into reliability etc, see Hersbach paper.

        Args:
            threshs: An array of the levels to apply to the data
            mean (bool): If `True`, return a mean score averaged over the domain
            QC (bool): If `True`, set anything less than zero to zero.

        Discretising using:
        https://www.kaggle.com/c/how-much-did-it-rain#evaluation
        """
        # Number of ensemble members
        nens, ny, nx = self.xfs.shape
        # Ensemble is ranked smallest to largest at each grid pt.
        xfs_sort = N.sort(self.xfs,axis=0)

        if QC:
            xfs_sort[xfs_sort < 0] = 0
            # xfs_sort[xfs_sort == N.nan] = 0

        c = N.zeros_like(threshs)
        for thidx,th in enumerate(threshs):
            # Sum of 2D field of probs/heaviside
            Px = N.sum(N.greater(xfs_sort,th),axis=0)/nens
            Hx = (self.heaviside(th-self.xa)).astype(int)
            c[thidx] = N.sum((Px - Hx) **2)
        # Average over all thresholds and grid points
        crps = (N.sum(c))/(nens*Px.size)

        return crps

        # Decompose? What about uncertainty

    def compute_crps_old(self,mean=True,debug=False):
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
            if (not count % 1000) and debug:
                print("{} of {}".format(count,ss))

            # Px is the % that self.xa will be smaller than x.
            # Sort ensemble forecasts
            xs = N.sort(self.xfs[:,x,y])
            # Generate CDF
            # Px = scipy.stats.norm.pdf(x)
            Px = norm.cdf(xs)

            # CDF of observation (0 or 1)
            Hv = self.heaviside(xs-self.xa[x,y])


            # dx = N.zeros(len(xs)+1)
            dx = N.zeros(len(xs))

            # WHICH:   ???
            # dx[1:-1] = N.diff(xs)
            # dx[:-1] = N.diff(xs)
            dx[1:] = N.diff(xs)

            # Px[0] should be 0
            # Px[-1] should be 1.

            method = 3
            if method == 1: # METHOD 1
                integrand = dx*((Px-Hv)**2)
                crps[x,y] = N.sum(integrand)

            elif method == 2:
                # METHOD 2
                from scipy.integrate import quad
                def integrand(argin):
                # def integrand(Px0,Hv0):
                    Px0, Hv0 = argin
                    return (Px0 - Hv0)**2

                # def get_crps(Px0,Hv0):
                    # return quad(integrand,-N.inf,N.inf,args=(Px0,Hv0))
                # get_crps_vec = N.vectorize(get_crps)
                # crps[x,y] = get_crps_vec(Px,Hv)
                crps[x,y] = quad(integrand,-N.inf,N.inf,args=(Px,Hv))



            elif method == 3:
                #https://www.mathworks.com/matlabcentral/fileexchange/47807-continuous-rank-probability-score?requestedDomain=www.mathworks.com
                problist = N.arange(xs.size)/(xs.size)
                problist2 = problist ** 2
                problistm2 = (1-problist)**2

                ob = self.xa[x,y]
                fcst = xs
                idxs = N.where(fcst <= ob)[0]
                if len(idxs) > 0:  # when obs > fcst
                    idx = idxs[-1]
                    crps_left = 0
                    if idx > 0:
                        fcst_left = fcst[:idx]
                        dx_left = N.diff(fcst_left)
                        p_left = problist2[:idx-1]
                        # crps_left = p_left * dx_left
                        crps_left = N.dot(p_left,dx_left)
                    if ob < fcst[-1]:
                        fcst_right = fcst[idx:]
                        dx_right = N.diff(fcst_right)
                        if len(dx_right) == 0:
                            crps_right = 0
                        else:
                            p_right = problistm2[idx:-1]
                            # crps_right = p_right * dx_right
                            crps_right = N.dot(p_right, dx_right)
                        # CDF crossing obs
                        crps_centreleft = problist2[idx]**(ob-fcst[idx])
                        crps_centreright = problistm2[idx]**(fcst[idx+1]-ob)
                        crps_vals = crps_left + crps_right + crps_centreleft + crps_centreright
                    else: # if obs are > all fcst
                        crps_right_outside = 1**(2*(ob-fcst[-1]))
                        crps_vals = crps_right_outside + crps_left
                else: # obs is < all fcst
                    dx_right = N.diff(fcst)
                    p_right = problistm2[:-1]
                    crps_right =N.dot(p_right,dx_right)
                    # crps_right = p_right*dx_right
                    crps_left_outside = 1**(2*(fcst[0]-ob))
                    crps_vals = crps_left_outside + crps_right
                crps[x,y] = crps_vals

        if mean:
            return crps.mean()
        else:
            return crps

    def heaviside(self,x):
        # if x < 0:
            # return 0
        # else:
            # return 1
        return x >= 0

    def cdf_(self,):
        pass


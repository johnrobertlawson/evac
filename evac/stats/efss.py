import os
import pdb
import itertools
import multiprocessing

import numpy as N
from scipy import signal

class EFSS:
    """ Ensemble fractions skill score from Duc et al (2013, Tellus).

    Optimisation helped by Faggian et al (2015, Mausam).

    Todo:
        * Implement faster algorithms than FFT from Faggian paper.
        * Round the fraction fields off, convert to fractions, bound to 0-1.
        * Address scores that aren't realistic - bug somewhere

    Notes:
        * The time, lat, and lon dimensions of fcst4d and obs3d must match.
    Args:
        fcst4d (numpy.ndarray): ensemble x time x lat x lon array of fcst
        obs3d (numpy.ndarray): time x lat x lon array of obs
        threshs (tuple, list): thresholds to evaluate FSS
        spatial_ns, temporal_ms (tuple, list): neighbourhood sizes

    Returns:
        A dictionary of FSS scores per threshold and neighbourhood
        eFSS[spatial_n][temporal_m][threshold] = float
    """
    def __init__(self,fcst4d,obs3d,threshs,spatial_ns,temporal_ms,
                    ncpus=1):
        """

        self.Im is threshold x ensemble member x time x lat x lon
        """
        # INSTANCE VARIABLE ASSIGNMENT
        self.fcst4d = fcst4d
        self.obs3d = obs3d
        self.threshs = threshs
        self.spatial_ns = spatial_ns
        self.temporal_ms = temporal_ms
        self.ncpus = ncpus

        self.nns = len(self.spatial_ns)
        self.nms = len(self.temporal_ms)
        self.nthreshs = len(self.threshs)

        # CHECK SANITY
        self.nmems, self.ntimes, self.nlats, self.nlons = self.fcst4d.shape
        assert self.ntimes == self.obs3d.shape[0]
        assert self.nlats == self.obs3d.shape[1]
        assert self.nlons == self.obs3d.shape[2]
        
        # COMPUTE BINARY FIELDS
        self.Im = self.compute_Im()
        self.Io = self.compute_Io()

        self.eFSS = self.parallelise()

    def generate_loop(self):
        for n,m,th in itertools.product(self.spatial_ns,self.temporal_ms,
                                    self.threshs):
            yield n,m,th

    def parallelise(self):
        eFSS = {n: {m: {th: {} for th in self.threshs} for m in self.temporal_ms} 
                    for n in self.spatial_ns}
        itr = self.generate_loop()
        if self.ncpus == 1:
            for i in itr:
                n,m,th = i
                eFSS[n][m][th] = self.compute_efss(i)
        else:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.compute_efss,itr)
                # eFSS = results
            # pdb.set_trace()
            for tup in results:
                score, n, m, th = tup
                eFSS[n][m][th] = score
        return eFSS

    def compute_efss(self,itr,m=None,n=None,th=None):
        """
        Args:
            function_of: if None, return FBS left as a function
                of spatial/time windows and threshold.

        Note:
            * MO_diff has size: n x m x thresh
        FBS has size: thresh
        """
        if itr is not None:
            n,m,th = itr
        M = self.compute_M(m=m,n=n,th=th)
        O = self.compute_O(m=m,n=n,th=th)
        MO_diff = (M - O)**2
        FBS_ref = M**2 + O**2
        eFSS = 1 - (N.nanmean(MO_diff) / N.nanmean(FBS_ref))
        return eFSS, n, m, th

    def compute_M(self,m,n,th):
        """
        n = spatial window size
        m = temporal window size
        th = threshold
        """
        conv = N.ones([self.nmems, m, n, n])
        M = signal.fftconvolve(self.Im[th],conv)
        return M

    def compute_O(self,m,n,th):
        conv = N.ones([m, n, n])
        O = signal.fftconvolve(self.Io[th],conv)
        return O

    def compute_Im(self):
        """ 
        """
        Im = dict()
        for th in self.threshs:
            Im[th] = self.fcst4d > th
        return Im

    def compute_Io(self):
        """
        Returns:
            Io is a binary 4d array: threshold x time x lat x lon
        """
        Io = dict()
        for th in self.threshs:
            Io[th] = self.obs3d > th
        return Io


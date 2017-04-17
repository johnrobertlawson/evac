"""Relative Operative Statistics.

Mixin of deterministic and probabilistic scores. Can be used to generate a
2x2 contingency table for ensembles.
"""
import pdb

import numpy as N

from .detscores import DetScores
from .probscores import ProbScores

class ROC(DetScores,ProbScores):
    def __init__(self,nens=None,pf=None,ob=None,ensemble=None,):
        """
        Args:
            nens        :   number of ensemble members
            pf          :   probability of an event (1d array)
            ob (bool)   :   whether it happened (1d array)
            thresh      :   probability threshold(s)
                                These must be 1/nens, 2/nens... nens/nens.
                                Turned off for now.
                                Default is all thresholds.
            ensemble    :   Ensemble instance.

            obs_array   :   2D array (lats,lons) for verification
        """
        if nens is None and pf is None and ob is None:
            # Logic here to compute the arrays
            pass
        #5D array (members, times, lvs, lats, lons)
        #Should be (nens,1,1,nlats,nlons) in shape


        self.nens = nens
        self.pf = pf
        self.ob = ob

        self.a,self.b,self.c,self.d = self.compute_roc()

    def compute_roc(self):
        thresholds = []
        X = N.zeros([self.nens])
        Y = N.zeros([self.nens])
        for nidx,n in enumerate(range(1,self.nens+1)):
            thresholds.append(n/self.nens)
            if self.pf[nidx] < thresholds[-1]:
                if self.ob[nidx]:
                    X[nidx] =+ 1
                else:
                    Y[nidx] =+ 1

        a = N.zeros([len(thresholds)])
        b = N.copy(a)
        c = N.copy(a)
        d = N.copy(a)
        for thidx,th in enumerate(thresholds):
            a[thidx] = N.sum(X[thidx+1:])
            b[thidx] = N.sum(Y[thidx+1:])
            c[thidx] = N.sum(X[:thidx+1])
            d[thidx] = N.sum(Y[:thidx+1])

        self.thresholds = thresholds
        return a,b,c,d

    def compute_roca(self):
        """Area under ROC curve.

        TODO."""
        FAR = self.compute_falsealarm()
        HR = self.compute_hitrate()
        ROCA = N.trapz(HR,x=FAR)
        return ROCA

    def compute_rocas(self):
        """ROC area skill score.
        """
        pass

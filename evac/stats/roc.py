"""Relative Operative Statistics.

Mixin of deterministic and probabilistic scores. Can be used to generate a
2x2 contingency table for ensembles.
"""
import pdb

import numpy as N

from .detscores import DetScores
from .probscores import ProbScores
from evac.stats.misc_stats import compute_contingency

class ROC(DetScores,ProbScores):
    def __init__(self,nens=None,pf=None,ob=None,ensemble=None,
                    obs_array=None,thresholds=None,overunder=None):
        """
        Args:
            nens        :   number of ensemble members
            pf          :   probability of an event (1d array)
            ob (bool)   :   whether it happened (1d array)
            ensemble    :   Ensemble instance.
            obs_array   :   2D array (lats,lons) for verification
            threshold   :   float to create binary condition
            overunder   :   operator for threshold.
        """
        if nens is None and pf is None and ob is None:
            # Logic here to compute the arrays
            pass
        #5D array (members, times, lvs, lats, lons)
        #Should be (nens,1,1,nlats,nlons) in shape


        self.nens = nens
        self.pf = pf
        self.ob = ob

        # self.a,self.b,self.c,self.d = self.compute_roc()
        self.a,self.b,self.c,self.d = self.compute_roc_contingency()
        # self.a,self.b,self.c,self.d = compute_contingency(forecast,obs,thresh,overunder)

    def compute_roc_contingency(self):
        # ws contains all prob thresholds
        ws = N.linspace(0,1,self.nens+1)
        a = N.zeros([len(ws)])
        b = N.copy(a)
        c = N.copy(a)
        d = N.copy(a)

        for i,w in enumerate(ws):
            a[i],b[i],c[i],d[i] = compute_contingency(self.pf,self.ob,w,'over')
        return a,b,c,d

    def _compute_roc(self):
        """
        w is the probability threshold list.
        """
        w = []
        X = N.zeros([self.nens]) 
        Y = N.zeros([self.nens])
        for nidx,n in enumerate(range(1,self.nens+1)):
            w.append(n/self.nens)
            if self.pf[nidx] < w[-1]:
                if self.ob[nidx]:
                    X[nidx] =+ 1
                else:
                    Y[nidx] =+ 1

        a = N.zeros([len(w)])
        b = N.copy(a)
        c = N.copy(a)
        d = N.copy(a)
        for thidx,th in enumerate(w):
            a[thidx] = N.sum(X[thidx+1:])
            b[thidx] = N.sum(Y[thidx+1:])
            c[thidx] = N.sum(X[:thidx+1])
            d[thidx] = N.sum(Y[:thidx+1])

        self.w = w
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

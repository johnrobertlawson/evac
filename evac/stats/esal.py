""" SAL extended to ensembles.

Logic taken from Radanovics et al 2018, WAF.
"""

import os
import pdb

import numpy as N

from evac.stats.sal import SAL
from evac.stats.probscores import ProbScores

class ESAL(SAL):
    """ Ensemble version of SAL.

    Note:
        * We use "S", "A", etc for the component names, neglecting the "e" prefix
            used in the Radanovics paper.
        * We assume the observation is just a single object, not an ensemble.
            Because I'm confused why you'd need an ensemble of verification.
    Args:
        OBC : ObjectBased instance for the "Control", usually observations.
        OBMd : dictionary of ObjectBased instances for the "model", usually
            ensemble forecasts. The key is the ensemble name; value is
            an ObjectBased instance.
    """
    def __init__(self,OBC,OBMd,dx=1,dy=1,footprint=None,thresh=None):
        self.OBC = OBC
        self.OBMd = OBMd
        self.footprint = footprint
        self.thresh = thresh

        self.member_names = sorted(self.OBMd.keys())
        self.nmems = len(self.member_names)
        self.nlats, self.nlons = self.OBMd[self.member_names[0]].raw_data.shape
        self.dx = dx
        self.dy = dy

        self.d = self.compute_d()
        self.A = self.compute_amplitude()
        self.L, self.L1, self.L2 = self.compute_location()
        self.S = self.compute_structure()

    def compute_amplitude(self,):
        """ Overridden from parent.
        """
        # First, compute ensemble average stats
        rr_all = N.ones([self.nmems])
        for n, (memname, memOB) in enumerate(self.OBMd.items()):
            rr_all[n] = N.mean(memOB.raw_data)
        rr_mod = N.mean(rr_all,axis=0)
        rr_obs = N.mean(self.OBC.raw_data)

        A = (rr_mod - rr_obs)/(0.5 * (rr_mod + rr_obs))
        return A

    def compute_L1(self):
        # A list of vectors
        xrr_all = N.ones([self.nmems,2])
        for n, (memname, memOB) in enumerate(self.OBMd.items()):
            # xrr_all[n,:] = memOB.objects['x_CoM']
            xrr_all[n,:] = memOB.x_CoM
        # Vector mean?
        xrr_mod = N.mean(xrr_all,axis=0)

        # vector subtraction
        dist_km = self.vector_diff_km(xrr_mod,self.OBC.x_CoM)
        # dist_km = self.vector_diff_km(xrr_mod,self.OBC.objects['x_CoM'])
        L1 = dist_km/self.d
        print(("L1 = {0}".format(L1)))
        return L1

    def compute_L2(self):
        """Defined as double the CRPS of weighted average distances between
        the centre of mass of an object and the total centre of mass divided by d.
        """
        rmodd = N.ones([self.nmems])
        for n, (memname, memOB) in enumerate(self.OBMd.items()):
            rmodd[n] = self.compute_r(memOB)/self.d
        
        robs = self.compute_r(self.OBC)
        robsd = robs/self.d

        PS = ProbScores(xa=robsd,xfs=rmodd,allow_1d=True)
        print("Using CRPS thresholds for dBZ!")
        crps = PS.compute_crps(threshs=N.arange(0,95,0.1))
        return 2*crps

    def compute_location(self):
        L1 = self.compute_L1()
        L2 = self.compute_L2()
        L = L1 + L2
        return L, L1, L2

    def compute_structure(self):
        V_all = N.ones([self.nmems])
        for n, (memname, memOB) in enumerate(self.OBMd.items()):
            V_all[n] = self.compute_V(memOB)
        V_mean = N.mean(V_all)

        V_obs = self.compute_V(self.OBC)

        S = (V_mean - V_obs)/(0.5 * (V_mean + V_obs))
        return S

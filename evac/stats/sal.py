import pdb

import numpy as N
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

from evac.datafiles.wrfout import WRFOut
from evac.datafiles.radar import Radar
from evac.stats.objectbased import ObjectBased

class SAL:
    """ Structure-Amplitude-Location verification method.

    From Wernli et al. with modifications from Lawson & Gallus.

    Plotting can be done by `evac.plot.salgraph` by passing in
    each SAL instance (Will this be too memory intensive?)
    or maybe just lists of S, A, L1, L2 components. Then there's
    eSAL.

    Args:
        OBC: ObjectBased instance for the "Control", usually observation.
        OBM: ObjectBased instance for the "Model", usually a forecast
    """
    def __init__(self,OBC,OBM,thresh,vrbl=None,utc=None,lv=None,
                    accum_hr=None,radar_datadir=None,
                    footprint=500,dx=1,dy=1,
                    autofactor=1/15):
        self.OBC = OBC
        self.OBM = OBM
        # Check grid are identical before proceeding
        assert self.OBC.raw_data.ndim == 2
        assert self.OBC.raw_data.shape == self.OBM.raw_data.shape

        self.thresh = thresh
        self.vrbl = vrbl
        self.utc = utc
        self.footprint = footprint
        self.autofactor = autofactor

        # if dx and dy are not specified, we just use unit distance
        self.dx = dx
        self.dy = dy
        self.nlats, self.nlons = self.OBC.raw_data.shape

        # Create some aliases for object data?

        # TODO: QC so all negative values set to 0?
        # self.C['data'][self.C['data']<0] = 0
       
        self.d = self.compute_d()
        self.A = self.compute_amplitude()
        self.L, self.L1, self.L2 = self.compute_location()
        self.S = self.compute_structure()

        print(("S = {0}    A = {1}     L = {2}".format(self.S,self.A,self.L)))


    def compute_d(self):
        xside = self.dx * self.nlons
        yside = self.dy * self.nlats
        d = N.sqrt(xside**2 + yside**2)
        return d

    def compute_amplitude(self,):
        A = (N.mean(self.OBM.raw_data) - N.mean(self.OBC.raw_data))/(
                0.5*(N.mean(self.OBM.raw_data) + N.mean(self.OBC.raw_data)))
        return A

    def compute_location(self,):
        L1 = self.compute_L1()
        L2 = self.compute_L2()
        L = L1 + L2
        return L, L1, L2

    def compute_L1(self,):
        # vector subtraction
        dist_km = self.vector_diff_km(self.OBM.x_CoM,self.OBC.x_CoM)
        L1 = dist_km/self.d
        print(("L1 = {0}".format(L1)))
        return L1

    def vector_diff_km(self,v1,v2):
        # From grid coords to km difference
        dist_gp = N.subtract(v1,v2)
        dist_km = self.dx * N.sqrt(dist_gp[0]**2 + dist_gp[1]**2)
        return dist_km

    def compute_L2(self,):
        rc = self.compute_r(self.OBC)
        rm = self.compute_r(self.OBM)
        L2 = 2*(N.abs(rc-rm)/self.d)
        print(("L2 = {0}".format(L2)))
        return L2
    
    def compute_r(self,OB):
        Rn_sum = 0
        for k,v in list(OB.objects.items()):
            if N.isnan(v['Rn']) or N.isnan(v['CoM'][0]):
                print("NaN detected in r computation. Ignoring.")
            else:
                Rn_sum += v['Rn'] * self.vector_diff_km(OB.x_CoM,v['CoM'])
        try:
            r = Rn_sum / OB.R_tot
        except ZeroDivisionError:
            r = 0
        return r

    def compute_structure(self):
        Vm = self.compute_V(self.OBM)
        Vc = self.compute_V(self.OBC)
        try:
            S = (Vm - Vc)/(0.5*(Vm+Vc))
        except ZeroDivisionError:
            S = 0
        return S

    def compute_V(self,OB):
        Vn_sum = 0
        for k,v in list(OB.objects.items()):
            Vn_sum += v['Rn'] * v['Vn']
        try:
            V = Vn_sum / OB.R_tot
        except ZeroDivisionError:
            V = 0
        return V

    def __str__(self):
        return ("This SAL object scores the follows: \n"
                    "S = {:.3f} \n A = {:.3f} \n L = {:.3f}".format(
                    self.S, self.A, self.L))

""" Identify objects in an array.
"""
import pdb
import os

import numpy as N
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

from evac.datafiles.wrfout import WRFOut
from evac.datafiles.radar import Radar

class ObjectBased:
    """ Object-based verification superclass.

    Todo:
        * Incorporate SAL logic
        * Be able to look up properties of objects
        * Link with other scores, e.g. CRPS
        * Link with other ObjectBased() instances (e.g. for SAL)

    Args:
        arr: 2D numpy array to identify objects in
    """
    def __init__(self,arr,thresh='auto',footprint=1/15):
        self.raw_data = arr
        # self.metadata = {}
        self.thresh = thresh
        self.footprint = footprint
        self.identify_objects()


    def identify_objects(self,):
        # self.metadata['Rmax'] = N.max(self.raw_data)
        self.Rmax = N.max(self.raw_data)

        if thresh == 'auto':
            self.Rstar = self.footprint * self.Rmax
            # self.metadata['Rstar'] = self.footprint * self.metadata['Rmax']
        else:
            self.Rstar = float(self.thresh)

        self.object_operators()

    def object_operators(self):
        """ This needs some explaining.
        """
        mask = N.copy(self.raw_data)
        mask[self.raw_data < self.Rstar] = False
        mask[self.raw_data >= self.Rstar] = True
        labeled, num_objects = ndimage.label(mask)

        sizes = ndimage.sum(mask, labeled, list(range(num_objects+1)))

        masksize = sizes < nsize
        remove_pixel = masksize[labeled]
        labeled[remove_pixel] = 0

        labels = N.unique(labeled)
        label_im = N.searchsorted(labels, labeled)

        # dic['objects'] = {}
        self.objects = {}

        # Total R for objects
        R_objs_count = 0

        for ln,l in enumerate(labels):
            cy, cx = ndimage.measurements.center_of_mass(self.raw_data,labeled,l)
            if ln == 0:
                # This is the centre of mass for the entire domain.
                self.x_CoM = (cx,cy)
            else:
                self.objects[l] = {}
                self.objects[l]['CoM'] = (cx,cy)
                self.objects[l]['Rn'] = ndimage.sum(self.raw_data,labeled,l)
                self.objects[l]['RnMax'] = ndimage.maximum(self.raw_data,labeled,l)
                self.objects[l]['Vn'] = self.objects[l]['Rn']/self.objects[l]['RnMax']
                R_objs_count += self.objects[l]['Rn']

        self.R_tot = R_objs_count
        self.obj_array = labeled
        return

    def active_px(self,fmt='pc'):
        """ Return number of pixels included in objects.
        Args:
            fmt (bool): if True, returns the active pixel count
                expressed as percentage.
        """
        active_px = N.count_nonzero(self.obj_array)
        tot_px = self.obj_array.size

        if fmt == 'pc':
            return (active_px/(tot_px*1.0))*100.0
        else:
            return active_px, tot_px

    def get_diagnostics(self):
        """ Compute various characteristics for each of the objects identified.
        """
        pass

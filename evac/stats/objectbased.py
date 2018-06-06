""" Identify objects in an array.
"""
import pdb
import os

import numpy as N
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

from evac.datafiles.wrfout import WRFOut
from evac.plot.scales import Scales
from evac.datafiles.radar import Radar
from evac.plot.figure import Figure
from evac.plot.birdseye import BirdsEye

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
    def __init__(self,arr,thresh='auto',footprint=500,f=1/15):
        assert arr.ndim == 2

        self.raw_data = arr
        # self.metadata = {}
        self.thresh = thresh
        self.footprint = int(footprint)
        self.f = f
        self.identify_objects()


    def identify_objects(self,):
        # self.metadata['Rmax'] = N.max(self.raw_data)
        self.Rmax = N.max(self.raw_data)

        if self.thresh == 'auto':
            self.Rstar = self.f * self.Rmax
            # self.metadata['Rstar'] = self.f * self.metadata['Rmax']
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

        masksize = sizes < self.footprint
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

    def get_updraughts(self,W_data,):
        """ Returns stats on object updraughts.

        Max updraught is computed at all levels at all grid points
        in the object.

        Args:
            W_data: 3-D array of W, where the lats/lons correspond directly
                to the raw_data passed into self originally.
        """
        if W_data.ndim == 3:
            W_data = N.max(W_data,axis=0)
        assert W_data.ndim == 2
        nobjs = len(self.objects)
        updraught_dict = dict()
        updraught_arr = N.zeros([nobjs])

        for nidx,(n,o) in enumerate(self.objects.items()):
            # The grid points of each object is labelled
            # in self.obj_array.
            grid_2d_idx = N.where(self.obj_array == n)
            arr_obj = W_data[grid_2d_idx[0],grid_2d_idx[1]]
            max_w = N.max(arr_obj)

            # Data in a dictionary and array
            updraught_dict[n] = max_w
            updraught_arr[nidx] = max_w
        
        return updraught_arr, updraught_dict

    def plot(self,fpath,fmt='default',W=None,vrbl='REFL_comp',
                # Nlim=None,Elim=None,Slim=None,Wlim=None):
                ld=None,lats=None,lons=None):
        """ Plot basic quicklook images.

        Setting fmt to 'default' will plot raw data,
        plus objects identified.
        """
        if ld is None:
            ld = dict()
        nobjs = len(self.objects)

        if fmt == 'default':
            F = Figure(ncols=2,nrows=1,figsize=(8,4),
                        fpath=fpath)
            # F.W = W
            with F:
                ax = F.ax[0]
                # Plot raw array
                BE = BirdsEye(ax=ax,fig=F.fig)

                # Discrete colormap
                import matplotlib as M
                cmap_og = M.cm.get_cmap('tab20')
                # cmap_colors = [cmap_og(i) for i in range(cmap_og.N)]
                color_list = cmap_og(N.linspace(0,1,nobjs))
                # cmap = M.colors.ListedColormap(M.cm.tab20,N=len(self.objects))
                cmap = M.colors.LinearSegmentedColormap.from_list('discrete_objects',color_list,nobjs)
                # bounds = N.linspace(0,nobjs,nobjs+1)
                # norm = M.colors.BoundaryNorm(bounds,cmap_og.N)
                masked_objs = N.ma.masked_less(self.obj_array,1)
                
                BE.plot2D(plottype='pcolormesh',data=masked_objs,save=False,
                            cb='horizontal',
                            #clvs=N.arange(1,nobjs),
                            W=W,
                            cmap=cmap,mplkwargs={'vmin':1},**ld,lats=lats,lons=lons)

                ax = F.ax[1]
                S = Scales(vrbl)
                BE = BirdsEye(ax=ax,fig=F.fig)
                BE.plot2D(data=self.raw_data,save=False,
                            W=W,
                            cb='horizontal',lats=lats,lons=lons,
                            cmap=S.cm,clvs=S.clvs,**ld)
        return

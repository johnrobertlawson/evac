""" Postage stamps or thumbnails.
"""
import os
import pdb

import numpy as N

from evac.plot.birdseye import BirdsEye
from evac.plot.figure import Figure
from evac.datafiles.radar import Radar
from evac.plot.scales import Scales

class ThumbNails(Figure):
    """ Plot multiple subplots.
    """
    def __init__(self,rowscols,fpath=None,outdir=None,fname='test_thumbs.png',
                    verif=True,proj='merc'):
        self.proj = proj
        if not fpath:
            fpath = os.path.join(outdir,fname)
        super().__init__(fpath,proj=proj,**self.setup_figure(rowscols))
        self.naxes = self.ax.size

    def setup_figure(self,rowscols):

        init_dict = dict(
            nrows = rowscols[0],
            ncols = rowscols[1],
            figsize = (10,8),
            )
        return init_dict

    def plot_verif(self,data,grid,lats=None,lons=None,subplotn=0,
                    save=False,utc=None,vrbl=None,scales=None,):
        """ Plot verification postage stamp.
        """
        ax = self.ax.flatten()[subplotn]
        cmap,clvs = self.get_scales(vrbl=vrbl)

        # Set up panel
        BE = BirdsEye(fpath=None,ax=ax,fig=self.fig,grid=grid,proj=self.proj)

        # mplkwargs = dict('cmap':cmap, 'levels':clvs)
        BE.plot2D(data,save=False,drawcounties=True,
                    cb=False,
                    # mplkwargs=mplkwargs,
                    # cmap=cmap,cen_lat=cen_lat,cen_lon=cen_lon,W=W,**ld,
                    cmap=cmap,levels=clvs,
                    )
        ax.set_title("Verif")
        if save:
            self.save()
        print("Plotted verification panel.")
        return

    @staticmethod
    def get_scales(scales=None,vrbl=None):
        if scales:
            S = scales
            cmap = S.cm
            clvs = S.clvs
        else:
            if vrbl:
                S = Scales(vrbl)
                cmap = S.cm
                clvs = S.clvs
            else:
                cmap = None
                clvs = None
        return cmap,clvs

    def plot_fcst(self,data,grid,titles=None,scales=None,
                    # W=None,cen_lon=None,cen_lat=None,
                    vrbl=None,ld=None,save=True):
        """
        Args:
            data: 3D numpy array, (members,lat,lon).
        """ 
        self.nmems = data.shape[0]

        cmap, clvs = self.get_scales(scales=scales,vrbl=vrbl)

        if ld == None:
            ld = {}

        # if W == None:
            # raise Exception("WRFOut object needed until refactoring complete.")

        # if cen_lat == None:
            # raise Exception("Need cen_lat and cen_lon until refactoring complete.")

        if titles is not None:
            title_itr = iter(titles)

        assert data.ndim == 3
        
        memlist = N.arange(self.naxes) - (self.naxes-self.nmems)

        for naxis, nmem, ax in zip(range(self.naxes),memlist,self.ax.flat):
            if naxis == 0:
                # This is for verification
                continue
            elif naxis == 1:
                # This is just for space
                # TODO: logic that skips this when there are 19 or whatever
                ax.grid("off")
                ax.axis("off")
                # continue
            else:
                assert nmem > -1
                BE = BirdsEye(fpath=None,ax=ax,fig=self.fig,grid=grid,proj=self.proj)
                if titles is not None:
                    title = ax.set_title(next(title_itr))
                print("Plotting forecast member #{} on subplot #{}".format(nmem,naxis))
                BE.plot2D(data[nmem,:,:],save=False,drawcounties=True,cb=True,
                            # clvs=clvs,cmap=cmap,cb=False,
                            # cen_lat=cen_lat,cen_lon=cen_lon,W=W,
                            # **ld,
                            levels=clvs,cmap=cmap)
        if save:
            self.save()

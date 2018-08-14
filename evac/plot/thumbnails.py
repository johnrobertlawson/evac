""" Postage stamps or thumbnails.
"""
import os
import pdb

import numpy as N
from mpl_toolkits.basemap import Basemap

from evac.plot.birdseye import BirdsEye
from evac.plot.figure import Figure
from evac.datafiles.radar import Radar
from evac.plot.scales import Scales

class ThumbNails(Figure):
    """ Plot multiple subplots.
    """
    def __init__(self,rowscols,fpath=None,outdir=None,fname='test_thumbs.png',
                    verif=True,proj='merc',use_basemap=False,
                    figsize=(10,8)):
        self.proj = proj
        self.use_basemap = use_basemap
        if not fpath:
            fpath = os.path.join(outdir,fname)
        super().__init__(fpath,proj=proj,use_basemap=use_basemap,
                        **self.setup_figure(rowscols,figsize))
        self.naxes = self.ax.size

    def setup_figure(self,rowscols,figsize):

        init_dict = dict(
            nrows = rowscols[0],
            ncols = rowscols[1],
            figsize = figsize,
            )
        return init_dict

    def plot_verif(self,data,grid,lats=None,lons=None,subplotn=0,
                    save=False,utc=None,vrbl=None,scales=None,
                    cb=True,levels=None):
        """ Plot verification postage stamp.
        """
        ax = self.ax.flatten()[subplotn]
        cmap,clvs = self.get_scales(vrbl=vrbl,lvs=levels)

        # Set up panel
        BE = BirdsEye(fpath=None,ax=ax,fig=self.fig,grid=grid,proj=self.proj,
                    use_basemap=self.use_basemap)

        # mplkwargs = dict('cmap':cmap, 'levels':clvs)
        BE.plot2D(data,save=False,cb=cb,
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
    def get_scales(scales=None,vrbl=None,lvs=None):
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
        if lvs is not None:
            clvs = lvs
        return cmap,clvs

    def plot_fcst(self,data,grid,titles=None,scales=None,
                    # W=None,cen_lon=None,cen_lat=None,
                    vrbl=None,ld=None,save=True,cb=False,
                    levels=None):
        """
        Args:
            data: 3D numpy array, (members,lat,lon).
        """ 
        self.nmems = data.shape[0]

        cmap, clvs = self.get_scales(scales=scales,vrbl=vrbl,lvs=levels)

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
                BE = BirdsEye(fpath=None,ax=ax,fig=self.fig,grid=grid,proj=self.proj,
                            use_basemap=self.use_basemap)
                if titles is not None:
                    title = ax.set_title(next(title_itr))
                print("Plotting forecast member #{} on subplot #{}".format(nmem,naxis))
                BE.plot2D(data[nmem,:,:],save=False,cb=cb,
                            # clvs=clvs,cmap=cmap,cb=False,
                            # cen_lat=cen_lat,cen_lon=cen_lon,W=W,
                            # **ld,
                            levels=clvs,cmap=cmap)
        if save:
            self.save()

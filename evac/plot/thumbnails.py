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
    def __init__(self,outdir,fname='test_thumbs.png',verif=True):
        self.outdir = outdir
        self.fname = fname

        super().__init__(self,**self.setup_figure())
        self.naxes = self.ax.size

    def setup_figure(self,):
        self.nmems = 18

        init_dict = dict(
            nrows = 4,
            ncols = 5,
            figsize = (10,8),
            )
        return init_dict

    def plot_verif(self,data,lats=None,lons=None,subplotn=0,
                    save=False,utc=None):
        """ Plot verification postage stamp.
        """
        ax = self.ax.flatten()[subplotn]
        # if isinstance(data,Radar):
            # data.plot_radar(fig=self.fig,ax=ax,
                    # drawcounties=True,cb=False)
        # else:
            # raise Exception("Not implemented")
        data.plot(fig=self.fig,ax=ax,
                drawcounties=True,cb=False,utc=utc)
        ax.set_title("Verif")
        if save:
            self.fig.save()
        return

    def plot_fcst(self,data,titles=None,scales=None,
                    W=None,cen_lon=None,cen_lat=None,
                    vrbl=None,
                    ld=None,save=True):
        """
        Args:
            data: 3D numpy array, (members,lat,lon).
        """ 
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

        if ld == None:
            ld = {}

        if W == None:
            raise Exception("WRFOut object needed until refactoring complete.")

        if cen_lat == None:
            raise Exception("Need cen_lat and cen_lon until refactoring complete.")

        if titles is not None:
            title_itr = iter(titles)

        assert data.ndim == 3
        assert data.shape[0] == self.nmems
        
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
                BE = BirdsEye(ax=ax,fig=self.fig,proj='merc')
                if titles is not None:
                    title = ax.set_title(next(title_itr))
                print("Plotting forecast member #{} on subplot #{}".format(nmem,naxis))
                BE.plot2D(data[nmem,:,:],save=False,drawcounties=True,
                            clvs=clvs,cmap=cmap,cb=False,
                            cen_lat=cen_lat,cen_lon=cen_lon,W=W,
                            **ld,)
        if save:
            self.save()

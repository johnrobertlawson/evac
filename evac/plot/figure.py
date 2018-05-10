import pdb
import os

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import evac.utils as utils
from evac.plot.scales import Scales

class Figure:
    """ Parent class for creating a matplotlib figure.

    This will use basemap or cartopy for images.

    Args:
        ax: optional matplotlib.axes instance
        fig: optional matplotlib.figure instance
        fpath: absolute path to figure file (.png)
        nrows, ncols (int): number of rows and columns, 
            if subplotting
        figsize (tuple): width x height of figure canvas
        initkwargs (dict): keyword arguments to pass to
            initiation of figure canvas, e.g., DPI.
    """
    def __init__(self,fpath,ax=None,fig=None,nrow=1,ncol=1,
                    figsize=(8,6),initkwargs=None):
        self.fpath = fpath

        if initkwargs is None:
            initkwargs = dict()
        
        if not ax and not fig:
            self.fig, self.ax = plt.subplots(nrow=nrow,ncol=ncol,
                                    figsize=figsize,**initkwargs)
        elif ax:
            assert fig
            self.ax = ax
            self.fig = fig

    @classmethod
    def create_fname(cls,*naming,joiner='_',append_ext=True,
                        use_time=True,filetype='.png'):
        """Creates file name from list of arguments.

        Optional:
            append_ext (bool)   :   If True, add .png extension
        """
        if len(naming) == 0:
            if use_time:
                utc = utils.generate_timestamp_fname(filetype)
            else:
                raise Exception("What are we creating the filename from?")
        fname = joiner.join([str(a) for a in naming])
        if append_ext:
            fname = cls.enforce_png_ext(fname)
        return fname

    @staticmethod
    def enforce_png_ext(fname):
        """Make sure file name or path ends with png.
        """
        if not fname.endswith('.png'):
            fname = fname + '.png'
        return fname

    @staticmethod
    def get_title_time(t):
        """Create pretty formatted date/time.
        """
        return utils.padded_times(t)
        
    def save(self,tight=True,close=True):
        utils.trycreate(self.fpath)
        if tight:
            self.fig.tight_layout()
        self.fig.savefig(self.fpath)
        print("Saving figure to {}".format(self.fpath))
        if close:
            plt.close(self.fig)
        return

    def just_one_colorbar(self,cb_fpath,cf,label=False,ticks=False):
        """ Create one colorbar separately.

        Useful for when manually creating subplots for papers etc.
        
        Todo:
            * Move to separate class, ColorBar(Figure).

        Args:
            cb_fpath: absolute path to look for, or create,
                the colorbar figure.
            cf: plot object from matplotlib plot function.            
            label: if True, apply label to colorbar axis
            ticks: if True, set custom ticks.
        """
        if not os,path.exists(cb_fpath):
            self.create_colorbar(cb_fpath,cf,labels=labels,ticks=ticks)
        return

    def create_colorbar(self,cb_fpath,cf,label=None,ticks=False):
        """ Create colorbar.

        Todo:
            * Move to ColorBar(Figure)

        Inputs:
        fpath   :   path to file
        fname   :   filename
        cf      :   contour filling for legend
        label   :   colorbar label

        """
        fig = plt.figure()
        CBax = fig.add_axes([0.15,0.05,0.7,0.02])
        CB = plt.colorbar(cf,cax=CBax,orientation='horizontal',ticks=tix)
        CB.set_label(label)
        fig.savefig(cb_fpath,tight=False)
        return

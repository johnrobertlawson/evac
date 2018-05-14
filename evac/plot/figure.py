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

    Allows geographical plotting with cartopy.

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
                    figsize=(8,6),initkwargs=None,proj=None,
                    grid=None):
        self.fpath = fpath
        self.grid = grid

        if initkwargs is None:
            initkwargs = dict()
        
        if proj:
            initkwargs['subplot_kw'] = dict(projection=proj)

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

    @staticmethod
    def just_one_colorbar(cb_fpath,cf,label=False,ticks=False):
        """ Create one colorbar separately.

        Useful for when manually creating subplots for papers etc.

        Args:
            cb_fpath: absolute path to look for, or create,
                the colorbar figure.
            cf: plot object from matplotlib plot function.            
            label: if True, apply label to colorbar axis
            ticks: if True, set custom ticks.
        """
        from evac.plot.colorbar import ColorBar
        if not os.path.exists(cb_fpath):
            CB = ColorBar(cb_fpath)
            CB.plot(cf=cf,labels=labels,ticks=ticks)
        return

    
    def __enter__(self):
        """ Set up the `with` block compatibility.
        """
        self.hold_opt = True
        self.save_opt = False
        print("Hold is on. Start plotting.")
        return self

    def __exit__(self,*args):
        """ 
        Todo:
            * Don't save if error?
        """
        self.hold_opt = False
        print("Figure is saved.")
        self.save()
        self.plt.close(fig)
        pass

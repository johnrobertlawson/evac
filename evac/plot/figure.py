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
    def __init__(self,fpath,ax=None,fig=None,nrows=1,ncols=1,
                    figsize=(8,6),proj=None,initkwargs=None,
                    grid=None,mplkwargs=None,**kwargs):
        self.fpath = fpath
        self.grid = grid

        # Option dictionaries
        self.update_kwargs(kwargs,cls=True,init=True)
        
        if initkwargs is None:
            initkwargs = dict()
        if proj:
            initkwargs['subplot_kw'] = dict(projection=proj)

        self.update_kwargs(mplkwargs,mpl=True,init=True)

        # Setting up figure canvas
        if (not ax) and (not fig):
            self.fig, self.ax = plt.subplots(nrows=nrows,ncols=ncols,
                                    figsize=figsize,**initkwargs)
        elif ax:
            assert fig
            self.ax = ax
            self.fig = fig


    def update_kwargs(self,kwdict,cls=False,mpl=False,init=False):
        if cls:
            assert not mpl
            if init:
                self.clskwargs = dict()
            skw = self.clskwargs
        elif mpl:
            if init:
                self.mplkwargs = dict()
            skw = self.mplkwargs
        else:
            raise Exception("Specify which kwargs to update")

        if isinstance(kwdict,dict):
            skw.update(kwdict)
        return

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
        if self.clskwargs.get('legend',False):
            self.ax.legend()
        self.hold_opt = False
        print("Figure is saved.")
        self.save()
        plt.close(self.fig)
        pass

    def _get_plot_options(self,mplkwargs=None,plotkwargs=None,**kwargs):
        """ Plot options common to all plotting methods.

        kwargs become clskwargs, which are options used by the class

        Args:
            mplkwargs (dict): options to be set after plot method (e.g., setting
                titles, xticks, etc)
            plotkwargs (dict): options passed to matplotlib when plotting (e.g.,
                x,y,z data, line colours, etc)
        """
        if mplkwargs == None:
            mplkwargs = {}
        if plotkwargs == None:
            plotkwargs = {}
        clskwargs = {}

        clskwargs['hold'] = kwargs.get('hold',False)
        try:
            clskwargs['save'] = kwargs['save']
        except:
            if clskwargs['hold']:
                clskwargs['save'] = False
            else:
                clskwargs['save'] = True
            print("Setting default save setting to {}".format(clskwargs['save']))
        # else:
        return clskwargs, mplkwargs, plotkwargs


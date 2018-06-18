import pdb
import os

import numpy as N
import matplotlib as M
# TODO backend should be set in matplotlibrc
# M.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from cartopy.feature import NaturalEarthFeature
from cartopy.io import shapereader
import cartopy

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
                    grid=None,mplkwargs=None,use_basemap=False,
                    **kwargs):
        self.fpath = fpath
        self.grid = grid
        self.proj = proj
        self.use_basemap = use_basemap
        self.basemap_done = False
        # Option dictionaries
        self.update_kwargs(kwargs,cls=True,init=True)
        
        if initkwargs is None:
            initkwargs = dict()
        if proj and not use_basemap:
            ccrsproj = utils.get_ccrs_proj(proj)
            initkwargs['subplot_kw'] = dict(projection=ccrsproj)
            print("Setting figure to {} projection.".format(proj))

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
        if self.fpath is not None:
            self.save()
            print("Figure is saved.")
        else:
            print("Figure not saved, as fpath is not set.")
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

    def basemap_setup(self,proj,resolution='i',draw_geography=True):
        G = self.grid

        if proj == 'merc':
            self.m = Basemap(projection=proj,llcrnrlat=G.llcrnrlat,
                    llcrnrlon=G.llcrnrlon,urcrnrlat=G.urcrnrlat,
                    urcrnrlon=G.urcrnrlon,lat_ts=G.lat_ts,
                    resolution=resolution,ax=self.ax)
        elif proj == 'lcc':
            raise Exception
        if draw_geography:
            self.draw_counties()
        self.basemap_done = True
        return

    def draw_counties(self,method=3,use_basemap=True):
        print("About to draw geographical info.")

        if use_basemap:
            self.m.drawcountries()
            self.m.drawcoastlines()
            self.m.drawstates()
            self.m.drawcounties()
        else:
            scale = '50m'
            # resolution = '10m'
            category = 'cultural'
            # name1 = 'admin_1_countries'
            name2 = 'admin_1_states_provinces_shp'
            # names = (name1, name2)
            names = (name2,)

            if method == 1:

                shp_fname = shapereader.natural_earth(resolution,category,name)
                df = geopandas.read_file(shp_fname)
                border = df.loc[df['ADMIN'] == 'USA']['geometry'].values[0]
                self.ax.add_geometries(border)
            elif method == 2:
                for name in names:
                    states = NaturalEarthFeature(category=category,
                                scale=scale,facecolor=None,name=name)
                    self.ax.add_feature(states,edgecolor='gray')
            elif method == 3:
                self.ax.add_feature(cartopy.feature.LAND)
                self.ax.add_feature(cartopy.feature.OCEAN)
                self.ax.add_feature(cartopy.feature.COASTLINE)
                self.ax.add_feature(cartopy.feature.BORDERS)
                self.ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
                self.ax.add_feature(cartopy.feature.RIVERS)

                states_provinces = cartopy.feature.NaturalEarthFeature(
                                        category='cultural',
                                        name='admin_1_states_provinces_lakes_shp',
                                        scale='50m',)
                self.ax.add_feature(states_provinces, edgecolor='gray')
                self.ax.background_patch.set_visible(False)
                ax.outline_patch.set_visible(False)
                # self.ax.coastlines()

        return

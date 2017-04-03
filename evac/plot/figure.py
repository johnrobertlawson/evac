"""All figures are subclasses of Figure.

"""
import pdb
import os

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from evac.utils.defaults import Defaults
import evac.utils.reproject as reproject
import evac.utils.unix_tools as unix_tools
import evac.utils.gis_tools as gis_tools

class Figure:
    def __init__(self,ax=None,fig=None,layout='normal',
                    mplargs=[],mplkwargs={},use_defaults=False):
        """ Parent class for creating a figure in matplotlib.

        Arguments:
            (No mandatory)

        Optional:
            ax                  :   matplotlib.axes instance to plot to
            fig                 :   matplotlib.figure instance to plot to
            layout (str)        :   if "insetv" or "inseth", uses GridSpec to create subplots
                                    of different sizes.
            mplargs (tuple,list):   other arguments to pass to matplotlib figure creation
            mplkwargs (dict)    :   other keyword arguments to pass to matplotlib
            use_defaults (bool) :   whether to load Defaults class. Set to be False
                                    if e.g. using your own matplotlibrc.
        """

        self.D = Defaults()
        self.save_figure = True

        # Create main figure
        if ax is not None and fig is not None:
            self.ax = ax
            self.fig = fig
            self.save_figure = False
        elif layout == 'insetv':
            self.fig = plt.figure(figsize=(8,6))
            self.gs = M.gridspec.GridSpec(1,2,width_ratios=[1,3])
            self.ax0 = plt.subplot(self.gs[0])
            self.ax1 = plt.subplot(self.gs[1])
        elif layout == 'inseth':
            self.fig = plt.figure(figsize=(6,8))
            self.gs = M.gridspec.GridSpec(2,1,height_ratios=[1,3])
            self.ax0 = plt.subplot(self.gs[0])
            self.ax1 = plt.subplot(self.gs[1])
        else:
            self.fig, self.ax = plt.subplots(*mplargs,**mplkwargs)

        if self.use_defaults:
            self.fig.set_dpi(self.D.dpi)

    def create_fname(self,*naming,joiner='_',append_ext=True):
        """Creates file name from list of arguments.

        Optional:
            append_ext (bool)   :   If True, add .png extension
        """
        fname = joiner.join([str(a) for a in naming])
        if append_ext:
            fname = self.enforce_png_ext(fname)
        return fname

    def enforce_png_ext(self,fname):
        """Make sure file name or path ends with png.
        """
        if not fname.endswith('.png'):
            fname = fname + '.png'
        return fname

    def get_title_time(self,t):
        """Create pretty formatted date/time.
        """
        return utils.padded_times(t)

    def set_figsize(self,width,height,fig):
        fig.set_size_inches(width,height)
        return

    def get_meshgrid_coords(self):
        self.mx,self.my = N.meshgrid(self.xx,self.yy)
        return

    def save(self,outpath,fname,tight=True):
        fname = self.enforce_png_ext(fname)
        unix_tools.trycreate(outpath)
        fpath = os.path.join(outpath,fname)
        if tight:
            self.fig.tight_layout()
        self.fig.savefig(fpath,bbox_inches='tight')
        print(("Saving figure {0}".format(fpath)))
        plt.close(self.fig)

#    def add_colorbar(self,datacb):
#        self.fig.colorbar(datacb)
#        return

    def just_one_colorbar(self,fpath,fname,cf,label=False,tix=False):
        try:
            with open(fpath): pass
        except IOError:
            self.create_colorbar(fpath,fname,cf,label=label,tix=tix)

    def create_colorbar(self,fpath,fname,cf,label='',tix=False):
        """
        Create colorbar.

        Inputs:
        fpath   :   path to file
        fname   :   filename
        cf      :   contour filling for legend
        label   :   colorbar label

        """
        self.fig = plt.figure()
        CBax = self.fig.add_axes([0.15,0.05,0.7,0.02])
        CB = plt.colorbar(cf,cax=CBax,orientation='horizontal',ticks=tix)
        CB.set_label(label)
        self.save(fpath,fname,tight=False)

    def basemap_from_newgrid(self,Nlim=None,Elim=None,Slim=None,Wlim=None,proj='merc',
                    lat_ts=None,resolution='i',nx=None,ny=None,
                    tlat1=30.0,tlat2=60.0,cen_lat=None,cen_lon=None,
                    lllon=None,lllat=None,urlat=None,urlon=None,
                    drawcounties=False):
        """Wrapper for utility method.
        """
        # m,lons,lats,xx,yy = reproject.create_new_grid(*args,**kwargs)
        # return m, lons, lats, xx[0,:], yy[:,0]
        self.m, self.lons, self.lats, self.xx, self.yy = reproject.create_new_grid(
                    Nlim=Nlim,Elim=Elim,Slim=Slim,Wlim=Wlim,proj=proj,
                    lat_ts=lat_ts,resolution=resolution,nx=nx,ny=ny,
                    tlat1=tlat1,tlat2=tlat2,cen_lat=cen_lat,cen_lon=cen_lon,
                    lllon=lllon,lllat=lllat,urlat=urlat,urlon=urlon)
        self.m.drawcoastlines()
        self.m.drawstates()
        self.m.drawcountries()
        if isinstance(drawcounties,str):
            self.m.readshapefile(drawcounties,'counties')
        return

    def basemap_setup(self,smooth=1,lats=False,lons=False,proj='merc',
                        Nlim=False,Elim=False,Slim=False,Wlim=False,
                        drawcounties=False):
        """
        TODO:
        Merge with above method.
        Needs rewriting to include limited domains based on lats/lons.
        Currently, assuming whole domain is plotted.
        """

        # Fetch settings
        if self.use_defaults:
            basemap_res = self.D.basemap_res

        if proj=='lcc':
            width_m = self.W.dx*(self.W.x_dim-1)
            height_m = self.W.dy*(self.W.y_dim-1)

            m = Basemap(
                projection=proj,width=width_m,height=height_m,
                lon_0=self.W.cen_lon,lat_0=self.W.cen_lat,lat_1=self.W.truelat1,
                lat_2=self.W.truelat2,resolution=basemap_res,area_thresh=500,
                ax=self.ax)

        elif proj=='merc':
            if self.W and not Nlim and not isinstance(lats,N.ndarray):
                Nlim,Elim,Slim,Wlim = self.W.get_limits()
            elif Nlim==False:
                Nlim = lats.max()
                Slim = lats.min()
                Elim = lons.max()
                Wlim = lons.min()

            m = Basemap(projection=proj,
                        llcrnrlat=Slim,
                        llcrnrlon=Wlim,
                        urcrnrlat=Nlim,
                        urcrnrlon=Elim,
                        lat_ts=(Nlim-Slim)/2.0,
                        resolution='l',
                        ax=self.ax)

        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        if isinstance(drawcounties,str):
            # if not drawcounties.endswith('/'):
                # drawcounties = drawcounties + '/'
            m.readshapefile(drawcounties,'counties',linewidth=0.3,color='grey')
        # import pdb; pdb.set_trace()
        # Draw meridians etc with wrff.lat/lon spacing
        # Default should be a tenth of width of plot, rounded to sig fig

        # s = slice(None,None,smooth)
        if self.W and not isinstance(lats,N.ndarray):
            x,y = m(self.W.lons,self.W.lats)
        else:
            x,y = m(*N.meshgrid(lons,lats))
        # pdb.set_trace()
        return m, x, y
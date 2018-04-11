""" Figure superclass.

Todo:
	* Merge __get_plot_options.
"""

import pdb
import os

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from evac.utils.defaults import Defaults
import evac.utils as utils
# import evac.utils.reproject as reproject
# import evac.utils.unix_tools as unix_tools
# import evac.utils.gis_tools as gis_tools

class Figure:
    """ Parent class for creating a figure in matplotlib.

    Todo:
        * Rewrite (along with other plotting scripts) for cartopy, due to 
            the depreciation of basemap.

    Args:
        ax                  :   matplotlib.axes instance to plot to
        fig                 :   matplotlib.figure instance to plot to
        layout (str)        :   if "insetv" or "inseth", uses GridSpec to create subplots
                                of different sizes.
        mplargs (tuple,list):   other arguments to pass to matplotlib figure creation
        mplkwargs (dict)    :   other keyword arguments to pass to matplotlib
        use_defaults (bool) :   whether to load Defaults class. Set to be False
                                if e.g. using your own matplotlibrc.
    """

    def __init__(self,ax=None,fig=None,layout='normal',
                    mplargs=[],mplkwargs={},use_defaults=False):
        # self.D = Defaults()
        self.save_figure = True
        self.use_defaults = use_defaults

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

    def create_fname(self,*naming,joiner='_',append_ext=True,
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

    def save(self,outpath=None,fname=None,tight=True):
        if (outpath is None) and (fname is None):
            try:
                outpath = self.outpath
                fname = self.fname
            except:
                raise Exception("Provide both output directory and filename.")
        fname = self.enforce_png_ext(fname)
        utils.trycreate(outpath)
        fpath = os.path.join(outpath,fname)
        if tight:
            self.fig.tight_layout()
            self.fig.savefig(fpath,bbox_inches='tight')
        else:
            self.fig.savefig(fpath)#,bbox_inches='tight')
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

    def basemap_from_newgrid(self,Nlim=None,Elim=None,Slim=None,Wlim=None,proj=None,
                    lat_ts=None,resolution='i',nx=None,ny=None,
                    tlat1=30.0,tlat2=60.0,cen_lat=None,cen_lon=None,
                    lllon=None,lllat=None,urlat=None,urlon=None,
                    drawcounties=False,xx=None,yy=None,
                    lats=None,lons=None):
        """Wrapper for utility method.
        """
        # m,lons,lats,xx,yy = utils.create_new_grid(*args,**kwargs)
        # return m, lons, lats, xx[0,:], yy[:,0]
        if isinstance(proj,str):
            self.proj = proj
        elif not hasattr(self,'proj'):
            self.proj = 'merc'

        if nx is None:
            if xx is None:
                if lons is None:
                    raise Exception("Need to give either nx/ny or xx/yy")
                else:
                    if lons.ndim == 2:
                        ny,nx = lons.shape
                    else:
                        nx = len(lons)
                        ny = len(lats)
            
            else:
                ny,nx = xx.shape

        # for merc
        # if None in (Nlim,Elim,Slim,Wlim,lat_ts,nx,ny):

        # def create_new_grid(Nlim=None,Elim=None,Slim=None,Wlim=None,proj='merc',
                    # lat_ts=None,resolution='i',nx=None,ny=None,
                    # tlat1=30.0,tlat2=60.0,cen_lat=None,cen_lon=None,):
        self.m, self.lons, self.lats, self.xx, self.yy = utils.create_new_grid(
                    Nlim=Nlim,Elim=Elim,Slim=Slim,Wlim=Wlim,proj=self.proj,
                    lat_ts=lat_ts,resolution=resolution,nx=nx,ny=ny,
                    tlat1=tlat1,tlat2=tlat2,cen_lat=cen_lat,cen_lon=cen_lon,)
                    # lllon=lllon,lllat=lllat,urlat=urlat,urlon=urlon)
        self.m.drawcoastlines()
        self.m.drawstates()
        self.m.drawcountries()
        if isinstance(drawcounties,str):
            self.m.readshapefile(drawcounties,'counties')
        return

    def basemap_setup(self,smooth=1,lats=False,lons=False,proj=None,
                        Nlim=False,Elim=False,Slim=False,Wlim=False,
                        drawcounties=False,res='i'):
        """
        TODO:
        Merge with above method.
        Needs rewriting to include limited domains based on lats/lons.
        Currently, assuming whole domain is plotted.
        """
        # if hasattr(self,'proj'):
            # assert proj == self.proj

        if self.proj=='lcc':
            width_m = self.W.dx*(self.W.x_dim-1)
            height_m = self.W.dy*(self.W.y_dim-1)

            m = Basemap(
                projection=self.proj,width=width_m,height=height_m,
                lon_0=self.W.cen_lon,lat_0=self.W.cen_lat,lat_1=self.W.truelat1,
                lat_2=self.W.truelat2,resolution=res,area_thresh=500,
                ax=self.ax)

        elif self.proj=='merc':
            if hasattr(self,'W') and not Nlim and not isinstance(lats,N.ndarray):
                Nlim,Elim,Slim,Wlim = self.W.get_limits()
            elif Nlim==False:
                Nlim = lats.max()
                Slim = lats.min()
                Elim = lons.max()
                Wlim = lons.min()
            else:
                # raise Exception
                pass

            m = Basemap(projection=self.proj,
                        llcrnrlat=Slim,
                        llcrnrlon=Wlim,
                        urcrnrlat=Nlim,
                        urcrnrlon=Elim,
                        lat_ts=(Nlim-Slim)/2.0,
                        resolution='l',
                        ax=self.ax)
        else:
            raise Exception

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
        # if hasattr(self,'W') and not isinstance(lats,N.ndarray):
        if hasattr(self,'W'):
            if self.W is not None:
                x,y = m(self.W.lons,self.W.lats)
            else:
                x,y = m(*N.meshgrid(lons,lats))
                # x,y = m(N.meshgrid(lons,lats))
        # pdb.set_trace()
        return m, x, y


    def __get_plot_options2(self,*args,**kwargs):
        """ Plot options common to all plotting methods.

        Note:
            * hold and save are mutually exclusive and must be 
                opposite booleans.

        Args:
            hold (bool): If True, the figure is not saved after the 
                call, and further plots can be made. If False,
                it is saved immediately and `save()` is not needed.
        """
        clskwargs = {}
        mplkwargs = {}
        plotkwargs = {}

        clskwargs['hold'] = kwargs.get('hold',False)
        clskwargs['save'] = kwargs.get('save',True)
        assert clskwargs['hold'] != clskwargs['save']

        return clskwargs, mplkwargs, plotkwargs

      def __get_plot_options1(self,vrbl,*args,**kwargs):
          """ Filter arguments and key-word arguments for plotting methods.

          Whatever is in dictionary will overwrite defaults in the plotting
          method.

          These may be
              * fhrs (forecast hour plotting times - or all)
              * ensmems (ensemble members, or all)


          """
          # Get plotting levels if not already given
          # TODO: rewrite this using hasattr() or something.
          S = Scales(vrbl)
          if not 'levels' in kwargs:
              kwargs['levels'] = S.clvs
          if not 'cmap' in kwargs:
              kwargs['cmap'] = S.cm

          # Specific things for certain variables
          if vrbl in ('REFL_10CM',"REFL_comp"):
              pass

          # Save all figures to a subdirectory
          if subdir in kwargs:
              utils.trycreate(subdir,is_folder=True)

          # What times are being plotted?
          # If all, return list of all times
          if 'utc' in kwargs:
              pass
          elif ('fchr' not in kwargs) or (kwargs['fchr'] == 'all'):
              kwargs['utc'] = E.list_of_times
          # Does this pick up on numpy arange?
          elif isinstance(kwargs['fchr'], (list,tuple)):
              kwargs['utc'] = []
              for f in kwargs['fchr']:
                  utc = self.inittime + datetime.timedelta(seconds=3600*f)
                  kwargs['fchr'].append(utc)

          # Make domain smaller if requested

          # Save data before plotting
          clskwargs['save_data'] = kwargs.get('save_data',False)
          return clskwargs,plotkwargs,mplkwargs

    def 

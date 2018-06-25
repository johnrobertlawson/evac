""" Composite radar archive data from mesonet.agron.iastate.edu.
"""
import os
import pdb
import calendar
import glob
import datetime

import scipy
import numpy as N
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as N

from evac.datafiles.wrfout import WRFOut
from evac.utils.defaults import Defaults

from evac.datafiles.pngfile import PNGFile
import evac.utils as utils
import evac.plot.colourtables as ct
from evac.plot.birdseye import BirdsEye
from evac.datafiles.obs import Obs
from evac.utils.grid import Grid

class Radar(PNGFile,Obs):
    """ Radar data download, manipulation, and plotting via :any:`evac.plot.birdseye`.

    This uses composite reflectivity data stored on Iowa State Univ.
    servers. Radar refers to one .wld and one .png file. For a group, use ObsGroup.
    User should specify either fpath to the .png file (when the data has 
    already been downloaded), or fdir and utc (Radar will look up the files, and if
    not present, the data is downloaded).

    Example:
        Acquire and plot radar data on a limited domain::

            from evac.datafiles.radar import Radar
            import datetime

            utc = datetime.datetime(2017,4,4,21,30,0)
            limdict = dict(Nlim=45.0,Elim=-90.0,Slim=35.0,Wlim=-100.0)
            radardir = '/path/to/radar/data'

            fig,ax = plt.subplots(1)
            R = Radar(radardir,utc=utc)
            R.get_subdomain(**limdict,overwrite=True)
            R.plot_radar(fig=fig,ax=ax,drawcounties=True)
            fig.savefig('/path/to/save.png')

    Todos:
        * Remove redundant basemap generation method and integrate with
            any new cartopy plotting scripts in :any:`evac.plot.birdseye`.
        * Downloading appears to be broken

    Args:
        fpath: Absolute path to (a) folder with radar data or
            (b) .png file containing data (in same directory
            as .wld file, with same naming scheme)
        utc (datetime.datetime): Time of radar data in UTC.
        proj: Projection for plotting.
    """
    def __init__(self,fpath,utc=None,proj='merc'):
        self.proj = proj

        # The _root attributes do not have an extension (.png, .wld)
        if fpath.endswith('.png'):
            utils.wowprint("Looking for radar **data file**.")
            self.fpath_png = fpath
            self.fpath_root = fpath[:-4]
            self.fname_root = os.path.basename(fpath)[:-4]
            self.fdir = os.path.dirname(self.fpath_png)
            self.utc = self.date_from_fname(self.fpath_png)
            self.fmt = self.get_fmt()
        else:
            utils.wowprint("About to **download** radar data file.")
            self.fdir = fpath
            assert utc is not None
            self.utc = utc
            self.fmt = self.get_fmt()
            self.fname_root = self.get_radar_fname()
            self.fname_png = self.fname_root + '.png'
            self.fpath_root = os.path.join(self.fdir,self.fname_root)
            self.fpath_png = self.fpath_root + '.png'

            self.download_data()

        # This is to maintain API with ObsGroup etc
        self.fpath = self.fpath_png

        self.fpath_wld = self.fpath_root + '.wld'

        data = scipy.ndimage.imread(self.fpath_png,mode='P')
        self.data = N.flipud(data)

        self.xlen, self.ylen, = self.data.shape

        # Metadata
        f = open(self.fpath_wld,'r').readlines()

        # pixel size in the x-direction in map units/pixel
        self.xpixel = float(f[0])

        # rotation about y-axis
        self.roty = float(f[1])

        # rotation about x-axis
        self.rotx = float(f[2])

        # pixel size in the y-direction in map units,
        self.ypixel = float(f[3])

        # x-coordinate of the center of the upper left pixel
        self.ulx = float(f[4])

        # y-coordinate of the center of the upper left pixel
        self.uly = float(f[5])

        # lower right corner
        self.lrx = self.ulx + self.ylen*self.xpixel
        self.lry = self.uly + self.xlen*self.ypixel

        if self.fmt == 'n0r':
            self.clvs = N.arange(0,80,5)
        elif self.fmt == 'n0q':
            self.clvs = N.arange(0,90.0,0.5)
        # import pdb; pdb.set_trace()

        # self.lats1D = N.linspace(self.lry,self.uly,self.xlen)[::-1]
        self.lats1D = N.linspace(self.lry,self.uly,self.xlen)
        self.lons1D = N.linspace(self.ulx,self.lrx,self.ylen)
        self.lons, self.lats = N.meshgrid(self.lons1D,self.lats1D)
        # self.lats, self.lons = N.meshgrid(self.lats1D,self.lons1D)
        assert self.data.shape == self.lats.shape

        self.grid = Grid(self)
        self.dBZ = self.get_dBZ()

    def download_data(self,ext='both'):
        if ext == 'both':
            exts = ('.png','.wld')
        else:
            exts = (ext,)
        for ex in exts:
            url = self.get_radar_url(utc=self.utc)
            fname_root = self.get_radar_fname(utc=self.utc)
            urlf = os.path.join(url,fname_root+ex)
            dest_fpath = os.path.join(self.fdir,self.fname_root+ex)
            cmd1 = 'wget {0} -O {1}'.format(urlf,dest_fpath)
            os.system(cmd1)
        return

    def get_fmt(self):
        tt = utils.ensure_timetuple(self.utc)
        self.tt = tt

        # Date for format change
        change = (2010,10,25,0,0,0)

        if tt[0] < 1995:
            print("Too early in database.")
            raise Exception
        elif calendar.timegm(tt) < calendar.timegm(change):
            return 'n0r'
        else:
            return 'n0q'

    def get_radar_fname(self,utc=None):
        if utc is None:
            utc = self.utc

        # Assume get_fmt has been run, and self.fmt exists
        fname = '{0}_{1:04d}{2:02d}{3:02d}{4:02d}{5:02d}'.format(
                    self.fmt,*self.tt)
        return fname

    @staticmethod
    def get_radar_url(utc):
        dt = utils.ensure_timetuple(utc)
        return ('mesonet.agron.iastate.edu/archive/data/'
                '{0:04d}/{1:02d}/{2:02d}/GIS/uscomp/'.format(*dt))

    def generate_basemap(self,fig=None,ax=None,Nlim=False,Elim=False,Slim=False,
                            Wlim=False):
        """
        Generate basemap.

        """
        if fig or ax:
            self.fig = fig
            self.ax = ax
        if not isinstance(Nlim,float):
            Nlim = self.uly
            Elim = self.lrx
            Slim = self.lry
            Wlim = self.ulx

        self.m = Basemap(projection='merc',
                    llcrnrlat=Slim,
                    llcrnrlon=Wlim,
                    urcrnrlat=Nlim,
                    urcrnrlon=Elim,
                    lat_ts=(Nlim-Slim)/2.0,
                    resolution='l')
        # self.m.drawcoastlines()
        # self.m.drawstates()
        # self.m.drawcountries()
        return self.m

    def __get_subdomain(self,Nlim,Elim,Slim,Wlim,overwrite=False):
        """
        This has now been moved to superclass Obs.
        Return data array between bounds

        If overwrite is True, replace class data with new subdomain
        """
        data,lats,lons = utils.return_subdomain(self.data,self.lats,
                                self.lons,Nlim,Elim,Slim,Wlim)

        if overwrite:
            self.lats = lats
            self.lons = lons
            self.data = data
            return
        else:
            return data,lats,lons

    def get_dBZ(self,data='self'):
        if data == 'self':
            data = self.data

        # pdb.set_trace()
        if self.fmt == 'n0q':
            dBZ = (data*0.5)-32
            # dBZ = (data*0.25)-32
        elif self.fmt == 'n0r':
            dBZ = (data*5.0)-30
        return dBZ

    def get(self,*args,**kwargs):
        """Wrapper for get_dBZ.
        """
        return self.get_dBZ(*args,**kwargs)

    def plot(self,*args,**kwargs):
        """ Wrapper for plot_radar.
        """
        ret = self.plot_radar(*args,**kwargs)
        return ret
    
    @staticmethod
    def date_from_fname(f,fullpath=False):
        if fullpath or f.startswith('/'):
            f = os.path.basename(f)
        fn, _ext = f.split('.')
        n0q,d = fn.split('_')
        fmt = '%Y%m%d%H%M'
        utc = datetime.datetime.strptime(d,fmt)
        return utc

    def plot_radar(self,outdir=False,fig=False,ax=False,fname=False,Nlim=False,
                    Elim=False, Slim=False,Wlim=False,cb=True,
                    drawcounties=False,save='auto'):
        """
        Plot radar data.

        Args:

        save        :   (str,bool) - If 'auto', saving will only occur if
                        fig and ax are not specified.
                        If True, the figure is saved; if False, not.
        """
        # if not fig:
            # fig, ax = plt.subplots()
        # self.generate_basemap(fig,ax,Nlim,Elim,Slim,Wlim)
        #lons, lats = self.m.makegrid(self.xlen,self.ylen)
        if isinstance(Nlim,float):
            data, lats, lons = self.get_subdomain(Nlim,Elim,Slim,Wlim)
        else:
            data = self.data
            lats = self.lats #flip lats upside down?
            lons = self.lons
            # x,y = self.m(*N.meshgrid(lons,lats))

        # xx,uy = self.m(lons,lats)
        # x,y = self.m(*N.meshgrid(lons,lats))
        # x,y = self.m(*N.meshgrid(lons,lats[::-1]))

        # Custom colorbar
        radarcmap = ct.reflect_ncdc(self.clvs)
        # radarcmap = ct.ncdc_modified_ISU(self.clvs)

        # Convert pixel levels to dBZ
        dBZ = self.get_dBZ(data)


        # dBZ[dBZ<0] = 0

        # def plot2D(self,data,fname,outdir,plottype='contourf',
                    # save=1,smooth=1,lats=False,lons=False,
                    # clvs=False,cmap=False,title=False,colorbar=True,
                    # locations=False):
        if not fname:
            tstr = utils.string_from_time('output',self.utc)
            fname = 'verif_radar_{0}.png'.format(tstr)
        F = BirdsEye(fig=fig,ax=ax,proj=self.proj)
        if cb:
            cb = 'horizontal'
        if save is 'auto':
            if (fig is not False) and (ax is not False):
                save = False
            else:
                save = True
        F.plot2D(dBZ,fname,outdir=outdir,lats=lats,lons=lons,
                    cmap=radarcmap,clvs=N.arange(5,90,5),
                    cb=cb,cblabel='Composite reflectivity (dBZ)',
                    drawcounties=drawcounties,save=save,lat_ts=50.0,)
        # im = self.ax.contourf(x,y,dBZ,alpha=0.5,cmap=radarcmap,
                                # levels=N.arange(5.0,90.5,0.5))
        # outpath = os.path.join(outdir,fname)
        # self.fig.colorbar(im,ax=self.ax)
        # self.fig.savefig(outpath)



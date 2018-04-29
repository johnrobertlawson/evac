import os
import pdb
import calendar
import glob

import scipy
import numpy as N

from evac.datafiles.pngfile import PNGFile
import evac.utils as utils
import evac.plot.colourtables as ct
from evac.plot.birdseye import BirdsEye
from evac.datafiles.obs import Obs

class Radar(PNGFile,Obs):
    """ Radar data download, manipulation, and plotting via :any:`evac.plot.birdseye`.

    This uses composite reflectivity data stored on Iowa State Univ.
    servers.

    Example:
        Acquire and plot radar data on a limited domain::

            from evac.datafiles.radar import Radar
            import datetime

            utc = datetime.datetime(2017,4,4,21,30,0)
            limdict = dict(Nlim=45.0,Elim=-90.0,Slim=35.0,Wlim=-100.0)
            radardir = '/path/to/radar/data'

            fig,ax = plt.subplots(1)
            R = Radar(utc,radardir)
            R.get_subdomain(**limdict,overwrite=True)
            R.plot_radar(fig=fig,ax=ax,drawcounties=True)
            fig.savefig('/path/to/save.png')

    Note that, if the radar data does not exist in radardir, it will be 
    downloaded automatically.

    Todos:
        * Allow different inheritance depending on the data format.
            Probably using the __next__ built-in.
        * Remove redundant basemap generation method and integrate with
            any new cartopy plotting scripts in :any:`evac.plot.birdseye`.

    Args:
        utc (datetime.datetime): Time of radar data in UTC.
        datapath: Absolute path to folder with radar data.
        proj: Projection for plotting.
    """
    def __init__(self,utc,datapath,proj='merc'):
        """
        Composite radar archive data from mesonet.agron.iastate.edu.

        Args:
            datapath (str)      :   Absolute path to folder or .png file
            wldpath (str)       :   Absolute path to .wld file
            fmt (str)           :   format of data - N0Q or N0R
        """
        self.proj = proj
        self.utc = utc
        fname_root = self.get_radar_fname()

        # Check for file
        # Download if not available
        fpath = os.path.join(datapath,fname_root)
        for ex in ('.png','.wld'):
            scan = glob.glob(fpath+ex)
            print("Looking in",datapath)
            print("Contents:",scan)

            if len(scan) == 0:
                url = self.get_radar_url()
                urlf = os.path.join(url,fname_root+ex)
                cmd1 = 'wget {0} -P {1}'.format(urlf,datapath)
                os.system(cmd1)

        png = fpath+'.png'
        wld = fpath+'.wld'

        self.data = scipy.ndimage.imread(png,mode='P')
        # if len(self.data.shape) == 3:
            # self.data = self.data[:,:,0]

        self.xlen, self.ylen, = self.data.shape

        # Metadata
        f = open(wld,'r').readlines()

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

        self.lats = N.linspace(self.lry,self.uly,self.xlen)[::-1]
        self.lons = N.linspace(self.ulx,self.lrx,self.ylen)

    def get_radar_fname(self):
        tt = utils.ensure_timetuple(self.utc)

        # Date for format change
        change = (2010,10,25,0,0,0)

        if tt[0] < 1995:
            print("Too early in database.")
            raise Exception
        elif calendar.timegm(tt) < calendar.timegm(change):
            self.fmt = 'n0r'
        else:
            self.fmt = 'n0q'

        fname = '{0}_{1:04d}{2:02d}{3:02d}{4:02d}{5:02d}'.format(
                    self.fmt,*tt)
        return fname

    def get_radar_url(self):
        dt = utils.ensure_timetuple(self.utc)
        return ('mesonet.agron.iastate.edu/archive/data/'
                '{0:04d}/{1:02d}/{2:02d}/GIS/uscomp/'.format(*dt))

    def generate_basemap(self,fig,ax,Nlim=False,Elim=False,Slim=False,
                            Wlim=False):
        """
        Generate basemap.

        """
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
                    resolution='l',
                    ax=self.ax)
        self.m.drawcoastlines()
        self.m.drawstates()
        self.m.drawcountries()

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

    def get_dBZ(self,data):
        if data == 'self':
            data = self.data

        # pdb.set_trace()
        if self.fmt == 'n0q':
            dBZ = (data*0.5)-32
            # dBZ = (data*0.25)-32
        elif self.fmt == 'n0r':
            dBZ = (data*5.0)-30
        return dBZ

    def plot(self,*args,**kwargs):
        """ Wrapper for plot_radar.
        """
        ret = self.plot_radar(*args,**kwargs)
        return ret

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

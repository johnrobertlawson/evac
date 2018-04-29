"""
Todo:
    * Rework the _get/_set_subdomain logic. This needs to be done at the
        superclass level and have consistent API across all classes/methods.
"""
import os
import glob
import datetime
import pdb

import numpy as N
from mpl_toolkits.basemap import Basemap

from evac.plot.birdseye import BirdsEye
from evac.datafiles.gribfile import GribFile
from evac.datafiles.obs import Obs
import evac.utils as utils

class StageIV(GribFile,Obs):
    """ Class that holds Stage IV verification data.

    Args:
        fpath: Absolute path to Stage IV file. 
        loadobj (bool): if True, load the data as instances of 
            :class:`~evac.datafiles.gribfile.GribFile()`.

    Todos:
        * Check for 9999s (missing) or negative. Other codes
        * Consistency: are we loading one or a group of files?

    """
    def __init__(self,fpath):

        try:
            import pygrib
        except ImportError:
            print("You need to install pygrib for the StageIV class.")
            raise Exception
        else:
            self.pygrib = pygrib

        self.fpath = fpath
        self.G = self.pygrib.open(fpath)
        self.fname = os.path.basename(self.fpath)
        self.utc = self.date_from_fname(self.fname)

        self.projection()
        self.lats, self.lons = self.return_latlon()
        self.return_subdomain = False

    def get(self,utc=None,lv=None,return_masked=False,check=True):
        """
        Get a given time, in similar manner to WRFOut.

        Wrapper for return_array with reshape to 4D

        lv is always None (compatibility with other gets).

        Args:
            check: if True, raise Exception if the data is not
                a legitimate array (e.g., shape = (0,0)

        """
        data2D = self.return_array(utc,accum_hr=accum_hr,check=check)
        # pdb.set_trace()

        if self.return_subdomain is not False:
            assert return_masked is False
            data2D = self.return_subdomain(data2D.data)

        data4D = data2D[N.newaxis,N.newaxis,:,:]

        if return_masked:
            return data4D
        return data4D.data

    @staticmethod
    def date_from_fname(f):
        _1, d, _2 = f.split('.')
        fmt = '%Y%m%d%H'
        utc = datetime.datetime.strptime(d,fmt)
        return utc

    def load_gribpart(self):
        # self.G = self.load_data(self.G,loadobj=True)
        self.G.seek(0)
        gg = G.select(name='Total Precipitation')[0]
        return gg

    def load_accum(self):
        raise Exception("Now use return_array.")
        return

    def return_latlon(self):
        gg = self.load_gribpart()
        latlon = gg.latlons()
        lats, lons = latlon
        return lats,lons

    def check_quality(self,arr):
        """ Probably needs better logic than this. TODO.
        """
        if arr == arr[0,0]:
            return False
        return True

    def return_array(self,check=True):
        """ 
        Return the (masked) array.

        Args:
            check: if True, raise exception if no data.
        """
        gg = self.load_gribpart()
        arr = gg.values

        if check:
            OK = self.check_quality(arr)
            if not OK:
                raise Exception("These data are rubbish.")
        return arr

    def return_point(self,lat,lon):
        latidx,lonidx = utils.get_latlon_idx(self.lats,self.lons,lat,lon)
        arr = self.return_array()
        return arr[latidx,lonidx]

    def projection(self):
        self.m = Basemap(projection='npstere',lon_0=-105.0,#lat_1=60.0,
                # llcrnrlon=lllon,llcrnrlat=lllat,urcrnrlon=urlon,urcrnrlat=urlat,
                            boundinglat=24.701632)
        self.xx, self.yy = self.m(self.lons,self.lats)
        # pdb.set_trace()
        # self.mx, self.my = N.meshgrid(self.xx,self.yy)

        # lllon = -119.023
        self.lllon = self.lons[0,0]
        # lllat = 23.117
        self.lllat = self.lats[0,0]
        # urlon = -59.9044
        self.urlon = self.lons[-1,-1]
        # urlat = 45.6147234
        self.urlat = self.lats[-1,-1]

        self.shape = self.lats.shape
        assert self.lats.shape == self.lons.shape

    def __get_subdomain(self,Nlim,Elim,Slim,Wlim,overwrite=False):
        """ Overriden method from parent.

        This is because StageIV uses a catalogue of data arrays,
            rather than a single array. Future TODO is to 
            make ST4 just one file, and do another class
            that is a 'data ensemble'.
        """
        # This is just here for compatibility with API
        overwrite = False

        self.set_subdomain(Nlim,Elim,Slim,Wlim,enable=True)
        return

    def __set_subdomain(self,Nlim,Elim,Slim,Wlim,enable=True):
        """ A partner of get_subdomain().
        
        Call signature changed from Obs method, to
        enable a time to be chosen, and data to be
        returned (not overridden). These settings are
        kept for future get() commands.
        """
        self.Nlim = Nlim
        self.Elim = Elim
        self.Slim = Slim
        self.Wlim = Wlim

        if enable:
            self.return_subdomain = self._return_subdomain
        return 

    def _return_subdomain(self,data):
        data,lats,lons = super().get_subdomain(self.Nlim,self.Elim,self.Slim,
                                                self.Wlim,data=data)
        # Store the new lat/lon array just in case.
        self._subdom_lats = lats
        self._subdom_lons = lons
        return data

        
    def __plot(self,data=None,outdir=False,fig=False,ax=False,fname=False,
                    # Nlim=False, Elim=False, Slim=False,Wlim=False,
                    cb=True,drawcounties=False,save='auto',
                    lats=None,lons=None,proj='merc',
                    utc=None):
        """ Plot data.

        Todo:
            * This doesn't feel right. Plotting should be done elsewhere.

        Args:

        save        :   (str,bool) - If 'auto', saving will only occur if
                        fig and ax are not specified.
                        If True, the figure is saved; if False, not.
        """
        # if isinstance(Nlim,float):
            # data, lats, lons = self.get_subdomain(Nlim,Elim,Slim,Wlim)
        # else:
            # if not given, try to use self's

        if data == None:
            assert utc is not None
            data = self.get(utc=utc)
        # data = getattr(self,'data',data)
        if lats == None:
            lats = self.lats 
        # lats = getattr(self,'lats',lats)
        if lons == None:
            lons = self.lons
        # lons = getattr(self,'lons',lons)

        # Need to implement options for this
        cmap = None
        clvs = None

        if not fname:
            tstr = utils.string_from_time('output',utc)
            fname = 'verif_ST4_{0}.png'.format(tstr)
        F = BirdsEye(fig=fig,ax=ax,proj=proj)
        if cb:
            cb = 'horizontal'
        if save == 'auto':
            if (fig is not False) and (ax is not False):
                save = False
            else:
                save = True
        F.plot2D(data,fname,outdir=outdir,lats=lats,lons=lons,
                    cmap=cmap,clvs=clvs,
                    cb=cb,cblabel='Accumulated precipitation (mm)',
                    drawcounties=drawcounties,save=save,lat_ts=50.0,)

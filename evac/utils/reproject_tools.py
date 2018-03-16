import os
import pdb

from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
import numpy as N

from WEM.postWRF.postWRF.wrfout import WRFOut
from WEM.postWRF.postWRF.hrrr import HRRR

class VerifGrid:
    """
    The verification grid should be inside the 1 km domain.

    Roughly 5 km res?

    W is the 1 km WRFOut instance
    """
    def __init__(self,W,nx,ny,xy1D=True):
        self.proj = 'lcc'
        # These are urcrnrlon etc
        self.Nlim = W.kwargs['urcrnrlat']
        self.Elim = W.kwargs['urcrnrlon']
        self.Slim = W.kwargs['llcrnrlat']
        self.Wlim = W.kwargs['llcrnrlon']

        self.nx = nx
        self.ny = ny
        self.cen_lat = W.kwargs['lat_0']
        self.cen_lon = W.kwargs['lon_0']
        self.tlat1 = W.kwargs['lat_1']
        self.tlat2 = W.kwargs['lat_2']

        self.m, self.lons, self.lats, self.xx, self.yy = create_new_grid(Nlim=self.Nlim,
                        Elim=self.Elim,Slim=self.Slim,Wlim=self.Wlim,nx=self.nx,ny=self.ny,proj=self.proj,
                        tlat1=self.tlat1,tlat2=self.tlat2,cen_lat=self.cen_lat,cen_lon=self.cen_lon)
        # if xy1D:
            # self.xx = self.xx[0,:]
            # self.yy = self.yy[:,0]

        # Aliases
        self.lat_0 = self.cen_lat
        self.lon_0 = self.cen_lon
        return

class WRF_native_grid:
    def __init__(self,fpath,resolution='i',xy1D=True):
        """Generates a basemap object for a WRF file's domain.

        """
        # keyword arguments to basemap generation
        self.kwargs = {'resolution':resolution}

        self.load_wrfout(fpath)
        self.load_attrs()
        self.load_corners()
        self.generate_basemap()

        if xy1D:
            self.xx = self.xx[0,:]
            self.yy = self.yy[:,0]

        self.combine_scopes()

    def load_wrfout(self,fpath):
        self.W = WRFOut(fpath)
        return

    def load_attrs(self,):
        self.kwargs['lat_0'] = self.W.nc.CEN_LAT
        self.kwargs['lon_0'] = self.W.nc.CEN_LON
        self.kwargs['lat_1'] = self.W.nc.TRUELAT1
        self.kwargs['lat_2'] = self.W.nc.TRUELAT2

    def load_corners(self):
        self.kwargs['llcrnrlon'] = self.W.lons[0,0]
        self.kwargs['llcrnrlat'] = self.W.lats[0,0]
        self.kwargs['urcrnrlon'] = self.W.lons[-1,-1]
        self.kwargs['urcrnrlat'] = self.W.lats[-1,-1]

    def generate_basemap(self):
        self.m = Basemap(projection='lcc',**self.kwargs)
        self.lons, self.lats, self.xx, self.yy = self.m.makegrid(
                        self.W.lons.shape[1],self.W.lons.shape[0],returnxy=True)

    def combine_scopes(self):
        """
        Put kwargs dictionary into self's namespace
        """
        self.__dict__ = dict(self.__dict__,**self.kwargs)

        # Alias of domain limits
        self.Wlim = self.llcrnrlon
        self.Slim = self.llcrnrlat
        self.Nlim = self.urcrnrlat
        self.Elim = self.urcrnrlon


class HRRR_native_grid(WRF_native_grid):
    """AS WRF_native_grid, but with some info about operational HRRR.
    """
    def __init__(self,fpath):
        super(HRRR_native_grid,self).__init__(fpath)

    def load_wrfout(self,fpath):
        """Over-riden from parent.
        """
        self.W = HRRR(fpath)

    def load_attrs(self,):
        """Over-riden from parent.

        TODO: get rid of this and just load from HRRR object.
        """
        self.kwargs['lat_0'] = 38.5
        self.kwargs['lon_0'] = -97.5
        self.kwargs['lat_1'] = 38.5
        self.kwargs['lat_2'] = 38.5


def create_new_grid(Nlim=None,Elim=None,Slim=None,Wlim=None,proj='merc',
                    lat_ts=None,resolution='i',nx=None,ny=None,
                    tlat1=30.0,tlat2=60.0,cen_lat=None,cen_lon=None,):
                    # lllon=None,lllat=None,urlat=None,urlon=None):
    """Create new domain for interpolating to, for instance.

    The following are mandatory arguments for mercator ('merc'):
    Nlim,Elim,Slim,Wlim = lat/lon/lat/lon limits for north/east/south/west respectively
    lat_ts = latitude of true scale?
    nx,ny = number of points in the x/y direction respectively
    """
    if proj == 'merc':
        if None in (Nlim,Elim,Slim,Wlim,lat_ts,nx,ny):
            print("Check non-optional arguments.")
            raise Exception
        m = Basemap(projection=proj,llcrnrlat=Slim,llcrnrlon=Wlim,
                    urcrnrlat=Nlim,urcrnrlon=Elim,lat_ts=lat_ts,resolution='h')
    elif proj == 'lcc':
        # if None in (tlat1,tlat2,cen_lat,cen_lon,lllon,lllat,urlon,urlat,nx,ny):
        if None in (tlat1,tlat2,cen_lat,cen_lon,Nlim,Elim,Slim,Wlim,nx,ny):
            print("Check non-optional arguments.")
            raise Exception
        m = Basemap(projection='lcc',lat_1=tlat1,lat_2=tlat2,lat_0=cen_lat,
                            lon_0=cen_lon,llcrnrlon=Wlim,llcrnrlat=Slim,
                            urcrnrlon=Elim,urcrnrlat=Nlim,resolution='h')

    lons, lats, xx, yy = m.makegrid(nx,ny,returnxy=True)
    return m, lons, lats, xx[0,:], yy[:,0]

def create_WRFgrid(f):
    """Constructs an instance of WRF_native_grid to return the vitals.

    f   -   absolute path to wrfout file.
    """
    W = WRF_native_grid(f)
    xx = W.xx
    yy = W.yy
    lons = W.lons
    lats = W.lats
    m = W.m
    return xx,yy,lats,lons,m

def reproject(data_orig,xx_orig=False,yy_orig=False,lats_orig=False,lons_orig=False,newgrid=False,
                    xx_new=False,yy_new=False,method='linear'):
    """
    data_orig               -   N.ndarray of original data (2D?)
    xx_orig,yy_orig         -   x-/y-axis indices of original data straight out a m(lons,lats)
    lats_orig,lons_orig     -   lats/lons of original data (shape?)
    newgrid                 -   basemap class classobject
    """
    # xx_new,yy_new = newgrid(lons_orig,lats_orig)
    assert len(xx_orig.shape) == 2
    if len(xx_new.shape) == 1:
        xx_new_dim = len(xx_new)
        yy_new_dim = len(yy_new)
        mx,my = N.meshgrid(xx_new,yy_new)
    elif len(xx_new.shape)== 2:
        xx_new_dim = len(xx_new[0,:])
        yy_new_dim = len(yy_new[:,0])
        mx,my = xx_new, yy_new
    else:
        raise Exception("Dimensions of xx_new and yy_new are too large.")

    # data_new = griddata((xx_orig.flat,yy_orig.flat),data_orig.flat,(xx_new.flat,
        # yy_new.flat)).reshape(xx_new_dim,yy_new_dim)
    # pdb.set_trace()
    print("Starting reprojection...")
    data_new = griddata((xx_orig.flat,yy_orig.flat),
                            data_orig.flat,(mx.flat,my.flat),
                            method=method).reshape(xx_new_dim,yy_new_dim)
                        # (mx.flat,my.flat),method=method)
    print("Completed reprojection.")
    return data_new

import os
import pdb

import numpy as N 
import cartopy.crs as ccrs
import pyproj
import shapely
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as N

from evac.datafiles.wrfout import WRFOut
from evac.utils.defaults import Defaults

import evac.utils.reproject_tools as reproject_tools
import evac.utils as utils

class Grid:
    """ Geographical grid information for plotting etc.

    Can take numerous datafiles or observation objects and return
    a common format to pass into BirdsEye, etc.

    A list of potential geographical variables:
        * nx,ny: number of grid points in x/y directions
        * Nlim, Elim, Slim, Wlim: depreciated. Use urcrnr/llcrnr.
        * urcrnr{lat,lon}, llcrnr{lat,lon}: lat/lon at extremes
        * cen_lat, cen_lon: central lat/lon of domain
        * lats, lons: 2D latitudes and longitudes

    Todo:
        * Redo with getter/setter methods or property decorators
            for attributes like cen_lat (and methods like get_cen_lat)

    Args:
        inst: instance of WRFOut, StageIV, etc. If not provided,
            other keyword arguments must be provided
        reproject_opts: dictionary of options for creating a new grid
    """
    def __init__(self,base=None,lats=None,lons=None,use_basemap=True):
        # This is here to avoid circular imports
        from evac.datafiles.wrfout import WRFOut
        from evac.datafiles.obs import Obs
        from evac.datafiles.stageiv import StageIV
        from evac.datafiles.obsgroup import ObsGroup
        from evac.datafiles.radar import Radar

        # If a group of obs are passed in, pick one randomly
        # so we can tell what class it is.
        if isinstance(base,ObsGroup):
            base_check = base.arbitrary_pick(return_instance=True)
        else:
            base_check = base

        # This is for debugging
        self.base_check = base_check

        if isinstance(base_check,dict):
            print("Creating new grid from reprojection options.")
            # self.xx, self.yy, self.lats, self.lons = self.create_grid(base)
            if use_basemap:
                self.create_grid_basemap(base)
            else:
                self.create_grid_cartopy(base)
        elif isinstance(base_check,WRFOut):
            print("Creating grid from WRFOut provided.")
            self.I = base
            self.load_wrfout_info()
        elif isinstance(base_check,StageIV):
            print("Creating grid from StageIV provided.")
            self.I = base
            self.load_stageiv_info()
        elif isinstance(base_check,Radar):
            self.I = base
            self.load_radar_info()
        elif (lats is not None) and (lons is not None):
            self.I = None
            self.load_user_info(lats,lons)
        else:
            raise NotImplementedError


        assert self.lats.ndim == 2
        assert self.lons.ndim == 2
        # self.nlats, self.nlons = self.lats.shape
        self.fill_attributes()

        if isinstance(self.lats,N.ma.core.MaskedArray):
            self.lats = self.lats.data
            self.lons = self.lons.data
        self.shape = self.lats.shape

    def get_corners(self):
        self.llcrnrlon = self.lons[0,0]
        self.llcrnrlat = self.lats[0,0]
        self.urcrnrlon = self.lons[-1,-1]
        self.urcrnrlat = self.lats[-1,-1]
        # self.llcrnrlat = self.lats.min()
        # self.llcrnrlon = self.lons.min()
        # self.urcrnrlat = self.lats.max()
        # self.urcrnrlon = self.lats.min()
        return

    def fill_attributes(self):
        if not hasattr(self,'llcrnrlat'):
            self.get_corners()
        if not hasattr(self,'lat_ts'):
            self.lat_ts = self.compute_lat_ts()
        if not hasattr(self,'xx'):
            self.xx, self.yy = self.m(self.lons,self.lats)
        if not hasattr(self,'nlons'):
            self.nlats, self.nlons = self.lats.shape
        if not hasattr(self,'Nlim'):
            self.Nlim = self.urcrnrlat
            self.Elim = self.urcrnrlon
            self.Slim = self.llcrnrlat
            self.Wlim = self.llcrnrlon
        return

    def load_user_info(self,lats,lons):
        if lats.ndim == 1:
            self.lons, self.lats = N.meshgrid(lons,lats)
        else:
            self.lats = lats
            self.lons = lons
        self.cen_lat, self.cen_lon = self.get_cen_latlon()
        self.get_corners()
        self.compute_lat_ts()
        self.nlats, self.nlons = self.lats.shape
        self.m, *_ = self.basemap_grid()
        return

    def compute_lat_ts(self):
        self.lat_ts = (self.urcrnrlat-self.llcrnrlat)/2
        return

    def basemap_grid(self):
        m, lons, lats, xx, yy = reproject_tools.create_new_grid(
                Nlim=self.urcrnrlat,
                Elim=self.urcrnrlon,Slim=self.llcrnrlat,
                Wlim=self.llcrnrlon,lat_ts=self.lat_ts,
                nx=self.nlons,ny=self.nlats,)
        return m, lons, lats, xx, yy

    def create_grid_basemap(self,opts):
        """ Basemap grid.
        """
        self.urcrnrlat = opts['urcrnrlat']
        self.urcrnrlon = opts['urcrnrlon']
        self.llcrnrlat = opts['llcrnrlat']
        self.llcrnrlon = opts['llcrnrlon']

        if 'nx' in opts:
            nx = opts['nx']
            ny = opts['ny']
            dx = None
            print("nx = {} and ny = {} requested".format(nx,ny))
        else:
            dx = opts['dx_km']
            nx = None
            ny = None
            print("dx = {} km requested.".format(dx))


        self.compute_lat_ts()

        self.m, self.lons, self.lats, self.xx, self.yy = reproject_tools.create_new_grid(
                Nlim=self.urcrnrlat,
                Elim=self.urcrnrlon,Slim=self.llcrnrlat,
                Wlim=self.llcrnrlon,lat_ts=self.lat_ts,
                nx=nx,ny=ny,dx=dx)

        # self.nx = self.xx.shape[1]
        # self.ny = self.yy.shape[0]
        self.nx = len(self.xx)
        self.ny = len(self.yy)
        self.nlats = self.ny
        self.nlons = self.nx

        self.dx = N.diff(self.xx).mean()
        self.dy = N.diff(self.yy).mean()
        print("Average dx = {:.1f}km and dy = {:.1f}km.".format(self.dx/1000,self.dy/1000))
        # pdb.set_trace()
        return

    def create_grid_cartopy(self,opts):
        """ Generate a new grid without basemap.

        For consistency with WRF and basemap, the upper right and
        lower left corners are specified. This reverses the logic
        in array creation, but the output arrays of lat, lon, x, y
        are all in the same format as WRF output.

        Note:
            The dictionary `opts` may have these key/values:
                * urcrnrlat, urcrnrlon, llcrnrlon, llcrnrlat (lat/lon extremes)
                * dx_km (grid point spacing)
        """
        # lat/lon projection
        proj_ll = pyproj.Proj(init='epsg:4326')

        # This one is used for Google Maps, AKA EPSG:900913
        proj = pyproj.Proj(init='epsg:3857')

        # Corners of new grid (approximately)
        urcrnr = shapely.geometry.Point((opts['urcrnrlon'],opts['urcrnrlat']))
        llcrnr = shapely.geometry.Point((opts['llcrnrlon'],opts['llcrnrlat']))

        # Project corners to target projection
        urx,ury = pyproj.transform(proj_ll, proj, urcrnr.x, urcrnr.y) 
        llx,lly = pyproj.transform(proj_ll, proj, llcrnr.x, llcrnr.y) 

        # Assume dx = dy
        dx = opts['dx_km'] * 1000
        dy = dx

        # Iterate over 2D area
        xx = []
        yy = []
        lats = []
        lons = []
        xlen = 0

        x = urx
        while x >= llx:
            y = ury
            while y >= lly:
                p = shapely.geometry.Point(pyproj.transform(proj, proj_ll, x, y))
                lons.insert(0,p.x)
                lats.insert(0,p.y)
                xx.insert(0,x)
                yy.insert(0,y)
                y -= dy
            xlen += 1
            x -= dx
        ylen = int(len(yy)/xlen)
        xarr = N.array(xx).reshape(xlen,ylen)
        yarr = N.array(yy).reshape(xlen,ylen)
        latarr = N.array(lats).reshape(xlen,ylen)
        lonarr = N.array(lons).reshape(xlen,ylen)
        print("Completed grid point computations.")

        # return xarr, yarr, latarr, lonarr 
        self.xx = xarr.T
        self.yy = yarr.T
        self.lats = latarr.T
        self.lons = lonarr.T
        return
        
    def get_cartopy_proj(self,proj):
        """ Set cartopy projection instance.
        """
        kwargs = dict()
        if proj == 'lcc':
            kwargs['central_longitude'] = self.cen_lon
            kwargs['central_latitude'] = self.cen_lat
        elif proj == 'merc':
            kwargs['central_longitude'] = self.cen_lon
            kwargs['min_latitude'] = self.lats.min()
            kwargs['max_latitude'] = self.lats.max()
        else:
            raise Exception
        return utils.get_ccrs_proj(proj,ckws=kwargs)

    def load_stageiv_info(self):
        self.lats = self.I.lats
        self.lons = self.I.lons
        self.cen_lat, self.cen_lon = self.get_cen_latlon()
        self.m = self.I.generate_basemap()
        return

    def load_radar_info(self):
        self.lats = self.I.lats
        self.lons = self.I.lons
        self.cen_lat, self.cen_lon = self.get_cen_latlon()
        self.m = self.I.generate_basemap()
        return

    def get_cen_latlon(self):
        nlat, nlon = self.lats.shape
        cen_lat = self.lats[int(nlat/2),int(nlon/2)]
        cen_lon = self.lons[int(nlat/2),int(nlon/2)]
        return cen_lat, cen_lon

    def load_wrfout_info(self):
        self.lat_0 = self.I.nc.CEN_LAT
        self.lon_0 = self.I.nc.CEN_LON
        self.lat_1 = self.I.nc.TRUELAT1
        self.lon_1 = self.I.nc.TRUELAT2
        self.llcrnrlon = self.I.lons[0,0]
        self.llcrnrlat = self.I.lats[0,0]
        self.urcrnrlon = self.I.lons[-1,-1]
        self.urcrnrlat = self.I.lats[-1,-1]
        self.cen_lat = self.I.cen_lat
        self.cen_lon = self.I.cen_lon
        # self.xx = self.I.xx
        # self.yy = self.I.yy

        # The lats/lons are masked array and these are awkward
        # self.lons = self.I.lons.data
        self.lons = self.I.lons
        # self.lats = self.I.lats.data
        self.lats = self.I.lats
        *_, self.m = reproject_tools.create_WRFgrid(self.I.fpath)
        return

    # def create_basemap(self,bmkwargs=None,proj='merc'):
        # if bmkwargs is None:
            # bmkwargs = {}
        # self.bmap = Basemap(projection=proj,**bmkwargs)
        # 
        # pass

    def convert_latlon_xy(self,lats,lons):
        if lats.ndim == 1:
            mlons, mlats = N.meshgrid(lons,lats)
        else:
            mlons = lons
            mlats = lats
        return self.m(mlons,mlats)

    def convert_latlon_xy_cartopy(self,lats,lons):
        """ Convert lats/lons to x and y.
        """
        xx = N.ones_like(lats)
        yy = N.ones_like(lons)

        # lat/lon projection
        # proj_ll = pyproj.Proj(init='epsg:4326')
        proj_ll = pyproj.Proj(proj='latlong',datum='WGS84')

        # cylindrical projection
        proj_xy = pyproj.Proj(init='epsg:3857')

        for i,j in N.ndindex(*xx.shape):
            # p.x is x, p.y is y

            # Pyproj suggests this is needed but it's not?
            # la = proj_ll(lats[i,j])
            # lo = proj_ll(lons[i,j])

            p = shapely.geometry.Point(pyproj.transform(proj_ll,proj_xy,lons[i,j],
                                        lats[i,j],))
            xx[i,j] = p.x
            yy[i,j] = p.y
        return xx,yy

    def interpolate(self,data,lats=None,lons=None,grid=None):
        """ Interpolate data and lat/lon grid to current grid.

        Note:
            User specifies grid or lat/lons.
        """
        if lats is None:
            lats = grid.lats
            lons = grid.lons

        data = utils.enforce_2d(data)
        # First cut before interpolating
        # cut_data, cut_lats, cut_lons = self.cut(data=data,lats=lats,lons=lons)
        # xx,yy = self.convert_latlon_xy(cut_lats,cut_lons)
        xx,yy = self.convert_latlon_xy(lats,lons)
        # data_reproj = reproject_tools.reproject(data_orig=cut_data,xx_orig=xx,
        data_reproj = reproject_tools.reproject(data_orig=data,xx_orig=xx,
                        yy_orig=yy,xx_new=self.xx,yy_new=self.yy)
        # pdb.set_trace()
        return data_reproj
        
    def cut(self,data,grid=None,lats=None,lons=None,return_grid=False):
        """ Trim data provided to current grid.
        Note:
            User specifies grid or lat/lons.
        """
        data = utils.enforce_2d(data)
        ld = self.get_limits()
        if lats is None:
            old_lats = grid.lats
            old_lons = grid.lons
        else:
            old_lats = lats
            old_lons = lons
        cut_data, cut_lats, cut_lons = utils.return_subdomain(
                            data=data,lats=old_lats,lons=old_lons,**ld)
        # import mpl_toolkits.basemap
        # cut_data = mpl_toolkits.basemap.interp(datain=data,xin=old_lons,
                            # yin=old_lats,xout=self.lons,yout=self.lats,
                            # checkbounds=True,)
        # pdb.set_trace()
        if return_grid:
            newgrid = Grid(lats=cut_lats,lons=cut_lons)
            return cut_data, newgrid
        return cut_data, cut_lats, cut_lons

    def get_limits(self,):
        limits = dict(
                    Nlim = self.Nlim,
                    Elim = self.Elim,
                    Slim = self.Slim,
                    Wlim = self.Wlim)
        return limits

    @staticmethod
    def interp_common_grid(Grid1,Grid2):
        """ Interpolate two Grid instances to new common grid.
        """
        newgrid = False
        return newgrid

    def __str__(self):
        infostr = "Grid, based on {} instance type.".format(self.base_check.__class__)
        return infostr

    def __repr__(self):
        infostr = "{} by {} grid, based on {} instance.".format(self.nlons,self.nlats,
                            self.base_check.__class__)
        return infostr

    def __call__(self,*args):
        return (self.lat_RBS(*args)[0][0],self.lon_RBS(*args)[0][0])

    def __enter__(self):
        """ Set up the `with` block compatibility.
        """
        from scipy.interpolate import RectBivariateSpline as RBS
        #yy,xx = [sorted(zz) for zz in N.indices(self.lats.shape)]
        yidx,xidx = N.indices(self.lats.shape)
        xx = xidx[0,:]
        yy = yidx[:,0]
        self.lat_RBS = RBS(x=xx,y=yy,z=self.lats)
        self.lon_RBS = RBS(x=xx,y=yy,z=self.lons)
        print("Ready to receive x/y points for interpolation")
        return self

    def __exit__(self,*args):
        """ 
        Todo:
            * Don't save if error?
        """
        del self.lat_RBS
        del self.lon_RBS
        print("Interpolation RBS closed.")
        pass


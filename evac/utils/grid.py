import os
import pdb

import numpy as N 
import cartopy.crs as ccrs
import pyproj
import shapely

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
    def __init__(self,base=None,lats=None,lons=None):
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
            self.create_grid(base)
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
        # self.nlats, self.nlons = self.lats.shape
        self.nlats, self.nlons = self.lats.shape

    def load_user_info(self,lats,lons):
        if lats.ndim == 1:
            self.lons, self.lats = N.meshgrid(lons,lats)
        else:
            self.lats = lats
            self.lons = lons
        self.cen_lat, self.cen_lon = self.get_cen_latlon()
        return

    def create_grid(self,opts):
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
        return

    def load_radar_info(self):
        self.lats = self.I.lats
        self.lons = self.I.lons
        self.cen_lat, self.cen_lon = self.get_cen_latlon()
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
        self.lons = self.I.lons
        self.lats = self.I.lats
        return

    # def create_basemap(self,bmkwargs=None,proj='merc'):
        # if bmkwargs is None:
            # bmkwargs = {}
        # self.bmap = Basemap(projection=proj,**bmkwargs)
        # 
        # pass

    def convert_latlon_xy(self,lats,lons):
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

    def interpolate_from(self,data,lats,lons):
        """ Interpolate data and lat/lon grid to current instance.
        """
        xx,yy = self.convert_latlon_xy(lats,lons)
        # yy,xx = self.convert_latlon_xy(lats,lons)
        data_reproj = reproject_tools.reproject(data_orig=data,xx_orig=xx,
                        yy_orig=yy,xx_new=self.xx,yy_new=self.yy)
        # pdb.set_trace()
        return data_reproj
        

    def interpolate_to(self,Grid2):
        """ Interpolate this grid to a different instance.
        """
        pass

    def trim_to(self,Grid2):
        """ Trim this grid to the lat/lon box of another instance.
        """
        pass

    @staticmethod
    def interp_common_grid(Grid1,Grid2):
        """ Interpolate two Grid instances to new common grid.
        """
        newgrid = False
        return newgrid


import os
import pdb

import cartopy.crs as ccrs

from evac.datafiles.wrfout import WRFOut
from evac.datafiles.obs import Obs
from evac.datafiles.stageiv import StageIV

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

    Args:
        inst: instance of WRFOut, StageIV, etc. If not provided,
            other keyword arguments must be provided
        reproject_opts: dictionary of options for creating a new grid
    """
    def __init__(self,inst=None,reproject_opts=None):
        self.I = inst

        if isinstance(self.I,WRFOut):
            self.load_wrfout_info()
        elif (self.I is None) and (isinstance(reproject_opts,dict)):
            self.create_grid(reproject_opts)

    def create_grid(self,opts):
        """ Generate a new grid.

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

        # lats = N.ones([],dtype='f8')
        # lons = lats.copy()

        # Starting on upper row
        # This is the NE point, in northern hemisphere context
        # x = pyproj.transform(proj_ll, proj, urcrnr.x, urcrnr.y)
        # xx = N.linspace(

        # This method seems inefficient...
        # Project corners to target projection
        s = pyproj.transform(proj_ll, proj, urcrnr.x, urcrnr.y) # Transform NW point to 3857
        e = pyproj.transform(proj_ll, proj, llcrnr.x, llcrnr.y) # .. same for SE

        # Assume dx = dy
        dx = opts['dx_km'] * 1000
        dy = dx

        # Iterate over 2D area
        gridpoints = []
        x = s[0]
        # while x < e[0]:
        while x > e[0]:
            y = s[1]
            # while y < e[1]:
            while y > e[1]:
                p = shapely.geometry.Point(pyproj.transform(proj, proj_ll, x, y))
                gridpoints.append(p)
                y += dy
            x += dx
        pdb.set_trace()

        # return newgrid
        


    def get_cartopy_proj(self,proj):
        """ Set cartopy projection instance.
        """
        kwargs = dict()
        if proj == 'lcc':
            kwargs['central_longitude'] = self.cen_lon
            kwargs['central_latitude'] = self.cen_lat
            P = ccrs.LambertConformal(**kwargs)
        return P

    def load_wrfout_info(self):
        self.lat_0 = self.I.nc.CEN_LAT
        self.lon_0 = self.I.nc.CEN_LON
        self.lat_1 = self.I.nc.TRUELAT1
        self.lon_1 = self.I.nc.TRUELAT2
        self.llcrnrlon = self.I.lons[0,0]
        self.llcrnrlat = self.I.lats[0,0]
        self.urcrnrlon = self.I.lons[-1,-1]
        self.urcrnrlat = self.I.lats[-1,-1]
        self.xx = self.I.xx
        self.yy = self.I.yy
        self.lons = self.I.lons
        self.lats = self.I.lats
        return

    # def create_basemap(self,bmkwargs=None,proj='merc'):
        # if bmkwargs is None:
            # bmkwargs = {}
        # self.bmap = Basemap(projection=proj,**bmkwargs)
        # 
        # pass

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


import os
import pdb

import numpy as N 
import cartopy.crs as ccrs
import pyproj
import shapely

from evac.datafiles.wrfout import WRFOut
from evac.datafiles.obs import Obs
from evac.datafiles.stageiv import StageIV
import evac.utils.reproject_tools as reproject_tools

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
    def __init__(self,base=None):
        if isinstance(base,dict):
            print("Creating new grid from reprojection options.")
            # self.xx, self.yy, self.lats, self.lons = self.create_grid(base)
            self.create_grid(base)
        elif isinstance(base,WRFOut):
            print("Creating grid from WRFOut provided.")
            self.I = base
            self.load_wrfout_info()
        elif isinstance(base,StageIV):
            print("Creating grid from StageIV provided.")
            self.I = base
            self.load_stageiv_info()
        else:
            raise NotImplementedError

        assert self.lats.ndim == 2


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
        self.xx = N.swapaxes(xarr,0,1)
        self.yy = N.swapaxes(yarr,0,1)
        self.lats = N.swapaxes(latarr,0,1)
        self.lons = N.swapaxes(lonarr,0,1)
        return
        
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

    def convert_latlon_xy(self,lats,lons):
        """ Convert lats/lons to x and y.
        """
        xx = N.ones_like(lats)
        yy = N.ones_like(lons)

        # lat/lon projection
        proj_ll = pyproj.Proj(init='epsg:4326')

        # cylindrical projection
        proj_xy = pyproj.Proj(init='epsg:3857')

        for i,j in N.ndindex(*xx.shape):
            # p.x is x, p.y is y
            p = shapely.geometry.Point(pyproj.transform(proj_ll, proj_xy, 
                                lats[i,j], lons[i,j]))
            xx[i,j] = p.x
            yy[i,j] = p.y
        return xx,yy

    def interpolate_from(self,data,lats,lons):
        """ Interpolate data and lat/lon grid to current instance.
        """
        xx,yy = self.convert_latlon_xy(lats,lons)
        data_reproj = reproject_tools.reproject(data=data,newgrid_xx=xx,
                        newgrid_yy=yy,self.xx,self.yy)
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


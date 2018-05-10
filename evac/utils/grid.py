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
    """
    def __init__(self,inst=None):
        self.I = inst

        if isinstance(self.inst,WRFOut):
            self.load_wrfout_info()

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


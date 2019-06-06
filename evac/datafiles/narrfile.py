import os
import pdb

import numpy as N

from evac.datafiles.gribfile import GribFile

class NARRFile(GribFile):
    """
    corners:
        # ll:
        * 1.000001N, 145.5W
        # lr:
        * 0.897945N, 68.32005W
        # ur:
        * 46.3544N, 2.569891W
        * 46.63433N, 148.6418E
    """
    def __init__(self,fpath):
        super().__init__(fpath)
        # pdb.set_trace()

    def return_latlon(self):
        gg = self.arbitrary_pick()
        latlon = gg.latlons()
        lats, lons = latlon
        # pdb.set_trace()
        return lats,lons

    def idx_from_lv(self,vrblkey,lv=None):
        """ Return the index to load a requested level (in hPa).
        If level doesn't exist, throw an exception.
        """
        rows = N.where(self.available_fields_array == vrblkey)[0]
        if len(rows) == 1 and lv == None:
            return 0
        lv_choice = self.available_fields_array[rows,5]

        # JRL: below might need moving to a child class (NARR, etc)
        if lv is None:
            lv_str = lv_choice[0]
        elif isinstance(lv,str):
            lv_str = lv
        else:
            raise exception

        for nl,l in enumerate(lv_choice):
            if lv_str in l:
                return nl

    def lookup_vrbl(self,vrbl):
        """
        Going to use NCEP abbreviations from "table 131" below.

        # These two tables, though something is wrong with pygrib? Grib1v2?
        # https://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html#TABLE131
        # https://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html#TABLE2
        """
        LOOKUP = {
                "HGT":"Total precipitation anomaly of at least 10 mm",# geopot in gpm
                "TMP":"11", # drybulb in K
                "UGRD":"33", #u-wind in m/s
                "VGRD":34, #v-wind in m/s
                "VVEL":39, # vertical velocity is Pa/sec
                "SPFH":51, # specific humidity in kg/kg
                "MCONV":135, # is this moisture convergence? kg/kg/s
                "CLWMR":153, # cloud water in kg/kg
                "TKE":158, #TKE in J/kg
                "ICMR":178, # ice mixing ratio in kg/kg
                    }
        return LOOKUP[vrbl]

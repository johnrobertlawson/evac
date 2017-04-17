"""Generic superclass for data file. Could be netCDF, grib1,
grib2...

This is where derived functions should be attached (computing vorticity, etc)
"""

class DataFile:
    def __init__(self,fpath):
        self.fpath = fpath

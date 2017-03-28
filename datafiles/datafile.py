"""Generic superclass for data file. Could be netCDF, grib1,
grib2...

"""

class DataFile:
    def __init__(self,fpath):

        self.fpath = fpath

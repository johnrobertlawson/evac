from netCDF4 import Dataset

from .datafile import DataFile

class NCFile(DataFile):
    def __init__(self,fpath):
        """Generic netCDF import.

        Subclass of generic DataFile class.
        """
        super().__init__(fpath)
        self.nc = Dataset(fpath,'r')

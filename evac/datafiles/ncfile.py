import os
import pdb

from netCDF4 import Dataset

from evac.datafiles.datafile import DataFile

class NCFile(DataFile):
    """Generic netCDF class.

    Subclass of generic DataFile class.
    """
    def __init__(self,fpath):
        super().__init__(fpath)
        self.nc = Dataset(fpath,'r')

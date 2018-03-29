import pdb
import os

from evac.datafiles.datafile import DataFile

class PNGFile(DataFile):
    """ Generic superclass to read in png data.
    """
    def __init__(self,fpath):
        self.fpath = fpath

import os
import pdb

from evac.datafiles.gribfile import GribFile

class GFS(GribFile):
    """ Class to hold GFS data (in GRIB format).

    Args:
        fpath: Absolute path to data file.
    """
    def __init__(self,fpath):
        """Initialise GFS object, a child of GribFile,
        grandchild of DataFile.
        """
        super().__init__(fpath)

    def lookup_vrbl(self,vrbl):
        LOOKUP = {}
        #LOOKUP['accum_precip'] = {'key':'','idx':0}
        LOOKUP['Td2'] = {'key':'2 metre dewpoint temperature','idx':0}

        return LOOKUP[vrbl]['key'], LOOKUP[vrbl]['idx']

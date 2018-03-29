import os
import pdb

from evac.datafiles.gribfile import GribFile

class GEFS(GribFile):
    """ Class to hold GEFS data (grib format).

    Args:
        fpath: Absolute path to GEFS file.
    """
    def __init__(self,fpath):
        """Initialise GEFS object, a child of GribFile,
        grandchild of DataFile.
        """
        super().__init__(fpath)

    def lookup_vrbl(self,vrbl):
        """ Table to find GEFS variable key from WRF-style codes.
        """
        LOOKUP = {}
        #LOOKUP['accum_precip'] = {'key':'','idx':0}
        LOOKUP['Td2'] = {'key':'2 metre dewpoint temperature','idx':0}
        LOOKUP['T2'] = {'key':'Temperature','level':'2000','idx':0}
        LOOKUP['T'] = {'key':'Temperature','idx':0}

        return LOOKUP[vrbl]['key'], LOOKUP[vrbl]['idx']

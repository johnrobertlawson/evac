import os
import pdb

class Obs:
    """
    Analogue of :any:`evac.datafiles.ensemble` but for observation data across scales, 
    variables, and times.
    
    Args:
        fpath: Absolute path to a data set.

    Todo:
        * Everything!
    """
    def __init__(self,fpath):
        self.fpath = fpath

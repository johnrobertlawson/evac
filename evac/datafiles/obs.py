import os
import pdb

import evac.utils as utils

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


    def get_subdomain(self,Nlim,Elim,Slim,Wlim,data='self',overwrite=False):
        """
        Return data array between bounds

        If overwrite is True, replace class data with new subdomain

        
        """
        if data == 'self':
            data = self.data
        data,lats,lons = utils.return_subdomain(data,self.lats,
                                self.lons,Nlim,Elim,Slim,Wlim)

        if overwrite:
            self.lats = lats
            self.lons = lons
            self.data = data
            return
        else:
            return data,lats,lons

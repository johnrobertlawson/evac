import pdb

from evac.datafiles.wrfout import WRFOut

class DomainGrid:
    """ Create a common class for verification domains.
    
    Args:
        inherit: dynamic inheritance depending on the
            domain that is being passed in.
    """
    def __init__(self,inherit=None):
        if isinstance(inherit,WRFOut):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,WRF_native_grid):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,VerifGrid):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,str):
            inherit = WRFOut(inherit)
            self.__dict__ = inherit.__dict__.copy()

            # self.lat_0 = inherit.nc.CEN_LAT
            # self.lon_0 = inherit.nc.CEN_LON
        else:
            raise Exception("Inherited object is {}".format(type(inherit)))
        return

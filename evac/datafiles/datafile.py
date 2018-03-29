class DataFile:
    """Generic superclass for data file. 
    
    This could be netCDF, grib1, grib2, etc.

    Args:
        fpath: Absolute path to file.

    Todo:
        * This is where derived functions should be attached
            (computing vorticity, etc) that are in :func:`~evac.derived.derived`.
        * Also, see :func:`evac.derived.derived`.
        * Also, see :any:`~evac.derived.derived`.
        * Also, see :any:`evac.derived.derived`.
        * Also, see :method:`~evac.derived.derived`.
        * Also, see :method:`evac.derived.derived`.
    """
    def __init__(self,fpath):
        self.fpath = fpath

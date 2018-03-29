import pdb
import os

from evac.datafiels.binaryfile import BinaryFile

class MRMS(BinaryFile):
    """Multi-radar, multi-sensor data. Converts from their custom
    binary data format.
    
    Todo:
        * A helper function to automatically generate file name
            to load needs to be written.

    Args:
        fpath (str, optional): file path to data. If False, other
            information must be presented to search for data (below).
        rootdir (str, optional): if fpath is False, this is required
            to search for data based on product and time
        product (str, optional): if fpath is False, this is required
            to find correct product. Pick from PRECIPRATE. etc
        utc (datetime.datetime, optional): if fpath is False, this
            is required to find data file for valid time.
    """
    def __init__(self,fpath=False,rootdir=False,product=False,
                    utc=False):
        if not fpath:
            fpath = self.generate_fpath(rootdir,product,utc)

        super().__init__(fpath)

    def generate_fpath(self,rootdir,prod,utc):
        """Return an absolute path based on product and time of data
        that is required, in a root directory.
        """
        # Fill this out
        valid_products = ('PRECIPRATE',)
        pdb.set_trace()
        if prod not in valid_products:
            raise Exception("Product name not valid.")

        fname = '{0}.{1}{2}{3}.{4}{5}{6}'.format(prod,utc.year,utc.month,
                                utc.day,utc.hour,utc.minute,utc.second)
        fpath = os.path.join(rootdir,fname)
        return fpath

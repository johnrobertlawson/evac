from evac.datafiles.datafile import DataFile

class BinaryFile(DataFile):
    """ Superclass for binary files.

    Args:
        fpath: Absolute path to binary file. 

    Todo:
        * Nothing yet
        * Not using MRMS yet
    """
    def __init__(self,fpath):
        super().__init__(fpath)


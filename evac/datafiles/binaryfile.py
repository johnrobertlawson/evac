from .datafile import DataFile

class BinaryFile(DataFile):
    def __init__(self,fpath):
        """ Superclass for binary files.

        Args:
            fpath: Absolute path to binary file. 

        Todo:
            * Nothing yet
            * Not using MRMS yet
        """
        super().__init__(fpath)

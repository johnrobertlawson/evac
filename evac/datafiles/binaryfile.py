""" Any miscellaneous binary file.
"""
from .datafile import DataFile

class BinaryFile(DataFile):
    def __init__(self):
        super().__init__(fpath)

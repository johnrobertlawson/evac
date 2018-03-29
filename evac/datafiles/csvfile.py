from evac.datafiles.datafile import DataFile

class CSVFile(DataFile):
    """ Comma separated value data type.

    Args:
        fpath: Absolute file path to the CSV file.
    """
    def __init__(self,fpath):
        super().__init__(fpath)

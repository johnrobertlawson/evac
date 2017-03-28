"""
Analogue of ensemble.py but for observation data across scales, variables,
and times.
"""

class Obs:
    """
    An instance represents a data set.
    """
    def __init__(self,fpath):
        self.fpath = fpath

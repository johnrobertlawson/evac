""" Create a catalogue of objects, either in numpy array or dictionary format.

Catalogue may/may not allow for matching between different sets, so subcataloging?

May just be a function that returns a dataframe or labelled numpy array.
"""

class CombinedObjects:
    def __init__(self,utcs,objs,lats,lons,):
        self.utcs = utcs
        self.objs = objs
        # TI calculation got moved to catalogue...

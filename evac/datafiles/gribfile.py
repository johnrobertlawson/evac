import os
import pdb

import numpy as N

try:
    import pygrib
except ModuleNotFoundError:
    print("Not using pygrib")
    pygrib = None

from evac.datafiles.datafile import DataFile

class GribFile(DataFile):
    """Generic load for Grib 1/2 files.

    Args:
        fpath: Absolute path to grib file.

    Todo:
        * Make general, i.e. remove HRRR etc specifics
    """
    def __init__(self,fpath):
        super().__init__(fpath)
        # self.fpath = fpath
        self.G = pygrib.open(self.fpath)
        self.available_records = [gg for gg in self.G]
        self.G.seek(0) # Screw you GRIB!
        self.available_fields_list = [str(gg).split(':') for gg in self.G]
        # pdb.set_trace()
        self.available_fields_array = N.array(self.available_fields_list)
        # Remove weird variables that aren't variables
        self.available_fields_unfiltered = N.unique(self.available_fields_array[:,1])
        delrows = []
        for row in N.arange(self.available_fields_unfiltered.shape[0]):
            if len(self.available_fields_unfiltered[row]) < 5:
                delrows.append(row)
        self.available_fields = N.delete(self.available_fields_unfiltered,delrows)
        self.projection()

        self.y_dim, self.x_dim = self.lats.shape
        # Only one time per file, right?
        self.t_dim = 1
        self.get_units()

    def get_units(self,):
        """ The units used for variables within.
        This should be overridden in children.
        """
        self.UNITS = {'pressure':'Pa',}
        return

    def get_record(self,vrblkey,idx=None,level=None):
        self.G.seek(0)
        ngg = len(self.G.select(name=vrblkey))
        # print(vrbl)
        # pdb.set_trace()
        print("There are {0} entries for the key {1}.".format(ngg,vrblkey))
        # Not sure why this is needed - random wrong indices popping up
        if ngg == 1:
            idx = 0
        else:
            idx = self.idx_from_lv(vrblkey,lv=level)
        gg = self.G.select(name=vrblkey)[idx]
        print("Picking index {0}.".format(idx))
        return gg

    def idx_from_lv(self,vrblkey,lv=None):
        """ Return the index to load a requested level (in hPa).
        If level doesn't exist, throw an exception.
        """
        rows = N.where(self.available_fields_array == vrblkey)[0]
        # idxs = self.available_fields_array[rows,5]
        if len(rows) == 1 and lv == None:
            # return idxs[0]
            return 0
        lv_choice = self.available_fields_array[rows,5]
        # Just assume Pa
        lv_str = '{:d}'.format(lv*100)
        # pdb.set_trace()
        for nl,l in enumerate(lv_choice):
            if lv_str in l:
                return nl


    def get(self,vrbl,idx=0,utc=None,level=None,lons=None,lats=None):
        """
        I don't think we need utc/lv - aren't all grib files one time? - but
        it is here to enable compatibility with other get() APIs.
        """
        # TODO: look up level in available fields array, return index
        vrblkey, idx = self.lookup_vrbl(vrbl)
        print("Variable {0} has key {1} in Grib data.".format(vrbl,vrblkey))
        # gg = self.get_record(vrblkey,idx=idx)
        # idx lookup should happen in get_record
        gg = self.get_record(vrblkey,level=level)
        arr = gg.values
        return self.enforce_4d(arr)
        # pdb.set_trace()
        # return arr

    def enforce_4d(self,arr):
        """TODO: move higher up to super-superclass?
        """
        if len(arr.shape) == 2:
            arr4d = arr[N.newaxis,N.newaxis,:,:]
        else:
            raise Exception("Why is this data not 2d?")
        return arr4d

    def lookup_vrbl(self,vrbl):
        LOOKUP = {}
        # Children should use the right keys/indices here
        return (vrbl, 0)

    def return_latlon(self):
        gg = self.arbitrary_pick()
        latlon = gg.latlons()
        lats, lons = latlon
        return lats,lons

    def projection(self):
        # self.m = Basemap(projection='npstere',lon_0=-105.0,#lat_1=60.0,
                # llcrnrlon=lllon,llcrnrlat=lllat,urcrnrlon=urlon,urcrnrlat=urlat,
                            # boundinglat=24.701632)
        self.lats, self.lons = self.return_latlon()
        # self.xx, self.yy = self.m(self.lons,self.lats)
        # pdb.set_trace()
        # self.mx, self.my = N.meshgrid(self.xx,self.yy)

        # lllon = -119.023
        self.lllon = self.lons[0,0]
        # lllat = 23.117
        self.lllat = self.lats[0,0]
        # urlon = -59.9044
        self.urlon = self.lons[-1,-1]
        # urlat = 45.6147234
        self.urlat = self.lats[-1,-1]

        self.shape = self.lats.shape
        assert self.lats.shape == self.lons.shape

    def arbitrary_pick(self):
        vrbl= self.available_fields[0]
        gg = self.get_record(vrbl,idx=0)
        return gg

    def return_vrbl_records(self,vrbl):
        vidx = N.where(self.available_fields_array[:,1]==vrbl)[0]
        for idx in vidx:
            print(self.available_fields_array[idx,:])
        return

    def search_keyword(self,keyword):
        """ Look for record entries based on a keyword.
        """
        for entry in self.available_fields_list:
            for x in entry:
                if keyword in x:
                    print(entry)
                    break
        return

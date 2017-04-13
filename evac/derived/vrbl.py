"""Subclass of numpy ndarray. Vrbl holds a 4D array for a variable
in the dimensions (time, level, lats, lons), similarly to WRF output.

The subclassing allows attachment of methods that compute derived fields
(e.g. vorticity) and inspection of metadata (dimension names etc.).

Example:

arr = Dataset('wrfout.nc').variables['T'][0,0,:,:]

gives you a normal numpy array. To get a Vrbl() subclass, go via
a `datafiles/` module, for instance:

arr = WRFOut('wrfout.nc').get('T',utc=0,lv=0)

This, when completed, will give a Vrbl array with attributes such as
dimensions.
"""
import numpy as N

class Vrbl(N.ndarray):
        #def __new__(cls,*args,**kwargs):
        def __new__(cls,arr):
            #return N.ndarray.__new__(cls,*args,**kwargs)
            obj = N.asarray(arr).view(cls)

            # This is where we'd add attributes to our subclass instance
            # obj.attr = kwarg1

            return obj

        def __array_finalize__(self,obj):
            """ This is always called (by superclass).

            1. If constructed explicitly, obj type is None
            2. If constructed via casting (e.g., N.view), obj type is N.ndarray
            3. If constructed via slicing, obj type is the class Vrbl.

            In case #2, we return a superclass (subject to change...)

            """
            if obj is None:
                return

            # Add attributes from __new__ call
            # self.attr = getattr(obj,'attr',None)

            # All attributes go here.
            # self.dims = self.get_dims()

        def holder(self):
            pass

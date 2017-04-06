"""Derived fields application to severe/convective weather.

Todos:
    * Make methods general to all data, but also plug-ins to 
        datatypes.
    * Do all computations need to be 2D - 5D? What about ensemble members?
    * Could create data type that contains an array, the names of the fields,
        and inherits methods that act on it. (Slicing etc)
"""

import numpy as N

import evac.derived.general as general
import evac.utils as utils

def compute_updraught_helicity(u,v,w,dz=None,z=None,z0=2000,z1=5000):
    """Eq. 11 & 12 from Kain et al 2008, WAF.

    Requires either z or dz to be specified.

    Args:
        u,v,w (N.ndarray)           :   3-D fields of wind, all same size.
        dz (N.ndarray)              :   1-D vector of difference in z
        z (N.ndarray)               :   1-D vector height above ground level (m)
        z0 (int,float,optional)     :   bottom level in m
        z1 (int,float,optional)     :   top level in m
    """
    utils.enforce_same_dimensions(u,v,w)
    
    vort_z = general.compute_vorticity(u,v)
    # Need to do mid-point approximation here
    # If dz, don't need to do anything.
    # If z, need to do differences.
    UH = (vort_z * w) * dz

    # Logic for smoothing?

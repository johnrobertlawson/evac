"""Scripts to check the nature of numbers.
"""
import pdb

import numpy as N

def check_type(n,typ):
    """ Return true if input n is a given 'type'. If n is an array, it 
    checks if the dtype of the array matches the type. If the array is
    made up of many types, it raises an error.

    Arguments:
        n           :   Any object (N.ndarray, int, N.int32, N.float64...)
        typ (str)   :   'int' is any type or dtype resembling an integer
                        'float' is any type or dtype resembling a float

    TODO: raise error if array has more than one dtype.
    """
    if is_array(N.ndarray):
        ntyp = arr.dtype
    else:
        ntyp = type(n)

    if typ == 'int':
        return N.issubdtype(ntyp,N.int)
    elif typ == 'float':
        return N.issubdtype(ntyp,N.float)

def is_array(n):
    return isinstance(n,N.ndarray)

def is_integer(n):
    return check_type(n,'int')

def has_integers(n):
    return check_type(n,'int')

def has_floats(n):
    return check_type(n,'float')

def is_float(n):
    return check_type(n,'float')

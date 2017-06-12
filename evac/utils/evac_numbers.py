"""Scripts to check the nature of numbers.
"""
import pdb

import numpy as N

def check_type(n,typ,allow_array=True):
    """ Return true if input n is a given 'type'. If n is an array, it
    checks if the dtype of the array matches the type. If the array is
    made up of many types, it raises an error.

    Arguments:
        n           :   Any object (N.ndarray, int, N.int32, N.float64...)
        typ (str)   :   'int' is any type or dtype resembling an integer
                        'float' is any type or dtype resembling a float

    TODO: raise error if array has more than one dtype.
    """
    if is_array(N.ndarray) and allow_array:
        ntyp = arr.dtype
    else:
        if not allow_array:
            return False
        else:
            ntyp = type(n)

    if typ == 'int':
        return N.issubdtype(ntyp,N.int)
    elif typ == 'float':
        return N.issubdtype(ntyp,N.float)
    elif typ == 'int_or_float':
        checkint = N.issubdtype(ntyp,N.int)
        checkfloat = N.issubdtype(ntyp,N.float)
        return max(checkint,checkfloat)

def is_array(n):
    return isinstance(n,N.ndarray)

def is_integer(n,allow_array=False):
    return check_type(n,'int',allow_array=allow_array)

def has_integers(n):
    return check_type(n,'int')

def has_floats(n):
    return check_type(n,'float')

def is_float(n,allow_array=False):
    return check_type(n,'float',allow_array=allow_array)

def is_number(n,allow_array=False):
    return check_type(n,'int_or_float',allow_array=allow_array)

def has_numbers(n,):
    return check_type(n,'int_or_float')

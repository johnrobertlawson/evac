""" Custom type checks.
"""

import numpy as N

from evac.utils.exceptions import FormatError

def is_list(obj):
    return isinstance(obj,(tuple,list,N.ndarray))

def enforce_same_dimensions(*args):
    if not check_same_dimensions(*args):
        raise FormatError("Input arrays not all same size.")
    return

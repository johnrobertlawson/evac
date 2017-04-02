""" Custom type checks.
"""

import numpy as N

def is_list(obj):
    return isinstance(obj,(tuple,list,N.ndarray))

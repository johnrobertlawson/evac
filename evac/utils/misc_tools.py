""" Miscellaneous functions like decorators.
"""

import pdb
import time

def time_me(f):
    """ Decorator for timing functions.

    Usage:

    @time_me
    def function():
        pass
    """
    def timed(*args, **kwargs):
        t0 = time.time()
        result = f(*args, **kwargs)
        t1 = time.time()

        print('function {} took {:2.3f} sec'.format(f.__name__,t1-t0))
        return result

    return timed

def loop_me(*args):
    """
    TODO: a decorator function that loops an existing function
    over all arguments passed.

    args should be a list of 1D numpy arrays, tuples, or lists
    """
    pass

#import datafiles
#import derived
#import lazy
#import plot
#import stats
#import stoch
#import utils




__all__ = ["evac",]
from evac.evac import *
"""
import pkgutil
import inspect

for loader, name, is_pkg in pkgutil.walk_packages(__path__):
    module = loader.find_module(name).load_module(name)

    for name, value in inspect.getmembers(module):
        if name.startswith('__'):
            continue

        globals()[name] = value
        __all__.append(name)
"""


__all__ = ["datafiles","derived","lazy","plot","stats","stoch","utils"]

"""
# The following should skip failed imports (mainly for PyGrib)
import builtins
from types import ModuleType

class DummyModule(ModuleType):
    def __getattr__(self, key):
        return None
    __all__ = []   # support wildcard imports




#def tryimport(name, globals={}, locals={}, fromlist=[], level=5):
#def tryimport(name, globals=globals(), locals=locals(), fromlist=[], level=1):
def tryimport(*args,**kwargs):
#def tryimport(name, globals=None, locals=None, fromlist=[], level=10):
    #if kwargs['name'] == 'pygrib':
    if args[0] == 'pygrib':
        try:
            return realimport(*args,**kwargs)
            #return realimport(name, globals(), locals(), fromlist, level)
        #except ImportError:
        except:
            return DummyModule(args[0])
    else:
        #return realimport(name,globals,locals,fromlist,level=level)
        return realimport(*args,**kwargs)

realimport, builtins.__import__ = builtins.__import__, tryimport

#import .datafiles
#import derived
#import stats
#__all__ = []

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

import pdb

import evac.utils as utils

"""Custom exceptions.

Todo:
    * Extend the functionality
    * Make the wowprint() actually work
    * https://github.com/UltimateHackers/hue

"""

class PrettyException(Exception):
    """ Subclassed error that prints real pretty.
    """
    def __init__(self,message,color=False,bold=False,
                    underline=False,formatargs=False):
        wowmessage = utils.wowprint(message,color=color,bold=bold,
                        underline=underline,formatargs=formatargs)
        super().__init__(wowmessage)

class FormatError(PrettyException):
    """Data is the wrong format or shape."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

class QCError(PrettyException):
    """Data has failed a quality control check.
    """
    def __init__(self,error='Quality Control error.',pass_idx=False):
        """
        Optional:
        Error (str) -   Message to send to user
        pass_idx (tuple, list, N.ndarray)   -   Location in data of bad vals.
        """
        print(error)
        if pass_idx is not False:
            if isinstance(pass_idx,(list,tuple)):
                print("Bad data is here: \n")
                [print(p) for p in pass_idx]
            else:
                print(pass_idx)

class NonsenseError(PrettyException):
    """ Data or argument doesn't make any sense.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

class WRFError(PrettyException):
    """ WRF run has broken somewhere. """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

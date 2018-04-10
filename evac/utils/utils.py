""" General miscellaneous functions.

Todo:
    * Move many functions in `unix_tools` and `gis_tools` into this module.
"""
import pdb
import datetime
import pytz

from pathlib import PosixPath
from evac.utils.exceptions import NonsenseError

def bridge_multi(orders,fromlist=False,todir=False,
                    return_newlist=False):
    """
    Soft-link, copy, or move multiple items.

    Args:
    orders      :   dictionary with the following layout:
                    dict( (frompath,topath)=cmd)
                    If orders is string, treat it as command,
                    and require arguments fromlist and todir.

    where

    cmd is 'copy','softlink', or 'move'
    frompath/topath is absolute path of origin/destination

    return_newlist (bool): if True, return a list with the new absolute
                                path of all files.
    """
    # Get into correct format
    if isinstance(orders,str):
        cmd = orders
        orders = {}
        assert cmd in ('copy','move','softlink')
        for f in fromlist:
            orders[(f,todir)] = cmd

    newlist = []
    for (frompath,topath),cmd in orders.items():
        newf = bridge(cmd,frompath,topath,return_newpath=True)
        newlist.append(newf)

    if return_newlist:
        return newlist
    return

def bridge(command,frompath,topath,catch_overwrite=False,
                return_newpath=False):
    """ Copy, move, soft-link.

    Args:

    command     :   (str) - must be one of "copy", "softlink",
                    "move".
    frompath    :   (str,Path) - absolute path to origin directory or file
    topath      :   (str,Path) - absolute path to destination directory or
                    file. If topath is only a directory, the link will
                    be the same filename as the origin.
    catch_overwrite     :   (bool) - if True, raise Exception if
                            topath already exists
    return_newpath (bool): if True, return the new absolute path of the file
    """

    # Use path-like objects
    frompath = enforce_pathobj(frompath)
    topath = enforce_pathobj(topath)

    if topath.is_dir():
        topath = topath / frompath.name

    if catch_overwrite and topath.exists():
        raise Exception("{} already exists.".format(topath))

    trycreate(topath)

    if command == "copy":
        # This == ridiculous, but that's what I get for using sodding pathlib
        topath.write_bytes(frompath.read_bytes())
    elif command == "move":
        # This replaces an existing 'topath' silently
        frompath.rename(topath)
    elif command == "softlink":
        frompath.symlink_to(topath)
    else:
        raise NonsenseError("The command variable **{}** is invalid".format(
                        command),color='red')

    if not topath.exists():
        print("File not created at",topath)
    else:
        print("File created at",topath)

    if return_newpath:
        return topath

def wowprint(message,color='red',bold=True,underline=False,
                formatargs=False,marker='**'):
    """ Print with colour/bold emphasis on certain things!

    message      :   (str) - string with any emphasised
                    item surrounded by two asterisks on both
                    sides. These are replaced.
    color        :   (str,bool) - the colour of the emphasised
                    text. If False, don't emphasise.

    Usage:

    wowprint("The program **prog** is broken",color='red')
    wowprint("The program {} is broken",formatargs=('prog',),color='blue',
                            underline=True)
    """
    # If nothing is set, pass through without changes
    if not any([color,bold,underline,formatargs]):
        print(message)
        return

    def findindex(message,marker='**'):
        try:
            idx = message.index(marker)
        except ValueError:
            raise Exception("The wowprint string needs to have two"
                                " lots of double asterisks.")
        else:
            return idx

    colors = {}
    colors['red'] = "\033[91m" # fail
    colors['purple'] = "\033[95m" # header
    colors['blue'] = "\033[94m" # OK
    colors['yellow'] = "\033[93m" # warning
    colors['green'] = "\033[92m" # OK

    colors['bold'] = "\033[1m" # bold
    colors['underline'] = "\033[4m" # underline

    endcode = "\033[0m" # use at the end, always

    codes = []
    if bold:
        codes.append(colors['bold'])
    if underline:
        codes.append(colors['underline'])
    if color:
        codes.append(colors[color])

    colorcode = ''.join(codes)
    # CASE 1: asterisks round colorasised part
    if not formatargs:
        #idx0 = findindex(message,marker)
        message.replace(marker,colorcode,1)
        #idx1 = findindex(message,marker)
        message.replace(marker,endcode,1)
        print(message)
    else:
        # This needs to be implemented
        pass

    return


def trycreate(loc, parents=True,exist_ok=True,loc_is_dir=False):
    """
    Todo:
        * Implement exist_ok.

    Args:

        loc: path-like object or string to directory
            or file. If loc is a directory, this method will attempt to
            create a folder if it doesn't already exist. If loc is
            a file, its parent directory will be treated as loc.
        parents (bool): if True, then equivalent to mkdir -p </path/>
        exist_ok(bool): if True, then warnings are ignored about folder already
                        existing.
    """
    l = enforce_pathobj(loc)
    wowprint("Checking **{}** exists.".format(l),color='blue')

    # First, check if file or directory exists
    if l.exists():
        if l.is_dir():
            print("Directory already exists.")
        else:
            print("File already exists.")
        return

    # If not, create the directory enclosing the file, or
    # the directory.

    # Assume that l is a path to a file, so we need l.parent for the directory
    if not loc_is_dir:
        l = l.parent

    l.mkdir(parents=parents,exist_ok=True)

    # Does the location exist?
    assert l.exists()
    wowprint("The directory **{}** has been made.".format(l),color='red')
    # pdb.set_trace()
    return

def enforce_pathobj(obj,path_type='posix'):
    objs = {'posix':PosixPath}
    if isinstance(obj,str):
        return objs[path_type](obj)
    elif isinstance(obj,PosixPath):
        return obj
    else:
        raise Exception("Object is neither path nor path-like object.")

def generate_timestamp_fname(extension):
    """ Returns a file name based on current time.
    """
    nowutc = datetime.datetime.now(tz=pytz.utc)
    fname = string_from_time('output',nowutc,strlen='second')
    return '.'.join((fname,extension))

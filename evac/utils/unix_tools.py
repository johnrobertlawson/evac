"""A collection of *nix utility scripts.

These help with directory, file, etc issues.
"""

import pdb
import os
import base64
from pathlib import Path, PosixPath

try:
    import paramiko
except ImportError:
    print("Module paramiko unavailable. Ignoring import.")

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

    if command is "copy":
        # This is ridiculous, but that's what I get for using sodding pathlib
        topath.write_bytes(frompath.read_bytes())
    elif command is "move":
        # This replaces an existing 'topath' silently
        frompath.rename(topath)
    elif command is "softlink":
        frompath.symlink_to(topath)
    else:
        raise NonsenseError("The command variable **{}** is invalid".format(
                        command),color='red')
    assert topath.exists()

    if return_newpath:
        return topath

def wowprint(message,color='red',bold=True,underline=False,formatargs=False):
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
        print(s)
        return

    def findindex(s, lookup = '**'):
            try:
                idx = s.index('**')
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
        idx0 = findindex('**')
        message.replace('**',colorcode,1)
        idx1 = findindex('**')
        message.replace('**',endcode,1)
        print(message)
    else:
        # This needs to be implemented
        pass
    
    return
    

def _bridge(cmd,frompath,topath,mv=False,cp=False,ln=False):
    """
    Olde version.
    Soft-link, copy, or move item.

    Create folder if it doesn't exist
    """
    #if sum([mv,cp,ln]) != 1:
    #    raise Exception("Choose one of mv, ln, cp commands.")
    todir, tofile = topath.split()
    trycreate(todir)
    if cmd in ('mv','move'):
        c = 'mv'
    elif cmd in ('cp','copy'):
        c = 'cp'
    elif cmd in ('ln','link','softlink'):
        c = 'ln -sf'

    cmd = "{} {} {}".format(c,frompath,topath)
    return cmd

def trycreate(loc, parents=True,exist_ok=True,loc_is_dir=False):
    """
    Args:

    loc     :   (Path,str) - path-like object or string to directory
                or file. If loc is a directory, this method will attempt to
                create a folder if it doesn't already exist. If loc is
                a file, its parent directory will be treated as loc.
    parents  :   (bool) - if True, then equivalent to mkdir -p </path/>
    exist_ok    :   (bool) - if True, then warnings are ignored about folder already
                    existing.
    """
    l = enforce_pathobj(loc)
    # if not l.is_dir():
    # if not loc_is_dir:
        # l = l.parent

    wowprint("Checking **{}** exists.".format(l),color='blue')
    # Does the location exist?
    if not l.exists():
        l.mkdir(parents=parents,exist_ok=True)
        wowprint("The directory **{}** has been made.".format(l),color='red')
    # pdb.set_trace()
    return


def softlink(frompath,topath):
    """ Arguments must be path objects or strings.

    Args:

    frompath    :   (str,Path) - absolute path to origin directory or file
    topath      :   (str,Path) - absolute path to destination directory or
                    file. If topath is only a directory, the link will
                    be the same filename as the origin.
    """

    # Use path-like objects
    frompath = enforce_pathobj(frompath)
    topath = enforce_pathobj(topath)

    # Check whether paths are files or directories
    #frompath.isfile()?

    #if not fname:
        # todir is a path to the new link
        # assert # Check this is file name and not directory
        # topath = todir
    # else:
        #topath = enforce_pathobj(topath) / fname

    frompath.symlink_to(topath)

    return

def enforce_pathobj(obj,path_type='posix'):
    objs = {'posix':PosixPath}
    if isinstance(obj,str):
        return objs[path_type](obj)
    elif isinstance(obj,PosixPath):
        return obj
    else:
        raise Exception("Object is neither path nor path-like object.")


def dir_from_naming(root,*args):
    """
    Generate file path from arguments

    Inputs:
    root    :    file path base
    args     :    list of arguments to join as separate folders
    """
    l = [str(a) for a in args]
    path = os.path.join(root,*l)
    return path

def ssh_client(ky,domain,username,password):
    key = paramiko.RSAKey(data=base64.decodestring(ky))
    client = paramiko.SSHClient()
    client.get_host_keys().add(domain, 'ssh-rsa', key)
    client.connect(domain, username=username, password=password)
    return client
    stdin, stdout, stderr = client.exec_command('ls')
    for line in stdout:
        print('... ' + line.strip('\n'))
    client.close()


def edit_namelist(f,sett,newval,absolute=False,doms=1,samealldoms=True):
    if doms < 1:
        raise ValueError
    f = enforce_pathobj(f)
    #fs = open(f,'r')
    fs = f.open('r')
    # with open(f,'r') as fopen:
    # flines = open(f,'r').readlines()
    flines = fs.readlines()
    # pdb.set_trace()
    done= False
    for idx, line in enumerate(flines):
        if sett in line:
            # print("Setting {0} is in line.".format(sett))
            fs.close()
            # spaces = 38 - len(sett)
            spaces = 36
            # print(sett,spaces)
            if not isinstance(newval,(tuple,list)):
                newval = '{0},     '.format(newval)*doms
            elif doms == 2:
                newval = '{0},     {1},'.format(newval[0],newval[1])
            elif doms == 3:
                newval = '{0},     {1},     {2},'.format(newval[0],newval[1],newval[2])

            flines[idx] = " {0: <{sp}}= {1}\n".format(sett,newval,sp=spaces)
            nameout = f.open('w',1)
            # nameout = open(f,'w',1)
            nameout.writelines(flines)
            nameout.close()
            done = True
            break
    fs.close()
    if absolute and (not done):
            # import pdb;pdb.set_trace()
        raise ValueError("Setting",sett,"not found in namelist.")
    return

"""
Utility scripts to help with directory, file, etc issues
"""

import pdb
import os
import base64
from pathlib import Path, PosixPath

import gis_tools as utils

try:
    import paramiko
except ImportError:
    print("Module paramiko unavailable. Ignoring import.")

def bridge_multi(D):
    """
    Soft-link, copy, or move multiple items.

    D   :   dictionary with the following layout:
    dict( (frompath,topath)=cmd)

    where

    cmd is 'copy','softlink', or 'move'
    frompath/topath is absolute path of origin/destination

    """
    for (frompath,topath),cmd in D.items():
        bridge(cmd,frompath,topath)
    return

def bridge(command,frompath,topath,catch_overwrite=True):
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
    """

    # Use path-like objects
    frompath = enforce_pathobj(frompath)
    topath = enforce_pathobj(topath)

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
    return

def wowprint(string,color='red',bold=True,underline=False,formatargs=False):
    """ Print with colour/bold colorasis on certain things!

    string      :   (str) - string with any colorasised
                    item surrounded by two asterisks on both
                    sides. These are replaced.
    color        :   (str,bool) - the colour of the colorasised
                    text. If False, don't colorasise.

    Usage:

    wowprint("The program **prog** is broken",color='red')
    wowprint("The program {} is broken",formatargs=('prog',),color='blue',
                            underline=True)
    """
    # If nothing is set, pass through without changes
    if not all(color,bold,underline,formatargs):
        return string

def bridge(frompath,topath,mv=False,cp=False,ln=False):
    """
    Soft-link, copy, or move item.

    Create folder if it doesn't exist
    """
    if (mv,cp,ln).sum() != 1:
        raise Exception("Choose one of mv, ln, cp commands.")
    todir, tofile = os.path.split(topath)
    utils.trycreate(todir)
    if mv:
        cmd = "mv {} {}".format(frompath,topath)

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
    fs = open(f,'r')
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
            nameout = open(f,'w',1)
            nameout.writelines(flines)
            nameout.close()
            done = True
            break
    fs.close()
    if absolute and (not done):
            # import pdb;pdb.set_trace()
        raise ValueError("Setting",sett,"not found in namelist.")
    return

    """
Utility scripts to help with directory, file, etc issues
"""

import pdb
import os
import base64
from pathlib import Path, PosixPath

import utils

try:
    import paramiko
except ImportError:
    print("Module paramiko unavailable. Ignoring import.")

def bridge_multi(D):
    """
    Soft-link, copy, or move items.

    D   :   dictionary with the following layout:
            '<'
    """
    pass

def bridge(frompath,topath):
    """
    Soft-link, copy, or move item.
    
    Create folder if it doesn't exist
    """
    tofile, todir = os.path.split(topath)
    utils.trycreate()

def softlink(frompath,todir,fname=False):
    """ Arguments must be path objects or strings
    """
    
    # Use path-like objects
    frompath = enforce_pathobj(frompath)

    # Check frompath is a file (or directory?)
    frompath.isfile()?

    if not fname:
        # todir is a path to the new link
        assert # Check this is file name and not directory
        topath = todir
    else:
        topath = enforce_pathobj(topath) / fname

    # Check topath
    frompath.symlink_to(topath)
    return

def enforce_pathobj(obj,path_type='posix'):
    objs = {'posix':PosixPath}
    if isinstance(obj,str):
        return objs[path_type](obj)
    elif issubclass(obj,Path):
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

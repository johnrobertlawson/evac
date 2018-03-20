#!/usr/bin/env pytho    n
#
# Create experiment directory and ensemble member directories, and copy/link executables and auxiliary files

import sys, os, glob
import datetime
from optparse import OptionParser
import numpy as N
from multiprocessing import Pool
import netCDF4 as ncdf

# BEGIN user-defined parameters

WRFDIR     = '/scratch/software/Odin/git/WRF/Intel/WRFV3.6.1_JC/main/'  # WRF executable dir
DATADIR    = '/scratch/software/Odin/git/WRF/Intel/WRFV3.6.1_JC/run/'

HOMEDIR    = '/lus/scratch/mflora/'  # dir containing your namelist.input
ICBCDIR    = '/scratch/mflora/Initial_ens_files_20160524_I0000' # Initial/Boundary Conditions dir (don't change)
MICRODIR   = '/scratch/mflora/PBL_WORK/ELRENO/100m_ELRENO/'    # Thompson microphysics tables dir (don't change)
wrfname    = 'wrfout_d01_2016-05-25_00:00:00'                  # wrfout filenames in ICBCDIR

debug      = True

# END user-defined parameters

myhelp = """\n
            Setup script for WRF-LETKF
            Usage: srun python create_exp_monte.py -d /scratch/mflora/20160524_I0000_reduced_50 -n 36
         """
#--------------------------------------------------------------------------------------
# FILE DICTIONARY
#            name                location                  target                 copy
#                                                                          (0=no, 1=link, 2=cp)
#--------------------------------------------------------------------------------------
FILE_DICT= {'model':            [WRFDIR,       'wrf.exe' ,                  1 ],
            'init_namelist':    [HOMEDIR,      'namelist.input' ,           2 ],
            'qr_acr_qg.bin':    [MICRODIR,     'qr_acr_qg.bin' ,            2 ],
            'freezeH2O.bin':    [MICRODIR,     'freezeH2O.bin' ,            2 ],
            'qr_acr_qs.bin':    [MICRODIR,     'qr_acr_qs.bin' ,            2 ],
            'init_state':       [ICBCDIR,      'wrfinput_d01',              2 ],
            'LANDUSE.TBL':      [DATADIR,       'LANDUSE.TBL',               1 ],
            'GENPARM.TBL':      [DATADIR,       'GENPARM.TBL',               1 ],
            'SOILPARM.TBL':     [DATADIR,       'SOILPARM.TBL',              1 ],
            'ETAMPNEW_DATA':    [DATADIR,       'ETAMPNEW_DATA',             1 ],
            'gribmap.txt':      [DATADIR,       'gribmap.txt',               1 ],
            'RRTM_DATA':        [DATADIR,       'RRTM_DATA',                 1 ],
            'RRTMG_LW_DATA':    [DATADIR,       'RRTMG_LW_DATA',             1 ],
            'RRTMG_SW_DATA':    [DATADIR,       'RRTMG_SW_DATA',             1 ],
            'tr49t67':          [DATADIR,       'tr49t67',                   1 ],
            'tr49t85':          [DATADIR,       'tr49t85',                   1 ],
            'tr67t85':          [DATADIR,       'tr67t85',                   1 ],
            'VEGPARM.TBL':      [DATADIR,       'VEGPARM.TBL',               1 ]
           }

#=======================================================================
#
#
# Function to do the link in the need scripts, info, etc.

def linkit(current_dir,experiment_dir,program_dir,target,link_option):

    makefile_dir = program_dir
    if debug: print "Change directory to " + makefile_dir

    if link_option == 0:  return

    if link_option == 1:

        LINK_CMD = 'ln -sf ' + os.path.join(makefile_dir,target) + ' ' + os.path.join(experiment_dir,target)

        if os.system(LINK_CMD) != 0:
            print
            print "ERROR!!! "
            print "ERROR!!!  Failed to EXECUTE: " + LINK_CMD
            print "ERROR!!! "
            print
            os.chdir(current_dir)
            sys.exit(1)

        if debug: print "Linked " + target + " to directory  " + experiment_dir

    if link_option == 2:

        if (target=='input.nml' or target=='qr_acr_qg.bin' or target=='freezeH2O.bin' or target=='qr_acr_qs.bin'):
            CP_CMD = 'cp ' + os.path.join(makefile_dir,target) + ' ' + os.path.join(experiment_dir,target)
        elif ('namelist.input' in target):
            CP_CMD = 'cp ' + os.path.join(makefile_dir,target) + ' ' + os.path.join(experiment_dir,"namelist.input")
        elif ('wrfbdy' in target):
	    CP_CMD = 'cp ' + os.path.join(makefile_dir,target) + ' ' + os.path.join(experiment_dir,"wrfbdy_d01")
	else:
            CP_CMD = 'cp ' + os.path.join(makefile_dir,target) + ' ' + os.path.join(experiment_dir,"wrfinput_d01")
        if os.system(CP_CMD) != 0:
            print
            print "ERROR!!! "
            print "ERROR!!!  Failed to EXECUTE: " + CP_CMD
            if not os.path.exists(os.path.join(makefile_dir,target)):
                print "\nCOMMAND FAILED because %s does not exist\n " % os.path.join(makefile_dir,target)
            elif not os.path.join(experiment_dir,target):
                print "\nCOMMAND FAILED because %s does not exist\n " % os.path.join(experiment_dir,target)
            else:
                print "\nBOTH FILE AND DIRECTORY EXISTS....something else f__ked up here....\n"
            print "ERROR!!! "
            print
            os.chdir(current_dir)
            sys.exit(1)

        if debug: print "Copied " + target + " to directory  " + experiment_dir

    return

def fnormal(prng=N.random, scale=1.0, size=(1,)):
  return prng.normal(scale=scale, size=size)
#
#=======================================================================
# main script

print
print "<<<<<===========================================================================================>>>>>>"
print
print

usage  = "usage: %prog [options] arg \n" + myhelp
parser = OptionParser(usage)

parser.add_option("-d",  "--dir",    dest="run_dir", default = None, type="string", help="Name of experiment directory")
parser.add_option("-n",  "--ne",     dest="ne",      default = 10,   type="string", help="Number of ensemble members")
parser.add_option("-f",  "--force",  dest="force",   default = True,                help="Force overwriting of existing experiment directory", action="store_false")

(options, args) = parser.parse_args()

cwd = os.getcwd()

#/////////////////////////////////////////////////////////////////
# Section to create (and move over current) experiment directories
    
if options.run_dir == None:
    link_directory = cwd + "/experiment"
else:
    link_directory =  options.run_dir

if options.force:

# If the requested experiment directory exists, we have to do something with the current one

    if os.path.exists(link_directory):

# Create timestamp string to append to directory name (to move it over)

        timestamp = datetime.datetime.fromtimestamp( os.path.getctime(link_directory) )
        newdir   = link_directory + "_" + timestamp.isoformat().replace('T', '_')
        print "\nERROR:  EXPERIMENT DIRECTORY ALREADY EXISTS, MOVING IT TO: %s \n" % (newdir)

        os.rename(link_directory, newdir)

# Okay, we moved stuff around, now create new working directory

os.mkdir(link_directory)
  
#//////////////////////////////////////////////////
# Section for creating run directories
# and for linking scripts and codes
#//////////////////////////////////////////////////

targets = FILE_DICT.keys()
member_dir = []

if options.ne == 10:
    ne = 10
else:
    ne = int(options.ne)

print "Creating directory structure and linking/copying WRF files"

# Ensemble members

for n in N.arange(1,ne+1):
    member = "%s/member%3.3i" % (link_directory, n)
    os.mkdir(member)

    for key in targets:
        target = FILE_DICT[key]
        if (key!='init_state' and key!='BCs' and key!= 'init_namelist'):
            ret = linkit(cwd,member,os.path.join(cwd,target[0]),target[1],target[2])
        elif (key=='init_state'):
            #targ=target[1]+"_%d" %(n)
            targ="%s/%s_%d" %(ICBCDIR, wrfname, n)
            ret = linkit(cwd,member,os.path.join(cwd,target[0]),targ,target[2])
        elif (key=='BCs'):
            targ=target[1]+"_%d" %(n)
            #targ="member%03d/wrfbdy_d01" %(n)
            ret = linkit(cwd,member,os.path.join(cwd,target[0]),targ,target[2])
        elif (key=='init_namelist'):
            targ=target[1]
            ret = linkit(cwd,member,os.path.join(cwd,target[0]),targ,target[2])
	else:
	    targ=target[1]+"*_%d" %(n)
            ret = linkit(cwd,member,os.path.join(cwd,target[0]),targ,target[2])
    member_dir.append(member)


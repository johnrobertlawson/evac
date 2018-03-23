""" This script generates namelist.input files for
the specified experiment. 

Change "settings and constants" to customise (e.g., number
of ensemble members). Script will return N namelist.input.n files,
where N is the number of ensembles, and n is from
1 to N inclusive. This file is placed in the folder where
the script runs from.

Usage from command line:

$ python3 generate_namelists.py [--sp/--mp] [--stoch/--nostoch] [--onlydo x]
                                [--overwrite] [--donotenforce] [--final]

Options --sp or --mp to select single or multi physics (pick one!)
Options --stoch or --nostoch to turn on/off stochasticity (pick one!)
Option --onlydo x, where 0 < x < N, only does one member (optional)
Option --overwrite will overwrite any existing files in the script directory
Option --donotenforce will ignore any missing parameters in the namelist.input
Option --final will use the namelist for the finished version of WRF.exe

For example:

$ python3 generate_namelists.py --sp --stoch --onlydo 4 --overwrite

This will generate a SP-Stoch namelist for ensemble member 4, ignoring existing namelists.

Written as a present for Nusrat Yussouf
John Lawson Feb 23 2017 CIMMS/NSSL
"""
######## IMPORTS ########
import os
import argparse
import random
import glob
import numpy as N
import itertools
import pdb
import time
import sys

# Dirty fix so Nusrat can see lhsmdu package
try:
    from lhsmdu import sample
except ImportError:
    raise Exception("lhsmdu.py needs to be in this directory.")

######## COMMAND LINE LOGIC ########
parser = argparse.ArgumentParser()
parser.add_argument('--sp', action='store_true')
parser.add_argument('--mp', action='store_true')
parser.add_argument('--stoch', action='store_true')
parser.add_argument('--nostoch', action='store_true')
parser.add_argument('--onlydo', type=int, default=False)
parser.add_argument('--overwrite',action='store_true')
parser.add_argument('--donotenforce',action='store_false')
parser.add_argument('--final',action='store_true')
NS = parser.parse_args()

if (not NS.mp) and (not NS.sp):
    raise Exception("Specify --sp or --mp in commmand.")
elif (not NS.stoch) and (not NS.nostoch):
    raise Exception("Specify --stoch or --nostoch in command.")
elif NS.mp and NS.sp:
    raise Exception("Choose either --sp or --mp, not both.")
elif NS.stoch and NS.nostoch:
    raise Exception("Choose either --stoch or --nostoch, not both.")
elif NS.sp:
    if NS.nostoch:
        experiment = "SP"
    elif NS.stoch:
        experiment = "SP-Stoch"
    else:
        raise Exception("Weird value for stochasticity?")
elif NS.mp:
    if NS.nostoch:
        experiment = "MP"
    elif NS.stoch:
        experiment = "MP-Stoch"
    else:
        raise Exception("Weird value for stochasticity?")
else:
    raise Exception("Illegal values somewhere: \n {0}.".format(NS))

# option to only generate one namelist
onlydo = NS.onlydo 
# Kills script if the template namelist is missing key parameters
enforce_all = NS.donotenforce 
if NS.final:
    path_to_namelist_template = '/home/john.lawson/ens_design/pycode/namelist.input.final'
else:
    path_to_namelist_template = '/home/john.lawson/ens_design/pycode/namelist.input.test.og'

######## CONSTANTS AND SETTINGS ########
#### USER-CONFIGURABLE HERE! ! ! ! #####
history_outname = False # Change this to absolute path of output file and location
debug = False # Print more output
nens = 18 # Number of ensemble members

######## LOGIC FOR SETTINGS ########
ensnums = range(1,nens+1) # names of each ensemble
nidxs = range(nens) # indices of each ensemble
doms = 2 # number of domains
# cwd = os.cwd() # Where to output namelist files
cwd = os.path.dirname(os.path.abspath(__file__)) # Where to output namelist files
namelists = ['namelist.input.{0}'.format(n) for n in ensnums]
paths_to_nl = {}
for nidx,n in zip(nidxs,ensnums):
    paths_to_nl[n] = {'old': path_to_namelist_template, 'new':os.path.join(cwd,namelists[nidx])}

######## FIRST CHECK TO SEE IF NAMELISTS MIGHT BE OVERWRITTEN ########
# nlfs = glob.glob(os.path.join(cwd,'namelist.*'))
# for k,v in paths_to_nl[n].items():
for nidx,n in zip(nidxs,ensnums):
    if onlydo and (n != onlydo):
        print("Skipping this one.")
        continue
    if os.path.isfile(paths_to_nl[n]['new']):
        if not NS.overwrite:
            raise Exception("Add --overwrite to command line to ignore existing namelists.")
        else:
            rmcmd = 'rm {0}'.format(paths_to_nl[n]['new'])
            os.system(rmcmd)


######## FUNCTIONS ########
def denormalise(normvals,xmin,xmax):
    x = xmin + normvals*(xmax-xmin)
    return x

def edit_namelist(f,sett,newval,enforce_all=False,doms=1,precision=3):
    """ A single integer or float will be repeated to all domains.
    A list/tuple of values, length dom, will be assigned the each domain.
    A list/tuple of one value (len 1) will be assigned only do the first domain
        (i.e. the rest are blank).
    """
    if doms < 1:
        raise ValueError

    # Create list of values for writing to line
    if not isinstance(newval,(tuple,list)):
        newvals = [newval,]*doms
    else:
        if len(newval) == 1:
            # Only one value despite more than one domain
            newvals = [newval[0],]
        else:
            newvals = list(newval)
    
    # Set format of number
    if isinstance(newvals[0],(float,N.float)):
        fm = "{:.3f}"
    elif isinstance(newvals[0],(int,N.int)):
        fm = "{}"
    else:
       raise Exception("Format of number is",type(newvals[0]))

    # Get length of string so we can make the namelist pretty
    pad = 6 - len(eval(''.join(("'",fm,"'",'.format(newvals[0])'))))
    if pad<1:
        pad = 1

    fs = open(f,'r')
    # with open(f,'r') as fopen:
    # flines = open(f,'r').readlines()
    flines = fs.readlines()
    # pdb.set_trace()
    done= False
    for idx, line in enumerate(flines):
        if sett in line:
            # Make sure the match occurs at start of line
            # OR begins with a space
            # (No partial matches!)
            if ' '+sett in line:
                pass
            elif (line.find(sett)) == 0:
                pass
            else:
                continue

            # Make sure the match ends in a space or equals
            if sett+' ' in line:
                pass
            elif sett+'=' in line:
                pass
            else:
                continue

            # print("Setting {0} is in line.".format(sett))
            fs.close()
            # spaces = 38 - len(sett)
            spaces = 36
            # print(sett,spaces)
            # if not isinstance(newval,(tuple,list)):
            newval = ''.join([''.join((fm,',',' '*pad)).format(v) for v in newvals])
            # newval = ''.join(fm,',',' '*5).format(newval)*doms
            # elif doms == 2:
                # newval = '{0},     {1},'.format(newval[0],newval[1])
            # elif doms == 3:
                # newval = '{0},     {1},     {2},'.format(newval[0],newval[1],newval[2])

            flines[idx] = " {0: <{sp}}= {1}\n".format(sett,newval,sp=spaces)
            nameout = open(f,'w',1)
            nameout.writelines(flines)
            nameout.close()
            done = True
            break
    fs.close()
    if not done:
        if enforce_all:
            raise ValueError("Setting",sett,"not found in namelist.")
        else:
            print("Warning: setting",sett,"not found in namelist.")
    return

def generate_seeds(num):
    seeds = random.sample(range(1,1000),num)
    return seeds

def generate_LHS_draw(seed):
    # Dictionary of stochastic settings
    ST = {}

    # hard code the number of stochastic settings here
    # Number depends on stochastic cumulus scheme selected
    ns = 20
    nss = iter(range(ns))
    LHS = sample(ns,nens,randomSeed=seed)

    #### MICROPHYSICS #### (5)
    # Some of these distrubtions could be beta like Hacker et al 2011
    ST['nssl_ehw0'] = denormalise(LHS[next(nss),:],0.4,1.0)
    ST['nssl_alphar'] = denormalise(LHS[next(nss),:],0.0,2.5)
    ST['nssl_alphah'] = denormalise(LHS[next(nss),:],0.0,3.0)
    ST['nssl_cccn'] = denormalise(LHS[next(nss),:],0.3e9,1.3e9)

    # Hail collect eff. should be higher than graupel
    # ehlw0 must be computed last because of constraint
    ehlw0 = N.zeros([nens])
    ehln = next(nss)
    for xn, x in enumerate(ST['nssl_ehw0'][0,:].A1):
        ehlw0[xn] = denormalise(LHS[ehln,xn],ST['nssl_ehw0'][0,xn],1.0)
    ST['nssl_ehlw0'] = N.matrix(ehlw0)

    #### SHORTWAVE #### (0)
    # No solar scattering any more
    # swscat = denormalise(LHS[4,:],0.2,2.0)

    ##### CUMULUS ##### ( )

    ##### LAND SURFACE ##### (8)
    # Radii - must be square root for uniform distribution round circle
    for key in ('morphr_crs','morphr_btr','morphr_rad','morphr_tbot'):
        ST[key] = N.sqrt(denormalise(LHS[next(nss),:],0.0,1.0))
    # Angles - TODO: Must be square root!
    for key in ('morphth_crs','morphth_btr','morphth_rad','morphth_tbot'):
        ST[key] = denormalise(LHS[next(nss),:],0.0,2*N.pi)

    return ST


######## CREATE DICTIONARY OF NAMELIST SETTINGS ########
NL = {}
NL['SP'] = {}
NL['SP-Stoch'] = {}
NL['MP'] = {}
NL['MP-Stoch'] = {}
experiments = NL.keys()

######## SETTINGS CONSTANT FOR ALL ENSEMBLES ########
ALLDICT = {}

# Set output location/name of wrfout
if history_outname:
    ALLDICT['history_outname'] = "'{0}'".format(history_outname)

ALLDICT['sf_surface_physics'] = 4 # Noah-MP
# sf_sfclay_physics = ?
#sf_urban_physics = 0 
ALLDICT['mp_physics'] = 17 # NSSL 2-mom with CCNs from namelist


######## SETTINGS THAT DEPEND ON THE EXPERIMENT ########
### NOTE: a tuple or list value for a namelist parameter
#         indicates the setting is different for each domain

# These are variable
bl_pbl_physics = None 
# YSU (1)
# MYJ (2) <----
# MYNN 2.5 (5) <----
# MYNN 3.0 (6)
# ACM2 (7)
# Shin-Hong (11)

cu_physics = None
# KF (1)
# Grell-D (93)
# Grell-3D (5)
# Tiedtke (6)
# Grell-F (3)
# KF CuP (10)
# New Tiedtke (16)

ra_lw_physics = None
# RRTM (1)
# RRTMG (4)

ra_sw_physics = None
# Dudhia (1)
# RRTMG (4)

######## SP and SP-Stoch SETTINGS ########
# SP and SP-Stoch are pretty simple!
SPDICT = {}
SPDICT['bl_pbl_physics'] = 5
SPDICT['cu_physics'] = (3,0)
SPDICT['ra_lw_physics'] = 4
SPDICT['ra_sw_physics'] = 4

for ex in ("SP","SP-Stoch"):
    for n in ensnums: 
        NL[ex][n] = {k:v for k,v in SPDICT.items()}

######## MP and MP-Stoch SETTINGS ########
for ex in ("MP","MP-Stoch"):
    for n in ensnums:
        NL[ex][n] = {}
        if n%2: # odd
            NL[ex][n]['ra_sw_physics'] = 1
            NL[ex][n]['ra_lw_physics'] = 1
        else: # even
            NL[ex][n]['ra_sw_physics'] = 4
            NL[ex][n]['ra_lw_physics'] = 4

        if (n%3) ==  0:
            NL[ex][n]['bl_pbl_physics'] = 5
        elif (n%3) == 1:
            NL[ex][n]['bl_pbl_physics'] = 11 # ShinHong
        else:
            NL[ex][n]['bl_pbl_physics'] = 2

        if n in range(1,7):
            NL[ex][n]['cu_physics'] = (1,0)
        elif n in range(7,13):
            NL[ex][n]['cu_physics'] = (5,0)
        else:
            assert n < 19 # Sanity check
            NL[ex][n]['cu_physics'] = (6,0)

######## WRITE SETTINGS TO MAIN DICTIONARY ########

# Add the constant settings for all!
for ex in experiments:
    for n in ensnums:
        for k,v in ALLDICT.items():
            NL[ex][n][k] = v


######## COMMON STOCHASTIC SETTINGS ########
for ex in experiments:
    for n in ensnums:
        ### SKEB SETTINGS ###
        NL[ex][n]['tot_backscat_psi'] = 0.01
        NL[ex][n]['tot_backscat_t'] = 0.001
        NL[ex][n]['ztau_psi'] = (86400,)
        NL[ex][n]['ztau_t'] = (86400,)
        NL[ex][n]['zsigma2_eps'] = (0.8,)
        NL[ex][n]['zsigma2_eta'] = (0.8,)

        ### Kober/Craig modified SPPT ###
        NL[ex][n]['aml'] = 150.0 # This is overridden in MYJ scheme?
        NL[ex][n]['ashc'] = 15.0 # trying to get field on par with -10 to 10 like Kober
        #NL[ex][n]['iseed'] = 0 # random
        NL[ex][n]['lengthscale_sppt'] = 15000 # 3dx like Kober
        NL[ex][n]['timescale_sppt'] = 600 # 10 min like Kober
        NL[ex][n]['stddev_cutoff_sppt'] = 3  # trial and error - this was 5 in some tests

# this adds a number to run # for multiple STCH/SKEB runs later
#modifier = 0

### Incoming solar ###
# 1.0 is 10%, 0.5 is 5%, etc
#swrad_scat = (0.2,2)

# Soil morphing turned on

######## VARIABLE STOCHASTIC SETTINGS ########
for ex in experiments:
    if ex.endswith('Stoch'):
        # Generate the Latin Hypercube samples
        seeds_sppt = generate_seeds(nens)
        seeds_skeb = generate_seeds(nens)
        seed_lhs = generate_seeds(1)
        # ehw0,ehlw0,alphar,alphah,cccn = generate_LHS_draw(seed_lhs)
        # Dictionary of stochastic values
        STOCHS = generate_LHS_draw(seed_lhs,)
        print("Inserting stochasticity for",ex)
    else:
        print("No stochasticity for",ex)
    for nidx,n in zip(nidxs,ensnums):
        if ex in ('SP-Stoch','MP-Stoch'):

            # Turn on Noah MP morphing
            NL[ex][n]['opt_btr'] = 11
            NL[ex][n]['opt_rad'] = 11
            NL[ex][n]['opt_tbot'] = 11
            NL[ex][n]['opt_crs'] = 11

            # Turn on SKEB
            NL[ex][n]['skebs'] = 1

            # Turn on PSP
            NL[ex][n]['sppt'] = 1

            # Set randomness seeds
            NL[ex][n]['nens'] = (seeds_skeb[nidx],)
            NL[ex][n]['ISEED_SPPT'] = (seeds_sppt[nidx],)

            # Draw from the LH samples
            # Assign ensemble members to draw different SKEB spectrum
            for k,v in STOCHS.items():
                # NL[ex][n][k] = v
                NL[ex][n][k] = STOCHS[k][0,nidx]
        elif ex in ('SP','MP'):
            # Turn off SKEB
            NL[ex][n]['skebs'] = 0
            # Turn off PSP
            NL[ex][n]['sppt'] = 0

            # Fix Noah defaults
            NL[ex][n]['opt_btr'] = 1
            NL[ex][n]['opt_rad'] = 3
            NL[ex][n]['opt_tbot'] = 2
            NL[ex][n]['opt_crs'] = 1

            # NSSL 2-moment defaults
            NL[ex][n]['nssl_cccn'] = 0.5e9
            NL[ex][n]['nssl_alphah'] = 0.0
            NL[ex][n]['nssl_alphar'] = 0.0
            NL[ex][n]['nssl_ehw0'] = 0.5
            NL[ex][n]['nssl_ehlw0'] = 0.75
        else:
            raise Exception

######## GENERATE NAMELISTS ########
for nidx,n in zip(nidxs,ensnums):
    print("Ensemble run ",n)
    if onlydo and (n != onlydo):
        print("Skipping this one.")
        continue

    # COPY NAMELIST FROM TEMPLATE
    copy_cmd = 'cp {0} {1}'.format(paths_to_nl[n]['old'],paths_to_nl[n]['new'])
    os.system(copy_cmd)

    # EDIT NAMELISTS
    for k,v in NL[experiment][n].items():
        # The 'enforce_all' keyword enforces all namelist changes
        # try:
        edit_namelist(paths_to_nl[n]['new'],k,v,enforce_all=enforce_all,doms=2)
        # except ValueError:
            # print("Check namelist template - it is missing the setting",k)
                # raise Exception
        # else:
            # if debug:
        if enforce_all:
            print("Changed setting",k,"to",v)
    print("Completed editing {0}".format(paths_to_nl[n]['new']))

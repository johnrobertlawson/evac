""" This should copy the namelist, wrf_input, wrf_bdy files
needed for a given ensemble member. Then move the wrfout file
into the correct folder.
"""

import os
import glob
import numpy as N
import itertools
import pdb
import time
import datetime

from netCDF4 import Dataset
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from WEM.lazyWRF.lazyWRF import main
from WEM.postWRF.postWRF.main import WRFEnviron
from WEM.postWRF.postWRF.wrfout import WRFOut
import WEM.utils as utils


pathtoWRFout = '/scratch/john.lawson/stoch_forecasts'
pathtoWRF = '/scratch/john.lawson/WRF/stochWRF_3.8.1_SINGLEDOM/WRFV3/run'
fractaldir = '/home/john.lawson/lovejoy/data'
wrfout_fname = 'wrfout.nc'
submit_fname = 'run_wrf.job'
sched = 'slurm'

# Cases - watch out if time is after 2400 UTC - day!
dt = datetime.datetime
## Format is CASES[casedate] = [list of inittimes]
CASES = {dt(2015,4,26,0,0,0):[dt(2015,4,27,3,0,0),],
            dt(2015,5,6,0,0,0):[dt(2015,5,7,0,0,0),],
            dt(2015,5,23,0,0,0):[dt(2015,5,23,21,0,0),],
            dt(2015,5,25,0,0,0):[dt(2015,5,26,2,0,0),],
            dt(2016,7,30,0,0,0):[dt(2016,7,31,0,0,0),]}

enstypes = ('SP','MP','SP-Stoch','MP-Stoch')

############# FUNCTIONS ##############
def sleep(message,sec=10):
    print(message)
    time.sleep(sec)
    return

def normalise(arr,minmax=0.1):
    arr2 = (arr-N.min(arr))/(N.max(arr) - N.min(arr))
    arr3 = (arr2-0.5)*(2*minmax)
    return arr3


def interpolate_to_wrfgrid(arr,new_shp):
    from scipy.interpolate import RectBivariateSpline as RBS

    old_shp = arr.shape
    xx = N.arange(old_shp[1])
    yy = N.arange(old_shp[0])
    rbs = RBS(xx,yy,arr)

    xx_new = N.arange(new_shp[1])
    yy_new = N.arange(new_shp[0])
    newarr = rbs(xx,yy)
    # pdb.set_trace()
    return newarr

def fractalise(fractals,n,fpath):
    """Pick nth member of a randomised list of perturbation
    fields, and apply it to soil moistire in wrfinput.
    """
    # Use different fractals for MP
    n = nnn-n
    # Load numpy array
    for f in fractals:
        if f.endswith('{:02d}.npy'.format(n)):
            frac_arr = N.load(f)
            f_out = f
    
    # Normalise to +/- 1%
    frac_arr_pc = normalise(frac_arr)

    # Load virgin data
    nc = Dataset(fpath,'r+')
    virgindata = nc.variables['SMOIS'][:,:,:,:]
    newdata = N.zeros_like(virgindata)
    shp = virgindata[0,0,:,:].shape

    # Interpolate fractal to wrfinput grid
    frac_arr_interp = interpolate_to_wrfgrid(frac_arr_pc,shp)
            
    # Perturb each soil level (broadcasting is scary)
    for lv in N.arange(virgindata.shape[1]):
        newdata[0,lv,:,:] = virgindata[0,lv,:,:] + (virgindata[0,lv,:,:]*frac_arr_interp)

    # Write back
    assert N.all(newdata > 0.0)
    # pdb.set_trace()
    nc.variables['SMOIS'][:,:,:,:] = newdata
    nc.close()
    return f_out

#### PROCEDURE ####

for case,it in CASES.items:
    for enstype in enstypes:
        # Initial settings
        cstr = "{:04d}{:02d}{:02d}".format(case.year,case.month,case.day)
        istr = "{:02d}{:02d}".format(it.hour,it.minute)
        pathtoWRFout_ci = os.path.join(pathtoWRFout,cstr,istr)
        utils.trycreate(pathtoWRFout_ci)
        namelistdir = os.path.join(pathtoWRF,'nls_{}'.format(enstype))
        dirstr = "eachmember_{:04d}{:02d}{:02d}{:02d}{:02d}".format(
                    it.year,it.month,it.day,it.hour,it.minute)
        # ICs file that needs renaming. Specific to member and time.
        wrfstr = "wrfout_".format()
        pathtoICs = "/work/nyussouf/WoF/{cstr}/NSSL-2M/EXP/{dirstr}/{wrfstr}".format(
                            cstr=cstr,dirstr=dirstr,wrfstr=wrfstr)
        pathtoLBCs = 
        

        # Copy wrfout ICs and rename to wrfinput, make sure python namelist generator knows

        # Run namelist generator - make sure it modifies the right wrfinputs
        # Needs to edit time/date/grid lat/lon?
        # Do 6 hour runs

        # This instance is for each ensemble run
        L = main.Lazy(sched,pathtoWRF,pathtoWRFout1,submit_fname,)

        exs = glob.glob(os.path.join(namelistdir,'namelist.inp*'))

        # Get rid of rogue rsl files to not confuse script
        cmd = "rm -f {0}".format(os.path.join(pathtoWRF,'rsl*'))
        os.system(cmd)

        # Begin ensemble
        for nidx,n in enumerate(range(1,19)):
            # Copy in fractalised ICs

            # Copy in LBCs

            # Copy in namelist

            # Get directory name for this run
            member = 'run_{:02d}'.format(n)
            tofolder = os.path.join(pathtoWRFout,member)

            if os.path.isfile(os.path.join(tofolder,'wrfout.nc')):
                print("Already done {}, moving on.".format(member))
                continue

            # Copy wrfinput
            oldinput_fpath = os.path.join(pathtoWRF,'wrfin_normal','wrfinput_d01')
            newinput_fpath = os.path.join(pathtoWRF,'wrfinput_d01')
            cmd = 'cp {0} {1}'.format(oldinput_fpath,newinput_fpath)
            os.system(cmd)

            # Generate fractal perturbations in soil moisture
            f_out = fractalise(fractals,n,newinput_fpath)

            # Copy namelist
            nl = os.path.join(namelistdir,'namelist.input.{}'.format(n))
            cmd = 'cp {0} {1}'.format(nl,os.path.join(pathtoWRF,'namelist.input'))
            os.system(cmd)

            print("Running {} with {} and {} fractal; to move to {}.".format(
                    member,os.path.basename(nl),os.path.basename(f_out),tofolder))
            # Submit job
            L.submit_job(waitmin=0.25,sleepmin=0.25,wrf_only=True)
            
            # Move wrfout file to new folder named after ensemble member
            sleep("Sleeping for 10 seconds to make sure everything has been logged to disk.")
            utils.trycreate(tofolder)
            L.copy_files(pathtoWRFout1,tofolder,{'wrfout.nc':'mv'})
            
            # Copy namelist and rsl.out.0000, rsl.error.0000
            L.copy_files(pathtoWRF,tofolder,{'namelist.input':'cp',
                                    'rsl.out.0000':'cp',
                                    'rsl.error.0000':'cp',
                                    'wrfinput_d01':'mv'})

            cmd = 'cp {} {}'.format(f_out,tofolder)
            os.system(cmd)

            sleep("Sleeping for 10 seconds to make sure everything has been copied.")
            cmd = "rm -f {0}".format(os.path.join(pathtoWRF,'rsl*'))
            os.system(cmd)


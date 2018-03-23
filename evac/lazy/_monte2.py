#!/usr/bin/env python
#   
# Integrate ensemble forward (post-DA)

import os
import time as cpu
import sys, glob
import string
from optparse import OptionParser
import numpy as np
from multiprocessing import Pool
import datetime
from subprocess import Popen

# path: full path to directory containing all exp directories
# prefix: exp directory (contains member001...member036)
# start/stop: forecast begin/end times
# ne: ensemble size

usage  = "usage: python ens_fcst_wrf_monte.py --path /scratch/mflora --prefix 20160524 --start 2016-05-25_00:00:00 --stop 2016-05-25_03:00:00 --ne 36"

batch_size = 36 # number of ensemble members to run at once
debug = False # debugging mode 

c0 = cpu.time()

def RunMember(mempath):
    """
    Updates namelist file with new start/end times then
    begins WRF integration.
    """
    os.chdir(mempath)

    cmd = "rm run_wrf*done"
    os.system(cmd)

    if options.start != None:

	### Commented out the following code and set dmin = 0 
	### DATE: 2017 Nov 28 
	dmin = 0 

        #if (end_hour == start_hour):
        #    dmin = end_minute - start_minute
        #else:
        #    dmin = 60 + end_minute - start_minute

        fin = open("namelist.input", "r")
        fout = open("temp", "wt")
	### Added "3*" to following lines 
	### DATE: 2017 Nov 28 
        for line in fin:
            if "run_minutes" in line:
                fout.write(" run_minutes                         = %d,\n" % (dmin))
            elif "start_year" in line:
                fout.write(" start_year                          = 3*%d,\n" % (start_year))
            elif "start_month" in line:
                fout.write(" start_month                         = 3*%d,\n" % (start_month))
            elif "start_day" in line:
                fout.write(" start_day                           = 3*%d,\n" % (start_day))
            elif "start_hour" in line:
                fout.write(" start_hour                          = 3*%d,\n" % (start_hour))
            elif "start_minute" in line:
                fout.write(" start_minute                        = 3*%d,\n" % (start_minute))
            elif "start_second" in line:
                fout.write(" start_second                        = 3*%d,\n" % (start_second))
            elif "end_year" in line:
                fout.write(" end_year                            = 3*%d,\n" % (end_year))
            elif "end_month" in line:
                fout.write(" end_month                           = 3*%d,\n" % (end_month))
            elif "end_day" in line:
                fout.write(" end_day                             = 3*%d,\n" % (end_day))
            elif "end_hour" in line:
                fout.write(" end_hour                            = 3*%d,\n" % (end_hour))
            elif "end_minute" in line:
                fout.write(" end_minute                          = 3*%d,\n" % (end_minute))
            elif "end_second" in line:
                fout.write(" end_second                          = 3*%d,\n" % (end_second))

            else:
                fout.write(line)
        fin.close()
        fout.close()

        os.system("mv temp namelist.input")

    if 1==1:

      exp_path = "%s/%s/" % (path, prefix)
      ### Removed -n 24 
      ### 2017 Nov 28
      model_exe = "time srun %s/wrf.exe" % (mempath)
      cmd = "%s >& %s.out" % (model_exe, mempath)
      mem = mempath[-3:]
      print cmd
      f = open('run_wrf%s.csh' %(mem), 'w')
      script_text = '#!/bin/csh\n\
\n\
#SBATCH -J run_wrf_%s\n\
#SBATCH -o %s/run_wrf_%s.out\n\
#SBATCH -e %s/run_wrf_%s.err\n\
##SBATCH --nodes=1\n\
#SBATCH -n 24\n\
#SBATCH --ntasks-per-node=24\n\
##SBATCH --cpus-per-task=1\n\
##SBATCH -A smallqueue\n\
##SBATCH -p workq\n\
#SBATCH -t 06:00:00\n\n\
\n\
#source ~/.tcshrc\n\
#set echo\n\
\n\
# cd to directory where job was submitted from\n\
cd $SLURM_SUBMIT_DIR\n\
\n\
setenv MPICH_VERSION_DISPLAY 1\n\
setenv MPICH_ENV_DISPLAY 1\n\
setenv MPICH_MPIIO_HINTS_DISPLAY 1\n\
setenv MALLOC_MMAP_MAX 0\n\
#setenv MALLOC_TRIM_THRESHOLD 536870912\n\
\n\
setenv MPICH_GNI_RDMA_THRESHOLD 2048\n\
setenv MPICH_GNI_DYNAMIC_CONN disabled\n\
\n\
setenv MPICH_MPIIO_HINTS "wrfout*:striping_factor=4,cb_nodes=4"\n\
setenv MPICH_CPUMASK_DISPLAY 1\n\
\n\
#setenv OMP_NUM_THREADS 2\n\
\n\
%s\n\
\n\
touch run_wrf%s_done\n' %(mem, mempath, mem, mempath, mem, cmd, mem)

      f.write(script_text)
      f.close()

      cmd = "chmod +x run_wrf%s.csh" %(mem)
      os.system(cmd)
      cmd = "sbatch run_wrf%s.csh" %(mem)
      if debug:
        print cmd
      else:
        os.system(cmd)

#-------------------------------------------------------------------------------
# Command line arguments
#-------------------------------------------------------------------------------

parser = OptionParser(usage)

parser.add_option("--prefix",     dest="prefix",    type="string",    help = "Experiment directory")
parser.add_option("--start",      dest="start",     type="string",    help = "Start time of model integration")
parser.add_option("--stop",       dest="stop",      type="string",    help = "Stop time of model integration")
parser.add_option("--ne ",        dest="ne",        type="int",       help = "Number of ensemble members")
parser.add_option("-p", "--path", dest="path",      type="string",    help = "Path for model runfiles")

(options, args) = parser.parse_args()

if options.path == None:
  path = "./"
else:
  path = options.path

if options.prefix != None:
    prefix = options.prefix
else:
    prefix = "experiment"
    
if options.start == None:
    parser.print_help()
    print "\nRun_Fcst Script: Start time not supplied...using start time in WRF namelist (this should be the first integration)\n"
else:
    start = options.start
    print "\nRun_Fcst Script: Start time supplied...modifying start/end times in WRF namelist.input\n"

if options.stop == None:
    parser.print_help()
    print "\nRun_Fcst Script: Stop time not supplied...using stop time in WRF namelist\n"
else:
    stop = options.stop

if options.ne == None:
    ne = 1
else:
    ne = options.ne

cpu_model = 0.0
c0 = cpu.time()

print "\n ens_fcst_wrf.py: Integrating ensemble members from a local time of  %s to  %s \n" % (start, stop)

end_year = int(stop.split("-")[0])
end_month = int(stop.split("-")[1])
end_day = int(stop.split("_")[0].split("-")[2])
end_hour = int(stop.split("_")[1].split(":")[0])
end_minute = int(stop.split("_")[1].split(":")[1])
end_second = int(stop.split("_")[1].split(":")[2])

if options.start != None:

    start_year = int(start.split("-")[0])
    start_month = int(start.split("-")[1])
    start_day = int(start.split("_")[0].split("-")[2])
    start_hour = int(start.split("_")[1].split(":")[0])
    start_minute = int(start.split("_")[1].split(":")[1])
    start_second = int(start.split("_")[1].split(":")[2])

if 1==1:

   dir = path + "/" + prefix
   os.chdir(dir)
   if debug: print "Changed to ", os.getcwd()

   mempaths = []
   for item in os.listdir(dir):
     if (os.path.isdir(os.path.join(dir, item)) and "QOL" not in item):
       mempaths.append(os.path.join(dir, item))
   mempaths = sorted(mempaths)
   #if debug: print mempaths

   if 1==1:

     num_batch=int(ne/batch_size) # number of batches
     cmd = "rm runwrf*done"
     os.system(cmd)

     for iter in range(0, num_batch):

       c1 = cpu.time()

       n1 = iter*batch_size
       n2 = n1+batch_size
       print mempaths[n1:n2]

       os.chdir(dir)

       pool = Pool(processes=batch_size)
       pool.map(RunMember, mempaths[n1:n2])
       pool.close()
       pool.join()

       count = 0
       while count != batch_size and not debug:
         count=0
         cpu.sleep(60)
         for n in range(n1, n2):
           fname = mempaths[n]+'/run_wrf%03d_done' % (n+1)
           if os.path.isfile(fname): count +=1
         if count>0: print count, " members finished"


   cpu.sleep(5)

cpu_model = cpu_model + cpu.time() - c0

print "TOTAL ENS_FCST_WRF time = ", (cpu.time() - c0)/60, " min"

time = [stop.split("-")[0], stop.split("-")[1], stop.split("_")[0].split("-")[2], stop.split("_")[1].split(":")[0], stop.split("_")[1].split(":")[1], stop.split("_")[1].split(":")[2]]

cmd = 'ls -l %s/%s/*/%s | wc -l' % (path, prefix, stop)
#files = ens.find_wrf_files(path, prefix, time, restart=True)
os.system(cmd)

#if len(files)!=ne:
#   print "%d member forecasts failed!\n" % (ne-len(files))
#   sys.exit(1)


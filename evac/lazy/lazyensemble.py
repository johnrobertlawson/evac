import os
import pdb
import datetime
from pathlib import PosixPath

from evac.lazy.lazywrf import LazyWRF

class LazyEnsemble:
    """
    Generate ensemble by running multiple WRF runs simulataneously.
    This is done by creating temporary folder and linking/copying.
    Efficiency improvements from Monte Flora.
    
    Will not work on Windows systems due to hardcoded PosixPath.

    Parent script should be run via batch submission.
    """
    def __init__(self,path_to_exedir,path_to_datadir, path_to_namelistdir,
                    path_to_icbcdir,path_to_outdir,path_to_batch,initutc,sched='slurm',
                    path_to_rundir=False,):
        """
        Args:
        
        wrfsrcdir       :   directory of compiled WRF executables
        wrfrundir       :   directory for running directory.
                            By default, it is the data directory.
                            By default,
                            (maybe better to have permanant duplicates,
                                rather than spinning up and deleting
                                new wrf.exe)
        namelistdir     :   directory containing namelist(s)
                            templates or full versions.
        icbcdir         :   directory with initial and boundary
                            condition files
        outdir          :   where to move wrfout files to
        path_to_batch   :   absolute path to *.job script (for slurm)
        initutc         :   initialisation time (datetime.datetime)
        """
        # PATH OBJECTS
        self.exedir = PosixPath(path_to_exedir)
        self.namelistdir = PosixPath(path_to_namelistdir)
        self.icbcdir = PosixPath(path_to_icbcdir)
        self.outdir = PosixPath(path_to_outdir)
        self.batchscript = PosixPath(path_to_batch)

        # By default, the run directory is where the data will end up
        if not isinstance(path_to_rundir,str):
            path_to_rundir = path_to_datadir
        self.wrfrundir = PosixPath(path_to_rundir)

        # Time of initialisation
        self.initutc = initutc

        # Lookup dictionary of all members
        self.members = {}
        
    def copy(self,*args,**kwargs):
        """ Wrapper for self.bridge()
        """
        kwargs['copy'] = True
        return self.bridge(*args,**kwargs)
    
    def softlink(self,frompath,topath):
        pass
    
    def move(self,frompath,topath):
        pass
    
    def run_wps(self):
        print("Not written yet.")
        return
    
    def run_wrf(self):
        """
        Inherit from parent?
        """
        pass
    
    def print_readme(self,outdir=self.outdir,):
        """
        Save a pretty-printed README file to a given directory.
        
        Optional argument will print namelist.input
        
        Maybe this could be a decorator method that logs this data for each
        ensemble member.
        """
        pass
    
    @staticmethod
    def edit_batchscript(fpath,linekey,newline):
        """ Replace a line from the submission script
        that matches a key (linekey) with newline.
        """
        fs = open(fpath,'r')
        flines = fs.readlines()
        for idx, line in enumerate(flines):
            if linekey in line:
                fs.close()
                flines[idx] = "{}\n".format(newline)
                nameout = open(f,'w',1)
                nameout.writelines(flines)
                nameout.close()
                return
        raise ValueError("Setting",linekey,"not found in script.")
        return

            
            
                

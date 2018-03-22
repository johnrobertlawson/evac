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
                    path_to_icbcdir,path_to_outdir,path_to_batch,initutc,
                    sched='slurm',path_to_rundir=False,delete_exe_copy=False,
                    ):
        """
        Args:
        
        path_to_exedir      :   (str) - directory of compiled WRF executables
        path_to_datadir     :   (str) - directory where wrfout and other data
                                will be saved.
        path_to_namelistdir :   (str) - directory containing namelist(s)
                                templates or full versions.
        path_to_icbcdir     :   (str) - directory with initial and boundary
                                condition files
        path_to_outdir      :   (str) - where to move wrfout files to
        path_to_batch       :   (str) - absolute path to *.job script
                                for slurm - not sure about rocks etc
        initutc             :   (datetime.datetime) - initialisation time
        
        path_to_rundir      :   (str, optional) - directory where wrf.exe will
                                be copied to,
                                and rsl.error files will be generated.
                                By default, this is the datadir.
                                If running multiple ensembles, it is
                                faster to have permanant folders for each
                                ensemble, rather than spinning up and
                                deleting wrf.exe)
        sched               :   (str, optional) - job scheduler.
                                Only implemented for slurm right now
        delete_exe_copy     :   (bool, optional) - whether the .exe files
                                will be deleted after the run is finished.
                                The default is no, as numerous ensembles
                                can use the same run folders, potentially.
        """
        # PATH OBJECTS
        self.exedir = PosixPath(path_to_exedir)
        self.datadir = PosixPath(path_to_datadir)
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

        self.sched = sched
        
        # Options
        self.delete_exe_copy = delete_exe_copy
        
        # INIT PROCEDURE
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

            
            
                

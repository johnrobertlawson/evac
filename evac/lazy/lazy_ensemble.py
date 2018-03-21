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
    """
    def __init__(self,path_to_exedir,path_to_datadir,path_to_rundir=False,namelistdir,
                    icbcdir,outdir,initutc,):
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
        initutc         :   initialisation time (datetime.datetime)
        """
        
        self.wrfsrcdir = PosixPath(wrfsrcdir)
        self.wrfrundir = wrfrundir
        self.namelistdir = namelistdir
        self.icbcdir = icbcdir
        self.outdir = outdir
        self.initutc = initutc
        
    def bridge(self,):
        """ Wrapper for utils.unix_tools.bridge. This copies, moves, or links
        files, depending on arguments.
        """
        
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
    
    
    
        
import os
import pdb
import datetime


from evac.lazy.lazywrf import LazyWRF

class LazyEnsemble:
    """
    Generate ensemble by running multiple WRF runs simulataneously.
    This is done by creating temporary folder and linking/copying.
    Efficiency improvements from Monte Flora.
    """
    def __init__(self,wrfsrcdir,wrfrundir,namelistdir,
                    icbcdir,outdir,initutc):
        """
        Args:
        
        wrfsrcdir       :   directory of compiled WRF executables
        wrfrundir       :   new directory for our WRF run(s)
        namelistdir     :   directory containing namelist(s)
                            templates or full versions.
        icbcdir         :   directory with initial and boundary
                            condition files
        outdir          :   where to move wrfout files to
        initutc         :   initialisation time (datetime.datetime)
        """
        
        self.wrfsrcdir = wrfsrcdir
        self.wrfrundir = wrfrundir
        self.namelistdir = namelistdir
        self.icbcdir = icbcdir
        self.outdir = outdir
        self.initutc = initutc
        
    def bridge(self,):
        """ Wrapper for utils.unix_tools.bridge. This copies, moves, or links
        files, depending on arguments.
        """
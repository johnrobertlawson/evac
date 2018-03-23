import os
import pdb
import datetime
from pathlib import PosixPath

import evac.utils as utils
from evac.lazy.lazywrf import LazyWRF
from evac.utils.exceptions import WRFError

class LazyEnsemble:
    """
    Generate ensemble by running multiple WRF runs simulataneously.
    This is done by creating temporary folder and linking/copying.
    Efficiency improvements from Monte Flora.

    Will not work on Windows systems due to hardcoded PosixPath.

    Parent script should be run via batch submission. The procedure is:

    L = LazyEnsemble(*args,**kwargs)
    L.run_all_members()

    """
    def __init__(self,path_to_exedir,path_to_datadir, path_to_namelistdir,
                    path_to_icbcdir,path_to_outdir,path_to_batch,initutc,
                    sched='slurm',path_to_rundir=False,delete_exe_copy=False,
                    ndoms=1,nmems=0,membernames=False):
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
        ndoms               :   (int) - number of domains
        nmems               :   (int) - number of members. If zero,
                                try counting length of membernames.
        membernames         :   (list,tuple) - list of strings for
                                member names. If False, use automatic
        """
        # Check - must have either number of members or a list of names
        assert isinstance(membernames,(list,tuple)) or (nmems > 0)

        # PATH OBJECTS
        self.exedir = PosixPath(path_to_exedir)
        self.datadir = PosixPath(path_to_datadir)
        self.namelistdir = PosixPath(path_to_namelistdir)
        self.icbcdir = PosixPath(path_to_icbcdir)
        self.outdir = PosixPath(path_to_outdir)
        self.batchscript = PosixPath(path_to_batch)

        # Shortcut for accessing the script name
        self.batchname = self.batchscript.name
        # By default, the run directory is where the data will end up
        if not isinstance(path_to_rundir,str):
            path_to_rundir = path_to_datadir
        self.wrfrundir = PosixPath(path_to_rundir)

        # Time of initialisation
        self.initutc = initutc

        self.sched = sched

        # Number of domains
        self.ndoms = ndoms

        # Options
        self.delete_exe_copy = delete_exe_copy

        # INIT PROCEDURE
        # Get member names
        if self.membernames is False:
            self.nmems = nmems
            self.membernames = ['mem{:02d}'.format(n) for n in
                                N.arange(1,self.nmems+1)]
        else:
            self.membernames = membernames
            self.nmems = len(self.membernames)

        # Lookup dictionary of all members
        self.members = catalog_members()

        # TO DO - how to automate nests too
        # For now, everything in same folder.
        utils.wowprint("Ensemble **created!**",color='purple',bold=True)
        utils.wowprint("Ensemble **created!**",color='purple',bold=True)

    def catalog_members(self,):
        """ Create the members dictionary.
        """
        members = {}
        for mem in N.arange(self.membernames):
            mem_datadir = self.datadir / mem
            utils.trycreate(mem_datadir)

            mem_rundir = self.wrfrundir / mem
            members[mem] = {'datadir':mem_datadir,
                                    'rundir':mem_rundir,}
        return members

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

        linekey must be unambiguous, otherwise it will change the
        first occurence in the file that matches.
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

    def batchscript_for_member(self,member):
        """ Edit parts of the batch submission script.
        """
        jobname = 'wrf_{}'.format(member)
        path_to_err = (self.wrfrundir / jobname).with_suffix('.err')
        path_to_out = (self.wrfrundir / jobname).with_suffix('.out')
        path_to_wrfexe = self.wrfrundir / wrf.exe
        command = f'time srun {path_to_wrfexe})'
        rundir = self.members[member]['rundir']
        batchpath = rundir / self.batchname
        cpu = self.cpus_per_job

        # Need to change tasks per node?
        #tn = self.tasks_per_job

        # First, changing the job name.
        # Note the print literals marked with "f" from Python 3.6)
        changes = dict("#SBATCH -J" = f"#SBATCH -J {jobname}",

        # Next, the output and error files
                "#SBATCH -o" = f"#SBATCH -o {path_to_err}",
                "#SBATCH -e" = f"#SBATCH -e {path_to_out}",
        # Next the cpu/node settings
                #f"#SBATCH --ntasks-per-node" = "#SBATCH --ntasks-per-node={tn}",
                "#SBATCH -n" = f"#SBATCH -n {cpu}",
        # Make sure we're in the right directory
                "cd" = f"cd {rundir}",
        # Next the wrf submit command
                -1 = command)

        for key,newline in changes.items():
            edit_batchscript(self,batchpath,key,newline)
        return


    def namelist_for_member(self,member):
        """ Edit parts of the namelist

        start/end year, month, day, hour, minute, second

        Use multiplied notation by number of domains (self.ndoms)
        """

    def run_all_members(self,cpus=0,nodes=1,):
        """ Automatically submit all jobs to the scheduler.
        """
        self.cpus_per_job = cpus
        self.nodes_per_job = nodes
        for member in self.members:
            self.run_wrf_member(member)

        # Generate README for each dir?
        return

    def run_wrf_member(self,member,prereqs,dryrun=False):
        """ Submit a wrf run to batch scheduler.

        member  :   (str) - name of member to run
        prereqs :   (dict) - a dictionary containing prerequisite files
                    needed in rundir. use this format:

                    dict( path_to_file = cmd)

                    where path_to_file is the absolute path to
                    the file in question (as str or Path obj) and
                    cmd is from ['copy','move','softlink'].
        dryrun  :   stop just before submission in debug mode.
        """
        # These two may well be identical.
        rundir = self.members[member]['rundir']
        datadir =  self.members[member]['datadir']

        # Get dictionary in correct format for bridge_multi
        PRQs = {}
        for k,v in prereqs.items():
            PRQs[(k,rundir)] = v

        # Copy, link, move everything needed to rundir
        utils.bridge_multi(PRQs)

        # Edit batch script
        # use self.cpus_per_job and self.nodes_per_job

        # Edit namelist?

        # Check output directory exists
        # and is written in namelist

        # Submit script
        batchloc = datadir / self.batchname
        cmd = "sbatch f{batchloc}"
        if dryrun:
            pdb.set_trace()
            raise Exception("Exiting - dry run.")
        os.system(cmd)

        # Monitor processes

        # Create README?

        # Clean up files

    def cleanup(self,folder,files):
        """ Deletes files in given directory that match a glob.
        """
        pass

    def check_complete(self,member,raise_error=False,sleep=10,firstwait=3600,
                        maxtime=(24*3600)):
        """ Check ensemble member to see if wrf run has finished.

        Args:
        raise_error     :   (bool) - if True, pass a WRFError as return
                            if WRF breaks.
        sleep           :   (int) - number of seconds to wait between checks
        firstwait       :   (int) - number of seconds to wait after submission
                            before checks
        maxtime         :   (int) - number of seconds after which the
                            method will assume the member's run has died
        """
        # Don't bother wasting resources etc checking for a run to finish
        time.sleep(firstwait)
        while True:
            rsl = self.members[member]['rundir'] / 'rsl.error.0000'vim


    def run_real_member(self,):
        """ Submit real.exe to batch scheduler.
        """
        pass

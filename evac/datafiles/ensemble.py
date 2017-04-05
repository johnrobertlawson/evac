""" Ensemble of forecasts.

Todos:
    * Change from ensemble of WRFOut files to any sort of datafile.
    * Consistency with method names (e.g., get() )
"""

import os
import pdb
import datetime

import numpy as N

import WEM.utils as utils
from .wrfout import WRFOut

# Dummy variable in place of proper subclass of WRFOut
AuxWRFOut = object

class Ensemble:
    def __init__(self,rootdir,initutc,doms=1,ctrl='ctrl',aux=False,
        model='wrf',fmt='em_real',f_prefix=None,loadobj=True,
        ncf=False,debug=False):
        """Class containing all ensemble members. Default is a
            deterministic forecast (i.e. ensemble of one control member).
            Each ensemble member needs to have a separate folder (named
            as the ensemble member's name), but different domains can be
            within a single member's folder.

        Args:
            rootdir (str):  Directory at root of datafiles
            initutc (datetime.datetime): Initialization time
            doms (int, optional): Number of domains
            ctrl (bool, optional): Whether ensemble has control member
            aux (bool, dict, optional): Dictionary lists, per domain, data
                files that contain additional variables and each
                file's prefix. Default is False (no auxiliary files).
                Not implemented yet.
            enstype (str, optional): Type of ensemble. Default is for
                Weather Research and Forecast (WRF) model.
            fmt (str, optional): The type of simulations. Default is
                real-world simulations (em_real from WRF).
            f_prefix (tuple, optional): Tuple of prefixes for each
                ensemble member's main data files. Must be length /doms/.
                Default is None, which then uses a method to determine
                the file name using default outputs from e.g. WRF.
        """
        self.model = model.lower()
        if self.model == 'wrf':
            print("File type is wrfout.")
            loadfunc = WRFOut
        else:
            raise NotImplementedError()

        self.debug = debug
        self.ctrl = ctrl
        self.rootdir = rootdir
        self.initutc = initutc
        self.doms = doms
        self.fmt = fmt
        self.loadobj = loadobj
        self.aux = aux
        self.ncf = ncf

        self.isaux = True if isinstance(self.aux,dict) else False
        if f_prefix is not None and len(f_prefix) is not doms:
            raise Exception("Length of main datafile prefixes must "
                                "match number of domains.")
        self.isctrl = True if ctrl else False
        # self.ndoms = len(self.doms)
        self.member_names = []
        self.members, self.fdt = self.get_members()
        self.nmems = len(self.member_names)
        self.nperts = self.nmems - self.isctrl

        # Get start and end of whole dataset (inclusive)
        self.filetimes = self.list_of_filetimes(arb=True)
        self.nt_per_file = self.compute_nt_per_file()
        self.itime = self.filetimes[0]
        self.hdt = self.compute_history_dt()
        # pdb.set_trace()
        if self.fdt is None:
            # self.ftime = self.filetimes[-1] + datetime.timedelta(seconds=self.hdt)
            self.timespan_sec = self.hdt*(self.nt_per_file-1)*len(self.filetimes)
            self.ftime = self.itime + datetime.timedelta(seconds=self.timespan_sec)
        else:
            self.ftime = self.filetimes[-1] + (
                    (self.nt_per_file-1)*datetime.timedelta(seconds=self.fdt))

        # Difference in output times across whole dataset
        # Might be split between history files


    def compute_fdt(self):
        """Compute the difference in time between each data file's first entry.

        Returns the difference in seconds.
        """
        f_diff = self.filetimes[1] - self.filetimes[0]
        return f_diff.seconds

    def get_members(self,):
        """Create a dictionary with all data.

        Format is:
        members[member][domain][time][data]

        Returns members (dict): Dictionary of ensemble members

        Also returns fdt (int): Seconds between output files.
        """
        members = {}
        fdt = None
        for dom in range(1,self.doms+1):
            # Get file name for initialisation time and domain
            # main_fname = self.get_data_fname(dom=dom,prod='main')
            if not self.ncf:
                main_fname  = utils.get_netcdf_naming(self.model,self.initutc,dom)
            else:
                main_fname = self.ncf
            # if dom in self.aux:
                # aux_fname = self.get_data_fname(dom=dom,prod='aux')
            # Each ensemble member has a domain
            for dirname,subdirs,files in os.walk(self.rootdir):
                # If ensemble size = 1, there will be no subdirs.
                # pdb.set_trace()
                if main_fname in files:
                    # pdb.set_trace()
                    dsp =  dirname.split('/')
                    rsp = self.rootdir.split('/')

                    # The following logic merges subdir names
                    # e.g. ensemble grouped twice by perturbation type?
                    if dsp[:-1] == rsp:
                        member = dirname.split('/')[-1]
                    elif dsp[:-2] == rsp:
                        member = '_'.join(dirname.split('/')[-2:])
                    else:
                        # pdb.set_trace()
                        # raise Exception("What is this folder structure?", dsp, rsp, dirname)
                        print("Skipping file in {}".format(dirname))
                        continue
                    if self.debug:
                        print("Looking at member {0}".format(member))
                    if member not in members:
                        members[member] = {dom:{}}
                    if dom==1:
                        self.member_names.append(member)
                    t = self.initutc
                    while True:
                        if not self.ncf:
                            t_fname = utils.get_netcdf_naming(self.model,t,dom)
                        else:
                            t_fname = self.ncf
                        # Check for history output time
                        fpath = os.path.join(self.rootdir,dirname,t_fname)
                        try:
                            dataobj = self.datafile_object(fpath,loadobj=self.loadobj)
                        except IOError:
                            # All wrfout files have been found
                            break
                        else:
                            # print("Assigning file path and maybe object.")
                            members[member][dom][t] = {'dataobj':dataobj,
                                                'fpath':fpath,
                                                'control': (member is self.ctrl)}
                            # print("Done.")
                        if (self.aux is not False) and (dom in self.aux):
                            # TODO: implement
                            fpath = os.path.join(self.rootdir,dirname,aux_fname)
                            dataobj = self.datafile_object(fpath,loadobj=self.loadobj)
                            members[member][dom][t]['auxdataobj'] = dataobj
                            members[member][dom][t]['auxfpath'] = fpath
                            members[member][dom][t]['control'] = member is self.ctrl
                        # Move to next time
                        if not self.ncf:
                            if fdt is None:
                                f1, f2 = sorted(files)[:2]
                                fdt = utils.dt_from_fnames(f1,f2,'wrf')

                                # Loop through files and estimate dt based on fname
                            else:
                                t = t + datetime.timedelta(seconds=fdt)
                        else:
                            break

        return members, fdt

    def datafile_object(self,fpath,loadobj=False,**kwargs):
        #Extend to include other files (GEFS, RUC etc)
        #TODO: Implement auxiliary wrfout files
        # print(fpath)
        if loadobj:
            ops = {'wrf':WRFOut,'aux':AuxWRFOut}
            answer = ops[self.model](fpath,**kwargs)
        else:
            os.stat(fpath)
            answer = False
        return answer

    def get_gefs_ensemble(self):
        """
        All gefs data files should be in the same folder.
        Each forecast time is a different file.
        """
        members = {}
        # This will break with subset of members, times, a/b grib files...
        allmembers = ['c00',] + ['p{0:02d}'.format(n) for n in range(1,21)]
        allfiles = glob.glob(os.path.join(self.ncroot,'ge*'))
        alltimes = N.arange(0,390,6)

        for ens in allmembers:
            self.memnames.append(ens)
            for t in alltimes:
                utc = self.initt + datetime.timedelta(hours=int(t))
                fname = 'ge{0}.t{1:02d}z.pgrb2f{2:02d}.nc'.format(
                            ens,self.initt.hour,t)
                if os.path.join(self.rootdir,fname) in allfiles:
                    if ens not in members.keys():
                        members[ens] = {}
                    members[ens][utc] = {'data':GEFS(os.path.join(self.rootdir,
                                fname)),'control':ens is 'c00'}
        return members

    def get_exceedance_probs(self,vrbl,overunder,threshold,
                            level=None,itime=None,
                            ftime=None,fcsttime=None,Nlim=None,
                            Elim=None,Slim=None,Wlim=None,
                            dom=1):
        """
        Return probability of exceeding or reaching a threshold.

        Args:
            vrbl (str,N.array)      : variable. If N.array, use provided data
                                        (i.e. override the loading).
                                        Must be 5D.
                                        [ensemble_member,time,level,lat,lon]
            overunder (str)         : threshold evaluation
                                        'over', 'overeq','under','undereq','equal'
            threshold (float,int)   : the threshold in same (SI) units as data
            itime (datetime.datetime)   : initial time
            ftime (datetime.datetime)   : final time
            fcsttime (datetime.datetime): single time (must be present if itime
                                            and ftime are None.)
            Nlim, Elim, Slim, Wlim (float,optional) : bounding box of lats/lons
        """
        if not overunder in OU.keys():
            raise Exception("Pick over or under for threshold comparison.")

        if isinstance(vrbl,N.ndarray):
            assert len(vrbl.shape) == 5
            all_ens_data = vrbl
        else:
            all_ens_data = self.get(vrbl,level=level,itime=itime,ftime=ftime,
                            fcsttime=fcsttime,Nlim=Nlim,Elim=Elim,
                            Slim=Slim,Wlim=Wlim,dom=dom)

        OU = {'over':__gt__,'under':__lt__,'overeq':__ge__,'undereq':__le__,
                'equal':__eq__}

        # True/False if member meets condition (5D)
        bool_arr = N.where(all_ens_data.OU[overunder](threshold),1,0)

        # Count members that exceed the threshold (4D)
        count_arr = N.sum(bool_arr,axis=0)

        # And convert to percentage (4D) for each time
        percent_arr = 100*(count_arr/self.nmems)

        return percent_arr # [times,levels,lats,lons]


    def closest_to_mean(self,vrbl,level,fcsttime,Nlim=False,Elim=False,
                            Slim=False,Wlim=False,):
        """
        Find closest member to the mean for given variable
        Passing latlon/box allows calculation over given area
        (Box is in grid spaces)
        """
        all_ens_data = self.get(vrbl,level=level,fcsttime=fcsttime,
                                          Nlim=Nlim,Elim=Elim,
                                          Slim=Slim,Wlim=Wlim)
        mean = N.mean(all_ens_data,axis=0)[0,0,:,:]
        diff = N.abs(N.sum((all_ens_data[:,0,0,:,:]-mean),axis=(1,2)))
        ensidx = N.argmin(diff,axis=0)

        return self.members_names[ensidx]


    def return_DF_for_t(self,utc,member,dom=1):
        t, tidx = self.find_file_for_t(utc,member,dom=dom)
        if self.debug:
            print("Loading data for time {0}".format(utc))
        # pdb.set_trace()
        fpath = self.members[member][dom][t]['fpath']
        DF = self.datafile_object(fpath,loadobj=True)
        return DF

    def get(self,vrbl,level=None,itime=False,ftime=False,
                        fcsttime=False,Nlim=None,Elim=None,
                        Slim=None,Wlim=None,inclusive=False,
                        lats=None,lons=None,dom=1,members=None):
        """
        Returns 5D array of data for ranges.

        Needs to load WRFOut files if self.loadobj is False.

        Ordered in descending order on pert. members
        First dimension is ensemble members.

        Args:
            inclusive (bool, optional): if True, included time specified
                at ftime in the time range. Default is False (like Python).

        Todos:
            * lat/lon box is in the correct projection?
            * Implement bounding lat/lon box.
        """

        ens_no = 0
        # pdb.set_trace()
        if vrbl is 'accum_precip':
            qpf = self.accumulated(vrbl='RAINNC',itime=itime,ftime=ftime,
                            level=level,Nlim=Nlim,Elim=Elim,
                            Slim=Slim,Wlim=Wlim,inclusive=inclusive,
                            lons=lons,lats=lats)
            return qpf
        if members is None:
            members = self.member_names
        elif isinstance(members,str):
            members = (members,)
        else:
            pass


        for nm,mem in enumerate(members):
            if self.debug:
                print("Working on member {0}".format(mem))
            if mem is self.ctrl:
                print("Skipping control member.")
                continue
            else:
                ens_no += 1

               # if itime and ftime:
                if isinstance(itime,datetime.datetime) and isinstance(
                            ftime,datetime.datetime):
                    # fts = N.arange(itime,ftime,self.hdt)
                    fts = utils.generate_times(itime,ftime,self.hdt,
                                inclusive=inclusive)
                else:
                    fts = [fcsttime,]

                for tn, ft in enumerate(fts):
                    t, tidx = self.find_file_for_t(ft,mem,dom=dom)
                    if self.debug:
                        print("Loading data for time {0}".format(ft))
                    fpath = self.members[mem][dom][t]['fpath']
                    DF = self.datafile_object(fpath,loadobj=True)
                    m_t_data = DF.get(vrbl,utc=tidx,level=level,lons=lons,lats=lats)[0,...]

                if ens_no == 1:
                    nz,nlats,nlons = m_t_data.shape
                    nt = len(fts)
                    all_ens_data = N.zeros((self.nperts,nt,nz,nlats,nlons))

                all_ens_data[ens_no-1,tn,:,:,:] = m_t_data

            return all_ens_data

    def accumulated(self,vrbl='RAINNC',itime=0,ftime=-1,level=False,Nlim=False,
                    Elim=False,Slim=False,Wlim=False,inclusive=False,
                    lons=None,lats=None):
        """Accumulate, for every ensemble member, at each grid point,
        the variable specified. Usually precipitation.

        Todos:
            * Logic to work out if values are for each history output
                timestep, from the start of the simulation, from the
                start of the data file...

        """
        if itime==0:
            itime = self.itime
        if ftime==-1:
            ftime = self.ftime

        if vrbl is 'RAINNC':
            itime_rainnc = self.get('RAINNC',fcsttime=itime,
                    lons=lons,lats=lats)
            ftime_rainnc = self.get('RAINNC',fcsttime=ftime,
                    lons=lons,lats=lats)
            accum = ftime_rainnc - itime_rainnc
        else:
            all_ens_data = self.get(vrbl,itime=itime,ftime=ftime,
                                        inclusive=inclusive,lats=lats,lons=lons)
            # time axis is 1
            accum = N.sum(all_ens_data,axis=1)

        # Resulting matrix is size (nperts,1,nz,nlats,nlons).
        return accum

    def mean(self,vrbl,fcsttime=False,level=False,Nlim=False,Elim=False,
             Slim=False,Wlim=False,itime=False,ftime=False):
        """
        Returns mean.
        """
        all_ens_data = self.get(vrbl,level=level,fcsttime=fcsttime,
                                    Nlim=Nlim,Elim=Elim,Slim=Slim,Wlim=Wlim,
                                    itime=itime,ftime=ftime)
        mean = N.mean(all_ens_data,axis=0)

        if Nlim:
            return mean, lats, lons
        else:
            return mean

    def std(self,vrbl,fcsttime=False,itime=False,ftime=False,level=False,
            Nlim=False,Elim=False,Slim=False,Wlim=False):
        """Return standard devation
        """
        all_ens_data = self.get(vrbl,level=level,fcsttime=fcsttime,
                                    Nlim=Nlim,Elim=Elim,Slim=Slim,Wlim=Wlim,
                                    itime=itime,ftime=ftime)
        std = N.std(all_ens_data,axis=0)

        if Nlim:
            return std, lats, lons
        else:
            return std

    def arbitrary_pick(self,dataobj=False,give_keys=False,give_path=False):
        """Arbitrary pick of a datafile entry in the members dictionary.

        Args:
            dataobj (bool, optional): if True, return the DataFile subclass.
                Otherwise, return filepath.
            give_keys (bool, optional): if True, return a list of
                member, domain, time keys to enter into a dictionary.

        """
        mem = self.member_names[0]
        dom = 1
        t = self.initutc
        arb = self.members[mem][dom][t]
        if dataobj:
            if give_keys:
                raise Exception("Pick only one of give_keys and dataobj.")
            return self.datafile_object(arb['fpath'],loadobj=True)
        elif give_keys:
            return mem, dom, t
        elif give_path:
            return arb['fpath']
        else:
            return arb

    def list_of_filetimes(self,arb=False,member=False,dom=1):
        """Return list of times for each data file's first time entry,
        for a member and domain.

        Args:
            arb (bool, optional): If true, arbitrarily pick a
                member and domain to build the time list, i.e.,
                assuming all members/domains have same list
            member (str, optional): Name of member to build times
                from. Needed if arb is False.
            dom (int, optional): Number of domain to build times
                from. Default is 1.
        """
        if (arb is False) and (member is False):
            raise Exception("Specify member name if not picking arbitrarily.")
        elif arb is True:
            member = self.member_names[0]

        alltimes = sorted(list(self.members[member][dom].keys()))
        return alltimes

    def compute_nt_per_file(self):
        DF = self.arbitrary_pick(dataobj=True)
        return DF.t_dim

    def compute_history_dt(self,):
        """Calculate time difference between each history output
        time. This could be across multiple files or in one.
        """
        if self.nt_per_file == 1:
            hdt = self.fdt
        else:
            # arbitrarily pick data file
            DF = self.arbitrary_pick(dataobj=True)
            hdt = DF.dt
        return hdt


    def find_file_for_t(self,simutc,member,dom=1):
        """Determine file to load given required time.

        Raises exception if history time doesn't exist.

        Args:
            utc (datetime.datetime): Desired time
            member (str): Name of member to look up. If "arb", pick
                a member arbitrarily (this assumes members have
                identical structure of times in the data).
            dom (int, optional): Domain number to look up

        Returns:
            t (datetime.datetime): members dictionary key for right file,
            index (int): Index in that file.

        Todos:
            * Give nearest time (and/or file object) if time doesn't exist.
        """
        if member == 'arb':
            member = self.member_names[0]

        if not ((simutc <= self.ftime) and (simutc >= self.itime)):
            raise Exception("Time outside range of data times.")

        # Returns index of file containing data
        ftidx, tdiff = utils.closest_datetime(self.filetimes,simutc,round='beforeinc')

        # Make positive
        tdiff = abs(tdiff)

        if tdiff == 0:
            assert self.filetimes[ftidx] == simutc
            t = simutc
            tidx = 0
        else:
            t = self.filetimes[ftidx]
            # tidx = int(self.hdt/(self.hdt + tdiff))
            tidx = int(tdiff/self.hdt)

        # pdb.set_trace()

        return t, tidx

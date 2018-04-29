import os
import glob
import pdb
import datetime

import numpy as N

import evac.utils as utils
from evac.datafiles.wrfout import WRFOut
from evac.datafiles.gefs import GEFS

# Dummy variable in place of proper subclass of WRFOut
AuxWRFOut = object

class Ensemble:
    """ Class containing an ensemble of (netCDF, WRF) forecasts.

    A deterministic forecast (i.e. ensemble of one control member) can
    be created by using one member.

    Each ensemble member needs to have a separate folder (named
    as the ensemble member's name), but different domains can be
    within a single member's folder.

    Todos:
        * Change from ensemble of WRFOut files to any sort of datafile.
        * Consistency with method names (e.g., get() )

    Example:
        To load a 5-D array of 10-m wind at the second time (and using defaults)::

            import datetime
            from evac.datafiles.ensemble import Ensemble

            initutc = datetime.datetime(2016,3,31,18,0,0)
            E = Ensemble(rootdir='/path/to/data',initutc=initutc)
            E.get('wind10',utc=1)

    Args:
        rootdir (str):  Directory at root of datafiles
        initutc (datetime.datetime): Initialization time (assumed UTC)
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
            The string 'default' will look for the wrfout_d0* string.
            Set None for a method to determine
            the file name using default outputs from e.g. WRF.
    """
    def __init__(self,rootdir,initutc,ndoms=1,ctrl='ctrl',aux=False,
                model='wrf',fmt='em_real',f_prefix='default',loadobj=True,
                ncf=False,debug=False,onefolder=False,fileformat='netcdf',
                gefsformat=False):
        self.model = model.lower()
        if self.model == 'wrf':
            print("File type is wrfout.")
            loadfunc = WRFOut
        elif self.model == 'gefs':
            self.gefsformat = gefsformat
            print("File type is GEFS.")
            loadfunc = GEFS
        else:
            raise NotImplementedError()

        self.ndoms = ndoms
        self.doms = list(range(1,ndoms+1))
        if f_prefix is 'default':
            f_prefix = ['wrfout_d{0:02d}'.format(d) for d in self.doms]

        self.debug = debug
        self.fileformat = fileformat
        self.ctrl = ctrl
        # Remove any trailing slash!
        self.rootdir = rootdir.rstrip("/")
        self.initutc = initutc
        self.fmt = fmt
        self.loadobj = loadobj
        self.aux = aux
        self.ncf = ncf

        self.isaux = True if isinstance(self.aux,dict) else False
        if f_prefix is not None and len(f_prefix) is not self.ndoms:
            raise Exception("Length of main datafile prefixes must "
                                "match number of domains.")
        self.isctrl = True if ctrl else False
        self.members, self.fdt = self.get_members(f_prefix=f_prefix)
        self.member_names = sorted(self.members.keys())
        self.nmems = len(self.member_names)
        self.nperts = self.nmems - self.isctrl

        # Check to see if the ensemble isn't just an empty list!
        if self.nmems == 0:
            print("Files not found.")
            raise Exception

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

    def get_gefs_members(self,):
        """ Load GEFS grib files that are all in one folder."""
        members = {}
        self.member_names = ['gec00',] + ['gep{:02d}'.format(n)
                                for n in range(1,21)]
        ftimes = ['anl',] + ['f{:02d}'.format(n) for n in range(6,180,6)]
        utcs = [self.initutc + datetime.timedelta(seconds=hr*3600) for
                        hr in range(0,180,6)]
        for member in self.member_names:
            members[member] = {1:{}}
            for t,ft in zip(utcs,ftimes):
                suffx = '.'.join((self.gefsformat,ft))
                try_fname = '{0}.t{1:02d}z.{2}'.format(
                # try_fname = '{0}.t{1:02d}z.{2}{3}'.format(
                                # member,self.initutc.hour,self.gefsformat,ft)
                                member,self.initutc.hour,suffx)
                try_fpath = os.path.join(self.rootdir,try_fname)
                if os.path.isfile(try_fpath):
                    dataobj = self.datafile_object(try_fpath,loadobj=self.loadobj)
                    members[member][1][t] = {'dataobj':dataobj,
                                            'fpath':try_fpath,
                                            'control': (member is self.ctrl)}
        fdt = 6 # fdt is 6 for gefs - at least for sensible ranges?
        # pdb.set_trace()
        return members, fdt


    def get_members(self,f_prefix=None):
        """Create a dictionary with all data.

        Format is:
        members[member][domain][time][data]

        Returns members (dict): Dictionary of ensemble members

        Also returns fdt (int): Seconds between output files.
        """
        members = {}
        if self.model == 'gefs':
            members,fdt = self.get_gefs_members(f_prefix=f_prefix)
        elif self.model == 'wrf':
            members,fdt = self.get_wrf_members(f_prefix=f_prefix)
        else:
            raise NotImplementedError
        return members, fdt

    def get_wrf_members(self,f_prefix=None):
        members = {}
        fdt = None
        # TODO: refactor to use self.
        doms = self.doms 
        for domn, dom in enumerate(doms):
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
                if isinstance(f_prefix,str):
                    files = [f for f in files if f.startswith(f_prefix)]
                elif isinstance(f_prefix,(list,tuple)):
                    files = [f for f in files if f.startswith(f_prefix[domn])]
                if main_fname in files:
                    print("Found",main_fname,"in",dirname)
                    dsp =  dirname.split('/')
                    rsp = self.rootdir.split('/')

                    # The following logic merges subdir names
                    # e.g. ensemble grouped twice by perturbation type?
                    # pdb.set_trace()
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
                        members[member] = {d:{} for d in doms}
                    # if dom==1:
                        # self.member_names.append(member)
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
                                if len(files) == 1:
                                    break
                                else:
                                    # Loop through files and estimate dt based on fname
                                    f1, f2 = sorted(files)[:2]
                                    fdt = utils.dt_from_fnames(f1,f2,'wrf')
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
            ops = {'wrf':WRFOut,'aux':AuxWRFOut,'gefs':GEFS}
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
        allfiles = glob.glob(os.path.join(self.rootdir,'ge*'))
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

    def get(self,vrbl,level=None,utc=None,itime=False,ftime=False,
                        fcsttime=False,Nlim=None,Elim=None,
                        Slim=None,Wlim=None,inclusive=False,
                        lats=None,lons=None,dom=1,members=None,
                        accum_hr=1,fcsthr=None):
        """
        Returns 5D array of data for ranges.

        Needs to load WRFOut files if self.loadobj is False.

        Ordered in descending order on pert. members
        First dimension is ensemble members.

        Args:
            inclusive (bool, optional): if True, included time specified
                at ftime in the time range. Default is False (like Python).
            fcsthr: in hours
            fcsttime: as datetime.datetime
            itime,ftime: datetime.datetime or integer indices
        Todos:
            * lat/lon box is in the correct projection?
            * Implement bounding lat/lon box.

            * fcsttime and utc are the same to maintain compatibility
            because get() APIs in different areas of evac/WEM
        """
        if fcsthr and (not fcsttime) and (not itime) and (not ftime) and (not utc):
            fcsttime = self.initutc + datetime.timedelta(seconds=3600*fcsthr)
        elif utc is not None:
            fcsttime = utc
        ens_no = 0
        # pdb.set_trace()
        if vrbl is 'accum_precip':
            qpf = self.accumulated(vrbl='RAINNC',itime=itime,ftime=ftime,
                            level=level,Nlim=Nlim,Elim=Elim,
                            Slim=Slim,Wlim=Wlim,inclusive=inclusive,
                            lons=lons,lats=lats,dom=dom,fcsttime=fcsttime,
                            accum_hr=accum_hr)
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

        # pdb.set_trace()
        return all_ens_data

    def accumulated(self,vrbl='RAINNC',itime=0,ftime=-1,level=False,Nlim=False,
                    Elim=False,Slim=False,Wlim=False,inclusive=False,
                    lons=None,lats=None,dom=1,fcsttime=None,accum_hr=1):
        """Accumulate, for every ensemble member, at each grid point,
        the variable specified. Usually precipitation.

        Todos:
            * Logic to work out if values are for each history output
                timestep, from the start of the simulation, from the
                start of the data file...

        """
        if (itime == False) and (ftime == False) and (fcsttime != None):
            ftime = fcsttime
            itime = ftime - datetime.timedelta(seconds=3600*accum_hr)
        elif itime==0:
            itime = self.itime
        elif ftime==-1:
            ftime = self.ftime

        if vrbl is 'RAINNC':
            itime_rainnc = self.get('RAINNC',fcsttime=itime,
                    lons=lons,lats=lats,dom=dom)
            ftime_rainnc = self.get('RAINNC',fcsttime=ftime,
                    lons=lons,lats=lats,dom=dom)
            accum = ftime_rainnc - itime_rainnc
        else:
            all_ens_data = self.get(vrbl,itime=itime,ftime=ftime,
                                        inclusive=inclusive,lats=lats,lons=lons,
                                        dom=dom)
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

    def arbitrary_pick(self,dataobj=False,give_keys=False,give_path=False,dom=1):
        """Arbitrary pick of a datafile entry in the members dictionary.

        Args:
            dataobj (bool, optional): if True, return the DataFile subclass.
                Otherwise, return filepath.
            give_keys (bool, optional): if True, return a list of
                member, domain, time keys to enter into a dictionary.

        """
        mem = self.member_names[0]
        t = self.initutc
        # pdb.set_trace()
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

    def get_limits(self,dom=1,fmt='dict'):
        """ Return the limits of the domain.

        Args:

        dom     :   (int) - domain number to return
        fmt     :   (str) - format to return data.
                    If 'dict', dictionary.
        """
        W = self.arbitrary_pick(dataobj=True,dom=dom)
        if fmt is 'dict':
            return {'Nlim':W.Nlim,
                    'Elim':W.Elim,
                    'Slim':W.Slim,
                    'Wlim':W.Wlim,}
        else:
            raise Exception

    def get_latlons(self,dom=1):
        W = self.arbitrary_pick(dataobj=True,dom=dom)
        return W.lats, W.lons

    # def get_xx_yy(self,dom=1):
        # W = self.arbitrary_pick(dataobj=True,dom=dom)
        # return W.xx, W.yy

    def get_nx_ny(self,dom=1):
        W = self.arbitrary_pick(dataobj=True,dom=dom)
        # return W.nx, W.ny
        return W.x_dim, W.y_dim

    def get_cenlatlons(self,dom=1):
        W = self.arbitrary_pick(dataobj=True,dom=dom)
        return W.cen_lat, W.cen_lon

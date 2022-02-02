"""Load ensemble data (one model or a superensemble).

TODO:
    * Write and update data classes in this directory
    * Add in RawArray and fix parallel processing via staticmethods?

JRL, Valparaiso Univ., 2022
"""

import os
import glob
import pdb
import datetime
import copy
import multiprocessing
from multiprocessing import RawArray
import collections
import functools
import operator

import numpy as np

from evac.utils.grid import Grid
import evac.utils as utils
from evac.datafiles.wrfout import WRFOut
from evac.datafiles.gefs import GEFS

class Ensemble:
    def __init__(self,data_inf,ctrl=False,aux=None,lazy_load=True):
        """ Class containing an ensemble of forecasts.
    
        A deterministic forecast (i.e. ensemble of one control member) can
        be created by using one member. Information on each member must first
        be passed during instantiation. Data is converted to the 
        multiprocessing.RawArray format to enable parallelisation and try to
        avoid running out of memory. 
    
        Todos:
            * Change from ensemble of WRFOut files to any sort of datafile.
            * Consistency with method names (e.g., get() )
            * Enable lagged ensembles - the initialisation data must be
                loaded from the metadata and attached to the member
                attributes herein.
    
        Example:
                To load a 5-D array [ensemble_member,time,level,lat,lon] of 
                10-m wind at the second time (and using other defaults)::
    
                import datetime
                from evac.datafiles.ensemble import Ensemble
                
                # Create dictionary list of data information
                # Only one member here, called ctrl
                datainf = [{"model":"GFS",
                        "initutc": datetime.datetime(2016,3,31,18,0,0),
                        "path":"/path_to/data",name="ctrl"},]
    
                E = Ensemble(datainf)
                E.get('wind10',utc=1)
    
        Args:
            data_inf (list,tuple): a collection of dictionaries for loading 
                ensemble data. Each dictionary is a member, and must contain
                the model ("GEFS","GFS","WRF","ECMWF") and the absolute
                filepath. Optional keys are ["name","ctrl","wrf_type",
                "np_array","utc_list","init_utc"]. Include "name" as the
                ensemble-member name, and each must be unique; if this is not
                specified, each member will be named in sequential order of the
                list, with the first as "ctrl" if ctrl is True. 
                If the model is generic WRF output, the type of simulation must
                be specified, passed as"ideal" for idealised WRF 
                and "em_real" for the real-world simulations ("em_real" 
                from WRF). The key {"ctrl":True} sets 
                that member as a control member. The key "np_array", with
                a value of a 4-D numpy array, adds a member that skips the
                data format identification. Of course, this data cannot
                be lazy-loaded, and a seequence of times corresponding to
                each of the time dimension must be assigned to the member with
                "utc_list". For other files, if multiple times are loaded
                (e.g., a separate netCDF file for each valid time), use
                "init_utc" to assign a datetime object as the initialisation
                for that file or data. 
            lazy_load (bool): if False, load the filetype upon instantiation.
                Otherwise, lazy load, i.e., identify the class but do not
                load data. The data is acquired later via get().
        """        
        self.data_inf = data_inf
        self.lazy_load = lazy_load
        
        
        self.mandatory_keys = ["model","fpath",]
        self.optional_keys = ["name","ctrl","wrf_type","np_array","utc_list",
                                "init_utc"]
        self.member_keys = self.mandatory_keys + self.optional_keys
        
        self.members = self.load_metadata()
        self.n_members = len(self.members)
        
    def load_metadata(self,):
        member_list = []
        member_names = []
        
        for mem_dict in self.data_inf:
            assert key in mem_dict.keys() for key in self.mandatory_keys
            assert key in self.optional_keys for key in mem_dict.keys()
            member_ntuple = self.get_member_namedtuple(mem_dict,len(
                                    member_list))
            # To double-check name uniqueness
            assert member_ntuple.name is not in member_names
            member_names.append(member_ntuple.name)
            
            member_list.append(member_ntuple)
        return member_list
        
    def get_member_namedtuple(self,mem_dict,n):
        EnsembleMember = collections.namedtuple("EnsembleMember",
            self.member_keys,defaults=["None",]*len(self.optional_keys))
        name = self.make_member_name(mem_dict,n)
        return ntuple
        
    def make_member_name(self,mem_dict,n):
        if "name" in mem_dict.keys():
            return mem_dict["name"]
        else:
            return f"p{n+1:03d}"
        
    def get_dims(self,member_name):
        """Identify the size of data to assign to metadata.
        """
        # get() always returns 5D, with first dimension being n_members.
        nt,nz,nlats,nlons = self.get(member_name,vrbl="T").shape[1:]
        return nt,nz,nlats,nlons

    def compute_fdt(self,member_name):
        """Compute the difference in time.

        Returns the difference in seconds.
        
        Use for models with numerous times (GFS?) and individual files for
        each time (ICON?)
        """
        return

    def datafile_object(self,fpath,model):
        data_classes = {'wrf':WRFOut,'gefs':GEFS}
        return data_classes["model"](fpath)
        
    def get_exceedance_probs(self,vrbl,overunder,threshold,
                            level=None,itime=None,
                            ftime=None,fcsttime=None,Nlim=None,
                            Elim=None,Slim=None,Wlim=None,
                            dom=1,fcsthr=None):
        """
        MOVE TO ANOTHER FOLDER OR SCRIPT
        
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
        OU = {'over':__gt__,'under':__lt__,'overeq':__ge__,'undereq':__le__,
                'equal':__eq__}

        if not overunder in OU.keys():
            raise Exception("Pick over or under for threshold comparison.")

        if isinstance(vrbl,N.ndarray):
            assert len(vrbl.shape) == 5
            all_ens_data = vrbl
        else:
            all_ens_data = self.get(vrbl,level=level,itime=itime,ftime=ftime,
                            fcsttime=fcsttime,Nlim=Nlim,Elim=Elim,
                            Slim=Slim,Wlim=Wlim,dom=dom,fcsthr=fcsthr)


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

    def return_datafile(self,utc,member):
        # Find file for the required time and member
        t, tidx = self.find_file_for_t(utc,member,dom=dom)
        return 

    def get(self,vrbl,level=None,utc=None,itime=False,ftime=False,
                        fcsttime=False,Nlim=None,Elim=None,
                        Slim=None,Wlim=None,inclusive=False,
                        lats=None,lons=None,dom=1,members=None,
                        accum_hr=1,fcsthr=None,fcstmin=None,
                        member=None):
        """
        Returns 5D array of data for ranges.

        Needs to load WRFOut files if self.load_dataobj is False.

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
        if members is None:
            members = member

        if fcsthr and (not fcsttime) and (not itime) and (not ftime) and (not utc):
            fcsttime = self.initutc + datetime.timedelta(seconds=3600*float(fcsthr))
        elif fcstmin is not None:
            fcsttime = self.initutc + datetime.timedelta(seconds=60*float(fcstmin))
        elif utc is not None:
            fcsttime = utc
        # pdb.set_trace()
        if (members is None) or (members == 'all'):
            members = self.member_names
        elif isinstance(members,(str,int)):
            members = [members,]
            assert len(members) == 1
        elif isinstance(members,(list,tuple)):
            pass

        if isinstance(members[0],int):
            name_list = []
            for i,member in enumerate(members):
                name_list.append(self.member_names[i])
            members = name_list

        if vrbl == 'accum_precip':
            qpf = self.accumulated(vrbl='RAINNC',itime=itime,ftime=ftime,
                            level=level,Nlim=Nlim,Elim=Elim,
                            Slim=Slim,Wlim=Wlim,inclusive=inclusive,
                            lons=lons,lats=lats,dom=dom,fcsttime=fcsttime,
                            accum_hr=accum_hr,members=members)
            return qpf

        ens_no = 0
        for nm,mem in enumerate(members):
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
                    DF = self.datafile_object(fpath,load_dataobj=True)
                    try:
                        m_t_data = DF.get(vrbl,utc=tidx,level=level,lons=lons,lats=lats)[0,...]
                    except:
                        print("Fail for ensemble member",nm,mem,"\nAlso,",self.initutc)
                        raise

                if ens_no == 1:
                    nz,nlats,nlons = m_t_data.shape
                    nt = len(fts)
                    nens = len(members)
                    all_ens_data = N.zeros((nens,nt,nz,nlats,nlons))

                all_ens_data[ens_no-1,tn,:,:,:] = m_t_data

        # pdb.set_trace()
        return all_ens_data

    def accumulated(self,vrbl='RAINNC',itime=None,ftime=None,level=False,
                    Nlim=False,
                    Elim=False,Slim=False,Wlim=False,inclusive=False,
                    lons=None,lats=None,dom=1,fcsttime=None,accum_hr=1,
                    members=None):
        """Accumulate, for every ensemble member, at each grid point,
        the variable specified. Usually precipitation.

        Todos:
            * Logic to work out if values are for each history output
                timestep, from the start of the simulation, from the
                start of the data file...

        """
        if (itime is False) and (ftime is False) and (fcsttime is not None):
            ftime = fcsttime
            itime = ftime - datetime.timedelta(seconds=3600*accum_hr)
        elif itime==0:
            itime = self.itime
        elif ftime==-1:
            ftime = self.ftime
        print("Accumulation computed from {} to {}".format(itime,ftime))

        if vrbl is 'RAINNC':
            itime_rainnc = self.get('RAINNC',fcsttime=itime,
                    lons=lons,lats=lats,dom=dom,members=members)
            ftime_rainnc = self.get('RAINNC',fcsttime=ftime,
                    lons=lons,lats=lats,dom=dom,members=members)
            accum = ftime_rainnc - itime_rainnc
        else:
            all_ens_data = self.get(vrbl,itime=itime,ftime=ftime,
                                        inclusive=inclusive,lats=lats,lons=lons,
                                        dom=dom,members=members)
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


    def get_corners(self,dom,chop_inside=False):
        """ Return the corners of a given domain.
        
        Use to find mutual data coverage?
        """
        return W.get_corners(chop_inside=chop_inside)
        
    def find_mutual_geocoverage(self):
        """Find all unique lat/lon grids. Then create a grid that covers the 
        mutual area. Then interpolate each member to that grid during load.
        """
        return

    def __str__(self):
        return ("This ensemble contains {} members, was initialised"
                " at {}, and has {} domains. Root directory is {}".format(
                    self.nmems,self.initutc,self.ndoms,self.rootdir))

    def __iter__(self):
        self.pick_member_idx = -1
        return self

    def __next__(self):
        self.pick_member_idx += 1
        if self.pick_member_idx > self.nmems-1:
            raise StopIteration
        mem = self.member_names[self.pick_member_idx]
        # return mem, self.members[mem]
        return mem
        
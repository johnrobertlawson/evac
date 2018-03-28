"""Compute or load data from netCDF file.

Dimensions of 4D variable X are X.dimensions:
(Time,bottom_top,south_north,west_east_stag)
Time, levels, latitude, longitude

Todos:
    * Move computation of derived variables to derive/ folder!
    * Make method assign the returns to instance instead of returning.
"""

from netCDF4 import Dataset
import sys
import os
import numpy as N
import calendar
import pdb
import scipy.ndimage
import collections
import scipy.interpolate
import datetime

import evac.utils as utils
import evac.utils.met_constants as mc
from WEM.postWRF.postWRF.ncfile import NC
import evac.derived.derived as derived

debug_get = 0

class WRFOut(NC):

    """
    An instance of WRFOut contains all the methods that are used to
    access and process netCDF data.
    """
    def __init__(self,fpath,fmt='em_real',ncks=False):
        """
        Initialisation fetches and computes basic user-friendly
        variables that are most often accessed.

        Args:
            fpath    :    absolute path to netCDF4 (wrfout) file
            fmt (str,optional)      : em_real or ideal or met_em
            ncks (bool,optional)    : whether WRFOut file has been stripped to a few variables.
                                        Hence check for KeyErrors for variables

        """
        super().__init__(fpath)

        self.fmt = fmt
        self.fpath = fpath
        self.ncks = ncks

        self.dx = self.nc.DX
        self.dy = self.nc.DY
        self.get_dimensions(fmt)

        # if ncks:
        try:
            self.wrf_times = self.nc.variables['Times'][:]
        except KeyError:
            self.wrf_times = N.arange(self.t_dim)
        else:
            # Get times in nicer format
            self.utcs = self.wrftime_to_datenum()
            if len(self.utcs) == 1:
                self.dt = None
            else:
                self.dt = self.utcs[1]-self.utcs[0]
                # Assuming constant time steps!

        if (ncks is False) and (fmt is 'em_real'):
            self.P_top = self.nc.variables['P_TOP'][0]

        # Loads variable lists
        self.fields = list(self.nc.variables.keys())
        self.computed_fields = list(self.return_tbl().keys())
        self.available_vrbls = self.fields + self.computed_fields

        # Now do specific loads for idealised or real runs
        if fmt is "em_real":
            self.em_real_init()

        elif fmt is "ideal":
            self.ideal_init()
        elif fmt is "met_em":
            self.ideal_init()

        self.Nlim, self.Elim, self.Slim, self.Wlim = self.get_limits()

    def get_dimensions(self,fmt='em_real'):
        self.t_dim = len(self.nc.dimensions['Time'])
        self.x_dim = len(self.nc.dimensions['west_east'])
        self.y_dim = len(self.nc.dimensions['south_north'])

        self.timekey = 'Time'
        self.lonkey = 'west'
        self.latkey = 'north'

        if fmt == 'met_em':
            self.z_dim = len(self.nc.dimensions['num_metgrid_levels'])
            self.lvkey = 'num_metgrid'
        else:
            self.z_dim = len(self.nc.dimensions['bottom_top'])
            self.lvkey = 'bottom'

        return

    def em_real_init(self):
        """
        Grab the 2D arrays of the wrfout file.

        self.lats and self.lons represent the 2D arrays.

        self.lats1D and self.lons1D represent the 1D arrays intersecting
        the middle of the domain. They shouldn't be used for plotting.

        """
        self.lats = self.nc.variables['XLAT'][0,...] # Might fail if only one time?
        self.lons = self.nc.variables['XLONG'][0,...]

        self.lats1D = self.lats[:,int(len(self.lats)/2)]
        self.lons1D = self.lons[int(len(self.lons)/2),:]

        self.cen_lat = float(self.nc.CEN_LAT)
        self.cen_lon = float(self.nc.CEN_LON)
        self.truelat1 = float(self.nc.TRUELAT1)
        self.truelat2 = float(self.nc.TRUELAT2)


    def ideal_init(self):
        pass


    def wrftime_to_datenum(self,fmt='timegm'):
        """
        Convert wrf's weird Times variable to datenum or datetime.

        """
        times = self.wrf_times
        wrf_times_epoch = N.zeros([times.shape[0]])

        for n,t in enumerate(times):
            tstr = ''.join(t.astype(str))

            yr = int(tstr[0:4])
            mth = int(tstr[5:7])
            day = int(tstr[8:10])
            hr = int(tstr[11:13])
            mins = int(tstr[14:16])
            sec = int(tstr[17:19])

            if fmt == 'timegm':
                wrf_times_epoch[n] = calendar.timegm([yr,mth,day,hr,mins,sec])
            elif fmt == 'datetime':
                wrf_time_epoch[n] = datetime.datetime(yr,mth,day,hr,mins,sec)
        return wrf_times_epoch


    def get_time_idx(self,utcs):
        """
        Get closest index to desired time

        Args:
            utcs (list,tuple, int)    :    times
        Returns:
            tidx (int): closest index to desired time

        """
        # import pdb; pdb.set_trace()
        dn = utils.ensure_datenum(utc)
        dns = utils.get_sequence(dn)
        tidx = []
        for t in dns:
            tidx.append(utils.closest(self.utc,t))
        return N.array(tidx)


    def check_compute(self,vrbl):
        """This method returns the required variables
        that need to be loaded from the netCDF file.

        Args:
            vrbl    :   WRF variable desired
        Returns:
            bool: True if variable exists in wrfout file.
                        False if the variable needs computing.
        """

        if vrbl in self.fields:
            return True
        else:
            return False

    def return_tidx_range(self,utc0,utc1):
        """
        Give a start and end time. Returns an array of
        all indices. Useful for self.get() to return an
        array of data with all times between utc0 and utc1.
        """
        idx0 = self.get_time_idx(utc0)
        idx1 = self.get_time_idx(utc1)
        return N.arange(idx0,idx1)


    def get(self,vrbl,utc=None,level=None,lats=None,lons=None,
                smooth=1,other=False):
        """
        Get data.

        Will interpolate onto pressure, height coordinates if needed.
        Will smooth if needed.

        Level:
        * indices: integer or N.ndarray of integers
        * pressure: string ending in 'hPa'
        * height: string ending in 'm' or 'km'
        * isentropic: string ending in 'K'

        Lats:
        * indices: integer or N.ndarray of integers
        * lats: float or N.ndarray of floats

        Lons:
        * indices: integer or N.ndarray of integers
        * lons: float or N.ndarray of floats

        Args:
        vrbl    :        WRF or computed variable required
        utc     :        indices (<500): integer/N.ndarray of integers.
                            time tuple: 6-item tuple or list/tuple of these.
                            datenum: integer >500 or list of them.
        """
        # import pdb; pdb.set_trace()
        # Time

        if utc is None:
            tidx = None
        elif isinstance(utc,int) and (utc<500):
        # elif isinstance(utc,(int,N.int64)) and utc<500:
            if utc < 0:
                # Logic to allow negative indices
                tidx = self.t_dim + utc
            else:
                tidx = utc
        elif isinstance(utc,(list,tuple,int,datetime.datetime)): # and len(utc[0])==6:
            tidx = self.get_time_idx(utc)
        elif isinstance(utc,N.ndarray): #and isinstance(utc[0],int):
            tidx = utc
        else:
            print("Invalid time selection.")
            raise Exception

        # Level
        # if not level:
        coords = utils.check_vertical_coordinate(level)
        # import pdb; pdb.set_trace()
        if (level is None) or (coords is 'eta'):
            lvidx = None
        elif coords == 'index':
            lvidx = level
        elif isinstance(coords,str):
            if coords == 'surface':
                lvidx = 0
            else:
                lvidx = coords
        else:
            print("Invalid level selection.")
            raise Exception

        # Lat/lon
        if not type(lats)==type(lons):
            # What about case where all lats with one lon?
            raise Exception
        if lats is None:
            lonidx = None
            latidx = None
        elif isinstance(lons,(list,tuple,N.ndarray)):
            if isinstance(lons[0],int):
                lonidx = lons
                latidx = lats
            elif isinstance(lons[0],float):
                # Interpolate to lat/lon
                lonidx = None
                latidx = None
        elif isinstance(lons,(int,N.int64)):
            lonidx = lons
            latidx = lats
        elif isinstance(lons,float):
            # Interpolate to lat/lon
            lonidx = utils.closest(self.lons1D,lons)
            latidx = utils.closest(self.lats1D,lats)
        else:
            print("Invalid lat/lon selection.")
            raise Exception
        # Check if computing required
        # When data is loaded from nc, it is destaggered

        if debug_get:
            print(("Computing {0} for level {1} of index {2}".format(vrbl,level,lvidx)))

        if self.check_compute(vrbl):
            if debug_get:
                print(("Variable {0} exists in dataset.".format(vrbl)))
            if lvidx is 'isobaric':
                data = self.get_p(vrbl,tidx,level,lonidx,latidx)
            elif isinstance(lvidx,(tuple,list,N.ndarray,int,type(None))):
                data = self.load(vrbl,tidx,lvidx,lonidx,latidx)
            else:
                raise Exception
        else:
            if debug_get:
                print(("Variable {0} needs to be computed.".format(vrbl)))
            if lvidx is 'isobaric':
                # data = self.get_p(vrbl,tidx,level,lonidx, latidx)[N.newaxis,N.newaxis,:,:]
                data = self.compute(vrbl,tidx,level,lonidx,latidx,other)
            else:
                data = self.compute(vrbl,tidx,lvidx,lonidx,latidx,other)

        # if len(data.shape) == 2:
            # data = data[N.newaxis,N.newaxis,:,:]
        # elif len(data.shape) == 3:
            # data = data[N.newaxis,:,:,:]
        # if len(data.shape) == 3:
            # data = N.expand_dims(data,axis=0)
        # import pdb; pdb.set_trace()
        data = self.make_4D(data,vrbl=vrbl)

        return data

    def load(self,vrbl,tidx,lvidx,lonidx,latidx):
        """
        Fetch netCDF data for a given variable, for given time, level,
        latitude, and longitude indices.

        Args:
            vrbl        :   WRF variable
            tidx        :   time index. False fetches all.
            lvidx       :   level index. False fetches all

        Todos:
            * Get rid of integer arguments earlier in the method chain, and
                make them single-element numpy arrays.
        """
        # import pdb; pdb.set_trace()
        # First, check dimension that is staggered (if any)
        destag_dim = self.check_destagger(vrbl)

        # Next, fetch dimension names
        dim_names = self.get_dims(vrbl)

        vrbldata = self.nc.variables[vrbl]
        sl = self.create_slice(vrbl,tidx,lvidx,lonidx,latidx,dim_names)
        # If that dimension has a slice of indices, it doesn't need staggering.
        if destag_dim and isinstance(sl[destag_dim],N.ndarray):
            destag_dim = None

        # import pdb; pdb.set_trace()
        data = self.destagger(vrbldata[sl],destag_dim)
        return data


    def create_slice(self,vrbl,tidx,lvidx,lonidx,latidx,dim_names):
        """
        Create slices from indices of level, time, lat, lon.
        False mean pick all indices.
        """
        # See which dimensions are present in netCDF file variable
        sl = []
        # if vrbl.startswith('RAINNC'):
            # pdb.set_trace()
        if any(self.timekey in p for p in dim_names):
            if tidx is None:
                sl.append(slice(None,None))
            elif isinstance(tidx,slice) or isinstance(tidx,N.ndarray):
                sl.append(tidx)
            else:
                sl.append(slice(tidx,tidx+1))

        if any(self.lvkey in p for p in dim_names):
            if lvidx is None:
                sl.append(slice(None,None))
            elif isinstance(lvidx,int):
                sl.append(slice(lvidx,lvidx+1))
            elif isinstance(lvidx,N.ndarray):
                sl.append(lvidx)
            else:
                sl.append(slice(None,None))

        if any(self.lonkey in p for p in dim_names):
            if lonidx is None:
                sl.append(slice(None,None))
            elif isinstance(lonidx,slice) or isinstance(lonidx,N.ndarray):
                sl.append(lonidx)
            elif isinstance(lonidx,(int,N.int64)):
                sl.append(slice(lonidx,lonidx+1))
            else:
                sl.append(slice(None,None))

        if any(self.latkey in p for p in dim_names):
            if latidx is None:
                sl.append(slice(None,None))
            elif isinstance(latidx,slice) or isinstance(latidx,N.ndarray):
                sl.append(latidx)
            elif isinstance(latidx,(int,N.int64)):
                sl.append(slice(latidx,latidx+1))
            else:
                sl.append(slice(None,None))

        return sl

    def check_destagger(self,var):
        """ Looks up dimensions of netCDF file without loading data.

        Returns dimension number that requires destaggering

        """
        stag_dim = None
        for n,dname in enumerate(self.nc.variables[var].dimensions):
            if 'stag' in dname:
                stag_dim = n

        return stag_dim

    def get_dims(self,var):
        dims = self.nc.variables[var].dimensions
        return dims

    def destagger(self,data,ax):
        """ Destagger data which needs it doing.

        Theta always has unstaggered points in all three spatial dimensions (axes=1,2,3).

        Data should be 4D but just the slice required to reduce unnecessary computation time.

        Don't destagger in x/y for columns

        Args:
            data    :   numpy array of data requiring destaggering
            ax      :   axis requiring destaggering


        """
        # Check for dimensions of 1.
        # If it exists, don't destagger it.

        shp = data.shape
        for n,size in enumerate(shp):
            if (size==1) and (n==ax):
                ax = None
                break

        #pdb.set_trace()

        if ax==None:
            return data
        else:
            nd = data.ndim
            sl0 = []     # Slices to take place on staggered axis
            sl1 = []

            for n in range(nd):
                if n is not ax:
                    sl0.append(slice(None))
                    sl1.append(slice(None))
                else:
                    sl0.append(slice(None,-1))
                    sl1.append(slice(1,None))

            data_unstag = 0.5*(data[sl0] + data[sl1])
            return data_unstag

    def return_tbl(self):
        """
        Returns a dictionary to look up method for computing a variable
        Todos:
            * merge with derived/_derived.py
        """
        tbl = {}
        tbl['shear'] = derived.compute_shear
        tbl['thetae'] = derived.compute_thetae
        tbl['cref'] = derived.compute_comp_ref
        tbl['wind10'] = derived.compute_wind10
        tbl['wind'] = derived.compute_wind
        tbl['CAPE'] = derived.compute_CAPE
        tbl['Td'] = derived.compute_Td
        tbl['pressure'] = derived.compute_pressure
        tbl['drybulb'] = derived.compute_drybulb
        tbl['theta'] = derived.compute_theta
        tbl['geopot'] = derived.compute_geopotential
        tbl['Z'] = derived.compute_geopotential_height
        tbl['dptp'] = derived.compute_dptp #density potential temperature pert.
        tbl['T2p'] = derived.compute_T2_pertub
        tbl['dpt'] = derived.compute_dpt #density potential temperature .
        tbl['buoyancy'] = derived.compute_buoyancy
        tbl['strongestwind'] = derived.compute_strongest_wind
        tbl['PMSL'] = derived.compute_pmsl
        tbl['RH'] = derived.compute_RH
        tbl['dryairmass'] = derived.compute_dryairmass
        tbl['QTOTAL'] = derived.compute_qtotal
        tbl['olr'] = derived.compute_olr
        tbl['es'] = derived.compute_satvappres
        tbl['e'] = derived.compute_vappres
        tbl['q'] = derived.compute_spechum
        tbl['fluidtrapping'] = derived.compute_fluid_trapping_diagnostic
        tbl['lyapunov'] = derived.compute_instantaneous_local_Lyapunov
        tbl['REFL_comp'] = derived.compute_REFL_comp
        tbl['temp_advection'] = derived.compute_temp_advection
        tbl['omega'] = derived.compute_omega
        tbl['density'] = derived.compute_density
        # tbl['accum_precip'] = derived.compute_accum_rain
        tbl['PMSL_gradient'] = derived.compute_PMSL_gradient
        tbl['T2_gradient'] = derived.compute_T2_gradient
        tbl['Q_pert'] = derived.compute_Q_pert
        tbl['vorticity'] = derived.return_vorticity

        return tbl

    # def compute(self,vrbl,tidx,lvidx,lonidx,latidx,other,lookup=0):
    def compute(self,vrbl,tidx,lvidx,lonidx,latidx,other=False,lookup=False):
        """ Look up method needed to return array of data
        for required variable.

        Keyword arguments include settings for computation
        e.g. top and bottom of shear computation

        Args:
            vrbl        :   variable name
            tidx        :   time index/indices
            lookup      :   enables a check to see if something can be
                                computed. Returns true or false.
        """
        tbl = self.return_tbl()
        if lookup:
            response = lookup in tbl
        else:
            # response = tbl[vrbl](tidx,lvidx,lonidx,latidx,other)
            # Edited to pass in self to derived functions.
            response = tbl[vrbl](self,tidx,lvidx,lonidx,latidx,other)
        return response




    def compute_ave(self,va,z1,z2):
        """
        Compute average values for variable in layer.

        Args:
            va      :   variable
            z1      :   height at bottom
            z2      :   height at top

        Returns:
            data    :   the averaged variable

        Todos:
            * Move or remove?
        """

        # Get coordinate system
        vc = self.check_vcs(z1,z2)


    def make_4D(self,datain,vrbl=False,missing_axis=False):
        """
        If vrbl, look up the wrfout file's variable dimensions and
        adjust accordingly to get into 4D structure. If not vrbl,
        for instance a computed variable, the user needs to specify
        which axes are missing in a tuple.
        """
        dataout = datain
        if len(datain.shape)<4:
            if vrbl in self.fields:
                dims = self.nc.variables[vrbl].dimensions
                missing = self.get_missing_axes(dims)
                for ax in missing:
                    dataout = N.expand_dims(dataout,axis=ax)
            else:
                while len(dataout.shape)<4:
                    dataout = N.expand_dims(dataout,axis=0)
        # import pdb; pdb.set_trace()
        return dataout

    def get_missing_axes(self,dims):
        axes = {0:"Time",1:"bottom",2:"south",3:"west"}
        missing = []

        for ax, axname in axes.items():
            present = bool([True for d in dims if axname in d])
            if not present:
                missing.append(ax)

        return missing

    def check_vcs(self,z1,z2,exception=1):
        """
        Check the vertical coordinate systems

        If identical, return the system
        If not, raise an exception.
        """

        vc = utils.level_type(z1)
        vc = utils.level_type(z2)

        if vc1 != vc2:
            print("Vertical coordinate systems not identical.")
            return False
            if exception:
                raise Exception
        else:
            return vc1

    def get_XY(self,lat,lon):
        """Return grid indices for lat/lon pair.
        """
        pass

    def get_limited_domain(self,da,skip=1,return_array='idx'):
        """
        Return smaller array of lats, lons depending on
        input dictionary of N, E, S, W limits.

        Args:
            skip            :   for creating thinned domains
            return_type     :   if idx, return array of idx
                                    if slice, return slice only
                                    if latlon, return lat/lon values
        """

        if isinstance(da,dict):
            N_idx = self.get_lat_idx(da['Nlim'])
            E_idx = self.get_lon_idx(da['Elim'])
            S_idx = self.get_lat_idx(da['Slim'])
            W_idx = self.get_lon_idx(da['Wlim'])
        else:
            N_idx = self.lats1D.shape[0]
            E_idx = self.lons1D.shape[0]
            S_idx = 0
            W_idx = 0

        if return_array=='latlon' or return_array=='slice':
            if isinstance(da,dict):
                lat_sl = slice(S_idx,N_idx,skip)
                lon_sl = slice(W_idx,E_idx,skip)
            else:
                lat_sl = slice(None,None,skip)
                lon_sl = slice(None,None,skip)

            if return_array == 'latlon':
                return self.lats1D[lat_sl], self.lons1D[lon_sl]
            else:
                return lat_sl, lon_sl
        elif return_array=='idx':
            return N.arange(S_idx,N_idx,skip), N.arange(W_idx,E_idx,skip)
        else:
            print("Invalid selection for return_array.")
            raise Exception

    def get_latlon_idx(self,lat,lon):
        latidx, lonidx = utils.get_latlon_idx(self.lats,self.lons,
                                    lat,lon)
        # coords = N.unravel_index(N.argmin((lat-self.lats)**2+
                    # (lon-self.lons)**2),self.lons.shape)
        # lon, lat
        # return [int(c) for c in coords]
        return latidx, lonidx

    def get_lat_idx(self,lat):
        lat_idx = N.where(abs(self.lats-lat) == abs(self.lats-lat).min())[0][0]
        return int(lat_idx)

    def get_lon_idx(self,lon):
        lon_idx = N.where(abs(self.lons-lon) == abs(self.lons-lon).min())[0][0]
        return int(lon_idx)

    def get_p(self,vrbl,tidx=None,level=None,lonidx=None,latidx=None):
        """
        Return an pressure level isosurface of given variable.
        Interpolation is linear so watch out.

        Dimensions returns as (height,lat,lon)
        Or is it (height,lon, lat!?)

        Todos:
            * Need to include limited domain functionality

        if vrbl=='pressure',create constant grid.
        """
        # print('GET_P:',vrbl,tidx,level)
        if isinstance(level,str) and level.endswith('hPa'):
            hPa = 100.0*int(level.split('h')[0])
            nlv = 1
        elif isinstance(level,(float,int)):
            hPa = level*100.0
            nlv = 1
        elif isinstance(level,(tuple,list)):
            hPa = [l*100.0 for l in level]
            nlv = len(hPa)
        else:
            print("Use XXXhPa, an integer, or list of integers for level.")
            raise Exception

        # If this breaks, user is requesting non-4D data
        # Duck-typing for the win

        P = self.get('pressure',utc=tidx,lons=lonidx,lats=latidx)[...]

        if vrbl=='pressure':
            dataout = N.ones([P.shape[0],nlv,P.shape[-2],P.shape[-1]])*hPa
        else:
            datain = self.get(vrbl,utc=tidx,lons=lonidx,lats=latidx)[...]
            # import pdb; pdb.set_trace()
            # What about RUC, pressure coords
            # dataout = N.zeros([nlv,P.shape[-1],P.shape[-2]])
            dataout = N.zeros([P.shape[0],nlv,P.shape[-2],P.shape[-1]])
            # pdb.set_trace()
            for t in range(P.shape[0]):
                for (i,j), p in N.ndenumerate(dataout[0,0,:,:]):
                    dataout[t,:,i,j] = N.interp(hPa,P[t,:,i,j][::-1],datain[t,:,i,j][::-1])
            # dataout = scipy.interpolate.griddata(P.flatten(),datain.flatten(),hPa)
        # if nlv == 1:
            # Return 2D if only one level requested
            # return self.make_4D(dataout[:,0,:,:])
        # else:
            # return self.make_4D(dataout)
        return dataout

    def interp_to_p_fortran(self,config,nc_path,var,lv):
        """ Uses p_interp fortran code to put data onto a pressure
        level specified.

        Args:
            config  :   contains directory of p_interp files
            nc_path :   path to original netCDF file data
            var     :   variable(s) to compute
            lv      :   pressure level(s) to compute

        Returns:
            fpath   :   path to new netCDF file with p co-ords
        """
        # Fetch paths
        p_interp_path = os.path.join(
                            config.p_interp_root,'p_interp')
        namelist_path = os.path.join(
                            config.p_interp_root,'namelist.pinterp')
        nc_root, nc_fname = os.path.split(nc_path)# Root directory of wrfout file
        output_root = nc_root # Directory to dump output file (same)

        """
        Can we add a suffix to the new netCDF file?
        Check to see if file already exists
        """
        # Copy old p_interp for backup
        command1 = ' '.join(('cp',p_interp_path,p_interp_path+'.bkup'))
        os.system(command1)

        # Edit p_interp's namelist settings
        edit_namelist(path_to_interp,'path_to_input',nc_root,col=18)
        edit_namelist(path_to_interp,'input_name',nc_fname,col=18)
        edit_namelist(path_to_interp,'path_to_output',output_root,col=18)
        edit_namelist(path_to_interp,'process','list',col=18)
        edit_namelist(path_to_interp,'fields',var,col=18)
        edit_namelist(path_to_interp,'met_em_output','.FALSE.',col=18)
        edit_namelist(path_to_interp,'fields',var,col=18)

        command2 = os.path.join('./',p_interp_path)
        os.system(command2) # This should execute the script

        return fpath


    def __edit_namelist(self,fpath,old,new,incolumn=1,col=23):
        """col=23 is default for wps namelists.
        Todos:
            * remove and put into lazy/
        """
        flines = open(fpath,'r').readlines()
        for idx, line in enumerate(flines):
            if old in line:
                # Prefix for soil intermediate data filename
                if incolumn==1:
                    flines[idx] = flines[idx][:col] + new + " \n"
                else:
                    flines[idx] = ' ' + old + ' = ' + new + "\n"
                nameout = open(fpath,'w')
                nameout.writelines(flines)
                nameout.close()
                break

    def get_limits(self):
        Nlim = float(self.lats1D[-1])
        Elim = float(self.lons1D[-1])
        Slim = float(self.lats1D[0])
        Wlim = float(self.lons1D[0])
        return Nlim, Elim, Slim, Wlim

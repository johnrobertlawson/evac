""" A collection of geographical-related functions.

Todo:
    * Make sure other utils functions are imported explicitly, so
        we can remove the star import in util's `__init__.py`.
"""
#import scipy.ndimage as nd
import calendar
import collections
import fnmatch
import math
import os
import pdb
import itertools
import sys
import time
import glob
import pickle
import datetime
import heapq
import random
import base64
from pathlib import Path, PosixPath
import pytz

# import cartopy.crs as ccrs
import xarray
import matplotlib as M
import numpy as N
from netCDF4 import Dataset

def decompose_wind(wspd,wdir,wdir_fmt='deg'):
    """Turn a wind speed and direction into u and v components.

    Uses meteorological convention - so a westerly wind is 270.

    TODO: Add option to use bearings (westerly wind is 90) or mathematical
    (i.e., zero begins at 3 o'clock)

    Arguments:
        wspd        :   Wind speed. Can be float, integer, or N.ndarray
        wdir        :   Wind direction. Can be float, integer, or N.ndarray.
    Optional:
        wdir_fmt    :   wdir is degrees by default ('deg'). Use 'rad' to use radians.
    """
    #if (type(wspd) == N.array) & (type(wdir) == N.array):
    #    uwind = N.array([-s * N.sin(N.radians(d)) if ((s>-1)&(d>-1)) else -9999
    #                for s,d in zip(wspd,wdir)])
    #    vwind = N.array([-s * N.cos(N.radians(d)) if ((s>-1)&(d>-1)) else -9999
    #                for s,d in zip(wspd,wdir)])
    #else:

    if wdir_fmt == 'deg':
        wdir = N.radians(wdir)
    elif wdir_fmt == 'rad':
        pass
    else:
        raise Exception("Choose 'deg' or 'rad' for wdir_fmt.")

    u = -wspd * N.sin(wdir)
    v = -wspd * N.cos(wdir)

    return u,v

def convert_velocity_units(V,conversion):
    """Convert a velocity from one unit to another.

    Args:
        V       :   a number (int, float, array)
        convert :   'ms_kt' converts the output from m/s to knots.
                    'ms_mph' converts from m/s to mph.
                    'kt_ms' converts from knots to m/s.
                    'kt_mph' converts from knots to mph.
                    'mph_kt' converts from mph to knots.
                    'mph_ms' converts from mph to m/s.
    Todo:
        * add km/h
    """
    if conversion == 'ms_kt':
        factor = 1.94384449
    elif conversion == 'ms_mph':
        factor = 2.23694
    elif conversion == 'kt_ms':
        factor = 0.51444444
    elif conversion == 'kt_mph':
        factor = 1.15078
    elif conversion == 'mph_ms':
        factor = 0.44704
    elif conversion == 'mph_kt':
        factor = 0.868976
    elif conversion in (False,None):
        factor = 1.0
    else:
        raise Exception("Conversion argument {} is not recognised.".format(
                        conversion))

    return V * factor

def combine_wind_components(u,v):
    wdir = N.degrees(N.arctan2(u,v)) + 180.0
    wspd = N.sqrt(u**2 + v**2)
    return wspd, wdir

def dewpoint(T,RH): # Lawrence 2005 BAMS?
    #T in C
    #RH in 0-100 format
    es = 6.11 * (10**((7.5*T)/(237.7+T)))
    e = es * RH/100.0
    alog = 0.43429*N.log(e) - 0.43429*N.log(6.11)
    Td = (237.7 * alog)/(7.5-alog)
    #pdb.set_trace()
    return Td

def csvprocess(data,names,convert=0):
    # Get stations
    stationlist = N.unique(data['stid']) # Stations in the record
    stnum = len(stationlist) # Number of them


    # Create dictionary from data
    D = {} # Initialise dictionary of data
    for s in stationlist:
        D[s] = {} # Initialise dictionary of station/s
        print('Loading station data for ' + s)
        those = N.where(data['stid']==s) # Indices of obs
        obNum = len(those[0]) # Number of obs

        # Entries from MesoWest data
        for n in names:
            D[s][n] = data[n][those]

        # Sort times
        year = [int(t[0:4]) for t in D[s]['tutc']]
        month = [int(t[4:6]) for t in D[s]['tutc']]
        day = [int(t[6:8]) for t in D[s]['tutc']]
        hour = [int(t[9:11]) for t in D[s]['tutc']]
        minute = [int(t[11:13]) for t in D[s]['tutc']]

        # Create python time
        D[s]['pytime'] = N.empty(obNum)
        for i in range(0,obNum-1):
            D[s]['pytime'][i] = cal.timegm([year[i],month[i],day[i],hour[i],minute[i],0])

        # Convert to S.I. units
        if convert == 1:
            D[s]['tmpc'] = (5/9.0) * (D[s]['tmpf']-32)
            D[s]['dptc'] = (5/9.0) * (D[s]['dptf']-32)
            D[s]['wsms'] = D[s]['wskn']*0.51444
            D[s]['wgms'] = D[s]['wgkn']*0.51444

        ### DATA PROCESSING
        # Find time of biggest wind increase (DSWS initiation)

        D[s]['dVdt'] = N.zeros(obNum) # This is our array of d(wind)/d(time)
        for i in range(0,obNum-1):
            if i == 0:
                pass # First record is zero as there's no gradient yet
            else:
                D[s]['dVdt'][i] = ((D[s]['wsms'][i] - D[s]['wsms'][i-1])/
                                    (D[s]['pytime'][i] - D[s]['pytime'][i-1]))
        if any(D[s]['dVdt']) == False:
            # Data were absent (-9999) or rubbish (0)
            D[s]['mwi'] = -9999
            D[s]['mwi_t'] = -9999
        else:
            D[s]['mwi'] = D[s]['dVdt'].max() # max wind increase
            loc = N.where(D[s]['dVdt']==D[s]['mwi'])
            D[s]['mwi_t'] = D[s]['pytime'][loc][0] # time of max wind increase

        # Find time of maximum wind gust
        if any(D[s]['wgms']) == False:
            D[s]['mg'] = -9999
            D[s]['mg_t'] = -9999
        else:
            D[s]['mg'] = D[s]['wgms'].max() # Maximum gust at the station
            loc = N.where(D[s]['wgms']==D[s]['mg'])
            D[s]['mg_t'] = D[s]['pytime'][loc][0] # time of max gust

        # Find time of maximum wind speed
        if any(D[s]['wsms']) == False:
            D[s]['ms'] = -9999
            D[s]['ms_t'] = -9999
        else:
            D[s]['ms'] = D[s]['wsms'].max() # Maximum gust at the station
            loc = N.where(D[s]['wsms']==D[s]['ms'])
            D[s]['ms_t'] = D[s]['pytime'][loc][0] # time of max gust


        # Frequency of observation in minutes
        try:
            D[s]['dt'] = (D[s]['pytime'][1] - D[s]['pytime'][0]) / 60
        except IndexError:
            D[s]['dt'] = -9999

        # Find lowest pressures
        try:
            D[s]['lowp'] = D[s]['alti'].min()
        except IndexError:
            D[s]['lowp'] = -9999
            D[s]['lowp_t'] = -9999
        finally:
            if D[s]['lowp'] != -9999:
                D[s]['lowp_t'] = N.where(D[s]['alti']==D[s]['lowp'])
            else:
                pass

    return D


# Constants
# TO DO: Move this to constants file
rE = 6378100 # radius of Earth in metres


def get_cross_section(Alat, Alon, Blat, Blon):
    # Now find cross-sections
    Ax, Ay = gridded_data.getXY(lats,lons,Alat,Alon)
    Bx, By = gridded_data.getXY(lats,lons,Blat,Blon)

    # Number of points along cross-section
    xspt = int(N.hypot(Bx-Ax,By-Ay))
    xx = N.linspace(Ay,By,xspt).astype(int)
    yy = N.linspace(Ax,Bx,xspt).astype(int)

    # Get terrain heights along this transect
    heights = topodata[xx,yy]

    # Work out distance along this cross-section in m
    xsdistance = gridded_data.xs_distance(Alat,Alon,Blat,Blon)

    xvals = N.linspace(0,xsdistance,xspt)
    xlabels = ['%3.1f' %(x/1000) for x in xvals]
    # Now plot cross-sections
    fig = plt.figure(figsize=(width,height))
    plt.plot(xvals,heights)
    delta = xspt/10
    plt.xticks(xvals[::delta],xlabels[::delta])
    plt.xlabel('Distance along cross-section (km)')
    fname = 'test1.png'
    plt.savefig(outdir+fname,bbox_inches='tight',pad_inches=0.3)
    plt.clf()


def generate_image_loop(dom='d03',fpath='./'):
    print("Creating a pretty loop...")
    os.system('convert -delay 50 '+fpath+dom+'*.png > -loop '+fpath+dom+'_windloop.gif')
    print("Finished!")

def p_interpol(dom,var,ncfolder,datestr):
    datestr2 = datestr[0:4]+'-'+datestr[4:6]+'-'+datestr[6:8]+'_'+datestr[8:10]+':00:00'
    fname = ncfolder+'wrfout_'+dom+'_'+datestr2
    nc = Dataset(fname,'r')
    var_data = nc.variables[var][:]
    if var == 'T': # pert. theta
        var_data += 300 # Add base state
    pert_pressure = nc.variables['P'][:] #This is perturbation pressure
    base_pressure = nc.variables['PB'][:] # This is base pressure
    pressure = pert_pressure + base_pressure # Get into absolute pressure in Pa.
    p_levels = N.arange(10000,100000,1000)
    var_interp = N.zeros((pressure.shape[0],len(p_levels),pressure.shape[2],pressure.shape[3]))
    #pdb.set_trace()
    for (t,y,x),v in N.ndenumerate(var_interp[:,0,:,:]):
        var_interp[t,:,y,x] = N.interp(p_levels,pressure[t,::-1,y,x],var_data[t,::-1,y,x])
    return var_interpi

def hgt_from_sigma(nc):
    vert_coord = (nc.variables['PH'][:] + nc.variables['PHB'][:]) / 9.81
    return vert_coord

# This converts wrf's time to python and human time
def find_time_index(wrftime,reqtimetuple,tupleformat=1):
    # wrftime = WRF Times in array
    # reqtime = desired time in six-tuple

    # Convert required time to Python time if required
    if tupleformat:
        reqtime = calendar.timegm(reqtimetuple)
    else:
        reqtime = reqtimetuple
    nt = wrftime.shape[0]
    pytime = N.zeros([nt,1])
    t = wrftime
    # Now convert WRF time to Python time
    for i in range(nt):
        yr = int(''.join(t[i,0:4]))
        mth = int(''.join(t[i,5:7]))
        day = int(''.join(t[i,8:10]))
        hr = int(''.join(t[i,11:13]))
        min = int(''.join(t[i,14:16]))
        sec = int(''.join(t[i,17:19]))
        pytime[i] = calendar.timegm([yr,mth,day,hr,min,sec])

    # Now find closest WRF time
    timeInd = N.where(abs(pytime-reqtime) == abs(pytime-reqtime).min())[0][0]
    return timeInd

# Find data for certain time, locations, level, variables
# For surface data
def TimeSfcLatLon(nc,varlist,times,latslons='all'):
    # Time (DateTime in string)
    if times == 'all':
        timeInds = list(range(nc.variables['Times'].shape[0]))
    elif len(times)==1: # If only one time is desired
        # Time is in 6-tuple format
        timeInds = find_time_index(nc.variables['Times'],times) # This function is from this module
    elif len(times)==2: # Find all times between A and B
        timeIndA = find_time_index(nc.variables['Times'],times[0])
        timeIndB = find_time_index(nc.variables['Times'],times[1])
        timeInds = list(range(timeIndA,timeIndB))
    # Lat/lon of interest and their grid pointd
    lats = nc.variables['XLAT'][:]
    lons = nc.variables['XLONG'][:]
    if latslons == 'all':
        latInds = list(range(lats.shape[-2]))
        lonInds = list(range(lons.shape[-1]))
    else:
        xmin,ymax = gridded_data.getXY(lats,lons,Nlim,Wlim)
        xmax,ymin = gridded_data.getXY(lats,lons,Slim,Elim)
        latInds = list(range(ymin,ymax))
        lonInds = list(range(xmin,xmax))

    # Return sliced data
    data = {}
    for v in varlist:
        data[v] = nc.variables[v][timeInds,latInds,lonInds]
        # Reshape if only one time
        if len(times)==1:
            data[v] = N.reshape(data[v],(len(latInds),len(lonInds)))
    return data

# Get lat/lon as 1D arrays from WRF
def latlon_1D(nc):
    Nx = nc.getncattr('WEST-EAST_GRID_DIMENSION')-1
    Ny = nc.getncattr('SOUTH-NORTH_GRID_DIMENSION')-1
    lats = nc.variables['XLAT'][0,:,Nx/2]
    lons = nc.variables['XLONG'][0,Ny/2,:]
    return lats, lons

def netcdf_files_in(folder,dom=1,init_time=0,model='auto',return_model=False):
    """ Hunt through given folder to find the right netcdf file for data.
    Return ncpath (absolute path to file) and model (RUC, WRF, etc).


    Args:
        folder      :   Absolute path to directory
        dom         :   specify domain. None specified if zero.
        init_time   :   initialisation time. Can be tuple or datenum.
                        If zero, then folder must contain one unambiguous file.
        model       :   Default: automatically detect the type of netcdf file
                        (RUC data, wrfout file, etc)

    Returns:
        ncpath      :   Absolute path to file

    """
    t = 'auto'
    if init_time:
        t = ensure_timetuple(init_time,fmt='single')

    # Set the model type to load.
    if model=='auto':
        # Get files, check prefix
        files = glob.glob(os.path.join(folder,'*'))
        matches = 0
        model_test = []
        for f in files:
            model_test.append(determine_model(f.split('/')[-1]))
            # print(model_test)
            model_set = set(model_test)
            # import pdb; pdb.set_trace()
            # if model_test:
                # matches += 1
                # model = model_test

        model_set.discard(False)
        matches = len(model_set)
        # import pdb; pdb.set_trace()
        if matches < 1:
            print("No netcdf files found.")
            raise Exception
        elif matches > 1 and isinstance(t,str):
            print("Ambiguous netcdf file selection. Specify model?")
            raise Exception
        else:
            model = list(model_set)[0]

    # import pdb; pdb.set_trace()
    if model=='wrfout':
        pfx = 'wrfout' # Assume the prefix
    elif model=='ruc':
        pfx = getdata.RUC_fname(t,filetype='netcdf')[:7]
        # TODO: logic that considers the four versions of RUC
    else:
        raise Exception

    # import pdb; pdb.set_trace()
    # Pick unambiguous
    if t=='auto':
        # We assume the user has wrfout files in different folders for different times
        f = glob.glob(os.path.join(folder,pfx+'*'))
        # import pdb; pdb.set_trace()
        if len(f) != 1:
            print("Ambiguous netCDF4 selection.")
            raise Exception
        else:
            if return_model:
                return f[0], model
            else:
                return f[0]
    else:
        if (dom > 8):
            print("Domain is out of range. Choose number between 1 and 8 inclusive.")
            raise IndexError

        fname = get_netcdf_naming(model,t,dom)
        f = glob.glob(os.path.join(folder,fname))

        if len(f) == 1:
            if return_model:
                return f[0], model
            else:
                return f[0]
        elif len(f) == 0:
            print("No netCDF4 file found.")
            raise Exception
        else:
            print("Ambiguous netCDF4 selection.")
            raise Exception


def wrfout_files_in(folders,dom=0,init_time='notset',descend=1,avoid=0,
                    unambiguous=0):
    """
    Hunt through given folder(s) to find all occurrences of wrfout
    files.

    Args:
        folders     :   list of absolute paths to directories
        dom         :   specify domain. None specified if zero.
        init_time   :   tuple of initialisation time
        descend     :   boolean: go into subfolders
        avoid       :   string of filenames. if a subfolder contains
                    the string, do not descend into this one.
        unambiguous :   only return a single absolute path, else throw
                    an Exception.

    Returns:
        wrfouts     :   list of absolute paths to wrfout files
    """

    folders = get_sequence(folders)
    avoids = []
    if 'avoid':
        # Avoid folder names with this string
        # or list of strings
        avoid = get_sequence(avoid)
        for a in avoid:
            avoids.append('/{0}/'.format(a))


    w = 'wrfout' # Assume the prefix
    if init_time=='notset':
        suffix = '*0'
        # We assume the user has wrfout files in different folders for different times
    else:
        try:
            it = utils.string_from_time('wrfout',init_time)
        except:
            print("Not a valid wrfout initialisation time; try again.")
            raise Error
        suffix = '*' + it

    if not dom:
    # Presume all domains are desired.
        prefix = w + '_d'
    elif (dom > 8):
        print("Domain is out of range. Choose number between 1 and 8 inclusive.")
        raise IndexError
    else:
        dom = 'd{0:02d}'.format(dom)
        prefix = w + '_' + dom

    wrfouts = []
    if descend:
        for folder in folders:
            for root,dirs,files in os.walk(folder):
                for fname in fnmatch.filter(files,prefix+suffix):
                    skip_me = 0
                    fpath = os.path.join(root,fname)
                    if avoids:
                        for a in avoids:
                            if a in fpath:
                                skip_me = 1
                    else:
                        pass
                    if not skip_me:
                        wrfouts.append(fpath)

    else:
        for folder in folders:
            findfile = os.path.join(folder,prefix+suffix)
            files = glob.glob(findfile)
            # pdb.set_trace()
            for f in files:
                wrfouts.append(os.path.join(folder,f))
    # pdb.set_trace()
    if unambiguous:
        if not len(wrfouts) == 1:
            print(("Found {0} wrfout files.".format(len(wrfouts))))
            raise Exception
        else:
            return wrfouts[0]
    else:
        return wrfouts

def getXY(lats,lons,ptlat,ptlon):
    """
    Output is lat, lon so y,x
    """
    # Find closest lat/lon in array
    minlat = abs(lats-ptlat).min()
    minlon = abs(lons-ptlon).min()
    # Find where these are in the grid
    wherelat = N.where(abs(lats-ptlat) == minlat)
    wherelon = N.where(abs(lons-ptlon) == minlon)
    # pdb.set_trace()
    lat_idx = N.where(lats==lats[wherelat])[0][0]
    lon_idx = N.where(lons==lons[wherelon])[0][0]
    exactlat = lats[wherelat]
    exactlon = lons[wherelon]
    return lat_idx,lon_idx, exactlat, exactlon

def gettopo():
    fname = '/uufs/chpc.utah.edu/common/home/u0737349/dsws/topodata/globe30.bin'
    f = open(fname,'r')
    fdata = N.fromfile(f,dtype='int16')
    # Transposes and reshapes to a lat-lon grid
    # Changes negative values to 0 (sea level)
    xnum = 43200.0
    ynum = 18000.0
    topodata = N.flipud(N.reshape(fdata,(ynum,xnum))).clip(0)
    #topodata = ((N.reshape(fdata,(xnum,ynum))).clip(0))
    f.close(); del fdata
    # Define size of pixels
    xpixel = 360/xnum
    ypixel = 150/ynum # Note only 150 degrees!
    # Create lat/lon grid
    lats = N.arange(-60,90,ypixel)#[::-1]
    lons = N.arange(-180,180,xpixel)#[::-1]
    print('Topographic data has been loaded. Everest is but a mere pixel.')
    return topodata, lats, lons

def xs_distance(Alat, Alon, Blat, Blon):
    phi1 = N.radians(90.0-Alat)
    phi2 = N.radians(90.0-Blat)
    theta1 = N.radians(Alon)
    theta2 = N.radians(Blon)
    arc = math.acos(math.sin(phi1)*math.sin(phi2)*math.cos(theta1-theta2) +
                    math.cos(phi1)*math.cos(phi2))
    xsdistance = rE * arc
    return xsdistance

# This dstacks arrays, unless it's the first time through, in which case it initialises the variable
def dstack_loop(data, Dict, Key):
    # Try evaluating dict[key]. If it doesn't exist, then initialise it
    # If it does exist, stack data
    try:
        Dict[Key]
    except KeyError:
        stack = data
        #Dict[Key] = data
    else:
        stack = N.dstack((Dict[Key],data))
    return stack
    pass

# Create thinned pressure levels for skew T barb plotting
def thinned_barbs(pres):
    levels = N.arange(20000.0,105000.0,5000.0)
    plocs = []
    for l in levels:
        ploc = N.where(abs(pres-l)==(abs(pres-l).min()))[0][0]
        plocs.append(ploc)
    thin_locs = N.array(plocs)
    return thin_locs # Locations of data at thinned levels


def padded_times(timeseq):
    padded = ['{0:04d}'.format(t) for t in timeseq]
    return padded

def string_from_time(usage,t,dom=0,strlen=None,convention=0,**kwargs):
    """ Generate a string from a given time.

    Example:
        Various formats for 2016/3/31/18/0/0::
        
            # usage == 'output':
            201603311800
            # usage == 'title':
            18:00Z on 31/03/2016
            # usage == 'skip':
            (2016, 3, 31, 18, 0, 0)
            # usage == 'wrfout':
            wrfout_d01_2016-03-31_18:00:00
            # usage == 'ruc':
            ruc2_252_{0:04d}{1:02d}{2:02d}_1800_000.nc
            # usage == 'dir':
            2016033118

    Args:
        usage (str): the format to output
        t: date/time, in tuple or datetime.datetime formats
        strlen: Minimum unit to output. "hour","minute" etc
        convention: convention of MM/DD versus DD/MM

    Returns:
        Formatted string.
    """
    t = ensure_timetuple(t)

    if isinstance(t,str):
        if usage == 'output':
            usage = 'skip' # Time is already a string
        elif usage == 'title':
            pass
        #    if kwargs['itime']: # For averages or maxima over time
        #        itime = kwargs['itime']
        #        ftime = kwargs['ftime']
        #    else:
        #        pass
        else:
            raise Exception
    elif isinstance(t,float) or isinstance(t,int):
        # In this case, time is in datenum. Get it into tuple format.
        t = time.gmtime(t)
    else:
        pass

    if usage == 'title':
        # Generates string for titles
        if not 'itime' in kwargs: # i.e. for specific times
        #if not hasattr(kwargs,'itime'): # i.e. for specific times
            strg = '{3:02d}:{4:02d}Z on {2:02d}/{1:02d}/{0:04d}'.format(*t)
        else: # i.e. for ranges (average over time)
            s1 = '{3:02d}:{4:02d}Z to '.format(*kwargs['itime'])
            s2 = '{3:02d}:{4:02d}Z'.format(*kwargs['ftime'])
            strg = s1 + s2
    elif usage == 'wrfout':
        # Generates string for wrfout file finding
        # Needs dom
        if not dom:
            print("No domain specified; using domain #1.")
            dom = 1
        strg = ('wrfout_d{:02d}_{:04d}-{:02d}-{:02d}_{:02d}:{:02d}:{:02d}'.format(
                                            dom,*t))
    elif usage == 'ruc':
        # This depends on the RUC version? Will break?
        strg = ('ruc2_252_{0:04d}{1:02d}{2:02d}_' +
                '{3:02d}{4:02d}_{5:02d}0.nc'.format(*t))
    elif usage == 'output':
        if not convention:
            # No convention set, assume DD/MM (I'm biased)
            convention = 'full'
        # Generates string for output file creation
        if convention == 'DM':
            strg = '{2:02d}{1:02d}_{3:02d}{4:02d}'.format(*t)
        elif convention == 'MD':
            strg = '{1:02d}{2:02d}_{3:02d}{4:02d}'.format(*t)
        elif convention == 'full':
            strg = '{0:04d}{1:02d}{2:02d}{3:02d}{4:02d}'.format(*t)
        else:
            print("Set convention for date format: DM or MD.")
    elif usage == 'dir':
        # Generates string for directory names
        # Needs strlen which sets smallest scope of time for string
        if not strlen:
             print("No timescope strlen set; using hour as minimum.")
             strlen = 'hour'
        n = lookup_time(strlen)
        strg = "{0:04d}".format(t[0]) + ''.join(
                ["{0:02d}".format(a) for a in t[1:n+1]])
    elif usage == 'skip':
        strg = t
    else:
        print("Usage for string not valid.")
        raise Exception
    return strg

def lookup_time(str):
    D = {'year':0, 'month':1, 'day':2, 'hour':3, 'minute':4, 'second':5}
    return D[str]

def get_level_naming(va,lv,**kwargs):
    #lv = kwargs['lv']

    if lv < 1500:
        return str(lv)+'hPa'
    elif lv == 2000:
        return 'sfc'
    elif lv.endswith('K'):
        return lv
    elif lv.endswith('PVU'):
        return lv
    elif lv.endswith('km'):
        return lv
    elif lv == 'all':
        if va == 'shear':
            name = '{0}to{1}'.format(kwargs['bottom'],kwargs['top'])
            return name
        else:
            return 'all_lev'


def check_vertical_coordinate(level):
    """ Check to see what type of level is requested by user.

    """
    if isinstance(level,(str,int,type(None))):
        lv = level
    elif isinstance(level,(list,tuple,N.ndarray)):
        lv = level[0]
    else:
        print(("What have you given me here? Level is"
                "{0}".format(type(level))))
        raise Exception

    # import pdb; pdb.set_trace()
    if isinstance(lv,int):
        if lv<100:
            return 'index'
        else:
            return 'isobaric'
    elif (lv is 'all') or (lv is None):
        return 'eta'

    elif lv.endswith('hPa'):
        # import pdb; pdb.set_trace()
        if lv[:4] == '2000':
            return 'surface'
        elif int(lv.split('h')[0]) < 2000:
            return 'isobaric'
        else:
            print("Pressure is in hPa. Requested value too large.")
            raise Exception

    elif lv.endswith('K'):
        return 'isentropic'

    elif lv.endswith('PVU'):
        return 'PV-surface'

    elif lv.endswith('km'):
        return 'geometric'

    else:
        print('Unknown vertical coordinate.')
        raise Exception

def closest(arr,val):
    """
    Find index of closest value.
    Only working on 1D array right now.

    Inputs:
    val     :   required value
    arr     :   array of values

    Output:

    idx     :   index of closest value

    """
    # pdb.set_trace()
    if isinstance(arr,N.ndarray):
        idx = N.argmin(N.abs(arr - val))
    elif isinstance(arr,(list,tuple)):
        # diffs = [N.argmin(N.abs(a - val)) for a in arr]
        # idx = diffs
        arr2 = N.array(arr)
        idx = N.argmin(N.abs(arr2 - val))
    else:
        raise Exception
    return idx

def closest_datetime(times,t,round=False):
    """Find closest value in list of datetimes.
    Return index of closest.

    Args:
        times (list,tuple)      : collection of datetimes.
        t (datetime.datetime)   : required time
        round (bool,str)        : If False, return closest index only.
            If 'afterinc', return index of first time after t.
            (If closest time is identical to t, return that index)
            If 'afterexc', same, but if closest time = t, return one after.
            If 'beforeinc', return index of last time before t.
            (If closest time is identical to t, return that index)
            If 'beforeexc', same, but if closest time = t, return one before.

    Returns:
        idx (int): Index of times requests
        dtss[idx] (int): Number of seconds difference between the two.
    """
    stimes = N.array(sorted(times))
    dts = stimes-t
    dtss = [(d.days*86400)+d.seconds for d in dts]

    # Closest index
    cidx = N.argmin(N.abs(dtss))

    if round is False:
        idx = cidx
    else:
        if dtss[cidx] == 0:
            bidx_inc = cidx
            aidx_inc = cidx
            bidx_exc = cidx-1
            aidx_exc = cidx+1
        elif times[cidx] < t:
            bidx_inc = cidx
            bidx_exc = cidx
            aidx_inc = cidx+1
            aidx_exc = cidx+1
        else:
            bidx_exc = cidx-1
            bidx_inc = cidx-1
            aidx_exc = cidx
            aidx_nxc = cidx

        if round is 'afterinc':
            idx = aidx_inc
        elif round is 'afterexc':
            idx = aidx_exc
        elif round is 'beforeinc':
            idx = bidx_inc
        elif round is 'beforeexc':
            idx = bidx_exc
        else:
            raise Exception("Enter valid value for round.")

    return idx, dtss[idx]

def dstack_loop(data, obj):
    """
    Tries to stack numpy array (data) into 'stack' object (obj).
    If obj doesn't exist, then initialise it
    If obj does exist, stack data.
    """
    if isinstance(obj,N.ndarray):
        stack = N.dstack((obj,data))
    else:
        stack = data

    return stack

def dstack_loop_v2(data, obj):
    """
    Need to set obj = 0 at start of loop in master script

    Tries to stack numpy array (data) into 'stack' object (obj).
    If obj doesn't exist, then initialise it
    If obj does exist, stack data.
    """
    try:
        print(obj)
    except NameError:
        stack = data
    else:
        stack = N.dstack((obj,data))

    return stack

def vstack_loop(data, obj):
    """
    Need to set obj = 0 at start of loop in master script

    Tries to stack numpy array (data) into 'stack' object (obj).
    If obj doesn't exist, then initialise it
    If obj does exist, stack data.
    """

    if isinstance(obj,N.ndarray):
        stack = N.vstack((obj,data))
    else:
        stack = data

    return stack


def generate_times(idate,fdate,interval,fmt='timetuple',inclusive=False):
    """
    :param itime:       Start date/time. Format is
                        YYYY,MM,DD,HH,MM,SS (calendar.timegm).
    :type itime:        list,tuple
    :param ftime:       End date/time. Same format as itime.
    :type ftime:        list,tuple
    :param interval:    interval between output times, in seconds.
    :type interval:     int
    :returns:           list of times in datenum format.

    """
    if isinstance(idate,datetime.datetime):
        # idate = (idate.year,idate.month,idate,day,idate.hour,
                    # idate.minute,idate.second)
        # fdate = (fdate.year,fdate.month,fdate,day,fdate.hour,
                    # fdate.minute,fdate.second)
        idate = datetime_to_timetuple(idate)
        fdate = datetime_to_timetuple(fdate)
    it = calendar.timegm(idate)
    ft = calendar.timegm(fdate)
    if inclusive:
        ft = ft + interval
    times = N.arange(it,ft,interval,dtype=int)
    tttimes = [ensure_timetuple(t) for t in times]
    if fmt=='datetime':
        dttimes = [timetuple_to_datetime(t) for t in tttimes]
        return dttimes
    else:
        return times

def generate_colours(M,n):
    """
    M       :   Matplotlib instance
    n       :   number of colours you want

    Returns

    Usage: when cycling over n plots, the colour should
    be colourlist[n].
    """

    colourlist = [M.cm.spectral(i) for i in N.linspace(0.08,0.97,n)]
    return colourlist

def get_sequence(x,sos=0):
    """ Returns a sequence (tuple or list) for iteration.
    Avoids an error for strings/integers.
    SoS = 1 enables the check for a sequence of sequences (list of dates)
    """
    # If sos is True, then use its first element as an example.
    if sos:
        y = x[0]
    else:
        y = x

    if isinstance(y, collections.Sequence) and not isinstance(y, str):
        # i.e., if y is a list or tuple
        return x
    else:
        # make a one-element list
        return [x,]

def convert_tuple_to_dntimes(times):
    """
    Convert tuple or tuple of tuples to datenum date format.
    """
    timelist = get_sequence(times,sos=1)

    dntimes = []
    for t in timelist:
        dntimes.append(calendar.timegm(t))

    return dntimes

def ensure_datenum(times,fmt='int'):
    """
    Make sure times are in list-of-datenums format.
    If not, convert them.

    Possibilities:
    times = 123456                                      #1
    times = (123456,)                                   #2
    times = (123456,234567)                             #3
    times = (2011,12,1,18,0,0)                          #4
    times = ((2011,12,1,18,0,0),(2011,12,2,6,0,0))      #5

    fmt     :   whether to return list of integers or an integer
                'int' or 'list'

    Output:
    dntimes = (123456,) or (123456,234567)
    """
    if isinstance(times,int):
        dntimes = [times,] #1
    elif isinstance(times,str):
        print("Don't give me strings...")
        raise Exception
    elif isinstance(times,datetime.datetime):
        dntimes = convert_tuple_to_dntimes(datetime_to_timetuple(times))
    elif isinstance(times,(list,tuple)): #2,3,4,5
        if not isinstance(times[0],int): #5
            dntimes = convert_tuple_to_dntimes(times)
        elif times[0]<3000: #4
            dntimes = convert_tuple_to_dntimes(times)
        elif isinstance(times[0],datetime.datetime):
            dntimes = [convert_tuple_to_dntimes(datetime_to_timetuple(t))
                        for t in times]
        else: #2,3
            dntimes = times

    if (fmt == 'list') or (len(dntimes)>1):
        return dntimes
    elif (fmt == 'int') or (len(dntimes)==1):
        return dntimes[0]
    else:
        print("Nonsense format choice.")
        raise Exception

    import pdb; pdb.set_trace()

def ensure_timetuple(times,fmt='single'):
    """
    MAke sure time(s) are in six-item tuple format
    (YYYY,MM,DD,HH,MM,SS)

    fmt     :   whether to return a list of tuples or single tuple.
                'list' or 'single'

    Possibilities:
    times = 123456                                      #1
    times = (123456,)                                   #2
    times = (123456,234567)                             #3
    times = (2011,12,1,18,0,0)                          #4
    times = ((2011,12,1,18,0,0),(2011,12,2,6,0,0))      #5
    """
    if isinstance(times,(int,N.int64)):
        tttimes = [list(time.gmtime(times)),] #1
    elif isinstance(times,str):
        print("Don't give me strings...")
        raise Exception
    elif isinstance(times,datetime.datetime):
        tttimes = datetime_to_timetuple(times)
    # elif isinstance(times[0],datetime.datetime):
        # tttimes = datetime_to_timetuple(times)
    elif isinstance(times,(list,tuple)): #2,3,4,5
        if not isinstance(times[0],int): #5
            tttimes = times
        elif times[0]<3000: #4
            tttimes = [times,]
        elif isinstance(times[0],datetime.datetime):
            tttimes = [datetime_to_timetuple(t) for t in times]
        elif isinstance(times[0]>3000): #2,3
            tttimes = [list(time.gmtime(t)) for t in times]

    # import pdb; pdb.set_trace()

    if (fmt == 'list') or (len(tttimes)>1):
        return tttimes
    elif (fmt == 'single') or (len(tttimes)==1):
        return tttimes[0]
    else:
        print("Nonsense format choice.")
        raise Exception

def ensure_datetime(t):
    """
    Possibilities:
    times = 123456                                      #1
    times = (123456,)                                   #2
    times = (123456,234567)                             #3
    times = (2011,12,1,18,0,0)                          #4
    times = ((2011,12,1,18,0,0),(2011,12,2,6,0,0))      #5
    times = datetime.datetime                           #6
    times = (datetime.datetime, datetime.datetime)      #7
    """
    if isinstance(t,(int,N.int64)):  #1
        utc = timetuple_to_datetime(list(time.gmtime(t)))
    elif isinstance(t,datetime.datetime): #6
        utc = t
    elif isinstance(t,(list,tuple)):
        if isinstance(t[0],datetime.datetime): #7
            utc = t
        elif isinstance(t[0],(int,N.int64)):
            if t[0] < 3000: # 4
                utc = timetuple_to_datetime(t)
            else: #2, 3
                utc = [timetuple_to_datetime(list(time.gmtime(t))) for
                        t in t]
        elif isinstance(t[0],(list,tuple)): # 5
            utc = [timetuple_to_datetime(t) for t in t]
    else:
        raise Exception("Unidentified format.")
    return utc

def datetime_to_timetuple(utc):
    tttime = (utc.year,utc.month,utc.day,
                utc.hour,utc.minute,utc.second)
    return tttime

def timetuple_to_datetime(utc):
    dttime = datetime.datetime(*utc[:6])
    return dttime

def dt_from_fnames(f1,f2,model):
    """Work out time difference between two data files
    from their naming scheme.


    Arguments:
        f1,f2 (str): filename, with or without extension
        model (str): model used to generate data
    Returns:
        Difference between files, in seconds.
    """
    if f1.endswith('.nc'):
        f1 = f1.split('.')[0]
        f2 = f2.split('.')[0]
    if (model=='wrfout') or (model=='wrf'):
        # We assume default naming
        t = []
        for f in (f1,f2):
            _1, _2, tstr = f.split('_',2)
            fmt = '%Y-%m-%d_%H:%M:%S'
            t.append(datetime.datetime.strptime(tstr,fmt))
        dt = t[1] - t[0]
        return dt.seconds

def get_netcdf_naming(model,t,dom=0):
    """
    By default:
    wrfout files don't have an extension
    other files have .nc extension (convert first)
    """

    t = ensure_datetime(t)
    # import pdb; pdb.set_trace()
    if (model=='wrfout') or (model=='wrf'):
        if not dom:
            print("No domain specified; using domain #1.")
            dom = 1
        # fname = ('wrfout_d{0:02d}_{1:04d}-{2:02d}-{3:02d}_{4:02d}:{5:02d}:{6:02d}'.format(dom,*t))
        fname = ('wrfout_d{0:02d}_{1:04d}-{2:02d}-{3:02d}_{4:02d}:{5:02d}:{6:02d}'.format(dom,
                                     t.year,t.month,t.day,t.hour,t.minute,t.second))
    elif model == 'ruc':
        # This depends on the RUC version? Will break?


        # prefix = ruc_naming_prefix(t)

        # fname = (prefix+'{0:04d}{1:r2d}{2:02d}_{3:02d}{4:02d}_{5:02d}0.nc'.format(*t))
        fname = getdata.RUC_fname(t,filetype='netcdf')
        # import pdb; pdb.set_trace()
    else:
        print("Model format not supported yet.")
        raise Exception
    return fname

def determine_model(fname):
    """
    Return model depending on naming convention.

    If no model exists, return false.
    """

    models = {'wrfout_d':'wrfout','ruc':'ruc','rap':'ruc'}

    for k,v in models.items():
        if k in fname[:10]:
            return v

    return False

def save_data(data,folder,fname,format='pickle'):
    """
    Save array to file.
    """

    # Strip file extension given
    fname_base = os.path.splitext(fname)[0]
    # Check for folder, create if necessary
    trycreate(folder)
    # Create absolute path
    fpath = os.path.join(folder,fname_base)

    if format=='pickle':
        with open(fpath+'.pickle','wb') as f:
            pickle.dump(data,f)
    elif format=='numpy':
        N.save(fpath,data)
    elif format=='json':
        j = json.dumps(data)
        with open(fpath+'.json','w') as f:
            print(j, file=f)
    else:
        print("Give suitable saving format.")
        raise Exception

    print(("Saved file {0} to {1}.".format(fname,folder)))

def load_data(folder,fname,format='pickle'):
    """
    Load array from file.
    """

    fname2 = os.path.splitext(fname)[0]
    fpath = os.path.join(folder,fname2)
    if format=='pickle':
        with open(fpath+'.pickle','rb') as f:
            data = pickle.load(f)
    elif format=='numpy':
        data = N.load(fpath+'.npy')
    elif format=='json':
        print("JSON stuff not coded yet.")
        raise Exception
    else:
        print("Give suitable loading format.")
        raise Exception

    print(("Loaded file {0} from {1}.".format(fname,folder)))
    return data


def return_subdomain(data,lats,lons,Nlim,Elim,Slim,Wlim,
                        fmt='latlon',axes=False):
    """
    Returns smaller domain of data and lats/lons based
    on specified limits.
    """
    # TODO - might need changing
    if lats.ndim == 2:
        lats = lats[:,0]
    if lons.ndim == 2:
        lons = lons[0,:]

    Nidx = closest(lats,Nlim)
    Eidx = closest(lons,Elim)
    Sidx = closest(lats,Slim)
    Widx = closest(lons,Wlim)

    # Assuming [lats,lons]
    if fmt=='latlon':
        if Nidx<Sidx:
            xmin,xmax = Nidx,Sidx
        else:
            xmin,xmax = Sidx,Nidx

        if Widx<Eidx:
            ymin,ymax = Widx,Eidx
        else:
            ymin,ymax = Eidx,Widx
    elif fmt=='lonlat':
        if Nidx<Sidx:
            ymin,ymax = Nidx,Sidx
        else:
            ymin,ymax = Sidx,Nidx

        if Widx<Eidx:
            xmin,xmax = Widx,Eidx
        else:
            xmin,xmax = Eidx,Widx
    else:
        print("Need right format")
        raise Exception
    # Radar data: latlon, N<S, W<E
    # data = data[Nidx:Sidx+1,Widx:Eidx+1]
    # pdb.set_trace()

    data = data[...,xmin:xmax+1,ymin:ymax+1]

    if Nidx<Sidx:
        lats = lats[Nidx:Sidx+1]
    else:
        # flipud for RUC data - does this break WRF?
        lats = lats[Sidx:Nidx+1]
        # lats = N.flipud(lats[Sidx:Nidx+1])
    if Widx<Eidx:
        lons = lons[Widx:Eidx+1]
    else:
        lons = lons[Eidx:Widx+1]

    # pdb.set_trace()
    return data,lats,lons

def interp_latlon(data,lat,lon,lats,lons):
    ntimes = data.shape[0]
    nlvs = data.shape[1]
    dataout = N.zeros([ntimes,nlvs,1,1])
    for lv in range(nlvs):
    # for t in range(ntimes):
        dataout[:,lv,0,0] = interp2point(data[:,lv:lv+1,:,:],lat,lon,lats,lons)
    return dataout


def interp2point(data,lat_loc,lon_loc,lat,lon,lvidx=0,xyidx=False):
        er = 6370000
        if xyidx:
            # Don't need data, ignore
            field = None
            data = None
        else:
            field = data[:,lvidx,:,:]
            # field = self.make_4D(data)[:,lvidx,:,:]
        templat = lat.ravel()
        templon = lon.ravel()

        #calculates the great circle distance from the inquired point
        delta = haversine_baker(lon_loc,lat_loc,templon,templat,earth_rad=er)
        #grid distance
        dxdy = haversine_baker(templon[0],templat[0],templon[1],templat[0],earth_rad=er)

        #9 smallest values to find the index of
        smallest = heapq.nsmallest(9,delta.ravel())

        wtf = N.in1d(delta.ravel(),smallest)
        tf2d = N.in1d(delta.ravel(),smallest).reshape(lat.shape)
        ix,iy = N.where(tf2d == True)
        if xyidx:
            xidx = N.median(ix)
            yidx = N.median(iy)
            return xidx, yidx

        weights =  1.0 - delta.ravel()[wtf]/(dxdy.ravel() * 2)
        weighted_mean = N.average(field[:,ix,iy],axis=1,weights=weights.ravel())

        return weighted_mean

def haversine_baker(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
    """
    Allows to calculate geographical distance
    using the haversine formula.
    :param lon1: longitude of the first set of locations
    :type lon1: numpy.ndarray
    :param lat1: latitude of the frist set of locations
    :type lat1: numpy.ndarray
    :param lon2: longitude of the second set of locations
    :type lon2: numpy.float64
    :param lat2: latitude of the second set of locations
    :type lat2: numpy.float64
    :keyword radians: states if locations are given in terms of radians
    :type radians: bool
    :keyword earth_rad: radius of the earth in km
    :type earth_rad: float
    :returns: geographical distance in km
    :rtype: numpy.ndarray
    """

    if radians == False:
        cfact = N.pi / 180.0
        lon1 = cfact * lon1
        lat1 = cfact * lat1
        lon2 = cfact * lon2
        lat2 = cfact * lat2

    # Number of locations in each set of points
    if not N.shape(lon1):
        nlocs1 = 1
        lon1 = N.array([lon1])
        lat1 = N.array([lat1])
    else:
        nlocs1 = N.max(N.shape(lon1))
    if not N.shape(lon2):
        nlocs2 = 1
        lon2 = N.array([lon2])
        lat2 = N.array([lat2])
    else:
        nlocs2 = N.max(N.shape(lon2))
    # Pre-allocate array
    distance = N.zeros((nlocs1, nlocs2))
    i = 0
    while i < nlocs2:
        # Perform distance calculation
        dlat = lat1 - lat2[i]
        dlon = lon1 - lon2[i]
        aval = (N.sin(dlat / 2.) ** 2.) + (N.cos(lat1) * N.cos(lat2[i]) * (N.sin(dlon / 2.) ** 2.))
        distance[:, i] = (2. * earth_rad * N.arctan2(N.sqrt(aval), N.sqrt(1 - aval))).T
        i += 1
    return distance.ravel()

def get_latlon_idx(lats,lons,lat,lon):
    coords = N.unravel_index(N.argmin((lat-lats)**2+
                (lon-lons)**2),lons.shape)
    # lon, lat
    return [int(c) for c in coords]

def make_subplot_label(ax,label):
    if not label.endswith(')'):
        label = label + ')'
    ax.text(0.1,0.15,label,transform=ax.transAxes,
        bbox={'facecolor':'white'},fontsize=13,zorder=1000)
    return ax

def check_same_dimensions(*args):
    """
    Args:
        ``*args``       :   A number of N.ndarray types

    Returns:
        True if all arrays have the same size.
    """
    arr = []
    for arr in args:
        shps.append(arr.shape)
    return all([a == arr[0] for a in arr])

def merge_netcdfs(outpath,filelist=None,globpath=None):
    """ Merge numerous netcdf files.

    Args:
        filelist (list,tuple): Absolute paths to all netcdf files.
        globpath (str): Absolute path with glob search
        outpath (str): Absolute path to merged file.
    """
    print("Merging netCDF files...")

    if filelist is None:
        assert globpath is not None
        gp = unix_tools.enforce_pathobj(globpath)
        globdir = gp.parent
        filelist = globdir.glob(gp.name)

    # Make sure paths are strings
    fl = [str(f) for f in filelist]

    xnc = [xarray.open_dataset(nc) for nc in fl]
    merged = xarray.concat(xnc, 'Time')
    # merged = xarray.concat(xnc, 'forecast_time')
<<<<<<< HEAD
    merged.to_netcdf(str(outpath),format="NETCDF4")
    print("Completed merge with NETCDF4 format; saved to {}".format(outpath))
=======
    merged.to_netcdf(str(outpath),format='NETCDF4')
    print("Completed merge; saved to {}".format(outpath))
>>>>>>> 996b97f8a996124bea953099fe3aeada369ae50c
    return

def pretty_fhr(fhr,in_fmt='hours',out_fmt=1):
    """ Returns a pretty time stamp suitable for looping plots
    over numerous forecast times.

    Todo:
        * Merge with :class:`~evac.utils.timetool`?

    Args:
        fhr: offset from initialisation time (units given by in_fmt).
        in_fmt (str): 'hours', 'minutes', or 'seconds'
        out_fmt (int): selects the format of the output string.
                        out_fmt=1 is "{h}h_{mm}m"
    """
    if in_fmt is 'hours':
        fhr = fhr * 60
    else:
        raise Exception("Not implemented yet.")

    if out_fmt == 1:
        h,m = divmod(fhr,60)
        outstr = '{:d}h_{:02d}min'.format(int(h),int(m))
    return outstr

def hr_to_sec(hr,fmt=int):
    """ Converts hours to seconds.

    Args:
        hr (int,float): number of hours
        fmt (type): file type of output. Integer, by default.
    """
    sec = hr * 3600
    return fmt(sec)

def get_random(seq,unique=False):
    """ Pick random member from a collection or sequence.

    Note:
        As a set is unique, there is an equal chance of picking each.
            The rest are not, so the chance of picking a given value
            is proportional to its frequency, unless unique is set.
        Also note this method is a lot slower than just converting to
            a list and choosing the first element, if an arbirary
            pick is needed rather than random each time.

    Args:
        seq (list, tuple, set, N.ndarray): pick from this.
        unique (bool): if True, pick random from unique values only.
    """
    og_seq_type = type(seq)
    if isinstance(seq,N.ndarray):
        seq = list(seq.flatten())
    if unique:
        seq = set(seq)
    if len(seq) == 0:
        raise Exception("The {} specified for seq is empty.".format(og_seq_type))
    return random.sample(seq,1)[0]

def get_ccrs_proj(proj,ckws=None):
    """ Convert basemap style strings to cartopy
    projections. 

    Args:
        proj (str): name of projection
        ckws (dict): dictionary of keyword arguments to pass
            to the crs instantiation.
    """
    if ckws is None:
        ckws = dict()

    PROJS = dict(
                lcc=ccrs.LambertConformal,
                merc=ccrs.Mercator,
                plate=ccrs.PlateCarree,
                )

    return PROJS[proj](**ckws)

def colorblind_friendly(keys,fmt='RGB'):
    colors = collections.OrderedDict(
                blue=(0,114,178),
                vermillion=(213,94,0),
                reddishpurple=(204,121,167),
                black=(0,0,0),
                orange=(230,159,0),
                skyblue=(86,180,233),
                bluishgreen=(0,158,115),
                yellow=(240,228,66),
                )

    NEWDICT = collections.OrderedDict()
    colorloop = itertools.cycle(list(colors.keys()))
    for key,color in zip(keys,colorloop):
        if fmt == 'RGB':
            NEWDICT[key] = [n/255 for n in colors[color]]
        else:
            raise Exception
    return NEWDICT



def time_me(f):
    """ Decorator for timing functions.

    Usage:

    @time_me
    def function():
        pass
    """
    def timed(*args, **kwargs):
        t0 = time.time()
        result = f(*args, **kwargs)
        t1 = time.time()

        print('function {} took {:2.3f} sec'.format(f.__name__,t1-t0))
        return result

    return timed

def loop_me(*args):
    """
    TODO: a decorator function that loops an existing function
    over all arguments passed.

    args should be a list of 1D numpy arrays, tuples, or lists
    """
    pass


def bridge_multi(orders,fromlist=False,todir=False,
                    return_newlist=False):
    """
    Soft-link, copy, or move multiple items.

    Args:
    orders      :   dictionary with the following layout:
                    dict( (frompath,topath)=cmd)
                    If orders is string, treat it as command,
                    and require arguments fromlist and todir.

    where

    cmd is 'copy','softlink', or 'move'
    frompath/topath is absolute path of origin/destination

    return_newlist (bool): if True, return a list with the new absolute
                                path of all files.
    """
    # Get into correct format
    if isinstance(orders,str):
        cmd = orders
        orders = {}
        assert cmd in ('copy','move','softlink')
        for f in fromlist:
            orders[(f,todir)] = cmd

    newlist = []
    for (frompath,topath),cmd in orders.items():
        newf = bridge(cmd,frompath,topath,return_newpath=True)
        newlist.append(newf)

    if return_newlist:
        return newlist

def bridge(command,frompath,topath,catch_overwrite=False,
                return_newpath=False):
    """ Copy, move, soft-link.

    Args:

    command     :   (str) - must be one of "copy", "softlink",
                    "move".
    frompath    :   (str,Path) - absolute path to origin directory or file
    topath      :   (str,Path) - absolute path to destination directory or
                    file. If topath is only a directory, the link will
                    be the same filename as the origin.
    catch_overwrite     :   (bool) - if True, raise Exception if
                            topath already exists
    return_newpath (bool): if True, return the new absolute path of the file
    """

    # Use path-like objects
    frompath = enforce_pathobj(frompath)
    topath = enforce_pathobj(topath)
   
    if topath.is_dir():
        topath = topath / frompath.name

    if catch_overwrite and topath.exists():
        raise Exception("{} already exists.".format(topath))

    trycreate(topath)

    if command is "copy":
        # This is ridiculous, but that's what I get for using sodding pathlib
        topath.write_bytes(frompath.read_bytes())
    elif command is "move":
        # This replaces an existing 'topath' silently
        frompath.rename(topath)
    elif command is "softlink":
        frompath.symlink_to(topath)
    else:
        raise NonsenseError("The command variable **{}** is invalid".format(
                        command),color='red')

    if not topath.exists():
        print("File not created at",topath)
    else:
        print("File created at",topath)

    if return_newpath:
        return topath

def wowprint(message,color='red',bold=True,underline=False,formatargs=False):
    """ Print with colour/bold emphasis on certain things!
   
    message      :   (str) - string with any emphasised
                    item surrounded by two asterisks on both
                    sides. These are replaced.
    color        :   (str,bool) - the colour of the emphasised
                    text. If False, don't emphasise.
    
    Usage:
    
    wowprint("The program **prog** is broken",color='red')
    wowprint("The program {} is broken",formatargs=('prog',),color='blue',
                            underline=True)
    """
    # If nothing is set, pass through without changes
    if not any([color,bold,underline,formatargs]):
        print(message)
        return

    def findindex(s, lookup = '**'):
            try:
                idx = s.index('**')
            except ValueError:
                raise Exception("The wowprint string needs to have two"
                                    " lots of double asterisks.")
            else:
                return idx
            
    colors = {}
    colors['red'] = "\033[91m" # fail
    colors['purple'] = "\033[95m" # header
    colors['blue'] = "\033[94m" # OK
    colors['yellow'] = "\033[93m" # warning
    colors['green'] = "\033[92m" # OK
    
    colors['bold'] = "\033[1m" # bold
    colors['underline'] = "\033[4m" # underline
    
    endcode = "\033[0m" # use at the end, always
    
    codes = []
    if bold:
        codes.append(colors['bold'])
    if underline:
        codes.append(colors['underline'])
    if color:
        codes.append(colors[color])
    
    colorcode = ''.join(codes)
    # CASE 1: asterisks round coloured part
    if not formatargs:
        idx0 = findindex('**')
        message = message.replace('**',colorcode,1)
        idx1 = findindex('**')
        message = message.replace('**',endcode,1)
        print(message)
    else:
        # This needs to be implemented
        raise Exception
    # pdb.set_trace()
    
    return
    

def _bridge(cmd,frompath,topath,mv=False,cp=False,ln=False):
    """
    Olde version.
    Soft-link, copy, or move item.

    Create folder if it doesn't exist
    """
    #if sum([mv,cp,ln]) != 1:
    #    raise Exception("Choose one of mv, ln, cp commands.")
    todir, tofile = topath.split()
    trycreate(todir)
    if cmd in ('mv','move'):
        c = 'mv'
    elif cmd in ('cp','copy'):
        c = 'cp'
    elif cmd in ('ln','link','softlink'):
        c = 'ln -sf'

    cmd = "{} {} {}".format(c,frompath,topath)
    return cmd

def trycreate(loc, parents=True,exist_ok=True,isdir=False):
    """
    Args:

    loc     :   (Path,str) - path-like object or string to directory
                or file. If loc is a directory, this method will attempt to
                create a folder if it doesn't already exist. If loc is
                a file, its parent directory will be treated as loc.
    parents  :   (bool) - if True, then equivalent to mkdir -p </path/>
    exist_ok    :   (bool) - if True, then warnings are ignored about folder already
                    existing.
    """
    l = enforce_pathobj(loc)

    # First, check if file or directory exists
    if l.exists():
        if l.is_dir():
            wowprint("Checking **{}** exists.".format(l),color='blue')
            # print("Directory already exists.")
        else:
            wowprint("Checking **{}** exists.".format(l.parent),color='blue')
            # print("File already exists.")
        print("Directory already exists.")
        return

    # If not, create the directory enclosing the file, or
    # the directory.

    # Assume that l is a path to a file, so we need l.parent for the directory
    if not isdir:
        l = l.parent

    l.mkdir(parents=parents,exist_ok=True)

    # Does the location exist?
    assert l.exists()
    wowprint("The directory **{}** has been made.".format(l),color='red')
    # pdb.set_trace()
    return


def softlink(frompath,topath):
    """ Arguments must be path objects or strings.

    Args:

    frompath    :   (str,Path) - absolute path to origin directory or file
    topath      :   (str,Path) - absolute path to destination directory or
                    file. If topath is only a directory, the link will
                    be the same filename as the origin.
    """

    # Use path-like objects
    frompath = enforce_pathobj(frompath)
    topath = enforce_pathobj(topath)

    # Check whether paths are files or directories
    #frompath.isfile()?

    #if not fname:
        # todir is a path to the new link
        # assert # Check this is file name and not directory
        # topath = todir
    # else:
        #topath = enforce_pathobj(topath) / fname

    frompath.symlink_to(topath)

    return

def enforce_pathobj(obj,path_type='posix'):
    objs = {'posix':PosixPath}
    if isinstance(obj,str):
        return objs[path_type](obj)
    elif isinstance(obj,PosixPath):
        return obj
    else:
        raise Exception("Object is neither path nor path-like object.")


def dir_from_naming(root,*args):
    """
    Generate file path from arguments

    Inputs:
    root    :    file path base
    args     :    list of arguments to join as separate folders
    """
    l = [str(a) for a in args]
    path = os.path.join(root,*l)
    return path

def ssh_client(ky,domain,username,password):
    try:
        import paramiko
    except ImportError:
        print("Module paramiko unavailable. Ignoring import.")
        raise

    key = paramiko.RSAKey(data=base64.decodestring(ky))
    client = paramiko.SSHClient()
    client.get_host_keys().add(domain, 'ssh-rsa', key)
    client.connect(domain, username=username, password=password)
    return client
    stdin, stdout, stderr = client.exec_command('ls')
    for line in stdout:
        print('... ' + line.strip('\n'))
    client.close()


def edit_namelist(f,sett,newval,absolute=False,doms=1,samealldoms=True):
    if doms < 1:
        raise ValueError
    f = enforce_pathobj(f)
    #fs = open(f,'r')
    fs = f.open('r')
    # with open(f,'r') as fopen:
    # flines = open(f,'r').readlines()
    flines = fs.readlines()
    # pdb.set_trace()
    done= False
    for idx, line in enumerate(flines):
        if sett in line:
            # print("Setting {0} is in line.".format(sett))
            fs.close()
            # spaces = 38 - len(sett)
            spaces = 36
            # print(sett,spaces)
            if not isinstance(newval,(tuple,list)):
                newval = '{0},     '.format(newval)*doms
            elif doms == 2:
                newval = '{0},     {1},'.format(newval[0],newval[1])
            elif doms == 3:
                newval = '{0},     {1},     {2},'.format(newval[0],newval[1],newval[2])

            flines[idx] = " {0: <{sp}}= {1}\n".format(sett,newval,sp=spaces)
            nameout = f.open('w',1)
            # nameout = open(f,'w',1)
            nameout.writelines(flines)
            nameout.close()
            done = True
            break
    fs.close()
    if absolute and (not done):
            # import pdb;pdb.set_trace()
        raise ValueError("Setting",sett,"not found in namelist.")
    return

def generate_timestamp_fname(extension):
    """ Returns a file name based on current time.
    """
    nowutc = datetime.datetime.now(tz=pytz.utc)
    fname = gis_tools.string_from_time('output',nowutc,strlen='second')
    return '.'.join(fname,extension)

def enforce_2d(data):
    if data.ndim == 2:
        return data
    elif data.ndim == 3:
        return data[0,:,:]
    elif data.ndim == 4:
        return data[0,0,:,:]
    elif data.ndim == 5:
        return data[0,0,0,:,:]
    else:
        raise Exception


def enforce_same_dimensions(*args):
    theshape = args[0].shape
    for arg in args:
        assert arg.shape == theshape
    return

def exceed_probs_2d(arr3D,val,overunder='over',fmt='pc'):
    """ Calculates the exceedence probability lat/lon field.
    """
    assert arr3D.ndim == 3

    nmems = arr3D.shape[0]

    comparefunc = dict(over = N.greater,under = N.less)

    # True/False if member meets condition (5D)
    bool_arr = N.where(comparefunc[overunder](arr3D,val),1,0)

    # Count members that exceed the threshold (4D)
    count_arr = N.sum(bool_arr,axis=0)

    if fmt == 'pc':
        # And convert to percentage (4D) for each time
        percent_arr = 100*(count_arr/nmems)

        return percent_arr # [times,levels,lats,lons]
    elif fmt == 'decimal':
        dec_arr = count_arr/nmems
        return dec_arr
    else:
        raise Exception

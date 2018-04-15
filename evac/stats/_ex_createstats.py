import pdb
import pickle
import random
import datetime
import glob
import os
import sys
import copy
import operator
import itertools
import time

from scipy.interpolate import griddata
import scipy.ndimage.filters
import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import RectBivariateSpline as RBS
# from scipy.interpolate import griddata

import WEM
import WEM.utils as utils
from WEM.postWRF.postWRF.wrfout import WRFOut
from WEM.postWRF.postWRF.ensemble import Ensemble
from WEM.postWRF.postWRF.obs import MRMS, StageIV
from WEM.postWRF.postWRF.hrrr import HRRR
from WEM.postWRF.postWRF.sal import SAL
from WEM.postWRF.postWRF.casati import Casati

from evac.stats.detscores import DetScores
from evac.stats.roc import ROC

# print(WEM.__path__)

##### SETTINGS #####
serial = False
testplots = False
debug = False
ncpu = 30
# stats = ('FSS',)
# newstype = 'all'
stats = ('SAL',)
newstype = 'mean'
# stats = ('CASATI',)
salthresh_d = 24
threshs = [(1/64),(1/32),(1/16),(1/8),(1/4),(1/2),1,2,4,8,16,32]
thresholds = threshs
models = ('NEWSe',)
# models = ('HRRRx','HRRRo','NEWSe')

outdir = '/home/john.lawson/HRRRvNEWSe/pyoutput'
cases = [datetime.datetime(2016,5,n,0,0,0) for n in (7,8,9,10,16,17,22,23,24,25,26,27)]
# cases = [datetime.datetime(2016,5,17,0,0,0),]

latidxx = N.arange(5,260,15)
lonidxx = N.arange(5,260,15)

statsdir = '/work/john.lawson/HRRRvNEWSe/stats'
# statsdir = '/scratch/john.lawson/HRRRvNEWSe/stats'
rootdir = {}
rootdir['NEWS'] = '/work/wof/realtime/FCST'
# rootdir['myNEWS'] = '/work/john.lawson/NEWS_postproc'
# rootdir['myHRRR'] = '/scratch/john.lawson/HRRR_postproc'
# rootdir['myNEWS'] = '/scratch/john.lawson/NEWS_postproc'
rootdir['HRRRx'] = '/work/john.lawson/HRRR_data/experimental'
rootdir['HRRRo'] = '/work/john.lawson/HRRR_data/operational'
rootdir['STAGEIV'] = '/work/john.lawson/STAGEIV_data'
rootdir['MRMS'] = '/work/john.lawson/MRMS_data/May2016'

ensnames = ["ENS_MEM_{0}".format(n+1) for n in range(18)]

process_HRRRo = 1
process_HRRRx = 0

vrbl = 'accum_precip'

### FUNCTIONS ###
def dt_from_wrfout(fpath):
    if fpath.startswith('/'):
        f = os.path.basename(fpath)
    else:
        f = fpath
    fmt = 'wrfout_d01_%Y-%m-%d_%H:%M:%S'
    dt = datetime.datetime.strptime(f,fmt)
    return dt

def append_missing(missingdata,rootdir,initutc,itime,ftime):
    missingdata.append((rootdir,initutc,itime,ftime))
    return missingdata
          

def tryget(E,prod,vrbl,itime,ftime,lv,missingdata,member=None,latidx=None,
                    lonidx=None,thresholds=None):
    try:
        if prod == 'mean':
            fcstdata = E.mean(vrbl,itime=itime,ftime=ftime,level=lv)[0,0,:,:]
        elif prod == 'member':
            fcstdata = E.get(vrbl,members=member,itime=itime,ftime=ftime,level=lv)[0,0,:,:]
            # pdb.set_trace()
        elif prod == 'allpoints':
            fcstdata = E.get(vrbl,ftime=ftime,itime=itime,lats=latidx,lons=lonidx,level=lv)
        elif prod == 'exceedprob':
            fcstdata = []
            for thresh in thresholds:
                fcstdata.append(E.get_prob_threshold(vrbl,'over',thresh,
                        itime=itime,ftime=ftime,level=lv))#[:,:])
        else:
            raise Exception("Specify a product.")
    except: #Change this to actual missing time error (in ensemble class)
        fcstdata = None
        missingdata = append_missing(missingdata,E.rootdir,E.initutc,itime,ftime)
    return fcstdata, missingdata


### DATASETS ###

# Stage Iv verification
print("Loading STAGE IV data")
ST4 = StageIV(rootdir['STAGEIV'])
print("Done.") 

print("Loading NEWS-e catalogue.")
# NEWS-e data
NEWS_CAT = {}
for case in cases:
    NEWS_CAT[case] = {}
    init_time = case + datetime.timedelta(seconds=60*60*15)
    tlim = case + datetime.timedelta(seconds=60*60*37)
    casestr = "{0:04d}{1:02d}{2:02d}".format(case.year,case.month,case.day)
    casedir = os.path.join(rootdir['NEWS'],casestr)
    while init_time < tlim:
        ls = glob.glob(os.path.join(casedir,'*'))
        itime_str = "{0:02d}{1:02d}".format(init_time.hour,init_time.minute)
        if os.path.join(casedir,itime_str) in ls:
            NEWS_CAT[case][itime_str] = {}
            NEWS_CAT[case][itime_str]['init_time'] = init_time
            NEWS_CAT[case][itime_str]['fdir'] = os.path.join(casedir,itime_str)
            for en in ensnames:
                fdir = os.path.join(casedir,itime_str,en)
                NEWS_CAT[case][itime_str][en] = {}
                NEWS_CAT[case][itime_str][en]['fdir'] = fdir
                wrfouts = glob.glob(os.path.join(fdir,'wrfout*'))
                for w in wrfouts:
                    validtime = dt_from_wrfout(w)
                    NEWS_CAT[case][itime_str][en][validtime] = w
        init_time = init_time + datetime.timedelta(seconds=60*30)
print("Done.")

# HRRR operational and/or experimental data
for vers in ('HRRRo','HRRRx'):
    HRRR_CAT = {}
    for case in cases:
        HRRR_CAT[case] = {}
        init_time = case + datetime.timedelta(seconds=60*60*15)
        tlim = case + datetime.timedelta(seconds=60*60*37)
        casestr = "{0:04d}{1:02d}{2:02d}".format(case.year,case.month,case.day)
        casedir = os.path.join(rootdir[vers],casestr)
        HRRR_CAT[case]['fdir'] = casedir
        while init_time < tlim:
            if init_time.day > case.day:
                case2 = case + datetime.timedelta(seconds=60*60*24)
                case2str = "{0:04d}{1:02d}{2:02d}".format(case2.year,case2.month,case2.day)
                case2dir = os.path.join(rootdir[vers],case2str)
            else:
                case2 = case
                case2str = casestr
                case2dir = casedir

            ls = glob.glob(os.path.join(case2dir,'*'))
            lsdir = [f for f in ls if not os.path.basename(f).endswith('zip')]
            itime_str = "{0}{1:02d}{2:02d}".format(case2str,init_time.hour,init_time.minute)
            fdir = os.path.join(case2dir,itime_str)
            itime_key = itime_str[-4:]

            # print(init_time,itime_str,fdir,case2)
            if fdir in ls:
                HRRR_CAT[case][itime_key] = {}
                HRRR_CAT[case][itime_key]['fdir'] = fdir
                HRRR_CAT[case][itime_key]['init_time'] = init_time
                fcsts = glob.glob(os.path.join(fdir,'*'))
                fcst_count = 0
                # pdb.set_trace()
                for f in fcsts:
                    if f.endswith('npy'):
                        continue
                    fcsthr = int(f[-2:])
                    fdir2 = os.path.join(fdir,f)
                    HRRR_CAT[case][itime_key][fcsthr] = fdir2

            init_time = init_time + datetime.timedelta(seconds=60*60)
    if vers == 'HRRRo':
        HRRRo_CAT = HRRR_CAT
    else:
        HRRRx_CAT = HRRR_CAT

#### PROCEDURE ####
def RMSE(E,ST4,itime,ftime,missingdata):
    pass

def get_closest(E,fcsttime):
    try:
        mem = E.closest_to_mean('REFL_10CM',level=0,fcsttime=fcsttime)
    except:
        print("Missing data for this time?")
        mem = None
    return mem

def get_contingency(fc,ob):
    # Loop through thresholds
    ABCDs = {}
    for th in threshs:
        ABCD[th] = compute_contingency(fc,ob,th,'over')
        # ABCDs[th] = [DS.a,DS.b,DS.c,DS.d]
    return ABCDs

def compute_roc(E,ST4,itime,ftime,missingdata,thresholds,
                xx,yy,newgrid):
    obsdata = ST4.return_array(ftime)
    print("Starting the obs data reprojection")
    reproj_ob = reproj('ST4',obsdata.data,xx=xx,yy=yy,newgrid=newgrid)
    print("Done.")

    ABCDs = {}
    fcst,missingdata = tryget(E,'exceedprob',vrbl,itime,ftime,None,
                                    missingdata,thresholds=thresholds)
    # pdb.set_trace()
    if fcst is not None:
        print("About to loop over thresholds.")
        for thidx,th in enumerate(thresholds):
            roc = ROC(pf=fcst[thidx],ob=reproj_ob,nens=E.nmems)
            ABCDs[th] = [roc.a,roc.b,roc.c,roc.d]
            # ABCD[th] = compute_contingency(fcst[thidx],reproj_ob,th,'over')
        print("Done threshold loop.")
    else:
        ABCDs = None
    return ABCDs, missingdata


def rankhist(freq,E,ST4,itime,ftime,missingdata):
    random.shuffle(latidxx)
    random.shuffle(lonidxx)
    lv = None
    fcsts, missingdata = tryget(E,'allpoints',vrbl,itime,ftime,lv,missingdata)
    if fcsts is None:
        print("Data missing for this date.")
    else:
        lzip = list(zip(latidxx,lonidxx))
        for latidx, lonidx in lzip:
            # Randomise a little bit
            latidx += int(N.random.uniform(-3,3))
            lonidx += int(N.random.uniform(-3,3))
            # Find point
            fcst = fcsts[:,0,0,latidx,lonidx]
            
            lon = E.arbitrary_pick(dataobj=True).get('XLONG',lons=lonidx,lats=latidx)
            lat = E.arbitrary_pick(dataobj=True).get('XLAT',lons=lonidx,lats=latidx)
            lat = lat.flatten()[0]
            lon = lon.flatten()[0]
            if debug:
                print('lat/lon idx:',latidx,lonidx)
                print('lat/lon',lat,lon)
        
            ob = ST4.return_point(ftime,lat,lon,accum_hr='01h')

            # fcsts, missingdata = tryget(E,'allpoints',vrbl,itime,ftime,lv,missingdata,
                            # lonidx=lonidx,latidx=latidx)
            fcst = fcst.flatten()

            # Check for ties between observed and ensemble members
            # uniq,uniq_n = N.unique(fcsts,return_counts=True)
            # tieidx = N.where(uniq_n>2)
            # tievals = uniq[tieidx]
            # for val in tievals:
            for idx,val in enumerate(fcst):
                if ob == val:
                    if debug:
                        print("Tie between members and observed value found. Dithering.")
                    # fcstties = N.where(fcsts == val)
                    # for tie in fcstties:
                    fcst[idx] += N.random.uniform(-1E-8,1E-8)

            binidx = N.searchsorted(N.sort(fcst),ob)
            freq[binidx] += 1
            if debug:
                print(binidx)
    return freq, missingdata

def compute_FSS(indata,ST4,itime,ftime,lv,missingdata,
            # threshs=(4,),ns=(5,),
            threshs=(0.5,1,2,4,8,16,32),ns=None,ns_step=4,
            # threshs=(0.5,1,2,4,8,16,32),ns=(21,),
            member=None,xx=False,yy=False,newgrid=None):
    # Get data for this period
    obsdata = ST4.return_array(ftime)
    # pdb.set_trace()
    if "ensemble" in indata.__module__:
        fmt = 'NEWSe'
        fcstdata, missingdata = tryget(indata,'member',vrbl,itime,ftime,lv,missingdata,member=member)
    elif "hrrr" in indata.__module__:
        fmt = 'HRRR'
        # fcstdata = indata.get(vrbl,utc=ftime,level=lv)
        # Is this for HRRRo?
        # fcstdata = indata.get(vrbl,idx=1)
        # Is this for HRRRx?
        fcstdata = indata.get(vrbl,idx=0)
    else:
        raise Exception("Need valid input data")

    if fcstdata is None:
        return None
    # Reprojections
    print("Starting the forecast data reprojection")
    reproj_fc = reproj(fmt,fcstdata,xx=xx,yy=yy,newgrid=newgrid)
    print("Done.")
    # Make masked value 0.
    N.ma.masked_values(obsdata,0.0)

    print("Starting the obs data reprojection")
    reproj_ob = reproj('ST4',obsdata.data,xx=xx,yy=yy,newgrid=newgrid)
    print("Done.")

    lenx = xx.shape[1]
    leny = xx.shape[0]
    # maxlen = (2*max(lenx,leny))-1
    maxlen = max(lenx,leny)
    if ns is None:
        ns = N.arange(1,maxlen,ns_step)

    # pdb.set_trace()
    MSE = {}
    FSS = {}

    if fcstdata is not None:
        for th in threshs:
            MSE[th] = {}
            FSS[th] = {}
            # Convert to binary using thresholds
            A
            fc = N.copy(reproj_fc)
            ob = N.copy(reproj_ob)
            fc[fc < th] = False
            fc[fc >= th] = True
            ob[ob < th] = False
            ob[ob >= th] = True
            for n in ns:
                MSE[th][n] = {}
                FSS[th][n] = {}
                # print("FSS for threshold {0} mm and n={1}.".format(th,n))
                # FSS computation w/ fractions

                method=2
                if method == 1:
                    time0 = time.time()
                    On = N.zeros([lenx,leny])
                    Mn = N.zeros([lenx,leny])
                    for i in N.arange(lenx):
                        for j in N.arange(leny):
                            IO = N.zeros([n,n])
                            IM = N.zeros([n,n])
                            for k in N.arange(n):
                                for l in N.arange(n):
                                    v1 = int((i+k-1-((n-1)/2)))-1
                                    v2 = int((j+l-1-((n-1)/2)))-1
                                    if (lenx > v1 > -1) and (leny > v2 > -1):
                                    # if (v1 > -1) and (v2 > -1):
                                        IO[k,l] = ob[v1,v2]
                                        IM[k,l] = fc[v1,v2]
                                    else:
                                        IO[k,l] = N.nan
                                        IM[k,l] = N.nan
                            # IOsum = N.sum(IO)
                            # On[i,j] = 1/(n**2)*IOsum
                            On[i,j] = N.nanmean(IO)
                            Mn[i,j] = N.nanmean(IM)
                            # print(On[i,j],Mn[i,j])

                    time1 = time.time()
                    print("Method #1 took {0}.".format(time1-time0))
                    fig,axes = plt.subplots(nrows=1,ncols=2)
                    for ax,data in zip(axes.flat,(On,Mn)):
                        ax.pcolormesh(data,cmap=M.cm.Greys)
                        ax.set_aspect('equal')
                    fig.savefig(os.path.join(outdir,'MSE_test_1.png'))

                elif method == 2:
                    pad = int((n-1)/2)
                    On = scipy.ndimage.filters.uniform_filter(ob,size=n,
                                    mode='constant',cval=0)
                    Mn = scipy.ndimage.filters.uniform_filter(fc,size=n,
                                    mode='constant',cval=0)

                    # Delete meaningless smoothed data
                    cutrangex = list(range(0,pad)) + list(range(lenx-pad,lenx))
                    cutrangey = list(range(0,pad)) + list(range(leny-pad,lenx))
                    On = N.delete(On,cutrangey,axis=0)
                    Mn = N.delete(Mn,cutrangey,axis=0)
                    On = N.delete(On,cutrangex,axis=1)
                    Mn = N.delete(Mn,cutrangex,axis=1)
                    cutlenx = On.shape[1]
                    cutleny = On.shape[0]

                testplot1 = False
                if testplot1 and (n in (1,37,77,109,141,185,229,249)):
                    fig,axes = plt.subplots(nrows=1,ncols=2)
                    for ax,data in zip(axes.flat,(On,Mn)):
                        ax.pcolormesh(data,cmap=M.cm.Greys)
                        ax.set_aspect('equal')
                    fig.savefig(os.path.join(outdir,'MSE_test_{}.png'.format(n)))
                # pdb.set_trace()

                # MSE
                sqdif = (On-Mn)**2
                MSE[th][n]['score'] = (1/(cutlenx*cutleny))*N.sum(sqdif)

                # Reference MSE
                MSE[th][n]['ref'] = (1/(cutlenx*cutleny))*(N.sum(On**2)+N.sum(Mn**2))
                # FSS
                FSS[th][n] = 1 - (MSE[th][n]['score'] / MSE[th][n]['ref'])
            # END n
        # END th

        testplot2 = False
        if testplot2:
            fig,ax = plt.subplots(1)
            for th in threshs:
                list_n = []
                for n in ns:
                    val = FSS[th][n]
                    list_n.append([0.0 if N.isnan(val) else val][0])
                ax.plot(ns,list_n,label="threshold = {} mm/hr".format(th))
            ax.legend()
            fig.savefig(os.path.join(outdir,'FSS_test.png'))
    # pdb.set_trace()
    return FSS
    
def reproj(model,data,xx=False,yy=False,lats=False,lons=False,newgrid=False):
    if model.startswith('HRRR'):
        HRRRexf = '/work/john.lawson/HRRR_data/operational/20160507/201605070200/1612802000003'
        HRRR_W = HRRR(HRRRexf)
        hrlons = HRRR_W.lons
        hrlats = HRRR_W.lats
        hrxx,hryy = newgrid(hrlons,hrlats)
        newdata = do_reproj(hrxx,hryy,data,xx,yy)
    elif model == 'ST4':
        st4lats = ST4.lats
        st4lons = ST4.lons
        stxx,styy = newgrid(st4lons,st4lats)
        newdata = do_reproj(stxx,styy,data,xx,yy)

    elif model == 'NEWSe':
        # newdata = data[::5,::5]
        newdata = data

    if len(newdata.shape) == 3:
        newdata = newdata[0,:,:]
    elif len(newdata.shape) == 4:
        newdata = newdata[0,0,:,:]

    # QC
    newdata[newdata>175] = 0
    newdata[newdata<0] = 0
    return newdata

def do_reproj(oldx,oldy,data,newx,newy,newxdim=False,newydim=False):
    newxdim = len(newx)
    newydim = len(newy)
    # pdb.set_trace()
    data_rpj = griddata((oldx.flat,oldy.flat),data.flat,(newx.flat,
        newy.flat)).reshape(newxdim,newydim)
    return data_rpj

def generate_fnames(newstype='mean'):
    models = ('NEWSe','HRRRo','HRRRx')
    STAT = {'RANKHIST':'rankhistdata',
                'FSS':'fssdata',
                'SAL':'saldata',
                'CASATI':'casatidata',
                'closest':'closest',
                'ABCD':'contingency',
                'ROC':'rocdata',
                }
    fnames = {}
    for m in models:
        fnames[m] = {}
        for s in STAT.keys():
            fnames[m][s] = {}
            for c in cases:
                fnames[m][s][c] = {}
                for i in NEWS_CAT[c].keys():
                    if m == 'NEWSe':
                        member = '_'+newstype
                    else:
                        member = ''
                    if s in ('FSS','CASATI','closest','ABCD','ROC'):
                        ext = 'pickle'
                    else:
                        ext = 'npy'
                    fname = '{0}_{1:%Y%m%d}_{2}{3}.{4}'.format(m,c,i,member,ext)
                    if s == 'SAL':
                        sp = 'SAL_{}'.format(salthresh_d)
                        fpath = os.path.join(statsdir,sp,fname)
                    else:
                        fpath = os.path.join(statsdir,s,fname)
                    fnames[m][s][c][i] = fpath
    return fnames
            
def testplot(data,grid,fname,xx,yy,salplot=False):
    mm = copy.copy(grid)
    if len(data.shape)== 3:
        assert data.shape[0] == 1
        data = data[0,...]
    if salplot:
        f1 = mm.pcolormesh(xx,yy,data)
    else:
        f1 = mm.contourf(xx,yy,data)
        # f1 = mm.contourf(xx,yy,data,levels=N.arange(0.1,25,0.1))
    mm.drawstates()
    plt.colorbar(f1)
    fpath = os.path.join(outdir,fname)
    plt.gcf().savefig(fpath)
    plt.close(plt.gcf())

    return

class get_HRRR_native_grid(object):
    def __init__(self,fpath):
        W = HRRR(fpath)
        cen_lat = 38.5
        cen_lon = -97.5
        tlat1 = 38.5
        tlat2 = 38.5
        # lllon = -105.43488 # [0,0]
        lllon = W.lons[0,0]
        #lllat = 35.835026 # [0,0]
        lllat = W.lats[0,0]
        #urlon = -96.506653 # [-1,-1]
        urlon = W.lons[-1,-1]
        #urlat = 42.708714 # [-1,-1]
        urlat = W.lats[-1,-1]
        self.m = Basemap(projection='lcc',lat_1=tlat1,lat_2=tlat2,lat_0=cen_lat,
                            lon_0=cen_lon,llcrnrlon=lllon,llcrnrlat=lllat,
                            urcrnrlon=urlon,urcrnrlat=urlat,resolution='i')
        self.lons, self.lats, self.xx, self.yy = self.m.makegrid(W.lons.shape[1],W.lons.shape[0],returnxy=True)
        # return m, lons, lats, xx[0,:], yy[:,0]

class get_NEWS_native_grid(object):
    def __init__(self,fpath):
        W = WRFOut(fpath)
        # cen_lat = 39.3569
        cen_lat = W.nc.CEN_LAT
        # cen_lon = -101.219
        cen_lon = W.nc.CEN_LON
        # tlat1 = 30.0
        tlat1 = W.nc.TRUELAT1
        # tlat2 = 60.0
        tlat2 = W.nc.TRUELAT2
        # [lats,lons]
        # lllon = -105.43488 # [0,0]
        lllon = W.lons[0,0]
        #lllat = 35.835026 # [0,0]
        lllat = W.lats[0,0]
        #urlon = -96.506653 # [-1,-1]
        urlon = W.lons[-1,-1]
        #urlat = 42.708714 # [-1,-1]
        urlat = W.lats[-1,-1]
        self.m = Basemap(projection='lcc',lat_1=tlat1,lat_2=tlat2,lat_0=cen_lat,
                            lon_0=cen_lon,llcrnrlon=lllon,llcrnrlat=lllat,
                            urcrnrlon=urlon,urcrnrlat=urlat,resolution='i')
        self.lons, self.lats, self.xx, self.yy = self.m.makegrid(W.lons.shape[1],W.lons.shape[0],returnxy=True)
        # return m, lons, lats, xx[0,:], yy[:,0]

def create_newgrid(f):
    # print(f)
    NEWS = get_NEWS_native_grid(f)
    # Create 15 km grid
    # xx = NEWS.xx[::5,::5]
    xx = NEWS.xx
    # yy = NEWS.yy[::5,::5]
    yy = NEWS.yy
    # lons = NEWS.lons[::5,::5]
    lons = NEWS.lons
    # lats = NEWS.lats[::5,::5]
    lats = NEWS.lats
    m = NEWS.m
    return xx,yy,lats,lons,m
    
def compute_casati(modeldata,obsdata,return_casati=False,threshs=threshs):
    # interpolate to finer grid
    reproj_data = {}
    for n,data in enumerate((modeldata,obsdata)):
        xx = N.arange(data.shape[-1])
        yy = N.arange(data.shape[-2])
        rbs = RBS(xx,yy,data)
        ix = N.linspace(xx.min(),xx.max(),256)
        iy = N.linspace(yy.min(),yy.max(),256)
        reproj_data[n] = rbs(ix,iy)

    casati = Casati(reproj_data[0],reproj_data[1],thresholds=threshs)
    return casati.MSE, casati.SS

def compute_sal(modeldata,obsdata,return_SAL=False,thresh=salthresh_d):
    # f should be lower than 1/15 default, as rain rates are high
    sal = SAL(obsdata,modeldata,mod_fmt='array',ctrl_fmt='array',
                footprint=100,f=(1/thresh),dx=3.0,dy=3.0)
    statsarr = N.array([sal.S,sal.A,sal.L])
    if return_SAL:
        return statsarr,sal
    else:
        return statsarr

# def compute_stats(c,i,lock=None,stats=False,models=False,vrbl='accum_precip',
def compute_stats(it,lock=None,vrbl='accum_precip',lv=None,):
    """
    m = model
    s = stat
    c = case date
    i = init times
    newstype = member (ignored if model is not NEWSe)
    """
    m,s,c,i,newstype = it

    # Get domain for this case
    for k,v in NEWS_CAT[c][i]['ENS_MEM_1'].items():
        if v.endswith('00'):
            dom_f = v
            break
    
    xx,yy,lats,lons,newgrid = create_newgrid(dom_f)
    missingdata = []
    count = 0
    print("======== Case",c,"========")
    print("-------- Init. time",i,"--------")

    print("Model",m)
    print("Stat/score",s)
    if m == 'NEWSe':
        print("NEWS-E type is",newstype)

    fnames = generate_fnames(newstype=newstype)
    npyf = fnames[m][s][c][i]

    # Check if computation was already done
    checkfiles = True
    if os.path.isfile(npyf) and checkfiles:
        print("File exists. Skipping.")
        return
    

    # List of all file times for this case/ inittime
    # filetimes = []

    if m == 'NEWSe':
        fdir = NEWS_CAT[c][i]['fdir']
        init_time = NEWS_CAT[c][i]['init_time']
        E = Ensemble(fdir,init_time,ctrl=False,loadobj=False)
        # filetimes = E.filetimes
    elif m == 'HRRRo' or 'HRRRx':
        if i.endswith('30'):
            # No HRRR for half-hour init times - no comparison
            return missingdata
        if m == 'HRRRo':
            HRRR_CAT = HRRRo_CAT
        else:
            HRRR_CAT = HRRRx_CAT
        # pdb.set_trace()
        try:
            fdir = HRRR_CAT[c][i]['fdir']
        except KeyError:
            # This file is missing
            missingdata = append_missing(missingdata,HRRR_CAT[c]['fdir'],i,'start','finish')
            print("We are missing the following {} directory: \n {}".format(
                        m,missingdata[-1]))
            return missingdata
        init_time = HRRR_CAT[c][i]['init_time']
        H = {}
        for v,w in HRRR_CAT[c][i].items():
            H[v] = w
            # filetimes.append(v)
        
    itime = init_time
    # print("******** Init hour",init_time,"********")
    # pdb.set_trace()
    if itime.minute == 30:
        itime = itime + datetime.timedelta(seconds=1800)
        endtime = init_time + datetime.timedelta(seconds=int(1.5*60*60))
        ntimes = 1
        print("Skipping the half-hour.")
    else:
        endtime = init_time + datetime.timedelta(seconds=3*60*60)
        ntimes = 3
    ftime = itime + datetime.timedelta(seconds=3600)
    print("End time is",endtime)
    # while ftime <= N.array(filetimes).max():

    # Initialise any stats arrays
    if s == 'RANKHIST':
        data = N.zeros([ntimes,19],dtype=N.int8)
    else:
        data = []

    ntime = -1
    while ftime <= endtime:
        ntime += 1
        print("++++++++ Forecast time ",ftime,"++++++++")
        dt = ftime-init_time
        fchr = int(dt.seconds/3600)

        ##### STATS METHODS #####
        if s == 'closest':
            if ntime == 0:
                print("Also adding IC closest-to-mean.")
                data.append(get_closest(E,itime))
            data.append(get_closest(E,ftime))

        if s == 'RANKHIST':
            data[ntime,:], missingdata = rankhist(data[ntime,:],E,ST4,itime,ftime,missingdata)
        elif s == 'ROC':
            data_t, missingdata = compute_roc(E,ST4,itime,ftime,missingdata, thresholds,
                                    xx,yy,newgrid)
        elif s == 'FSS':
            # pdb.set_trace()
            if m == 'NEWSe':
                indata = E
            elif m == 'HRRRo' or 'HRRRx':
                indata = HRRR(H[fchr])
            else:
                raise Exception
            data_t = compute_FSS(indata,ST4,itime,ftime,lv,missingdata,member=newstype,
                        xx=xx,yy=yy,newgrid=newgrid)
                # fcstdata = E.mean(vrbl,itime=itime,ftime=ftime,level=lv)
                # fcstdata, missingdata = tryget(E,'mean',vrbl,itime,ftime,lv,missingdata)
                # print("NEWSe mean computed from",fdir)
                # fcstdata = HRRR(H[fchr]).get(vrbl,utc=ftime,level=lv)
                # print("HRRR file is",H[fchr])
            # data_t = compute_FSS(fcstdata,ST4,itime,ftime,mean=True,closestmean=True)

        elif s in ("CASATI","SAL","ABCD"):
            if m.startswith('HRRR'):
                fcstdata = HRRR(H[fchr]).get(vrbl)
            elif m == 'NEWSe':
                if newstype == 'mean':
                    # fcstdata = E.mean(vrbl,itime=itime,ftime=ftime,level=lv)[0,0,:,:]
                    fcstdata, missingdata = tryget(E,'mean',vrbl,itime,ftime,lv,missingdata)
                    mstr = 'NEWSe_mean'
                elif isinstance(newstype,str):
                    # fcstdata = E.get(vrbl,member=newstype,itime=itime,ftime=ftime,level=lv)[0,0,:,:]
                    fcstdata, missingdata = tryget(E,'member',vrbl,itime,ftime,lv,missingdata,member=newstype)
                    mstr = 'NEWSe_{0}'.format(newstype)
                    
            if fcstdata is not None:
                # reproj_fc = reproj(m,fcstdata,xx,yy,lats,lons,newgrid)
                print("Starting the forecast data reprojection")
                reproj_fc = reproj(m,fcstdata,xx=xx,yy=yy,newgrid=newgrid)
                print("Done.")
                obsdata = ST4.return_array(ftime)
                # Make masked value 0.
                N.ma.masked_values(obsdata,0.0)

                print("Starting the obs data reprojection")
                reproj_ob = reproj('ST4',obsdata.data,xx=xx,yy=yy,newgrid=newgrid)
                print("Done.")
                if s == 'SAL':
                    data_t, sal = compute_sal(reproj_fc,reproj_ob,return_SAL=True)
                elif s == 'CASATI':
                    data_t = compute_casati(reproj_fc,reproj_ob)
                elif s == 'ABCD': # a dictionary of a,b,c,d with thresholds
                    data_t = get_contingency(reproj_fc,reproj_ob)

                # print(data_t)
                if testplots and (fchr == 1):
                    print("Doing test plots")
                    # pdb.set_trace()
                    datestr = '{0:%Y%m%d}'.format(c)
                    inits = '{0:02d}{1:02d}'.format(init_time.hour,init_time.minute)
                    fts = '{0:02d}00'.format(ftime.hour)
                    if m == 'NEWSe':
                        NEWS = get_NEWS_native_grid(dom_f)
                        oldgrid = NEWS.m
                        nexx,neyy = oldgrid(NEWS.lons,NEWS.lats)
                    elif m.startswith('HRRR'):
                        HRex = get_HRRR_native_grid(H[fchr])
                        # HRex = HRRR(H[fchr])
                        oldgrid = HRex.m
                        nexx,neyy = oldgrid(HRex.lons,HRex.lats)
                    # testplot(fcstdata,oldgrid,'test_{0}_{1}_{2}_{3}_orig.png'.format(datestr,inits,fts,mstr),nexx,neyy)
                    testplot(reproj_fc,newgrid,'test_{0}_{1}_{2}_{3}_reproj.png'.format(datestr,inits,fts,mstr),xx,yy)
                    if m == 'NEWSe':
                        stxx,styy = newgrid(ST4.lons,ST4.lats)
                        # testplot(obsdata.data,newgrid,'test_{0}_{1}_{2}_ST4_orig.png'.format(datestr,inits,fts),stxx,styy)
                        testplot(reproj_ob,newgrid,'test_{0}_{1}_{2}_ST4_reproj.png'.format(datestr,inits,fts),xx,yy)
                    # print("Done. Moving to next forecast time.")
                    
                    testplot(sal.M['obj_array'],newgrid,'test_SALobj_fcst.png',xx,yy,salplot=True)
                    testplot(sal.C['obj_array'],newgrid,'test_SALobj_obs.png',xx,yy,salplot=True)
                    assert True==False

        #########################
        if s in ("SAL",):
            if (fcstdata is None):
                print("Skipping missing data.")
                data.append(None)
            else:
                data.append(data_t)
        elif s in ("FSS","CASATI","ABCD","ROC"):
            data.append(data_t)

        print("Advancing time...")
        itime = itime + datetime.timedelta(seconds=3600)
        ftime = ftime + datetime.timedelta(seconds=3600)

    ##### SAVE OUTPUT #####
    print("~~~~~~~~~~~~~~~~~ SAVING DATA ~~~~~~~~~~~~~~~~~~")
    utils.trycreate(os.path.dirname(npyf))
    if s in ('SAL',):
        data = N.stack(data,axis=0)
        N.save(npyf,data)
    elif s in ('RANKHIST',):
        N.save(npyf,data)
    elif s in ("FSS","CASATI","closest",'ABCD','ROC'):
        with open(npyf, 'wb') as f:
            pickle.dump(data, f)
    else:
        raise Exception("No save logic specified for {}".format(s))
    print("Saved data to",npyf)
    del data
    return missingdata

def genit(stats=[],models=[],newstype=newstype,testplots=False):
    """
    m = model
    s = stat
    c = case date
    i = init times
    newstype = member (ignored if model is not NEWSe)
    m,s,c,i,newstype = it
    """
    assert isinstance(stats,(list,tuple))
    assert isinstance(models,(list,tuple))

    if testplots:
        yield ('NEWSe','SAL',datetime.datetime(2016,5,17,0,0,0),'0200','ENS_MEM_10')
    else:
        for m in models:
            if m != 'NEWSe':
                enslist = ('',)
            elif newstype == 'all':
                enslist = ensnames
            elif newstype == 'mean':
                enslist = ('mean',)
            else:
                raise Exception

            for s in stats:
                for c in cases:
                    for i in NEWS_CAT[c].keys():
                        for e in enslist:
                            yield (m,s,c,i,e)


itr = genit(stats=stats,models=models,newstype=newstype)

if serial:
    itr = genit(stats=stats,models=models,newstype=newstype,testplots=testplots)
    for r in range(10):
        compute_stats(next(itr))
    assert True == False
else:
    pass

from multiprocessing import Process, Lock, Pool
if __name__ == '__main__':
    print("Looping over cases.")
    pool = Pool(ncpu)
    # Pool.starmap()?
    # results = pool.map_async(compute_stats,itr,)#chunksize=1)
    results = pool.map(compute_stats,itr,chunksize=1)
    print("DONE #1")
    pool.close()
    print("DONE #2")
    pool.join()
    print("DONE #3")

writeout= True
if writeout:
    misf = open('/home/john.lawson/HRRRvNEWSe/missingdata.txt','w')
    for poolres in results:
        for res in poolres:
            if isinstance(res,list):
                for el in res:
                    misf.write('{0},'.format(el))
                misf.write('\n')
    misf.close()

# for c in cases:
# para_stats(c)

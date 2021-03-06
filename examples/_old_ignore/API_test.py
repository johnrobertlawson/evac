import sys
import numpy as N
import pdb
import matplotlib as M
import os
M.use('gtkagg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.append('/home/jrlawson/gitprojects/WEM')

from WEM.postWRF.postWRF import WRFEnviron

p = WRFEnviron()

outdir = ('/home/jrlawson/public_html/bowecho/'
            '20110419/GEFSR2/p04/WDM6_Hail/test')
ncdir = ('/chinook2/jrlawson/bowecho/'
            '20110419/GEFSR2/p04/WDM6_Hail')

plot2D = 1
streamlines = 0

# Create list of times
itime = (2013,8,16,1,0,0)
ftime = (2013,8,16,5,0,0)
hourly = 1
times = p.generate_times(itime,ftime,hourly*60*60)
levels = (2000,)
vrbls = ('cref',)

# import pdb; pdb.set_trace()
# Plot standard 2D plots
if plot2D:
    for t in times:
        for lv in levels:
            for v in vrbls:
                p.plot2D(v,t,lv,ncdir,outdir)

# Streamline 2D plots
if streamlines:
     for t in times:
         for lv in levels:
            p.plot_streamlines(t,lv,ncdir,outdir,smooth=5)

# Skew Ts
# skewT_time = (2013,8,16,3,0,0)
# skewT_latlon = (35.2435,-97.4708)
# p.plot_skewT(skewT_time,skewT_latlon)

# DKE
# p.compute_diff_energy('sum_z','total',path_to_wrfouts,times,
                    #   d_save=picklefolder, d_return=0,d_fname=pfname)

# latA, lonA, latB, lonB = (40.1,-89.0,38.1,-86.0)
# xstime = (2011,4,20,3,0,0)
# clvs = N.arange(-25,27.5,2.5)
#clvs = N.arange(0.0001,0.001,0.0001)
# p.plot_xs('parawind',xstime,latA=latA,lonA=lonA,latB=latB,lonB=lonB,
            # wrf_sd=wrf_sd, out_sd=out_sd,clvs=clvs,ztop=9)

# for t in times:
    # fig = plt.figure(figsize=(8,6))
    # gs = M.gridspec.GridSpec(1,2,width_ratios=[1,3])
    # ax0 = plt.subplot(gs[0])
    # ax1 = plt.subplot(gs[1])
    # cf0, cf1 = p.cold_pool_strength(besttime,wrf_sd=wrf_sd,out_sd=out_sd,swath_width=150,twoplot=0,fig=fig,axes=(ax0,ax1))

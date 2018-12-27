import os
import pdb
import pickle

import numpy as N
import matplotlib as M
M.use("agg")
import matplotlib.pyplot as plt

from evac.datafiles.ensemble import Ensemble
from evac.stats.fi import FI

compute = 0
plot = 1
fi = None
dx = 3

if compute:
    rootdir = '/home/nothijngrad/data'
    ncf = 'wrfout.nc'

    E = Ensemble(rootdir=rootdir,ncf=ncf)
    fcst = E.get('REFL_comp',fcsthr=1)[:,0,0,:,:]

    # Verify the ensemble against the ensemble mean, just for ease
    fracign = FI(xa=N.mean(fcst[:,N.newaxis,:,:],axis=0),xfs=fcst[:,N.newaxis,:,:],
                    # thresholds=N.arange(5,80,5),
                    thresholds = [10,30,50],
                    # neighborhoods=[3,10],
                    ncpus=6,
                    )
    fi = fracign.results
    with open("saved_fi.pickle", 'wb') as f:
        pickle.dump(obj=fi, file=f)

if plot:
    toplot = {}
    if fi is None:
        with open("saved_fi.pickle", 'rb') as f:
            fi = pickle.load(f)
    dBZs = sorted(fi.keys())
    neighs = sorted(fi[dBZs[0]].keys())

    for dBZ in dBZs:
        toplot[dBZ] = []
        for neigh in neighs:
            toplot[dBZ].append(fi[dBZ][neigh][1][0]['FI'])

    fig, ax = plt.subplots(1)
    for dBZ in dBZs:
        ax.plot(neighs,toplot[dBZ],label="{} dBZ".format(dBZ))
    ax.legend()
    ax.set_xticks(neighs)
    ax.set_xticklabels((0.5*dx*N.array(neighs)))
    ax.set_xlabel("Neighbourhood radius (km)")
    # ax.set_ylabel("FI and components (bits)")
    ax.set_ylabel("Fractional Ignorance (bits)")
    fig.savefig("fi_plot.png")



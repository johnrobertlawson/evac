import os
import pdb

import numpy as N
import matplotlib as M

from evac.plot.figure import Figure

class SALGraph(Figure):
    def __init__(self,fpath,medians=True,cb=True,
                    quartiles=True,figsize=(6,6)):
        super().__init__(fpath=fpath,figsize=figsize)

        self.medians = medians
        self.cb = cb
        self.quartiles = quartiles

        # Lists of all plotted values for later computation
        self.SS = []
        self.AA = []
        self.LL = []

        self.ax.axhline(0, color='k')
        self.ax.axvline(0, color='k')
        self.ax.set_xlim([-2,2])
        self.ax.set_ylim([-2,2])
        self.ax.set_xlabel("Structural component")
        self.ax.set_ylabel("Amplitude component")

        # self.ax.set_axisbgcolor('lightgrey')
        self.ax.set_facecolor('lightgrey')

    def plot(self,SAL=None,S=None,A=None,L=None,save=False,
                check_valid=True,marker='o',ignore_nans=True):
        """ Plot one scatter point for a given SAL object.

        Args:
            SAL : An instance of SAL or ESAL object.
        """
        if SAL is None:
            class SAL: pass
            SAL.S = S
            SAL.A = A
            SAL.L = L

        if check_valid:
            assert abs(SAL.S) <= 2.0
            assert abs(SAL.A) <= 2.0
            assert 0.0 <= SAL.L <= 2.0

        if ignore_nans:
            if N.isnan(SAL.S) or N.isnan(SAL.A) or N.isnan(SAL.L):
                print("Not plotting due to NaN(s).")
                return False
        
        cmap = M.cm.nipy_spectral_r
        self.sc = self.ax.scatter(SAL.S, SAL.A, c=SAL.L, vmin=0, vmax=2,
                    s=25,cmap=cmap,alpha=0.9,edgecolor='k',
                    linewidth=0.15,zorder=500,marker=marker)


        # Append info for medians computation
        self.SS.append(SAL.S)
        self.AA.append(SAL.A)
        self.LL.append(SAL.L)

        if save:
            self.save()

        return True

    def lists_to_arrays(self):
        self.SS = N.array(self.SS)
        self.AA = N.array(self.AA)
        return

    def plot_medians(self):
        medS = N.median(self.SS)
        medA = N.median(self.AA)

        medkwargs = dict(color='black',
                        zorder=300,
                        linewidth=0.8,
                        linestyle=":")
        self.ax.axvline(medS,**medkwargs)
        self.ax.axhline(medA,**medkwargs)

    def plot_quartiles(self):
        lbS = N.percentile(self.SS,25)
        ubS = N.percentile(self.SS,75)
        lbA = N.percentile(self.AA,25)
        ubA = N.percentile(self.AA,75)
        
        width = ubS - lbS
        height = ubA - lbA

        self.ax.add_patch(M.patches.Rectangle((lbS,lbA),width,height,
                            facecolor='white',alpha=1.0,linewidth=0.5,
                            zorder=100))

    def plot_colorbar(self):
        cbax = self.fig.add_axes([0.17,0.17,0.22,0.05])
        cblab = N.array([0.0,0.25,0.5,0.75,1.0,2.0])
        cb = self.fig.colorbar(self.sc,cax=cbax,
                            ticks=cblab,orientation="horizontal")
        cb.set_label("Location component",labelpad=-38)
        cbax.set_xticklabels(cblab)
        return

    def __exit__(self):
        """ Override Figure context exit by
        plotting medians, if requested.
        """
        if self.cb:
            self.plot_colorbar()
        if self.medians:
            self.lists_to_arrays()
            self.plot_medians()
        if self.quartiles:
            self.plot_quartiles()
        super().__exit__()
        return

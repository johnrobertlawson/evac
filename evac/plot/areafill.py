import pdb
import os

import matplotlib as M
import numpy as N

from evac.plot.figure import Figure

class AreaFill(Figure):
    def __init__(self,fpath,data,Ls,thresholds,ax=None,fig=None):
        """ Plot neighbourhood scores like Casati MSE.

        Args:
            data: numpy array 2d.
            Ls (tuple,list): wavelength multiples 
            thresholds (tuple,list): thresholds used
        """
        super().__init__(fpath=fpath,ax=ax,fig=fig,)
        self.data = data

        self.thresholds = sorted(thresholds)
        self.Ls = sorted(Ls)


    def plot(self,cb=True,plotkwargs=None,L_multiplier=1):
        """
        Handy plotkwargs include:
            * levels, to see contouring levels
            * cmap, for colormap
        """

        if plotkwargs is None:
            plotkwargs = dict(cmap=M.cm.plasma)

        xx = N.arange(len(self.thresholds))
        yy = N.arange(len(self.Ls))
        xxtx = xx
        yytx = yy
        xxtl = self.thresholds
        # if L_multiplier != 1:
            # yytl = [L_multiplier*L for L in self.Ls]
        # else:
        yytl = self.Ls

        cf = self.ax.contourf(xx,yy,self.data,
                        **plotkwargs)
        self.ax.set_xticks(xxtx)
        self.ax.set_xticklabels([str(x) for x in xxtl])
        self.ax.set_yticks(yytx)
        self.ax.set_yticklabels([str(L_multiplier*(2**y)) for y in yytl])
        self.ax.set_xlabel("Thresholds (mm/hr)")
        self.ax.set_ylabel("Spatial scale (km)")

        self.fig.colorbar(cf)
        return

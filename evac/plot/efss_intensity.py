import pdb
import os

import numpy as N
import matplotlib.pyplot as plt

from evac.plot.figure import Figure

class EFSS_Intensity(Figure):
    """ Heatmap-style plots.

    Note:
        * User should provide either efss_inst, or the rest.

    Args:
        fpath (str): absolute path to figure output
        efss_inst: an instance of EFSS class
        arr3d: a Numpy array of data, spatial x temporal x threshold
        spatial_windows (list, tuple): in same order as arr3d axis
        temporal_windows (list, tuple): ditto
        threshs (list, tuple): ditto
    """
    def __init__(self,fpath,efss_inst=None,arr3d=None,
                    threshs=None,spatial_windows=None,
                    temporal_windows=None,
                    figsize=(8,6),mplkwargs=None,**kwargs):
        
        super().__init__(fpath=fpath,figsize=figsize,mplkwargs=mplkwargs,**kwargs)

        if efss_inst is not None:
            self.threshs = sorted(efss_inst.threshs)
            self.spatial_windows = sorted(efss_inst.spatial_ns)
            self.temporal_windows = sorted(efss_inst.temporal_ms)
            self.data = efss_inst.eFSS

        else:
            raise NotImplementedError

        self.nthreshs = len(self.threshs)
        self.ntempwind = len(self.temporal_windows)
        self.nspatwind = len(self.spatial_windows)

        self.matrix = self.generate_matrix()

    def generate_matrix(self):
        nx = self.nthreshs * self.ntempwind
        ny = self.nspatwind
        matrix = N.zeros([nx,ny])

        for x1idx,th in enumerate(self.threshs):
            for x2idx,m in enumerate(self.temporal_windows):
                xidx = (x1idx * self.ntempwind) + x2idx
                for yidx, n in enumerate(self.spatial_windows):
                    matrix[xidx,yidx] = self.data[n][m][th]
        return matrix.T

    def plot_intensity(self,plotkwargs=None,invert_y=False,
                        tickx_top=True,
                        xticklabels=False,x2ticklabels=False,
                        yticklabels=False,):#**kwargs):
        """

        Args:

        Note:
            * cmap can be passed for colormap
            * colornorm (bool) can be passed. If True, the colour scheme 
                is normalised so the highest FSS value is max. 
                Otherwise, between 0 and 1.
        """
        if plotkwargs is None:
            plotkwargs = dict()
        heatmap = self.ax.pcolormesh(self.matrix,**plotkwargs)
        self.ax.set_aspect('equal')

        
        # Put ticks into centre of each row/column
        self.ax.set_yticks(N.arange(self.matrix.shape[0]) + 0.5, minor=False)
        self.ax.set_xticks(N.arange(self.matrix.shape[1]
                    )[::self.ntempwind] + (self.ntempwind/2), minor=False)
        if invert_y:
            self.ax.invert_yaxis()
        # if tickx_top:
            # self.ax.xaxis.tick_top()

        if xticklabels is not False:
            if xticklabels == 'auto':
                xticklabels = self.threshs
            # if isinstance(xticklabels[0],str):
            self.ax.set_xticklabels(xticklabels, minor=False)
        if yticklabels is not False:
            if yticklabels == 'auto':
                yticklabels = self.spatial_windows
            self.ax.set_yticklabels(yticklabels,minor=False)


        # Make grid prettier
        self.ax.grid(False)

        # for tk in self.ax.xaxis.get_major_ticks():
            # tk.tick10n = False
            # tk.tick20n = False
        # for tk in self.ax_top.xaxis.get_major_ticks():
            # tk.tick10n = False
            # tk.tick20n = False
        # for t in self.ax.yaxis.get_major_ticks():
            # tk.tick10n = False
            # tk.tick20n = False

        self.ax.tick_params(axis='both',which='both',bottom=False,
                    top=False,left=False,right=False,
                    labeltop=False,labelbottom=True)

        # Upper ticks
        self.ax_top = self.ax.twiny()
        self.ax_top.set_xticks(N.arange(self.matrix.shape[1]) + 0.5, minor=False)
        self.ax_top.set_xticklabels(x2ticklabels, minor=False)
        self.ax_top.tick_params(axis='both',which='both',bottom=False,
                    top=False,left=False,right=False,
                    labeltop=True,labelbottom=False)

        #vlines = N.arange(self.matrix.shape[1])[::self.ntempwind] #+ self.ntempwind
        #plt.vlines(x=vlines, ymin=0.0, ymax=self.matrix.shape[0])

        self.fig.subplots_adjust(right=0.85)
        cbar_ax = self.fig.add_axes([0.9,0.1,0.02,0.8])
        self.fig.colorbar(heatmap,cax=cbar_ax,orientation='vertical')
        # if self.clskwargs['save']:
            # self.save()
        return
        
    def save(self,tight=True,close=True):
        """ Overriden for optimal settings.
        """
        import evac.utils as utils

        utils.trycreate(self.fpath)
        self.fig.savefig(self.fpath)
        print("Saving figure to {}".format(self.fpath))
        plt.close(self.fig)
        return
        

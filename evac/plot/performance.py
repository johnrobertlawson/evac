import pdb
import os

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt

import evac.utils as utils
from evac.utils.evac_numbers import is_number, has_integers, is_integer, is_array
from evac.plot.figure import Figure

class Performance(Figure):
    """ Performance diagram, from Roebber 2009 WAF.

    Some checking with WFRT's verif package (Github)?

    If passing in stats, they should be:

    .. math::
        \\textrm{POD}     &=   A/A+C

        \\textrm{FAR}     &=   B/A+B

        \\textrm{BIAS}    &=   (A+B)/(A+C)

        \\textrm{CSI}     &=   A/(A+B+C)

        \\textrm{SR}      &=   1-\\textrm{FAR}

    Example:
        To maintain flexibility, the user should do the following
        outside the class once initialised::

            Performance.plot_data(...) # As many times as needed
            Performance.save() # When finished

    Some plotting settings, like bias contours, can be accessed::

        Performance.contours['bias'] = [0.5,1,1.5]

    Args:
        fpath     :   directory for saving figure
        fname       :   file name
        legendkwargs:   key-words for ax.legend().
    """
    def __init__(self,fpath=None,outdir=None,fname=None,legendkwargs=None,
                    legend=False):
        if fpath is None:
            fpath = os.path.join(outdir,fname)
        super().__init__(fpath=fpath)
        if legendkwargs is None:
            legendkwargs = dict(facecolor = 'white',
                                framealpha = 1.0)
        self.lk = legendkwargs
        self.legend = legend
        self.alpha = 0.2

        # Plotting settings
        self.contours = {}
        self.contours['bias'] = N.array([0.25,0.5,0.75,1.0,1.25,1.5,2.5,5,10])
        self.contours['csi'] = N.arange(0.1,1.1,0.1)
        # self.xticks = 201
        self.xticks = 601
        self.yticks = self.xticks
        self.sr_x = N.linspace(0,1,self.xticks)
        self.pod_y = N.linspace(0,1,self.yticks)

        self.create_axis()
        self.plot_bias_lines()
        self.plot_csi_lines()

        print("Ready to receive plotting data via plot_data().")

    def create_axis(self):
        self.fig, self.ax = plt.subplots(1,figsize=(5,5))
        self.ax.grid(False)

        self.ax.set_xlim([0,1])
        xlab = N.arange(0,1.1,0.1)
        self.ax.set_xticks(xlab)
        self.ax.set_xticklabels(xlab)
        self.ax.set_xlabel("Success Ratio")

        self.ax.set_ylim([0,1])
        ylab = xlab
        self.ax.set_yticks(ylab)
        self.ax.set_yticklabels(ylab)
        self.ax.set_ylabel("Probability of Detection")

        self.ax.set_aspect('equal')

    def plot_bias_lines(self):
        bias_arr = self.compute_bias_lines()
        for bn, b in enumerate(self.contours['bias']):
            # self.ax.plot(bias_arr[:,0,bn],bias_arr[:,1,bn],
            # Add a non-transparent line for clarity
            if 0.9 < b < 1.1:
                alpha = 1.0
            else:
                alpha = self.alpha
            self.ax.plot(bias_arr[:,1,bn],self.pod_y,
                            color='red',lw=1,linestyle='--',
                            alpha=alpha)
            bstr = '{:1.1f}'.format(b)
            # pdb.set_trace()
            if b < 1:
                xpos = 1.00
                ypos = b
            else:
                xpos = 1.00/b
                ypos = 1.00
            self.ax.annotate(bstr,xy=(xpos,ypos),xycoords='data',
                            # bbox=dict(fc='red'),color='white',
                            fontsize=7,color='red',
                            xytext=(2,3),textcoords='offset points',
                            )

    def plot_csi_lines(self):
        xpos = 0.945
        csi_arr = self.compute_csi_lines()
        # mid = N.floor(csi_arr.shape[0]/2)
        # mid = N.int(csi_arr.shape[0] - (0.975*csi_arr.shape[0]))
        nc = len(self.contours['csi'])
        for cn, c in enumerate(self.contours['csi']):
            # csi_c = self.ax.plot(csi_arr[:,0,cn],csi_arr[:,1,cn],

            self.ax.plot(self.sr_x,csi_arr[:,1,cn],
                                color='blue',lw=1,alpha=self.alpha)


            cstr = '{:1.1f}'.format(c)
            # self.ax.clabel(csi_c[1],fontsize=8,inline=True,fmt='%1.1f')
            # pdb.set_trace()
            # if not N.isnan(csi_arr[mid,1,cn]):
            if True:
                # negx = int(negx_factor*self.xticks)
                # yidx = N.abs(csi_arr[:,0,cn]-(1-negx_factor)).argmin()
                ypos = 1/((1/c) + 1 - (1/xpos))
                # self.fig.text(N.ones(nc)*0.93,csi_arr[mid,1,cn],
                                # cstr,color='blue',)#facecolor='white')
                self.ax.annotate(cstr,xy=(xpos-0.01,ypos),
                # self.ax.annotate(cstr,xy=(1-negx_factor,csi_arr[mid,1,cn]),
                                    xycoords='data', fontsize=7,
                                    #color='white', bbox=dict(fc='blue'),
                                    bbox=dict(fc='white',color='white',pad=0),
                                    xytext=(3,3),textcoords='offset points',
                                    color='blue')
        # pdb.set_trace()


    def compute_bias_lines(self):
        bias_arr = N.zeros([self.sr_x.size,2,self.contours['bias'].size])
        for bn, b in enumerate(self.contours['bias']):
            bias_arr[:,0,bn] = self.pod_y/b
            bias_arr[:,1,bn] = self.sr_x/b
        bias_arr[bias_arr>1.0] = 1.0
        # pdb.set_trace()
        return bias_arr

    def compute_csi_lines(self):
        # might need to switch around the axes
        csi_arr = N.zeros([self.sr_x.size,2,self.contours['csi'].size])
        for cn, c in enumerate(self.contours['csi']):
            # x coordinates for this CSI contour
            sr_x_c = 1/((1/c) + 1 - (1/self.pod_y))
            # and y coordinates
            pod_y_c = 1/((1/c) + 1 - (1/self.sr_x))

            # assign to array
            csi_arr[:,0,cn] = sr_x_c
            csi_arr[:,1,cn] = pod_y_c

        # pdb.set_trace()
        csi_arr[csi_arr < 0] = N.nan
        csi_arr[csi_arr > 1] = N.nan
        return csi_arr

    def plot_data(self,pod=None,far=None,sr=None,abcd=None,a=None,
                  b=None,c=None,d=None,label=None,
                    plotkwargs=None):
        """ Plot a datum point onto graph.

        Example plotkwargs include:
            s   :   size
            marker  :   marker, e.g.: v X + * D
            c   :   colour
            alpha   :   0-1
            edgecolors  :   colour


        Args:
            pod     :   probability of detection (one point)
            sr      :   success rate (one point)
            far     :   false alarm rate (one point)
            abcd    :   2x2 contingency table (to be implemented)
            a,b,c,d :   integers of contingency table
            plotkwargs  :   plotting system for matplotlib.
        """
        if plotkwargs is None:
            plotkwargs = dict()
        def get_scores(a,b,c,d):
            DS = DetScores(a=a,b=b,c=c,d=d)
            pod = DS.get('POD')
            sr = DS.get('SR')
            return pod,sr

        if is_array(abcd) and has_integers(abcd):
            # Compute scores here.
            a,b,c,d = abcd.flatten()
            pod,sr = get_scores(a,b,c,d)

        elif is_number(sr) or is_number(far):
                # x axis is present
                if not is_number(pod):
                    # y axis is not
                    raise Exception
                if sr is None:
                    sr = 1.0 - far
        else:
            for x in (a,b,c,d):
                if not is_integer(x):
                    raise Exception
            pod,sr = get_scores(a,b,c,d)

        self.ax.scatter(sr,pod,label=label,zorder=500,**plotkwargs)
        # print("Ready to save via save(), or plot more with plot_data().")
        from matplotlib.ticker import FormatStrFormatter
        self.ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    def save(self,close=True):
        """ Extends Figure.save() by completing tasks before saving.
        """
        if self.legend:
            utils.wowprint("**Drawing legend.**")
            self.ax.legend(**self.lk)
        # self.fig.tight_layout()
        # self.fig.subplots_adjust(top=0.95,right=0.95)
        # super(Performance,self).save(tight=False)
        super().save()
        if close:
            plt.close(self.fig)

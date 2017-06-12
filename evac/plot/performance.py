import pdb
import os

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt

from evac.numbers import is_number, has_integers, is_integer, is_array

class Performance(Figure):
    def __init__(self,outpath,fname):
        """ Performance diagram, from Roebber 2009 WAF.

        Some checking with WFRT's verif package (Github)

        To maintain flexibility, the user should do the following
        outside the class once initialised:

        Performance.plot_data(...) # As many times as needed
        Performance.save() # When finished

        If passing in stats, they should be:

        POD     =   A/A+C
        FAR     =   B/A+B
        BIAS    =   (A+B)/(A+C)
        CSI     =   A/(A+B+C)
        SR      =   1/FAR ...
                =   (A+B)/B

        Some plotting settings, like bias contours, can be accessed:

        Performance.contours['bias'] = [0.5,1,1.5]

        Args:
            outpath     :   directory for saving figure
            fname       :   file name
        """
        self.outpath = outpath
        self.fname = fname

        # Plotting settings
        self.contours = {}
        self.contours['bias'] = N.array([0.25,0.5,0.75,1.0,1.25,1.5,2.5,5,10])
        self.contours['csi'] = N.arange(0.1,1.1,0.1)
        self.xticks = 100
        self.yticks = self.xticks
        self.sr_x = N.arange(0,(1+1/self.xticks),self.xticks)
        self.pod_y = N.arange(0,(1+1/self.yticks),self.yticks)

        self.create_axis()
        self.plot_bias_lines()
        self.plot_csi_lines()

        print("Ready to receive plotting data via plot_data().")

    def create_axis(self):
        self.fig, self.ax = plt.subplots(1,figsize=(5,5))
        self.ax.grid(False)

        self.ax.set_xlim([0,1])
        xlab = N.arange(0,1.05,0.05)
        self.ax.set_xticks(xlab)
        self.ax.set_xticklabels(xlab)
        self.ax.set_xlabel("Success Ratio")

        self.ax.set_ylim([0,1])
        ylab = xlab
        self.ax.set_yticks(ylab)
        self.ax.set_yticklabels(ylab)
        self.ax.set_ylabel("Probability of Detection")

    def plot_bias_lines(self):
        bias_arr = self.compute_bias_lines()
        self.ax.plot(bias_arr,color='red',lw=1,linestyle='--')
        for bn, b in enumerate(self.contours['bias']):
            bstr = '{:1.1f}'.format(b)
            self.fig.text(bias_arr[-1,0,bn],bias_arr[-1,1,bn],bstr)

    def plot_csi_lines(self):
        csi_arr = self.compute_csi_lines()
        csi_c = self.ax.plot(csi_arr,color='blue',lw=1,)
        self.fig.clabel(csi_c,fontsize=8,inline=True,fmt='%1.1f')

    def compute_bias_lines(self):
        bias_arr = N.zeros(self.sr_x.size,2,self.contours['bias'].size)
        for bn, b in enumerate(self.contours['bias']):
            bias_arr[:,0,bn] = self.pod_y/b
            bias_arr[:,1,bn] = self.sr_x/b
        return bias_arr

    def compute_csi_lines(self):
        # might need to switch around the axes
        csi_arr = N.zeros(self.sr_x.size,2,self.contours['csi'].size)
        for cn, c in enumerate(self.contours['csi']):
            # x coordinates for this CSI contour
            sr_x_c = 1/((1/c) + 1 + (1/self.pod_y))
            # and y coordinates
            pod_y_c = 1/((1/c) + 1 + (1/self.sr_x))

            # assign to array
            csi_arr[:,0,cn] = sr_x_c
            csi_arr[:,1,cn] = pod_y_c

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
        else:
            for x in (a,b,c,d):
                if not is_integer(x):
                    raise Exception
            pod,sr = get_scores(a,b,c,d)

        self.ax.scatter(sr,pod,label=label,**plotkwargs)
        print("Ready to save via save(), or plot more with plot_data().")

    def save(self,):
        """ Extends Figure.save() by completing tasks before saving.
        """
        self.ax.legend()
        super(Performance,self).save()

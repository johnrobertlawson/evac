"""
Top-down 2D plots, which are so common in meteorology
that they get their own file here.

Subclass of Figure.
"""

import pdb
import matplotlib as M
M.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as N
import collections
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

from wrfout import WRFOut
from defaults import Defaults
from figure import Figure
import WEM.utils as utils
from scales import Scales
import stats

class BirdsEye(Figure):
    def __init__(self,nc=False,ax=0,fig=0):
        super(BirdsEye,self).__init__(nc=nc,ax=ax,fig=fig)

    def get_plot_arguments(self,cmap=False,clvs=False):
        """
        Returns colourmap and contouring levels

        Options keyword arguments:
        clvs    :   manually override contour levels
        """
        # import pdb; pdb.set_trace()
        data = self.data.reshape((self.la_n,self.lo_n))

        # List of args and dictionary of kwargs
        plotargs = [self.x,self.y,data]
        plotkwargs = {}

        # if self.mplcommand == 'contour':
            # multiplier = S.get_multiplier(vrbl,lv)
        if clvs is not False:
            plotkwargs['levels'] = clvs

        if cmap is not False:
            # cmap = eval('M.cm.{0}'.format(cmap))
            plotkwargs['cmap'] = cmap
        return plotargs, plotkwargs

    # Old plot_data
    def plot2D(self,data,fname,outdir,plottype='contourf',
                    save=1,smooth=1,lats=False,lons=False,
                    clvs=False,cmap=False,title=False,cb=True,
                    locations=False,m=False,x=False,y=False,
                    Nlim=False,Elim=False,Slim=False,Wlim=False):

        """
        Generic method that plots any matrix of data on a map

        Inputs:
        data        :   2D matrix of data
        outdir      :   path to plots
        outf        :   filename for output (with or without .png)

        Optional:
        plottype    :   matplotlib function for plotting
        smooth      :   Gaussian smooth by this many grid spaces
        clvs        :   scale for contours
        title       :   title on plot
        save        :   whether to save to file

        :param locations:       Locations to plot on the basemap.
                                Format: locations = {'label':(lat,lon),etc}
        :type locations:        dict
        """
        # INITIALISE
        self.data = data
        if x is False and y is False and m is False:
            if not Nlim:
                self.bmap,self.x,self.y = self.basemap_setup(smooth=smooth,lats=lats,
                                                    lons=lons,)#ax=self.ax)
            else:
                self.bmap,self.x,self.y = self.basemap_setup(smooth=smooth,lats=lats,
                                                    lons=lons,Nlim=Nlim,Elim=Elim,
                                                    Slim=Slim,Wlim=Wlim)

        else:
            self.bmap = m
            self.x = x
            self.y = y

        self.la_n = self.data.shape[-2]
        self.lo_n = self.data.shape[-1]

        plotargs, plotkwargs = self.get_plot_arguments(clvs=clvs,cmap=cmap)

        # import pdb; pdb.set_trace()
        if plottype == 'contour':
            f1 = self.bmap.contour(*plotargs,**plotkwargs)
        elif plottype == 'contourf':
            f1 = self.bmap.contourf(*plotargs,**plotkwargs)
        elif plottype == 'pcolor':
            f1 = self.bmap.pcolor(*plotargs,**plotkwargs)
        elif plottype == 'pcolormesh':
            f1 = self.bmap.pcolormesh(*plotargs,**plotkwargs)
        elif plottype == 'scatter':
            f1 = self.bmap.scatter(*plotargs,**plotkwargs)
        else:
            print("Specify correct plot type.")
            raise Exception

        if isinstance(locations,dict):
            for k,v in locations.iteritems():
                if isinstance(v,tuple) and len(v) == 2:
                    xpt, ypt = self.bmap(v[1],v[0])
                    # bbox_style = {'boxstyle':'square','fc':'white','alpha':0.5}
                    self.bmap.plot(xpt,ypt,'ko',markersize=3,zorder=100)
                    self.ax.text(xpt,ypt,k,ha='left',fontsize=7)
                    # self.ax.text(xpt-15000,ypt,k,bbox=bbox_style,ha='left',fontsize=7)
                else:
                    print("Not a valid location argument.")
                    raise Exception

        if isinstance(title,basestring):
            plt.title(title)
        if save:
            self.save(outdir,fname)
            plt.close(self.fig)
        if cb:
            self.fig.colorbar(f1,orientation='vertical')
        else:
            return f1

    def plot_streamlines(self,U,V,outdir,fname,lats=False,lons=False,smooth=1,
                            title=False,lw_speed=False):
        """
        Plot streamlines.

        U       :   U-component of wind (nx x ny)
        V       :   V-component of wind (same dimensions)

        lw_speed    :   linewidth is proportional to wind speed
        """
        m,x,y = self.basemap_setup()

        if lw_speed:
            wind = N.sqrt(U**2 + V**2)
            lw = 5*wind/wind.max()
        else:
            lw = 1

        if smooth>1:
            U = stats.gauss_smooth(U,smooth)
            V = stats.gauss_smooth(V,smooth)

        m.streamplot(x[self.W.x_dim/2,:],y[:,self.W.y_dim/2],U,V,
                        density=1.8,linewidth=lw,color='k',arrowsize=3)

        if isinstance(title,basestring):
            self.ax.set_title(title)

        self.save(outdir,fname)

    def spaghetti(self,t,lv,va,contour,wrfouts,outpath,da=0,dom=0):
        """
        wrfouts     :   list of wrfout files

        Only change dom if there are multiple domains.
        """
        m,x,y = self.basemap_setup()

        time_idx = self.W.get_time_idx(t)

        colours = utils.generate_colours(M,len(wrfouts))

        # import pdb; pdb.set_trace()
        if lv==2000:
            lv_idx = None
        else:
            print("Only support surface right now")
            raise Exception

        lat_sl, lon_sl = self.get_limited_domain(da)

        slices = {'t': time_idx, 'lv': lv_idx, 'la': lat_sl, 'lo': lon_sl}

        # self.ax.set_color_cycle(colours)
        ctlist = []
        for n,wrfout in enumerate(wrfouts):
            self.W = WRFOut(wrfout)
            data = self.W.get(va,slices)[0,...]
            # m.contour(x,y,data,levels=[contour,])
            ct = m.contour(x,y,data,colors=[colours[n],],levels=[contour,],label=wrfout.split('/')[-2])
            print("Plotting contour level {0} for {1} from file \n {2}".format(
                            contour,va,wrfout))
            # ctlist.append(ct)
            # self.ax.legend()

        # labels = [w.split('/')[-2] for w in wrfouts]
        # print labels
        # self.fig.legend(handles=ctlist)
        # plt.legend(handles=ctlist,labels=labels)
        #labels,ncol=3, loc=3,
        #                bbox_to_anchor=[0.5,1.5])

        datestr = utils.string_from_time('output',t,tupleformat=0)
        lv_na = utils.get_level_naming(va,lv)
        naming = ['spaghetti',va,lv_na,datestr]
        if dom:
            naming.append(dom)
        fname = self.create_fname(*naming)
        self.save(outpath,fname)


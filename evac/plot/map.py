import collections
import os
import pdb

import numpy as N
from mpl_toolkits.basemap import Basemap

from evac.plot.figure import Figure
from evac.utils.reproject_tools import WRF_native_grid, VerifGrid
from evac.datafiles.wrfout import WRFOut

class DomainGrid:
    def __init__(self,inherit=None):
        """
        Create a common class for verification domains, WRF out domains, etc
        """
        if isinstance(inherit,WRFOut):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,WRF_native_grid):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,VerifGrid):
            self.__dict__ = inherit.__dict__.copy()
        elif isinstance(inherit,str):
            inherit = WRFOut(inherit)
            self.__dict__ = inherit.__dict__.copy()

            # self.lat_0 = inherit.nc.CEN_LAT
            # self.lon_0 = inherit.nc.CEN_LON
        else:
            raise Exception("Inherited object is {}".format(type(inherit)))

class Map(Figure):
    def __init__(self,):
        super().__init__()

    def make_domain_obj(self,input_objs):
        """
        Turn WRFOut files or grids into domains.
        """
        # pdb.set_trace()
        DGs = []

        for obj in input_objs:
            DGs.append(DomainGrid(inherit=obj))
        
        self.domains = DGs
        return


    def plot_domains(self,domains,outdir,labels=None,latlons='auto',fname='domains.png',
                        colour=0,res='h',proj='auto',pad=0.2):
        """
        domains     :   list of WRFOut objects for each domain
                        Largest domain should be index 0, if using auto settings
        latlons     :   dictionary of Nlim,Elim,Slim,Wlim
                        for plot
        pad         :   multiply factor (0-1) for padding if latlons is 'auto'
        """
        if labels is not None:
            assert len(labels) == len(domains)
        else:
            labels = N.arange(len(domains))

        # Generate domain objects
        self.make_domain_obj(domains)

        # largest domain
        LD = self.domains[0]

        if proj == 'auto':
            proj = 'merc'
            # proj = W.proj
        elif proj is None:
            proj = 'merc'

        if latlons == 'auto':
            Nlim = LD.Nlim + abs(LD.Nlim-LD.Slim)*pad
            Elim = LD.Elim + abs(LD.Elim-LD.Wlim)*pad
            Slim = LD.Slim - abs(LD.Nlim-LD.Slim)*pad
            Wlim = LD.Wlim - abs(LD.Elim-LD.Wlim)*pad
            try:
                lat_0 = LD.lat_0
                lon_0 = LD.lon_0
            except AttributeError:
                lat_0 = LD.cen_lat
                lon_0 = LD.cen_lon

        else:
            Nlim = latlons['Nlim']
            Elim = latlons['Elim']
            Slim = latlons['Slim']
            Wlim = latlons['Wlim']
            lat_0 = latlons['lat_0']
            lon_0 = latlons['lon_0']

        self.m = Basemap(
            projection=proj,
            llcrnrlon=Wlim,llcrnrlat=Slim,
            urcrnrlon=Elim,urcrnrlat=Nlim,
            lat_0=lat_0,lon_0=lon_0,
            resolution=res,area_thresh=500,
            ax=self.ax
            )

        self.m.drawcoastlines()
        self.m.drawstates()
        self.m.drawcountries()

        if not isinstance(colour,collections.Sequence):
            colours = ['k',] * len(self.domains)
        else:
            colours = colour
        # Get corners of each domain
        for gridlabel,dom,colour in zip(labels,self.domains,colours):
            print(("Plotting domain {0} for {1}".format(gridlabel,dom)))
            x,y = self.m(dom.lons,dom.lats)
            xl = len(x[0,:])
            midpt = len(y[0,:])//2         
            self.ax.annotate(gridlabel,color=colour,fontsize=10,xy=(x[0,-int(0.12*xl)],y[0,midpt]),
                         bbox=dict(fc='white'),alpha=1,va='center',ha='left')    
            self.m.plot(x[0,:],y[0,:],colour,lw=2)
            self.m.plot(x[:,0],y[:,0],colour,lw=2) 
            self.m.plot(x[len(y)-1,:],y[len(y)-1,:],colour,lw=2)     
            self.m.plot(x[:,len(x)-1],y[:,len(x)-1],colour,lw=2)    

        fpath = os.path.join(outdir,fname)
        self.fig.savefig(fpath)

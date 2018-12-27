import pdb
import os

from mpl_toolkits.basemap import Basemap

import evac.utils as utils
from evac.plot.figure import Figure

class Domains(Figure):
    def __init__(self,fpath,*domdics,pad=0.07,figsize=None):
        """
        Args:
            domdics: pass in dictionaries for plotting.
            pad: if "auto", pad 10% of maxlat/minlat distance (and ditto for
                lons) outside of the domains.

        Note:
            Each dictionary should be done as follows:
            dom1 = {
                    "name":"fcst_domain",   # <-- for annotating or legend
                    "label":"fcst"          # <-- for labelling. Use name if missing.
                    "color":"black",        # <-- matplotlib colours
                    "lats":obj,          # <-- 2D numpy array of lats
                    "lons":obj,          # <-- 2D numpy array of lons
                    }
        """
        super().__init__(fpath=fpath,figsize=figsize)

        self.pad = pad
        self.domains = {}
        for dd in domdics:
            self.domains[dd['name']] = dd

    def add_domain(self,name,label,lats,lons,color="black"):
        self.domains[name] = {'name':name,
                                'label':label,
                                'color':color,
                                'lats':lats,
                                'lons':lons,}
        return

    def create_bmap(self,):
        maxlat = -9999
        maxlon = -9999
        minlat = 9999
        minlon = 9999

        for d,dd in self.domains.items():
            maxlat = max(maxlat,dd['lats'].max())
            maxlon = max(maxlon,dd['lons'].max())
            minlat = min(minlat,dd['lats'].min())
            minlon = min(minlon,dd['lons'].min())

        lat_offset = self.pad*abs(maxlat-minlat)
        lon_offset = self.pad*abs(maxlon-minlon)

        # pdb.set_trace()
        bmap = Basemap(
                urcrnrlat=maxlat+lat_offset,
                urcrnrlon=maxlon+lon_offset,
                llcrnrlat=minlat-lat_offset,
                llcrnrlon=minlon-lon_offset,
                ax=self.ax,
                resolution='h',
                projection='merc',
                area_thresh = 10000,
                lat_1=45.0,lat_2=55.0,lat_0=50.0,lon_0=-107.0)

        bmap.drawcoastlines()
        bmap.drawmapboundary(fill_color='lightgray')
        bmap.fillcontinents(color="gray",lake_color="lightgray")
        bmap.drawstates()
        bmap.drawcountries()

        return bmap

    def plot_domains(self, fontcolor='k'):
        """ Plot domain boxes for all in dictionary.
        """
        bmap = self.create_bmap()
        for d, dd in self.domains.items():
            self._outline_domain(bmap=bmap,label=dd['label'],lats=dd['lats'],
                                    lons=dd['lons'],color=dd['color'],
                                    fontcolor=fontcolor)
        return
        
    def _outline_domain(self,bmap,label,lats,lons,color='k',fontcolor='k'):
        # for gridlabel,dom,colour in zip(labels,self.domains,colours):
        x,y = bmap(lons,lats)
        xl = len(x[0,:])
        midpt = len(y[0,:])//2
        if label is not None:
            self.ax.annotate(label, color=fontcolor,fontsize=11, 
                        #color=color,
                        xy=(x[0,-int(0.12*xl)],y[0,midpt]),
                        bbox=dict(fc='white'),alpha=1,va='center',ha='left')
        bmap.plot(x[0,:],y[0,:],color,lw=2)
        bmap.plot(x[:,0],y[:,0],color,lw=2)
        bmap.plot(x[len(y)-1,:],y[len(y)-1,:],color,lw=2)
        bmap.plot(x[:,len(x)-1],y[:,len(x)-1],color,lw=2)
        return

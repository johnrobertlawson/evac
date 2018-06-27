import pdb
import collections
import os

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as M
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as N

from evac.datafiles.wrfout import WRFOut
from evac.utils.defaults import Defaults
from evac.plot.figure import Figure
from evac.plot.scales import Scales
import evac.stats as stats

class BirdsEye(Figure):
    """ Top-down 2D plots.

    Todo:
        * Massive rewrite with cartopy, now basemap is depreciated.
        * Move most options to mplkwargs (e.g. cmap, levels).

    Args:
        ax,fig (optional): If None (default), then a new figure/axis object
            will be created in the instance (via parent class 
            :class:`~evac.plot.figure.Figure`.
        ideal (bool): If True, options for an idealised plot will be used.
            For example, no lat/lon, geography data.
        grid: instance of Grid to assist cartopy setup.
        mplargs, mplkwargs: Arguments to pass to matplotlib.

    Returns:
        A figure is generated (if no fig/ax is specified), or plotting (onto an
            existing figure).
    """

    def __init__(self,fpath,ax=None,fig=None,ideal=False,grid=None,proj='merc',
        use_basemap=True):
        # ccrs_proj = grid.get_cartopy_proj(proj)
        super().__init__(fpath=fpath,ax=ax,fig=fig,grid=grid,use_basemap=use_basemap)
                            # proj=ccrs_proj,
        self.ideal = ideal
        self.use_basemap = use_basemap
        self.proj = proj

    def plot2D(self,data,plottype='contourf',locations=False,
                    save=True,smooth=1,title=False,cb=True,
                    # lats=False,lons=False,m=False,x=None,y=None,
                    # Nlim=False,Elim=False,Slim=False,Wlim=False,
                    # nx=None,ny=None,cen_lat=None,cen_lon=None,
                    # lat_ts=None,W=None,
                    inline=False,cblabel=False,ideal=False,
                    draw_geography=True,mplkwargs=None,
                    extend=False,hold=False,
                    lats=None,lons=None,proj=None,
                    # color='k',
                    *args,**kwargs):

        """
        Generic method that plots any matrix of data on a map

        Notes:
            mplkwargs settings include:
                * levels sets contour levels
                * cmap sets colormap
                * color sets plot line colours

        Todo:
            * just_one_colorbar integration

        Args:
            data: 2D matrix of data
            plottype (str): matplotlib function for plotting
            smooth (int): Gaussian smooth by this many grid spaces
            title (str): title on plot
            save (bool): whether to save to file
            locations: plot locations on map. Format is {'label':(lat,lon)}.
            cb (bool,str): if True, plot colorbar. If string, this
                should be an absolute path to plot one cb separately
            cblabel (str): label for colorbar
            ideal (bool): if True, plot without looking for geographical data
            drawcounties (bool): if True, plot US counties on map
            mplkwargs (dict): dictionary of keyword arguments to pass
                to plotting method
            kwargs: Alias for mplkwargs.

        """
        if (proj is None) and (not ideal):
            proj = getattr(self,'proj','merc')
        if mplkwargs is None:
            mplkwargs = {}
        mplkwargs.update(**kwargs)

        save = getattr(self,'save_opt',save)
        hold = getattr(self,'hold_opt',hold)

        self.data = data

        # if idealised, don't use geographical map
        if ideal:
            self.bmap = self.ax
            # pdb.set_trace()
            self.y = N.arange(len(data[:,0]))
            self.x = N.arange(len(data[0,:]))

        else:
            if lats is None:
                assert self.grid is not None
                # self.lats2d, self.lons2d = N.meshgrid(self.grid.lats,self.grid.lons)
                self.lats2d = self.grid.lats
                self.lons2d = self.grid.lons
            elif lats.ndim == 1:
                # self.lats2d, self.lons2d = N.meshgrid(lats,lons)
                self.lons2d, self.lats2d = N.meshgrid(lons,lats)
            elif lats.ndim == 2:
                self.lats2d = lats
                self.lons2d = lons
            else:
                raise Exception("lats and lons not valid.")
            assert self.lats2d.ndim == 2

            if self.use_basemap:
                if not self.basemap_done:
                    self.basemap_setup(proj=proj,draw_geography=draw_geography)
                    self.bmap = self.m
                    self.x, self.y = self.m(self.lons2d,self.lats2d)
                # self.y = self.grid.yy
            else:
                raise Exception("Need to make Cartopy work")
                self.x = self.lats2d
                self.y = self.lons2d

        # mplkwargs['X'] = self.x
        # mplkwargs['Y'] = self.y
        # mplkwargs['Z'] = self.data
        if plottype == 'scatter':
            mplargs = [*self.data]
        else:
            mplargs = [self.x, self.y, self.data]
        
        # pdb.set_trace()

        # TODO: move this logic to other methods
        if plottype == 'contour':
            cb = False
            f1 = self.bmap.contour(*mplargs,**mplkwargs)
            if inline:
                plt.clabel(f1,inline=True,fmt='%d',color='black',fontsize=9)
        elif plottype == 'contourf':
            f1 = self.bmap.contourf(*mplargs,**mplkwargs)
        elif plottype == 'pcolor':
            f1 = self.bmap.pcolor(*mplargs,**mplkwargs)
        elif plottype == 'pcolormesh':
            f1 = self.bmap.pcolormesh(*mplargs,**mplkwargs)
        elif plottype == 'scatter':
            cb = False
            f1 = self.bmap.scatter(*mplargs,**mplkwargs)
        elif plottype == 'quiver':
            f1 = self.bmap.quiver(*mplargs,**mplkwargs)
        elif plottype == 'plot':
            cb = False
            f1 = self.bmap.plot(*mplargs,**mplkwargs)
        else:
            print("Specify correct plot type.")
            raise Exception



        if isinstance(locations,dict):
            self.plot_locations(locations)
        if isinstance(title,str):
            self.ax.set_title(title)
        if cb != False:
            if cb==True:
                cb1 = plt.colorbar(f1,orientation='vertical',ax=self.ax)
            elif cb=='horizontal':
                cb1 = plt.colorbar(f1,orientation='horizontal',ax=self.ax)
            elif cb == 'only':
                save = False
                self.fig,self.ax = plt.subplots(figsize=(4,0.8))
                cb1 = plt.colorbar(f1,cax=self.ax,orientation='horizontal')
                if isinstance(cblabel,str):
                    cb1.set_label(cblabel)
                self.save(outdir,fname+'_cb')
            else:
                cb1 = plt.colorbar(f1,orientation='vertical',cax=cb)
        if cb and isinstance(cblabel,str):
            cb1.set_label(cblabel)
        if save:
            self.save()
            plt.close(self.fig)
        return

    def plot_locations(self,locations):
        for k,v in locations.items():
            if isinstance(v,tuple) and len(v) == 2:
                xpt, ypt = self.bmap(v[1],v[0])
                # bbox_style = {'boxstyle':'square','fc':'white','alpha':0.5}
                self.bmap.plot(xpt,ypt,'ko',markersize=3,zorder=100)
                self.ax.text(xpt,ypt,k,ha='left',fontsize=7)
                # self.ax.text(xpt-15000,ypt,k,bbox=bbox_style,ha='left',fontsize=7)
            else:
                print("Not a valid location argument.")
                raise Exception
        return

    def plot_streamlines(self,U,V,outdir,fname,lats=False,lons=False,smooth=1,
                            title=False,lw_speed=False,density=1.8,ideal=False):
        """
        Plot streamlines.

        U       :   U-component of wind (nx x ny)
        V       :   V-component of wind (same dimensions)

        lw_speed    :   linewidth is proportional to wind speed
        """
        if ideal:
            m = plt
            x = N.arange(len(U[:,0]))
            y = N.arange(len(U[0,:]))
        else:
            m,x,y = self.basemap_setup(lats=lats,lons=lons)

        if lw_speed:
            wind = N.sqrt(U**2 + V**2)
            lw = 5*wind/wind.max()
        else:
            lw = 1

        if smooth>1:
            U = stats.gauss_smooth(U,smooth)
            V = stats.gauss_smooth(V,smooth)

        # m.streamplot(x[self.W.x_dim/2,:],y[:,self.W.y_dim/2],U,V,
                        # density=density,linewidth=lw,color='k',arrowsize=3)
        m.streamplot(x,y,U,V,
                        density=density,linewidth=lw,color='k',arrowsize=3)

        if isinstance(title,str):
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
            print(("Plotting contour level {0} for {1} from file \n {2}".format(
                            contour,va,wrfout)))
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

    def make_subplot_label(ax,label):
        ax.text(0.05,0.85,label,transform=ax.transAxes,
            bbox={'facecolor':'white'},fontsize=15,zorder=1000)
        return

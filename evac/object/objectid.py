""" Identify objects from 2D (nlat x nlon) array of reflectivity, for
one time, from either observed or forecast data.

TODO:
    * Plug into object comparison, object attribute extraction objects.
"""

import pdb
import datetime
import os

import numpy as N
from mpl_toolkits.basemap import Basemap
import skimage.measure
import skimage.feature
#import skimage.preprocessing
import pandas
import scipy
from scipy.interpolate import RectBivariateSpline as RBS
import matplotlib as M
M.use("Agg")
import matplotlib.pyplot as plt

import evac.utils as utils
from evac.plot.birdseye import BirdsEye
from evac.plot.scales import Scales

class ObjectID:
    def __init__(self,data,lats,lons,dx=1,threshold=45.0,footprint=144,
                    prod_code=0,time_code=0,lead_time=0,classify=False,pca=None,
                    scaler=None,features=None):
        """
        Default footprint is multiple of 9 here, so that 108/(3km^2) = dx
            to compare domains.

        Args:
            dx: grid spacing, in km
            threshold: mimimum dBZ value to identify objects
            footprint: minimum area (in pixels) of objects
            time_code (datetime): Put in dataframe (otherwise 0)
            prod_code (float,int): put in dataframe
        """
        assert data.ndim == 2
        assert lats.ndim == 2

        self.classify = classify
        if classify:
            self.pca = pca
            self.scaler = scaler
            self.features = features

        #if isinstance(time_code,datetime.datetime):
            # time_code = time_code.toordinal()
        self.time_code = time_code
        self.prod_code = prod_code
        self.lead_time = lead_time

        self.dx = dx
        self.raw_data = data
        self.lats = lats
        self.lons = lons
        self.footprint = footprint
        self.threshold = threshold

        self.nlats, self.nlons = self.raw_data.shape
        # Create lat/lon RBS grids for fast interpolation of object location
        self.rbs = None
        # This is needed to remove storms on edge
        self.edgecoords = self.list_of_edgecoords()

        # For debugging... need to rename later as just 'maskout'
        self.radmask = N.zeros_like(self.raw_data)

        self.bmap = self.create_bmap(N.max(self.lats),N.max(self.lons),
                            N.min(self.lats),N.min(self.lons))
        self.xx, self.yy = self.bmap(self.lons,self.lats)

        self.object_field, self.object_props, self.OKidxs, self.idxarray, \
            self.OKidxarray = self.identify_objects()
        self.objects = self.create_structured_array(self.object_props,self.OKidxs)

    def create_bmap(self,urcrnrlat=None,urcrnrlon=None,llcrnrlat=None,
                        llcrnrlon=None,ax=None,pad=0.07,lats=None,lons=None):
        """ Basemap is suitable for both plotting and computation.

        Args:
            urcrnrlat... (float): extremes of the domain
            ax (optional): for plotting purposes
        """
        if ax is None:
            res = 'l'
        else:
            res = 'h'

        if urcrnrlat is None:
            if lats is None:
                lats = self.lats
                lons = self.lons

            maxlat = N.max(lats)
            maxlon = N.max(lons)#.max())
            minlat = N.min(lats)#.min())
            minlon = N.min(lons)#.min())

            lat_offset = pad*abs(maxlat-minlat)
            lon_offset = pad*abs(maxlon-minlon)

            urcrnrlat = maxlat+lat_offset
            urcrnrlon = maxlon+lon_offset
            llcrnrlat = minlat-lat_offset
            llcrnrlon = minlon-lon_offset

        meanlon = N.mean((urcrnrlon,llcrnrlon))
        meanlat = N.mean((urcrnrlat,llcrnrlat))
        bmap = Basemap(
                    # width=12000000,height=9000000,
                    urcrnrlat=urcrnrlat,urcrnrlon=urcrnrlon,
                    llcrnrlat=llcrnrlat,llcrnrlon=llcrnrlon,
                    rsphere=(6378137.00,6356752.3142),
                    resolution=res,projection='lcc',ax=ax,
                    lat_1=llcrnrlat,lat_2=urcrnrlat,lat_0=meanlat,lon_0=meanlon)

        if ax is not None:
            bmap.drawcoastlines()
            # bmap.drawmapboundary(fill_color='lightgray')
            # bmap.fillcontinents(color="white",lake_color="lightgray")
            bmap.drawstates()
            bmap.drawcountries()
        return bmap

    def identify_objects(self):
        data_masked = N.where(self.radmask,0.0,self.raw_data)
        obj_field_OK = N.zeros_like(data_masked)
        obj_init = N.where(data_masked >= self.threshold, data_masked, 0.0)
        obj_bool = obj_init > 0.0

        # 1, 2, 3...
        obj_labels = skimage.measure.label(obj_bool).astype(int)
        obj_labels_OK = N.zeros_like(obj_labels)

        # JRL: consider cache=False to lower memory but slow operation
        obj_props = skimage.measure.regionprops(obj_labels,obj_init)
        obj_OK = []
        for prop in obj_props:
            id = int(prop.label)
            size_OK = False
            not_edge_OK = False

            # Remove small storms
            if prop.area > (self.footprint/(self.dx**2)):
                size_OK = True

            # Remove storms touching perimeter
            # 1 km grid might have more storms that pass this check
            # This is 'fair', as that's the nature of the 1 km grid!
            # Also, a comparison might want to check if a storm was present
            # but removed - to avoid a false 'missed hit'. Instead, use a NaN.

            # print("Checking for an edge for object",id)
            # debug = True if id == 81 else False
            not_edge_OK = not self.check_for_edge(prop.coords,debug=False)

            # if all((size_OK,not_edge_OK)):
            if size_OK and not_edge_OK:
                obj_OK.append(id)

        # This replaces gridpoints within objects with original data
        for i,ol in enumerate(obj_OK):
            obj_field_OK = N.where(obj_labels==ol, obj_init, obj_field_OK)
            obj_labels_OK =  N.where(obj_labels==ol, obj_labels, obj_labels_OK)
        obj_props_OK = [obj_props[o-1] for o in obj_OK]
        # pdb.set_trace()
        return obj_field_OK, obj_props_OK, obj_OK, obj_labels, obj_labels_OK

    def list_of_edgecoords(self,pad=3):
        # TODO: consider doing outside 2 cells for 3km, and 6 cells for 1km
        pad = int(N.ceil(pad/self.dx))
        # These are all the values on the edge of the raw_data domain
        latidxs = N.indices(self.lats.shape)[0][:,0]
        lonidxs = N.indices(self.lats.shape)[1][0,:]

        ee0 = list(range(pad))
        elat0 = list(range(self.nlats-3,self.nlats))
        elon0 = list(range(self.nlons-3,self.nlons))
        elat = ee0 + elat0
        elon = ee0 + elon0

        _a = [(z,n) for n in lonidxs for z in elat]
        # _b = [(self.nlats,n) for n in lonidxs]
        _b = [(n,z) for n in latidxs for z in elon]
        # _d = [(n,self.nlons) for n in latidxs]

        edgecoords = _a + _b # + _c + _d
        # print("Edge coordinates have been calculated.")
        return edgecoords

    def check_for_edge(self,coords,debug=False):
        """ Returns True if any point is on the edge.
        """
        check_all = [(c[0],c[1]) in self.edgecoords for c in coords]
        if debug: pdb.set_trace()
        return any(check_all)

    def create_structured_array(self,obj_props,OKidx):
        """
        TODO:
            * A way to add to this array later, e.g., max updraught speed.
            * A way to compare 2+ of these arrays for numerous times, and
                match objects through it
            * Why is id one fewer than label? Which is used for
                calling up specific objects? It's related to the zeroth object
                being discarded somewhere in scikit logic.
        """
        # About 23?
        DTYPES = {
                    # "id":"i4",
                    "area":"i4",
                    "min_row":"i4",
                    "min_col":"i4",
                    "max_row":"i4",
                    "max_col":"i4",
                    "bbox_area":"i4", # why is this just the domain?
                    "centroid_row":"f4",
                    "centroid_col":"f4",
                    "centroid_lat":"f4",
                    "centroid_lon":"f4",
                    "convex_area":"i4",
                    "eccentricity":"f4",
                    "equivalent_diameter":"f4",
                    'extent':'f4',
                    "label":"i4",
                    "max_intensity":"f4",
                    "mean_intensity":"f4",
                    "min_intensity":"f4",
                    "perimeter":"f4",
                    "weighted_centroid_row":"f4",
                    "weighted_centroid_col":"f4",
                    "weighted_centroid_lat":"f4",
                    "weighted_centroid_lon":"f4",

                    "ratio":"f4",
                    "longaxis_km":"f4",
                    # "prod":"i4",
                    "prod_code":"object",
                    # "time":"int32",
                    "time":"datetime64",
                    "qlcsness":'f4',
                    # "qlcsness":'string_',
                    "lead_time":"f4",
                    "dx":"f4",
                    }

        dtypes = {'names':[], 'formats':[]}
        for k,v in DTYPES.items():
            dtypes['names'].append(k)
            dtypes['formats'].append(v)

        nobjs = len(obj_props)
        nprops = len(dtypes['names'])
        # N.dtype([(k,v) for k,v in DTYPES.items()])
        # objects = N.array(nobjs,dtype=N.dtype(dtypes))
        obj_arr = N.zeros([nprops,nobjs]).T
        # pdb.set_trace()
        objects = pandas.DataFrame(data=obj_arr,columns=dtypes['names'],)
                            # dtype='f4')

        for oidx,o in enumerate(obj_props):
            # ID number in regionprops (for looking up more stuff later)
            # JRL: this is misleading - always one lower than label. Use label.
            # objects['id'][oidx] = OKidx[oidx]

            # Number of pixels in object (footprint)
            # DTYPE: small int
            # objects['area'][oidx] = o.area
            objects.loc[oidx,'area'] =  o.area

            # Bounding box (min_row, min_col, max_row, max_col)
            # Pixels belonging to the bounding box are in
            # the half-open interval [min_row; max_row) and [min_col; max_col)
            # DTYPE: small ints
            min_row, min_col, max_row, max_col = o.bbox
            #objects['min_row'][oidx] = min_row
            objects.loc[oidx,'min_row'] = min_row
            #objects['max_row'][oidx] = max_row
            objects.loc[oidx,'max_row'] = max_row
            #objects['min_col'][oidx] = min_col
            objects.loc[oidx,'min_col'] = min_col
            #objects['max_col'][oidx] = max_col
            objects.loc[oidx,'max_col'] = max_col

            # Bounding box area, in pixels. The rectangle containing obj?
            # DTYPE: small int
            # JRL: this appears to be almost the whole domain?
            #objects['bbox_area'][oidx] = o.bbox_area
            objects.loc[oidx,'bbox_area'] = o.bbox_area

            # Centroid as (row,col). DTYLE: small floats
            cr, cc = o.centroid
            #objects['centroid_row'][oidx] = cr
            objects.loc[oidx,'centroid_row'] = cr
            #objects['centroid_col'][oidx] = cc
            objects.loc[oidx,'centroid_col'] = cc

            # Convert to lat/lon (DTYLE: small floats)
            clat, clon = self.interp_latlon(yrow=cr,xcol=cc)
            #objects['centroid_lat'][oidx] = clat
            objects.loc[oidx,'centroid_lat'] = clat[0][0]
            #objects['centroid_lon'][oidx] = clon
            objects.loc[oidx,'centroid_lon'] = clon[0][0]

            # Number of pixel in convex hull (polygon that encloses region)
            # DTYPE: small int
            #objects['convex_area'][oidx] = o.convex_area
            objects.loc[oidx,'convex_area'] = o.convex_area

            # A boolean array of the object's presence. Not used here.
            # Access via the number
            # convex_image

            # An array of coordinates contained in the object. Not used here
            # coords

            # Eccentricity of the ellipse with same second moments as object
            # A value of 0 is a circle.
            # DTYPE: small float
            #objects['eccentricity'][oidx] = o.eccentricity
            objects.loc[oidx,'eccentricity'] = o.eccentricity

            # Equivalent diameter: a circle with the same area
            # JRL: is this in grid points? Looks the same units as footprint
            # DTYPE: small float
            # objects['equivalent_diameter'][oidx] = o.equivalent_diameter
            objects.loc[oidx,'equivalent_diameter'] = o.equivalent_diameter

            # skipping euler_number (number of holes?)

            # skipping extent (area/ (rows x cols))
            # objects['extent'][oidx] = o.extent
            objects.loc[oidx,'extent'] = o.extent

            # skipped filled_area and filled_image (holes filled in)

            # skipped image (boundary box image)

            # skipped inertia_tensor - rotation around its mass?
            # skipping inertia_tensor_eigvals

            # skipped intensity_image, which is each object dBZ

            # Is this the same as ID? DTYPE: small int
            #objects['label'][oidx] = o.label
            objects.loc[oidx,'label'] = o.label

            # Centroid, relative to region bounding box... huh? Skipping

            # maximum/mean/minimum value in object
            # DTYPE: float
            #objects['max_intensity'][oidx] = o.max_intensity
            objects.loc[oidx,'max_intensity'] = o.max_intensity
            #objects['mean_intensity'][oidx] = o.mean_intensity
            objects.loc[oidx,'mean_intensity'] = o.mean_intensity
            #objects['min_intensity'][oidx] = o.min_intensity
            objects.loc[oidx,'min_intensity'] = o.min_intensity

            # perimeter - for looking at fractal dimension increase
            # DTYPE:float
            #objects['perimeter'][oidx] = o.perimeter
            objects.loc[oidx,'perimeter'] = o.perimeter

            # JRL: slice is missing from mine - to extract from original data

            # weighted centroid. DTYPE: floats
            wcr, wcc = o.weighted_centroid
            # objects['weighted_centroid_row'][oidx] = wcr
            objects.loc[oidx,'weighted_centroid_row'] = wcr
            # objects['weighted_centroid_col'][oidx] = wcc
            objects.loc[oidx,'weighted_centroid_col'] = wcc

            # Convert to lat/lon (DTYLE: small floats)
            wclat, wclon = self.interp_latlon(yrow=wcr,xcol=wcc)
            # objects['weighted_centroid_lat'][oidx] = wclat[0][0]
            # objects['weighted_centroid_lon'][oidx] = wclon[0][0]
            objects.loc[oidx,'weighted_centroid_lat'] = wclat[0][0]
            objects.loc[oidx,'weighted_centroid_lon'] = wclon[0][0]

            # JRL: custom ratio of longest to shortest sides
            #drow = abs(min_row-max_row)
            #dcol = abs(min_col-max_col)
            #if drow == dcol:
            #    ratio = 1
            #else:
            #    maxside = max((drow,dcol))
            #    minside = min((drow,dcol))
            #    ratio = maxside/minside
            ratio = o.minor_axis_length/o.major_axis_length
            #objects['ratio'][oidx] = ratio
            objects.loc[oidx,'ratio'] = ratio

            #objects['longaxis_km'][oidx] = o.major_axis_length * self.dx
            objects.loc[oidx,'longaxis_km'] = o.major_axis_length * self.dx

            #objects['prod'][oidx] = self.prod_code
            objects.loc[oidx,'prod_code'] = self.prod_code
            #objects['time'][oidx] = self.time_code
            objects.loc[oidx,'time'] = self.time_code

            if self.classify:
                series = objects.loc[oidx:oidx,self.features]
                qlcsness = self.pca.transform(self.scaler.transform(series))[0][0]
                # pdb.set_trace()
                #objects['qlcsness'][oidx] = qlcsness
                objects.loc[oidx,'qlcsness'] = qlcsness
                # aa = results[0].loc[0:0,features]
                # aat = pca.transform(scaler.transform(aa))

            # Lead time in minutes, usually.
            #objects['lead_time'][oidx] = self.lead_time
            objects.loc[oidx,'lead_time'] = self.lead_time
            objects.loc[oidx,"dx"] = self.dx
            # pdb.set_trace()
        return objects

    def interp_latlon(self,xcol,yrow):
        #points = (N.indices(self.xx.shape)[0].flat,N.indices(self.yy.shape)[0].flat)
        # clat = scipy.interpolate.griddata(points=(self.xx.flat,self.yy.flat),
        #clat = scipy.interpolate.griddata(points=points,
        #            values=self.lats.flat,xi=(xcol,yrow),method='linear')
        # clon = scipy.interpolate.griddata(points=(self.xx.flat,self.yy.flat),
        #clon = scipy.interpolate.griddata(points=points,
        #            values=self.lons.flat,xi=(xcol,yrow),method='linear')
        #pdb.set_trace()
        if self.rbs is None:
            self.rbs = dict()
            x = N.arange(self.raw_data.shape[0])
            y = N.arange(self.raw_data.shape[1])
            self.rbs['lats'] = RBS(x=x,y=y,z=self.lats)
            self.rbs['lons'] = RBS(x=x,y=y,z=self.lons)

        clat = self.rbs['lats'](yrow,xcol)
        clon = self.rbs['lons'](yrow,xcol)
        return clat, clon

    def determine_mode(self):
        """ Determine if a reflectivity object is likely cellular or linear.

        For each object, QLCS objects will have:
            1. High eccentricity (>0.75?)
            2. Major axis length
            3.
        """
        pass

    def plot_quicklook(self,outdir,what='all',fname=None,ecc=0.2):
        """ Plot quick images of objects identified.

        Args:
            what (str): the type of quickplot
            outdir (str): where to put the quickplot

            fname (str,optional): if None, automatically name.
            ecc (float): eccentricity of object, for discriminating (later)
        """
        assert what in ("all","qlcs","ecc","shapeindex","ratio","extent",
                            "4-panel","pca")

        def label_objs(bmap,ax,locs):
            for k,v in locs.items():
                xpt, ypt = bmap(v[1],v[0])
                # bbox_style = {'boxstyle':'square','fc':'white','alpha':0.5}
                # bmap.plot(xpt,ypt,'ko',)#markersize=3,zorder=100)
                # ax.text(xpt,ypt,k,ha='left',fontsize=15)
                ax.annotate(k,xy=(xpt,ypt),xycoords="data",zorder=1000,
                                fontsize=7,fontweight='bold')
                # pdb.set_trace()
            return

        def do_label_array(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                locs[str(int(o.label))] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.idxarray < 1, self.idxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object labels")
            return

        def __do_label(bmap,ax):
            locs = dict()
            for o, okidx in zip(self.objects.itertuples(),self.OKidxs):
                #data = skimage.preprocessing.normalize(
                raw_data = self.object_props[okidx].intensity_image
                pdb.set_trace()
                data = N.linalg.norm(axis=(0,1),x=raw_data)
                shidx = skimage.feature.shape_index(data)
                sistr = "{:0.3f}".format(shidx)
                locs[sistr] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.idxarray < 1, self.idxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object shape index")
            return

        def do_label_ecc(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                eccstr = "{:0.2f}".format(o.eccentricity)
                locs[eccstr] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.OKidxarray < 1, self.OKidxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,alpha=0.5)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object eccentricity")
            return

        def do_label_extent(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                lab = "{:1.2f}".format(o.extent)
                locs[lab] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.OKidxarray < 1, self.OKidxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,alpha=0.5)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object extent (fill pc. of bbox)")
            return

        def do_label_longest(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                lab = "{}".format(int(o.longaxis_km))
                locs[lab] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.OKidxarray < 1, self.OKidxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,alpha=0.5)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object longest-side length (km)")
            return

        def do_label_pca(bmap,ax,discrim_vals=(-0.2,0.5)):
            locs = dict()
            obj_discrim = N.zeros_like(self.OKidxarray)
            for o in self.objects.itertuples():
                lab = "{:0.2f}".format(o.qlcsness)
                locs[lab] = (o.centroid_lat,o.centroid_lon)
                id = int(o.label)
                if o.qlcsness < discrim_vals[0]:
                    discrim = 1
                elif o.qlcsness > discrim_vals[1]:
                    discrim = 3
                else:
                    discrim = 2
                obj_discrim = N.where(self.idxarray == id, discrim, obj_discrim)
            marr = N.ma.masked_where(obj_discrim < 1, obj_discrim)
            pcm = bmap.pcolormesh(data=marr,x=x,y=y,alpha=0.5,vmin=1,vmax=3,
                                    cmap=M.cm.get_cmap("magma",3),)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))

            mode_names = ["","Cellular","Ambiguous","Linear/Complex"]

            def format_func(x,y):
                return mode_names[x]


            # This function formatter will replace integers with target names


            formatter = plt.FuncFormatter(lambda val, loc: mode_names[val])

            # We must be sure to specify the ticks matching our target names
            cax = fig.add_axes([0.75, 0.1, 0.2, 0.04])
            #cax.set_xlim(0.5, 3.5)
            fig.colorbar(pcm,cax=cax,ticks=(1,2,3), format=formatter,
                            orientation='horizontal')

            # Set the clim so that labels are centered on each block

            label_objs(bmap,ax,locs)
            #plt.colorbar(pcm,cax=cax)
            ax.set_title("Object PC1 (QLCS-ness)")
            return

        def do_label_ratio(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                lab = "{:1.2f}".format(o.ratio)
                locs[lab] = (o.centroid_lat,o.centroid_lon)
            idxarray = N.ma.masked_where(self.OKidxarray < 1, self.OKidxarray)
            bmap.pcolormesh(data=idxarray,x=x,y=y,alpha=0.5)
                    #vmin=1,vmax=self.idxarray.max()+1)
            #        levels=N.arange(1,self.objects.shape[0]))
            label_objs(bmap,ax,locs)
            ax.set_title("Object side-ratio (short/long)")
            return

        def do_raw_array(bmap,ax):
            bmap.contourf(data=self.raw_data,x=x,y=y,
                    levels=S.clvs,cmap=S.cm)
            ax.set_title("Raw data")
            return

        def do_intensity_array(bmap,ax):
            locs = dict()
            for o in self.objects.itertuples():
                locs[str(int(o.label))] = (o.weighted_centroid_lat,o.weighted_centroid_lon)
            bmap.contourf(data=self.object_field,x=x,y=y,
                            levels=S.clvs,cmap=S.cm)
            label_objs(bmap,ax,locs)
            ax.set_title("Intensity array")
            return

        def do_table(bmap,ax):
            cell_text = []
            table = self.objects
            for row in range(len(table)):
                cell_text.append(table.iloc[row])

            tab = plt.table(cellText=cell_text, colLabels=table.columns, loc='center')
            ax.add_table(tab)
            plt.axis('off')
            return

        def do_highlow_ecc(bmap,ax,ecc,overunder):
            assert overunder in ("over","under")
            func = N.greater_equal if overunder == "over" else N.less
            locs = dict()
            qlcs_obj = []
            obj_field = N.zeros_like(self.object_field)
            for o in self.objects.itertuples():
                if func(o.eccentricity,ecc):
                    locs[str(int(o.label))] = (o.centroid_lat,o.centroid_lon)
                    qlcs_obj.append(o.label)
            for olab in qlcs_obj:
                obj_field = N.where(self.idxarray==olab, self.object_field, obj_field)
            bmap.contourf(data=obj_field,x=x,y=y,
                            levels=S.clvs,cmap=S.cm)
            label_objs(bmap,ax,locs)
            if overunder == "over":
                ax.set_title("QLCS objects (ecc > {:0.2f})".format(ecc))
            else:
                ax.set_title("Cellular objects (ecc < {:0.2f})".format(ecc))
            return

        # The above needs saving to pickle and loading/copying each time
        # a plot is made, if optimising.

        # 2x2:
        # ax1 is the raw field
        # ax2 is the labelled objects (ID)
        # ax3 is the object field, annotate row/col and lat/lon centroids
        # ax4 is the DataFrame?
        if what == "4-panel":
            fig,axes = plt.subplots(2,2,figsize=(9,7))
        else:
            fig,axes = plt.subplots(1,3,figsize=(9,4))

        if fname is None:
            fname = "quicklook.png"

        for n,ax in enumerate(axes.flat):
            # if n!=0:
            #     continue
            print("Plotting subplot #",n+1)
            bmap = self.create_bmap(ax=ax)
            x,y = bmap(self.lons,self.lats)
            S = Scales(vrbl='REFL_comp')

            if n == 0:
                do_raw_array(bmap,ax)

            elif n == 1:
                if what == "qlcs":
                    do_highlow_ecc(bmap,ax,ecc,"over")
                elif what in ("ecc","extent","4-panel"):
                    do_label_ecc(bmap,ax)
                else:
                    do_intensity_array(bmap,ax)
                #elif what == "ecc":
                #    do_intensity_array(bmap,ax)


            elif n==2:
                if what == "qlcs":
                    do_highlow_ecc(bmap,ax,ecc,"under")
                elif what == "all":
                    do_label_array(bmap,ax)
                elif what== "ecc":
                    #do_label_ecc(bmap,ax)
                    do_label_ratio(bmap,ax)
                elif what == "shapeindex":
                    do_label_shapeindex(bmap,ax)
                elif what == "ratio":
                    do_label_ratio(bmap,ax)
                elif what in ("4-panel","extent"):
                    do_label_extent(bmap,ax)
                elif what == "pca":
                    do_label_pca(bmap,ax)

            elif n==3:
                if what == "4-panel":
                    do_label_longest(bmap,ax)
                else:
                    do_table(bmap,ax)

        fpath = os.path.join(outdir,fname)
        utils.trycreate(fpath)
        fig.tight_layout()
        fig.savefig(fpath)
        print("Saved figure to",fpath)
        plt.close(fig)
        return

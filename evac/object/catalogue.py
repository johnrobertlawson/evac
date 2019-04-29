import os
import pdb
import multiprocessing
import pickle
import itertools
import datetime
import operator

import matplotlib as M
import matplotlib.pyplot as plt
import numpy as N
import pandas
import sklearn.preprocessing
import sklearn.decomposition
# from scipy.spatial import KDTree
import scipy.spatial.distance as ssd
from mpl_toolkits.basemap import Basemap

import evac.utils as utils

class Catalogue:
    def __init__(self,df_list,ncpus=1,tempdir='./'): #time_list=None,prod_list=None):
        """
        Args:
            df_list: list of DataFrames, or one big one.
            tempdir: optional absolute path to folder for saving temporary
                files to disc (to increase speed of script)

        Note:
            * The time and product/model/format of each object's dataframe object
                is present in the dataframe itself. This way, a mega-array is
                generated that can be sliced, or matching can be performed.

        TODO:
            * Pass in forecast data, and a new 'side' df is created that
                has megaframe_idx to sync with the main object data.
        """
        self.df_list = df_list
        self.ncpus = ncpus
        self.tempdir = tempdir

        if isinstance(self.df_list,pandas.DataFrame):
            self.megaframe = self.df_list
            self.ndfs = None
        else:
            self.ndfs = len(self.df_list)
            self.megaframe = pandas.concat(self.df_list,ignore_index=True)

        self.nobjs = len(self.megaframe)

        # TODO: make sure every object has a unique index, that is preserved
        # even when slicing

    def match_domains(self,members,initutcs,leadtimes,domains):
        """ Use of Total Interest to match objects in time/across models or obs.

        Maybe do CombinedObject class. This could score itself!

        Need a lat/lon grid for each product.
        """
        # JRL: not set up for more than 2 domains yet (plus verif)
        assert len(domains) == 2
        self.bmap = self.create_bmap()

        # Separate dataset into [domains and ensemble members] and verification
        def gen():
            for member, initutc, leadtime in itertools.product(members,
                                                    initutcs,leadtimes):
                yield member, initutc, leadtime, domains, None #verif_domain

        itr = gen()
        print("Submitting jobs to match objects...")

        if self.ncpus > 1:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.match,itr)
        else:
            for i in itr:
                self.match(i)

        print("Jobs done.")
        return results

    def match_verif(self,members,initutcs,leadtimes,domain,verif_domain):
        """ Use of Total Interest to match objects in time/across models or obs.

        Maybe do CombinedObject class. This could score itself!

        Need a lat/lon grid for each product.
        """
        # JRL: not set up for more than 2 domains yet (plus verif)
        assert isinstance(domain,str)
        self.bmap = self.create_bmap()

        # Separate dataset into [domains and ensemble members] and verification
        def gen():
            for member, initutc, leadtime in itertools.product(members,
                                                    initutcs,leadtimes):
                yield member, initutc, leadtime, domain, verif_domain

        itr = gen()
        print("Submitting jobs to match objects...")

        if self.ncpus > 1:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.match,itr)
        else:
            for i in itr:
                self.match(i)

        print("Jobs done.")
        return matches


    def match(self,i):
        def get_sub_df(mf,member,validtime,domain,leadtime):
            sub_df = mf[(mf['member'] == member) &
                        (mf['time'] == validtime) &
                        (mf['domain'] == domain) &
                        # This line is just for sanity
                        (mf['lead_time'] == leadtime)]
            # pdb.set_trace()
            return sub_df

        member, initutc, leadtime, domains, verif_domain = i
        mf = self.megaframe
        validtime = initutc + datetime.timedelta(seconds=60*int(leadtime))
        if False:
            print("=====MATCHING!=====\n DEBUG:",member, initutc, leadtime,
                        domains, verif_domain)

        both_df = dict()
        if verif_domain is None:
            for domain in domains:
                both_df[domain] = get_sub_df(mf,member,validtime,domain,leadtime)
        else:
            # To avoid confusion
            assert isinstance(domains,str)
            domain = domains
            del domains
            both_df[domains] = get_sub_df(mf,member,validtime,domain,leadtime)
            both_df[verif_domain] = mf[(mf['time']==validtime) & (mf['domain']==verif_domain)]

        # First, let's match the object between all forecast domains
        # Need to create a dataframe that slices all possible objects
        tups = {}
        for domain, df in both_df.items():
            tups[domain] = []
            for o in df.itertuples():
                tups[domain].append(o)

        listoflists = []
        for domain in domains:
            td = tups[domain]
            if len(td) == 0:
                return None
            listoflists.append(td)

        # Just for this member:
        member_TIs = {}
        print("=====MATCHING!=====\n DEBUG:",member, initutc, leadtime,
                                domains, verif_domain)

        for objA, objB in itertools.product(*listoflists):
            compare_tup = (int(objA.megaframe_idx),int(objB.megaframe_idx))
            # pdb.set_trace()
            TI = self.compute_total_interest(propA=objA,propB=objB)
            member_TIs[compare_tup] = TI
        # Find the maximum TI and the pair that produced it

        matched_pair = max(member_TIs.items(), key=operator.itemgetter(1))[0]

        # idx = N.argmax(N.array(ti_data))
        # the_combo = combo[idx]

        if member_TIs[matched_pair] < 0.2:
            return None

        # matches.append(matched_pair)
        return matched_pair

    def compute_total_interest(self,propA=None,propB=None,cd=None,md=None,td=None,
                            cd_max=40.0,md_max=40.0,td_max=25.0,use_bbox=False):
        """
        Returns the total interest between object pairs.

        Args:
            propA, propB: dictionary? struct array? panda dataframe? slice of
                one object's properties

            cd (float): the centroid distance between two objects (km)
            md (float): the minimum distance between two objects (km)
            td (float): the time difference between the objects (min)

            cd_max (float): the limit for centroid distance (km)
            md_max (float): the limit for minimum distance (km)
        """
        def gen_latlon(idx, lats, lons):
            pdb.set_trace()
            latsA = 0
            lonsA = 0
            latsB = 0
            lonsB = 0

            return latsA, lonsA, latsB, lonsB

        # Sanity checks
        if propA is not None:
            assert propB is not None
            assert not all((cd,md,td))
            method = 1
        else:
            assert propB is None
            assert all((cd,md,td))
            method = 2
            raise Exception
        assert all((cd_max,md_max,td_max))

        #if method == 1:
        # Compute distance between centroids (km)
        cd = utils.xs_distance(propA.centroid_lat, propA.centroid_lon,
                        propB.centroid_lat, propB.centroid_lon)/1000.0

        # Compute closest distance between objects

        if use_bbox:
            objidx = []
            # Assume bounding box is the whole object
            for obj in (propA,propB):
                minr = int(obj.min_row)
                maxr = int(obj.max_row)
                minc = int(obj.min_col)
                maxc = int(obj.max_col)

                # Open the right lat/lon file and find the lats/lons
                objidx.append(slice(minr,maxr),slice(minc,maxc))
        else:
            #lonpts = {0:[],1:[]}
            #latpts = {0:[],1:[]}
            coords = {0:[],1:[]}
            llpts = {0:[],1:[]}
            xypts = {0:[],1:[]}
            xx = {0:[],1:[]}
            yy = {0:[],1:[]}
            for n,obj in enumerate((propA,propB)):
                OID = utils.load_pickle(obj.fpath_save)
                # OID.lats, OID.lons
                # OID
                for o in OID.object_props:
                    if all((obj.mean_intensity == o.mean_intensity,
                                obj.centroid_row == o.centroid[0],
                                obj.centroid_col == o.centroid[1])):
                        coords[n] = o.coords

                lats = OID.lats[coords[n][:,0],coords[n][:,1]]
                lons = OID.lons[coords[n][:,0],coords[n][:,1]]
                xx[n], yy[n] = self.bmap(lons,lats)
                xypts[n] = N.array([(x,y) for x,y in zip(xx[n],yy[n])])


                #for cx,cy in coords[n]:
                    # latpts[n].append(OID.lats[cx,cy])
                    # lonpts[n].append(OID.lons[cx,cy])
                    #llpts[n].append((OID.lats[cx,cy],OID.lons[cx,cy]))
                # npts = coords.shape[0]
                # pdb.set_trace()

                #del coords


        # latpts[n] has all latitude values for points in object n
        #kdarr = N.zeros_like()
        #KDTree()
        #combo_lat = list(utils.combinate_lists(latpts[0],latpts[1]))
        #combo_lon = list(utils.combinate_lists(lonpts[0],lonpts[1]))
        #latsA, lonsA, latsB, lonsB = gen_latlon(combo_idx,lat,lon_slices)


        distances = ssd.cdist(N.array(xypts[0]),N.array(xypts[1]),'euclidean')


        #distances = haversine_baker(lonsA, latsA, lonsB, latsB)
        md = N.nanmin(distances)/1000
        # pdb.set_trace()

        td = ((propA.time - propB.time).total_seconds())/60

        TI = ((td_max-td)/td_max)*0.5*(((cd_max-cd)/cd_max)+((md_max-md)/md_max))
        return TI

    def compute_new_attributes(self,df_in,do_suite="W"):
        """ Adds new attributes from a forecast field to each object dataframe.

        data_in columns: fcst_vrbl, valid_time, fcst_min, prod_code, path_to_pickle
        (all that's needed to find unique matches)

        TODO:
            * do_suite for UH/AWS?
        """
        def get_commands(df_in):
            commands = []
            for fidx,fcst_df in enumerate(df_in.itertuples()):
                if (fcst_df.fcst_min == 0):
                    continue
                utils.print_progress(total=len(df_in),idx=fidx,every=500)
                obj_df = self.lookup_obj_df(self.megaframe,
                                            valid_time=fcst_df.valid_time,
                                            fcst_min=fcst_df.fcst_min,
                                            prod_code=fcst_df.prod_code)
                if obj_df.shape[0] == 0:
                    continue
                # pdb.set_trace()
                if do_suite == "W":
                    commands.append([obj_df, fcst_df.path_to_pickle, fidx])
                elif do_suite == "UH02":
                #    self.get_uh_attributes(lv='02')
                    raise Exception
                elif do_suite == "UH25":
                #    self.get_uh_attributes(lv='25')
                    raise Exception
            return commands

        if do_suite == "W":
            func = self.get_updraught_attributes

        fname = "commands.pickle"
        do_pickle = True
        fpath = os.path.join(self.tempdir,fname)
        if os.path.exists(fpath) and do_pickle:
            print("Loading pickle of commands")
            with open(fpath,"rb") as f:
                commands = pickle.load(file=f)
        else:
            print("Generating list of commands")
            commands = list(get_commands(df_in))
            # pdb.set_trace()
            with open(fpath,"wb") as f:
                pickle.dump(obj=commands,file=f)

        # Submit in parallel
        print("Submitting jobs to compute {} attributes...".format(do_suite))
        if self.ncpus > 1:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(func,commands)
        else:
            for i in commands:
                func(i)

        print("Jobs done.")

        # pdb.set_trace()
        df_out = pandas.concat(results,ignore_index=True)
        return df_out

    def lookup_obj_df(self,data,valid_time,fcst_min=0,prod_code=0):
        """
        Args:
            data: an itertuple from a data_in dataframe (appending new info)
        """

        little_slice = data[(data['prod_code'] == prod_code) &
                            (data['lead_time'] == fcst_min) &
                            (data['time'] == valid_time)]
        return little_slice

    # def get_updraught_attributes(self,obj_df,W_field):
    def get_updraught_attributes(self,i):
    #@classmethod
    #def get_updraught_attributes(cls,obj_df,W_field):
        """
        Returns dataframe of object label (same as in obj_df) and new
            attributes. This can be concatenated with an existing df.
        Args:
            obj_df : DataFrame with object information
            W_field : 2D array of W data
        Todo:
            * Updraught max level (hPa, km AGL), magnitude (m/s)
            * Updraught location along major axis or radius?
        """
        obj_df, W_field, fidx = i
        if isinstance(W_field,str):
            W_field = N.load(W_field)

        #pdb.set_trace()
        # This won't work because the index is the megaframe idx
        # prod = obj_df.loc[0,"prod_code"]
        #t = obj_df.loc[0,"time"]
        #lead_time = obj_df.loc[0,"lead_time"]

        mega_idx = obj_df.head(1).megaframe_idx.values[0]
        #if (mega_idx % 1000) == 0:
        if False:
            prod = obj_df.head(1).prod_code.values[0]
            t = obj_df.head(1).time.values[0]
            lead_time = obj_df.head(1).lead_time.values[0]
            print("About to compute updraught attributes for:",
                prod,t,"||| lead time =", lead_time,"min.")

        # utils.print_progress(total=self.nobjs,idx=fidx,every=300)

        DTYPES = {
                "megaframe_idx_test":"i4",

                "max_updraught":"f4",
                "max_updraught_row":"i4",
                "max_updraught_col":"i4",

                "min_updraught":"f4",
                "mean_updraught":"f4",
                "ud_distance_from_centroid":"f4",
                "ud_angle_from_centroid":"f4",
                "test_gridsize":"object",
                }

        nobjs = len(obj_df)
        new_df = utils.do_new_df(DTYPES,nobjs,)

        for oidx,obj in enumerate(obj_df.itertuples()):
            dx = obj.dx

            # Get index for object in megaframe
            # oidx = obj.label

            # Get index for object in new df
            #Ix = obj.Index
            Ix = obj.megaframe_idx
            new_df.loc[oidx,"megaframe_idx_test"] = Ix


            minr = int(obj.min_row)
            maxr = int(obj.max_row)
            minc = int(obj.min_col)
            maxc = int(obj.max_col)

            # Max/mean updraught in this object's bbox
            objidx = slice(minr,maxr),slice(minc,maxc)
            # objidx = obj['coords']
            W_slice = W_field[objidx]

            new_df.loc[oidx,'max_updraught'] = N.nanmax(W_slice)
            #N.mean(W_slice)

            # Location of max updraught in object
            maxmax = N.where(W_slice == N.nanmax(W_slice))
            maxur = maxmax[0][0]
            maxuc = maxmax[1][0]

            # Find this location back on the main grid
            new_df.loc[oidx,'max_updraught_row'] = maxur + minr
            new_df.loc[oidx,'max_updraught_col'] = maxuc + minc

            # Min/mean
            new_df.loc[oidx,'mean_updraught'] = N.nanmean(W_slice)
            new_df.loc[oidx,"min_updraught"] = N.nanmin(W_slice)

            # Relative to centroid? Could be radius, using equivalent circle
            # distance_from_centroid
            cr = obj.centroid_row
            cc = obj.centroid_col
            dist_km, angle = utils.distance_angle_from_coords(maxur+minr,
                                                    maxuc+minc,cr,cc,dx=dx)
            new_df.loc[oidx,'ud_distance_from_centroid'] = dist_km

            # angle_from_centroid
            new_df.loc[oidx,'ud_angle_from_centroid'] = angle



            new_df.loc[oidx,'test_gridsize'] = str(W_field.shape[0])

            # distance_from_wcentroid

            # Size of updraught, using W = ? m/s as threshold
            # updraught_area_km
            # updraught_width_km
            # pdb.set_trace()

        return new_df

    def do_pca(self,features=None):
        """
        Todo:
            * Move this to PCA stat class. Return that class.
        """
        if features is None:
            features = ['area','eccentricity','extent','max_intensity',
                        'mean_intensity','perimeter','longaxis_km']
        x = self.megaframe.loc[:,features].values
        y = self.megaframe.loc[:,['prod_code']].values
        scaler = sklearn.preprocessing.StandardScaler()
        x = scaler.fit_transform(x)

        pca = sklearn.decomposition.PCA(n_components=3)
        PCs = pca.fit_transform(x)
        cols = ['PC1','PC2','PC3']
        _PC_df = pandas.DataFrame(data=PCs,columns=cols)
        PC_df = pandas.concat([_PC_df, self.megaframe[['prod_code']]], axis=1)

        return pca, PC_df, features, scaler

    def do_lda(self,train_data,targets):
        """
        Args:
            train_data: this has the training data
            targets: the classes belonging to each data point
        """
        features = ['area','eccentricity','extent','max_intensity',
                        'mean_intensity','perimeter','longaxis_km']
        X = self.megaframe.loc[:,features].values

        lda = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()
        hmm = lda.transform(X)

    def create_bmap(self,):
        bmap = Basemap(width=12000000,height=9000000,
                    rsphere=(6378137.00,6356752.3142),
                    resolution='l',projection='lcc',
                    lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.0)
        return bmap

import os
import pdb
import multiprocessing
import pickle
import itertools
import datetime
import operator
import time
import math
import copy

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
from evac.stats.detscores import DetScores

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
        """ Generate 2x2 contingency for each domain.
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

    def match_two_groups(self,dictA,dictB,do_contingency=False,
                            td_max=20.0):
        """ Match object belonging to any subset of the megaframe -
        make them non-overlapping if you don't want weird results.

        Do dictionary as key = property and value = value, like:

        dictA = dict("member" = "m01")

        If do_contingency, also return 2x2 table (see that class module!)

        JRL TODO: must have a way to subset times to within the
        max time difference window, to limit calculation of permutations
        """
        self.bmap = self.create_bmap()
        rets = self.match_contingency(dictA,dictB,td_max=td_max)

        print("Jobs done.")
        return rets
        # iterator can be each of the dictA/dictB keys!

        # For each permutation of objects...
        # Need to pass back out the result:
        # a, b, c, d

    def match_contingency(self,dictA,dictB,td_max=20.0):
        def get_sub_df(mf,the_dict):
            sf = mf.copy()
            for k,v in the_dict.items():
                sf = sf[sf[k] == v]
            return sf

        mf = self.megaframe
        # validtime = initutc + datetime.timedelta(seconds=60*int(leadtime))

        both_df = dict()
        for nd,d in enumerate((dictA,dictB)):
            both_df[nd] = get_sub_df(mf,d)

        # At this point, we have a both_df[x].shape of (nobj, nprops) for
        # x = (0,1).

        # Need to create a dataframe that slices all possible objects
        # JRL TODO: do I need to do permutations with repetition,
        # especially for verification purposes?
        # Also need to only permute the list of times that is within the
        # time window

        # this_gen_len = 0
        method = 1
        if method == 1:
            def gen():
                # This will shuffle - faster?
                # itA = both_df[0].sample(frac=1).itertuples()
                # itB = both_df[1].sample(frac=1).itertuples()
                itA = both_df[0].itertuples()
                itB = both_df[1].itertuples()
                for objA, objB in itertools.product(itA,itB):
                    # Only yield if objA and objB are within td_max
                    if abs((objA.time-objB.time).total_seconds()) <= (60.0*td_max):
                        yield objA, objB, self.bmap

            def count_gen():
                this_gen_len = 0
                for objA, objB in itertools.product(both_df[0].itertuples(),
                                            both_df[1].itertuples()):
                    if abs((objA.time-objB.time).total_seconds()) <= (60.0*td_max):
                        this_gen_len += 1
                return this_gen_len

            gg = gen()
            # pdb.set_trace()
            this_gen_len = count_gen()

        elif method == 2:
            tups = {}
            # tups will be used to hold a list of each sub_df's objects
            for nd, df in both_df.items():
                # nd is indicator of dictA/B
                tups[nd] = []
                # Below, looping over all objects - doesn't parallelise
                # unless you use multiprocess, the fork that allows
                # dill instead of pickle. That's because itertuples() is NamedTuple.
                for o in df.itertuples():
                #for o in df.iterrows():
                # for o in df.iloc[]
                    tups[nd].append(o)
                    # pdb.set_trace()

            # Number of objects for this time in each field
            nobj = {}
            listoflists = []
            for nd in (0,1):
                td = tups[nd]
                nobj[nd] = len(td)
                listoflists.append(td)

            def gen():
                for objA, objB in itertools.product(*listoflists):
                    # Find valid time of objA
                    # Subset objB (is it a DF or Series?) for only those times
                    yield objA, objB, self.bmap
            gg = gen()

        if method != 3:
            time0 = time.time()

            # this_gen_len = sum(1 for x in (gg))
            # Estimate of chunk size
            cs = math.ceil(this_gen_len/self.ncpus)
            # pdb.set_trace()
            cs2 = math.ceil(this_gen_len/(10*self.ncpus))

            from multiprocess import Pool as mpPool
            print("About to parallelise!")
            if self.ncpus > 1:
                # with multiprocessing.Pool(self.ncpus) as pool:
                with mpPool(self.ncpus) as pool:
                # with mpPool(maxtasksperchild=1) as pool:
                    # results = pool.imap_unordered(self.parallel_match_verif,gg)
                    results = pool.map(self.parallel_match_verif,gg,chunksize=cs)
                    # results = pool.map(self.parallel_match_verif,gg)

            else:
                results = []
                for oo in gg:
                    results.append(self.parallel_match_verif(oo))
            print("Computation done.")

            time1 = time.time()
            dt = time1-time0
            dtm = int(dt//60)
            dts = int(dt%60)

            print(f"TI calculation for all times took {dtm:d} min  {dts:d} sec.")
        else:
            # results = utils.load_pickle("./debug_match.pickle")
            pass

        member_TIs = {}
        # pdb.set_trace()
        A_set = set()
        B_set = set()
        # if (len(A_set) == 0) or (len(B_set)==0):
            # return (None,None)
        for r in results:
            # ( (intA, intB), TI)
            # try:
            compare_tup, TI = r
            # except:
            #     return (None,None)
            member_TIs[compare_tup] = TI
            A_set.add(compare_tup[0])
            B_set.add(compare_tup[1])

        # A_set is all unique objects in dictA data, and vice versa for B_set

        # Can we just loop over A_set, or do we need to do union/intersection
        # of A_set and B_set?

        # Ready for verifying:
        _a = 0
        _b = 0
        _c = 0
        _d = 0 # Probably not going to use?

        # Find the maximum TI for each obj and the pair that produced it
        matched_pairs = {}
        for objID in A_set.union(B_set):
            all_TIs = []
            all_otherIDs = []
            # Find all comparisons that occurred with this object
            for comp, TI in member_TIs.items():
                # comp is the tuple of (objA,objB)
                if objID in comp:
                    all_TIs.append(TI)
                    otherID = comp[0] if comp[0] != objID else comp[1]
                    all_otherIDs.append(otherID)

            # matched_pairs[objID] = max(member_TIs.items(), key=operator.itemgetter(1))[0]
            maxTI = N.max(N.array(all_TIs))
            if maxTI < 0.2:
                matched_pairs[objID] = None
            else:
                maxidx = N.argmax(N.array(all_TIs))
                # otherID = all_otherIDs.index(maxidx)
                otherID = all_otherIDs[maxidx]
                matched_pairs[objID] = (otherID,maxTI)


            already_hit = []
            # Verify now
            # "fcst" is A and "obs" is B - convention for 2x2 table.
            if objID in A_set:
                assert objID not in B_set

                # If not match, this is a false alarm
                if matched_pairs[objID] is None:
                    _b += 1
                else:
                    already_hit.append((objID,otherID))
                    # If matched, this is a hit
                    _a += 1

            elif objID in B_set:
                assert objID not in A_set

                # If matched, this is a hit (already counted?)
                if (otherID,objID) in already_hit:
                    print("Already matched 1!")
                    pass
                elif (objID,otherID) in already_hit:
                    print("Already matched 2!")
                    pass
                elif matched_pairs[objID] is None:
                    # If not matched, this is a missed hit
                    _c += 1
                else:
                    _a += 1


        CONT = DetScores(a=_a,b=_b,c=_c,d=0)
        # pdb.set_trace()
        return matched_pairs, CONT

    def parallel_match_verif(self,oob,td_max=20.0):
        objA, objB, bmap = oob
        compare_tup = (int(objA.megaframe_idx),int(objB.megaframe_idx))
        TI = utils.compute_total_interest(bmap,propA=objA,propB=objB,
                                            td_max=td_max,)
        # pdb.set_trace()
        return (compare_tup, TI)

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

    def compute_new_attributes(self,df_in,do_suite="W",rot_exceed_vals=None,
                                suffix=None):
        """ Adds new attributes from a forecast field to each object dataframe.

        data_in columns: fcst_vrbl, valid_time, fcst_min, prod_code, path_to_pickle
        (all that's needed to find unique matches)

        TODO:
            * do_suite for UH/AWS?

        """
        if not suffix:
            suffix = "_" + do_suite

        def get_commands(df_in):
            commands = []
            for fidx,fcst_df in enumerate(df_in.itertuples()):
                # JRL: need to change this for obs verification?
                # if (fcst_df.fcst_min == 0):
                    # continue
                utils.print_progress(total=len(df_in),idx=fidx,every=500)
                obj_df = self.lookup_obj_df(self.megaframe,
                                            valid_time=fcst_df.valid_time,
                                            fcst_min=fcst_df.fcst_min,
                                            prod_code=fcst_df.prod_code)
                if obj_df.shape[0] == 0:
                    # print("Skipping: no object found.")
                    # pdb.set_trace()
                    continue
                if do_suite == "W":
                    commands.append([obj_df, fcst_df.path_to_pickle, fidx])
                elif do_suite.startswith("UH"):
                    commands.append([obj_df, fcst_df.path_to_pickle, fidx,
                                        rot_exceed_vals,do_suite])
                else:
                    raise Exception
            return commands

        assert do_suite in ("UH02","UH25","W")
        if do_suite == "W":
            func = self.get_updraught_attributes
        elif do_suite.startswith("UH"):
            assert rot_exceed_vals
            func = self.get_rotation_attributes

        do_pickle = True if suffix != do_suite else False
        fname = f"commands{suffix}.pickle"
        fpath = os.path.join(self.tempdir,fname)
        print("Looking for",fpath)
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
        #try:
        df_out = pandas.concat(results,ignore_index=True)
        #except ValueError: # When no objects are identified
        #    return
        #else:
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

    def get_data_slice(self,obj,field):
        minr = int(obj.min_row)
        maxr = int(obj.max_row)
        minc = int(obj.min_col)
        maxc = int(obj.max_col)

        # Max/mean updraught in this object's bbox
        objidx = slice(minr,maxr),slice(minc,maxc)
        # objidx = obj['coords']
        return field[objidx]

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

            W_slice = self.get_data_slice(obj,W_field)

            new_df.loc[oidx,'max_updraught'] = N.nanmax(W_slice)
            #N.mean(W_slice)

            # Location of max updraught in object
            maxmax = N.where(W_slice == N.nanmax(W_slice))
            maxur = maxmax[0][0]
            maxuc = maxmax[1][0]

            minc = int(obj.min_col)
            minr = int(obj.min_row)
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

    def get_rotation_attributes(self,i):
        """
        rot_exceed_vals is either

        (1) a list/tuple of uh/aws values
        to exceed. The property name is rot_exceed_ID_#, where the last
        character(s) is the index of rot_exceed_vals. This way, the user
        can e.g., specify percentiles in terms of absolute values.

        or

        (2) a Threshs object, which is currently BROKEN! Needs to be implemented
        from VSE_dx project TODO JRL
        """
        obj_df, rot_field, fidx, rot_exceed_vals, layer = i
        if layer.endswith("02"):
            layer = "02"
        elif layer.endswith("25"):
            layer = "25"
        else:
            raise Exception


        if isinstance(rot_field,str):
            rot_field = N.load(rot_field)

        mega_idx = obj_df.head(1).megaframe_idx.values[0]

        if False:
            prod = obj_df.head(1).prod_code.values[0]
            t = obj_df.head(1).time.values[0]
            lead_time = obj_df.head(1).lead_time.values[0]
            print("About to compute rotation attributes for:",
                prod,t,"||| lead time =", lead_time,"min.")

        DTYPES = {
                "megaframe_idx_test":"i4",
                }

        nobjs = len(obj_df)
        new_df = utils.do_new_df(DTYPES,nobjs,)

        for oidx,obj in enumerate(obj_df.itertuples()):
            dx = obj.dx
            if obj.member.startswith("m"):
                vrbl = "UH" + layer
            else:
                vrbl = "AWS" + layer

            if layer == "02":
                rot = "lowrot"
            elif layer == "25":
                rot = "midrot"
            else:
                raise Exception

            # Get index for object in megaframe
            # oidx = obj.label

            # Get index for object in new df
            #Ix = obj.Index
            Ix = obj.megaframe_idx
            new_df.loc[oidx,"megaframe_idx_test"] = Ix

            if "AWS" in vrbl:
                # Not matching
                return new_df

            rot_slice = self.get_data_slice(obj,rot_field)

            new_df.loc[oidx,f'max_{rot}'] = N.nanmax(rot_slice)
            #N.mean(rot_slice)

            # Location of max updraught in object
            maxmax = N.where(rot_slice == N.nanmax(rot_slice))
            maxur = maxmax[0][0]
            maxuc = maxmax[1][0]

            minc = int(obj.min_col)
            minr = int(obj.min_row)
            # Find this location back on the main grid
            new_df.loc[oidx,f'max_{rot}_row'] = maxur + minr
            new_df.loc[oidx,f'max_{rot}_col'] = maxuc + minc

            # Min/mean
            new_df.loc[oidx,f'mean_{rot}'] = N.nanmean(rot_slice)
            new_df.loc[oidx,f"min_{rot}"] = N.nanmin(rot_slice)

            # Relative to centroid? Could be radius, using equivalent circle
            # distance_from_centroid
            cr = obj.centroid_row
            cc = obj.centroid_col
            dist_km, angle = utils.distance_angle_from_coords(maxur+minr,
                                                    maxuc+minc,cr,cc,dx=dx)
            new_df.loc[oidx,f'{rot}_distance_from_centroid'] = dist_km

            # angle_from_centroid
            new_df.loc[oidx,f'{rot}_angle_from_centroid'] = angle

            # Exceedence of various percentiles
            # Would need to know if UH or AWS, and which level (0-2, 2-5 km etc)
            # Also, do for exceedence of magnitudes.
            if not isinstance(rot_exceed_vals,(list,tuple)):
                TH = rot_exceed_vals
                dx = int(obj.dx)

                evs = TH.get_threshs(vrbl=vrbl,fmt=dx)


            for n, v in enumerate(evs):
                prop_name = f"{rot}_exceed_ID_{n}"

                if N.where(rot_slice > v)[0].size > 0:
                    answer = True
                else:
                    answer = False

                new_df.loc[oidx,prop_name] = answer
                new_df.loc[oidx,prop_name+"_val"] = v
            # pdb.set_trace()

        # print("Generated vrbl stats.")
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

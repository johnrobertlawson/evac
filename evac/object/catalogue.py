import os
import pdb
import multiprocessing
import pickle

import matplotlib as M
import matplotlib.pyplot as plt

import numpy as N
import pandas
import sklearn.preprocessing
import sklearn.decomposition

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

    def do_matching(self,models,times):
        """ Use of Total Interest to match objects in time/across models or obs.

        Maybe do CombinedObject class. This could score itself!

        Need a lat/lon grid for each product.
        """
        # Loop over a given df (one time, one model, many objects)
        # Loop over 2x product of objects, within a "sphere" of times and models
        for inittime, validtime, model in itertools(inittimes,times,models):
            ti_data = []
            combo = []
            # Need to create a dataframe that slices all possible objects
            for objA, objB in itertools.combinations(objects,2):
                TI = self.compute_total_interest(propA=objA,propB=objB)
                combo.append((objA,objB))
                ti_data.append(TI)
            # Find the maximum TI and the pair that produced it
            idx = N.argmax(N.array(ti_data))
            the_combo = combo[idx]


    def compute_total_interest(self,propA=None,propB=None,cd=None,md=None,td=None,
                            cd_max=40.0,md_max=40.0,td_max=25.0):
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
        # Sanity checks
        if propA is not None:
            assert propB is not None
            assert not all(cd,md,td)
            method = 1
        else:
            assert propB is None
            assert all(cd,md,td)
            method = 2
        assert all(cd_max,md_max,td_max)

        if method == 1:
            # Compute cd, md, td
            cd = 0
            md = 0
            td = 0

        TI = ((td_max-td)/t_max)*0.5*(((cd_max-cd)/cd_max)+((md_max-md)/md_max))
        return TI

    def compute_new_attributes(self,df_in,do_suite="W"):
        """ Adds new attributes from a forecast field to each object dataframe.

        data_in columns: fcst_vrbl, valid_time, fcst_min, prod_code, path_to_pickle
        (all that's needed to find unique matches)

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
        fpath = os.path.join(self.tempdir,fname)
        if os.path.exists(fpath):
            print("Loading pickle of commands")
            with open(fpath,"rb") as f:
                commands = pickle.load(file=f)
        else:
            print("Generating list of commands")
            commands = list(get_commands(df_in))
            with open(fpath,"wb") as f:
                pickle.dump(obj=commands,file=f)

        # pdb.set_trace()
        # Submit in parallel
        print("Submitting jobs to compute {} attributes...".format(do_suite))
        if self.ncpus > 1:
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(func,commands)
        else:
            for i in commands:
                func(i)
        print("Jobs done.")

        df_out = pandas.concat(results,ignore_index=True)
        # pdb.set_trace()
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
                #"megaframe_idx":"i4",

                "max_updraught":"f4",
                "max_updraught_row":"f4",
                "max_updraught_col":"f4",

                "min_updraught":"f4",
                "mean_updraught":"f4",


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
            #new_df.loc[oidx,"megaframe_idx"] = Ix

            minr = int(obj.min_row)
            maxr = int(obj.max_row)
            minc = int(obj.min_col)
            maxc = int(obj.max_col)

            # Max/mean updraught in this object's bbox
            objidx = slice(minr,maxr),slice(minc,maxc)
            # objidx = obj['coords']
            W_slice = W_field[objidx]

            new_df.loc[oidx,'max_updraught'] = N.max(W_slice)
            #N.mean(W_slice)

            # Location of max updraught in object
            maxmax = N.where(W_slice == W_slice.max())
            maxur = maxmax[0][0]
            maxuc = maxmax[1][0]
            # Find this location back on the main grid, I think
            new_df.loc[oidx,'max_updraught_row'] = maxur
            new_df.loc[oidx,'max_updraught_col'] = maxuc
            #new_df.loc[oidx,'max_updraught_row'] = maxur + minr
            #new_df.loc[oidx,'max_updraught_col'] = maxuc + minc

            # Min/mean
            new_df.loc[oidx,'mean_updraught'] = N.mean(W_slice)
            new_df.loc[oidx,"min_updraught"] = N.min(W_slice)

            # Relative to centroid? Could be radius, using equivalent circle
            # distance_from_centroid
            cr = obj.centroid_row
            cc = obj.centroid_col
            dist_km, angle = utils.distance_angle_from_coords(maxur,maxuc,cr,
                                                                    cc,dx=dx)
            new_df.loc[oidx,'distance_from_centroid'] = dist_km

            # angle_from_centroid
            new_df.loc[oidx,'angle_from_centroid'] = angle

            # distance_from_wcentroid

            # Size of updraught, using W = ? m/s as threshold
            # updraught_area_km
            # updraught_width_km
            # pdb.set_trace()
            # pdb.set_trace()

        return new_df

    def do_pca(self):
        """
        Todo:
            * Move this to PCA stat class. Return that class.
        """

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

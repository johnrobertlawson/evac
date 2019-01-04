import os
import pdb

import matplotlib as M
import matplotlib.pyplot as plt

import numpy as N
import pandas
import sklearn.preprocessing
import sklearn.decomposition

import evac.utils as utils

class Catalogue:
    def __init__(self,df_list,): #time_list=None,prod_list=None):
        """
        Args:
            df_list: list of DataFrames, or one big one.

        Note:
            * The time and product/model/format of each object's dataframe object
                is present in the dataframe itself. This way, a mega-array is
                generated that can be sliced, or matching can be performed.
        """
        self.df_list = df_list
        self.ndfs = len(df_list)

    def do_matching(self,models,times):
        """ Use of Total Interest to match objects in time/across models or obs.

        Maybe do CombinedObject class. This could score itself!

        Need a lat/lon grid for each product.
        """
        # Loop over a given df (one time, one model, many objects)
        # Loop over 2x product of objects, within a "sphere" of times and models
        for inittime, validtime, model:
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

        TI = ((td_max-td)/t_max)*0.5*(((cd_max-cd)/cd_max)+((md_max-md)/md_max))
        return TI

    def add_new_attributes(self,field_list,time_list,model_list,do_suite="W"):
        """ Adds new attributes from a forecast field to each object dataframe.

        For instance, objects are identified in reflectivity. To look at
        attributes related to updrafts, add in W forecast data.

        Each field must be of identical lat/lon dimensions as its corresponding
        object dataframe, such that look-up of values is correct.

        Because the memory usage might be high, "field_list" could be a list
        of absolute paths to open numpy arrays.
        """
        for field, utc, model in zip(field_list,time_list,model_list):
            if isinstance(field,str):
                data = N.load(field)
            elif isinstance(field,N.ndarray):
                data = field
            else:
                raise Exception

            # Do stuff here
            if do_suite == "W":
                self.get_updraught_attributes()
            elif do_suite == "UH02":
                self.get_uh_attributes(lv='02')
            elif do_suite == "UH25":
                self.get_uh_attributes(lv='25')
        return

    def get_updraught_attributes(self):
        """
        Todo:
            * Updraught max level (hPa, km AGL), magnitude (m/s)
            * Updraught location along major axis or radius?
        """
        pass

    def do_pca(self):
        """
        Todo:
            * Move this to PCA stat class. Return that class.
        """
        if isinstance(self.df_list,pandas.DataFrame):
            self.megaframe = self.df_list
        else:
            self.megaframe = pandas.concat(self.df_list,ignore_index=True)
        features = ['area','eccentricity','extent','max_intensity',
                        'mean_intensity','perimeter','longaxis_km']
        x = self.megaframe.loc[:,features].values
        y = self.megaframe.loc[:,['prod']].values
        scaler = sklearn.preprocessing.StandardScaler()
        x = scaler.fit_transform(x)

        pca = sklearn.decomposition.PCA(n_components=3)
        PCs = pca.fit_transform(x)
        cols = ['PC1','PC2','PC3']
        _PC_df = pandas.DataFrame(data=PCs,columns=cols)
        PC_df = pandas.concat([_PC_df, self.megaframe[['prod']]], axis=1)

        return pca, PC_df, features, scaler

    def do_lda(self,train_data,targets):
        """
        Args:
            train_data: this has the training data
            targets: the classes belonging to each data point
        """
        self.megaframe = pandas.concat(self.df_list,ignore_index=True)
        features = ['area','eccentricity','extent','max_intensity',
                        'mean_intensity','perimeter','longaxis_km']
        X = self.megaframe.loc[:,features].values

        lda = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()
        hmm = lda.transform(X)

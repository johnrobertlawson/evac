import os
import operator
import pdb
import itertools
import multiprocessing
import functools
import random

import numpy as N
from scipy import signal

import evac.utils.utils as utils

class FI:
    def __init__(self,xa,xfs,thresholds,decompose=True,
                    neighborhoods='auto',temporal_window=1,
                    ncpus=1,efss=False):
        """Fractional Ignorance.

        The "continuous" part of this scheme is not over thresholds in the raw
            data (e.g., mm for rainfall), but for fractional values. There's
            a good reason for being able to choose thresholds: the user
            may only be interested in running FI for flash-flood cases, for
            example.

        Returns:
            dictionary of the form:

            dict[threshold][neighborhood][temporal_window][timestep][score] = x

            There is N.NaN returned for the first and last (W-1)/2 timesteps
            for temporal window W. The score keys are strings of FI, UNC,
            REL, RES, FISS (the skill score).

        Args:
            xa: 3-D observation array of size (ntimes,nlat, nlon).
            xfs: 4-D forecast array of size (N,ntimes,nlat,nlon) where N is
                ensemble size and ntimes is the number of time steps.
            thresholds (tuple,list): thresholding to apply to input data.
            neighborhoods: if 'auto', use smallest to largest spatial
                kernels. These can accessed after processing with the
                self.neighborhoods attribute, or the keys of the
                returned dictionary.
            temporal_window (int): by default, no temporal
                windowing is performed. This value
                must be integers of odd numbers only.
            decompose (bool): if True, return a named tuple of
                uncertainty (UNC), reliability (REL), and resolution (RES)
                dictionaries. The value of FI is the sum of all three.

        TODO:
            * The temporal window should be constant, such that it matches the
                number of times in the xfs/xa data.

        Returns:
            Either FI or the decomposition as a dictionary.
        """
        self.efss = efss

        assert xfs.shape[1:] == xa.shape
        if (xfs.ndim == 3) and (xa.ndim == 2):
            xfs = xfs[:,N.newaxis,:,:]
            xa = xa[N.newaxis,:,:]
        self.xa = xa
        self.xfs = xfs

        self.nens, self.nt, self.nlats, self.nlons = self.xfs.shape

        # 0 -> 0.01 and 1 -> 0.99 (to address underdispersion).
        self.probthreshs = self.no_binary(N.linspace(0,1,self.nens+1))

        check = self._assert_list(thresholds)
        self.thresholds = thresholds

        if neighborhoods == 'auto':
            # The biggest box we can fit in the domain
            # maxwindow = N.floor(min(self.nlats,self.nlons)/2)
            maxwindow = N.floor(min(self.nlats,self.nlons))
            # neighborhoods = N.arange(1,maxwindow,2,dtype=int)
            neighborhoods = N.linspace(1,maxwindow,num=15,dtype=int)
        check = self._assert_list(neighborhoods)
        self.neighborhoods = neighborhoods

        assert (temporal_window % 2 == 1)
        self.temporal_window = temporal_window

        # The "centre" of observations given the time window
        self.tidx = int((self.temporal_window-1)/2)

        self.decompose = decompose

        # Create results dictionary
        self.results = {th: {n: { tw: {} for tw in (self.temporal_window,)}
                        for n in self.neighborhoods} for th in self.thresholds}

        self.ncpus = ncpus
        namedtup = self.parallelise()
        # return self.results

    def parallelise(self):
        # Loop over all thresholds, neighbourhood sizes, and temporal windows
        loop = self.generate_loop()

        if self.ncpus == 1:
            for i in loop:
                thresh, neigh, tempwindow = i
                result = self.compute_fi(i)
                self.update_results(result)
        else:
            # Shuffling will optimise, because larger
            # neighbourhoods will take longer.
            with multiprocessing.Pool(self.ncpus) as pool:
                results = pool.map(self.compute_fi,loop)
            for result in results:
                self.update_results(result)
                # scores, thresh, neigh, tempwindow, tidx = tup
                # self.results[thresh][neigh][tempwindow][tidx] = fi

    def update_results(self,result):
        scores, thresh, neigh, tempwindow = result
        # for score in ("REL","RES","UNC","FI","FISS"):
        for score,val in scores.items():
            self.results[thresh][neigh][tempwindow][score] = val
                            # thresh=thresh,neigh=neigh,
                            # tempwindow=tempwindow,tidx=tidx,)
        return


    def compute_fi(self,i):
        """
        Consider coding a continuous version (loop over all possible
        fractions possible in ensemble members) or simply ranked
        (where thresholding is done earlier).
        """
        thresh, neigh, tempwindow = i
        print("Calculating FI for threshold = {}, neighborhood = {},"
                " tempwindow = {}, tidx = {}.".format(thresh,neigh,
                tempwindow, self.tidx))

        # First threshold the data into a binary array
        # Whether each gridpoint is over (1) or under (0) the threshold
        # Need to gather times either side of tidx for temporal window
        tidxs = slice(self.tidx-int((tempwindow-1)/2),
                        1+self.tidx+int((tempwindow-1)/2))
        Im = N.where(self.xfs[:,tidxs,:,:] > thresh, 1, 0)
        Io = N.where(self.xa[tidxs,:,:] > thresh, 1, 0)

        if self.efss:
            # "Bonus" score for efss
            M_kernel_efss = N.ones([self.nens,tempwindow,neigh,neigh])
            M_efss_raw = signal.fftconvolve(Im,M_kernel_efss,mode='same')
            total = (neigh**2) * tempwindow * self.nens
            M_efss = N.abs(N.around(M_efss_raw)/total)
            assert N.all(M_efss >= 0.0)
            assert N.all(M_efss <= 1.0)

        # Next, apply the 3D kernel
        M_kernel = N.ones([1,tempwindow,neigh,neigh])
        M = signal.fftconvolve(Im,M_kernel,mode='same')

        O_kernel = N.ones([tempwindow,neigh,neigh])
        O = signal.fftconvolve(Io,O_kernel,mode='same')
        # Note that M and O range from 0 to neigh**2

        # Round to integers then divide into fractional values
        # Get rid of any weird negatives, etc
        # TODO: is there a more elegant way to do this?
        total = (neigh**2) * tempwindow
        # M = N.abs(N.around(M)/(neigh**2))
        M = N.abs(N.around(M)/total)
        # O = N.abs(N.around(O)/(neigh**2))
        O = N.abs(N.around(O)/total)

        assert N.all(M >= 0.0)
        assert N.all(M <= 1.0)
        assert N.all(O >= 0.0)
        assert N.all(O <= 1.0)

        frac_get_method = 2

        if frac_get_method == 1:
            # As we're using fractional arrays, the distance is constant,
            # see eq. 27/28 in Toedter and Ahrens 2010
            # (e.g., spatial window = 3 means 0, 1/9, 2/9...1)
            if neigh == 1:
                distances_sw = set(N.array([0,1]))
            else:
                distance_sw = 1/(neigh**2)
                distances_sw = set(N.arange(0,1+distance_sw,distance_sw))

            # TODO: Needs amending for temporal windowing
            # distance would not be constant, for a start.
            # would have to merge sets of temporal fractions
            # (e.g., temp window = 3 means 0, 0.5, 1)

            if tempwindow == 1:
                distances_tw = set(N.array([0,1]))
            else:
                distance_tw = 1/(tempwindow-1)
                distances_tw = set(N.arange(0,1+distance_tw,distance_tw))

            fracset = sorted(functools.reduce(operator.or_, (distances_sw, distances_tw)))
        elif frac_get_method == 2:
            fracset = N.unique(N.hstack((N.unique(M),N.unique(O))))

        # Next, loop over all possible fractional values
        # Calculate REL, RES, UNC for each
        REL = N.zeros([len(fracset)])
        RES = N.copy(REL)
        UNC = N.copy(REL)

        # pdb.set_trace()
        # for fidx, fracval in enumerate(N.linspace(0,1,num=num)):
        for fidx, fracval in enumerate(fracset):
            REL[fidx] = self.compute_rel(M,O,fracval)
            RES[fidx] = self.compute_res(M,O,fracval)
            UNC[fidx] = self.compute_unc(O,fracval)
        FISS = (RES-REL)/UNC
        FI = REL - RES + UNC


        # TODO
        # Outside the domain (e.g. using a large neigh)
        # IS this nans or what? Bias in the score?
        # Does it affect uncertainty?

        frac_distances = N.diff(fracset)
        # pdb.set_trace()

        rel = N.nansum(frac_distances * REL[:-1])
        res = N.nansum(frac_distances * RES[:-1])
        unc = N.nansum(frac_distances * UNC[:-1])
        fiss = N.nansum(frac_distances * FISS[:-1])
        fi = N.nansum(frac_distances * FI[:-1])

        # fi = rel - res + unc
        # fiss = (res - rel)/unc

        scores = dict(REL=rel,RES=res,UNC=unc,FI=fi, FISS=fiss)
        print("Scores: REL = {:.3f}, RES = {:.3f}, UNC = {:.3f}; FI = {:.3f}; FISS = {:.3f}".format(
                scores['REL'],scores['RES'],scores['UNC'],scores['FI'],scores['FISS']))

        if self.efss:
            MO_diff = (M_efss - O)**2
            FBS_ref = M_efss**2 + O**2
            efss = 1 - (N.nanmean(MO_diff) / N.nanmean(FBS_ref))
            scores['eFSS'] = efss
            print("Bonus eFSS = ",efss)

        return scores, thresh, neigh, tempwindow

    def eq_14_16(self,M,O,f,func):
        """ Eqs 14 and 16 are very similar in Toedter and Ahrens
            (2012, MWR), so this combines the logic.
        """
        # enforce 2D obs
        O = O[0,:,:]
        M = M[:,0,:,:]
        Mf = utils.exceed_probs_2d(M,f,overunder='over',fmt='decimal')
        Mf = self.no_binary(Mf)
        results = N.zeros_like(self.probthreshs)
        for iidx,yi in enumerate(self.probthreshs):
            widx = N.where(Mf == yi)

            # pyi = N.count(Mf[widx] == yi)
            pyi = self.compute_pyi(Mf,yi)
            zi = self.compute_zi(O[widx],f,yi)
            z = self.compute_z(O,f)
            # pdb.set_trace()
            results[iidx] = func(pyi=pyi,yi=yi,zi=zi,z=z)
        # return N.sum(results)
        return N.nansum(results)

    def no_binary(self,arr):
        """ Avoids prob forecasts of 0% or 100%.

        This then avoids infinite values of Ignorance.
        """
        arr[arr<0.01] = 0.01
        arr[arr>0.99] = 0.99
        return arr

    def compute_pyi(self,probs,yi):
        """ Freq. of probs equal to yi.
        """
        pyi = (N.where(probs == yi)[0].size)/probs.size
        # pdb.set_trace()
        return pyi

    def compute_z(self,O,f):
        """
        z is the prob of observed occurrence in the
        sample (0-1 frequency).
        """
        z = (N.where(O > f)[0].size)/O.size
        # z = (N.where(O == f)[0].size)/O.size
        return z

    def compute_zi(self,o,f,yi):
        """
        Args:
            o: subset of fractional array O that were
                points with forecast prob yi.
            f: fraction value of interest (later, integrated)
            yi: forecast prob of interest (later, integrated)

        A well calibrated ensemble has zi == yi.
        It is the observed (0-1) frequency of a ob
        that occurred for prob forecasts for the ith
        probability (between 0 and 1).

        For a given yi (e.g. 70% or 0.7), what percentage
        of those points verified?

        Conditional freq of occurrence on all occasions where
        yi was forecasted.
        """
        try:
            # zi = (N.where(o == f)[0].size)/o.size
            zi = (N.where(o > f)[0].size)/o.size
        except ZeroDivisionError:
            zi = N.nan
        # pdb.set_trace()
        return zi

    def compute_rel(self,M,O,f):
        """
        Loop over all bins of probability
        For each, if the M grid point

        Generate probs of exceedence for the given f.
        bar zi = yi ... bar zi is frequency of
        obs for a given bin yi.
        """
        def rel_eq(pyi,yi,zi,**kwargs):
            """ Eq. 14 in Toedter and Ahrens 2012 MWR,
            without the summation.

            Shouldn't have divide-by-zero problems
            because yi is never 0 or 1 (we apply
            no_binary() to all data).
            """
            # try:
            rel = pyi * (zi*N.log2(zi/yi) +
                    (1-zi)*N.log2((1-zi)/(1-yi)))
            # except ZeroDivisionError:
                # rel = N.nan
            return rel
        REL = self.eq_14_16(M,O,f,rel_eq)
        # pdb.set_trace()
        return REL

    def compute_res(self,M,O,f):
        def res_eq(pyi,zi,z,**kwargs):
            """ Eq. 16 in Toedter and Ahrens 2012 MWR,
            without the summation.

            If z = 1, divide by zero error. This means all
            points have passed a threshold in the raw data.
            This is because there is no variability,
            so forecasts can't be rated for their ability
            to capture it.
            """
            try:
                res = pyi * (zi*N.log2(zi/z) +
                    (1-zi)*N.log2((1-zi)/(1-z)))
            except ZeroDivisionError:
                res = N.nan
            return res
        return self.eq_14_16(M,O,f,res_eq)

    def compute_unc(self,O,f):
        """
        z is the frequency of f in the sample O.
        """
        # z = N.count(O==f)/O.size
        z = self.compute_z(O,f)
        UNC = -z*N.log2(z) - (1-z)*N.log2(1-z)
        return UNC

    def generate_loop(self):
        """ The randomisation is to partly optimise the parallelisation.
        """
        th_rand = random.sample(list(self.thresholds),len(self.thresholds))
        n_rand = random.sample(list(self.neighborhoods),len(self.neighborhoods))
        for thresh, neigh in itertools.product(th_rand,n_rand):
            yield thresh, neigh, self.temporal_window

    def _assert_list(self,l):
        if isinstance(l,(N.ndarray,list,tuple,set)):
            return True
        else:
            raise Exception("This is not a valid list, tuple etc.")

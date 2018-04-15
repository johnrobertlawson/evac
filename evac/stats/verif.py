""" Helper scripts to run stats over ensembles for a given
initialisation time (looping over levels, ensemble members, valid times, and 
lat/lon if required)

The idea is that no actual plotting is done herein, but that this class
calls the relevant class methods outside of `Verif` to do this.

Note:
    Logic for parallelised stats, for both domains on an ensemble,
    for a range of times, and for QPF, goes as follows:

        * User needs to create Ensemble, StageIV etc instances
        * --- USER INSTANTIATES VERIF ---
        * init takes: directory info for data, images
        * Reprojection grids formed
        * --- USER COMPUTES STATS ---
        * Set up pool of CPU workers
        * Generate iterable from stats, times, domains (fixed vrbl)
        * Need to find correct verification object and pass it in
        * Compute stats by reprojecting and etc
        * Save output to disk with complex naming
        * Join/close pool
        * --- USER CALLS PLOT SCRIPTS ---
        * For list of stats/variable, load data from above, and plot
        * Call master script, process kwargs (user settings like hold),
            plotargs (e.g., posterior settings like title), and
            mplargs (e.g., plot-time settings like cmap). This sets default
            that aren't specified.
        * Plot might be simple linegraph at one time
        * Plot may also compare multiple times/domains, with averages
        * Hence having lots of little .npy data files is good

Todo:
    * Order methods so interface comes top.
    * Create Fields() class like create_newgrid(), and attach to all obs,
        ensemble, wrfout classes as an attribute.
    * Protect methods with _ and __ where users won't use it (find/replace!)
    * Create an Obs class (things like QC data, save/load). This should be
        inherited by all datatypes that are also classes (mix-in). Make
        sure any methods have exclusive names to not clash (use _/__?)
    * Create a `evac.plot.thumbnails` class
    * Ensure all figure subtypes has a "print" call to notify of path to figure
    * Save numpy files of computed data for CPU-intense calculations.
    * Decorate all plotting/computation processes with timeme().
    * Be aware if saved data already exists - load if that's the case.
    * Wrap __get_plot_options into figure().
    * Parallelise!

"""
import os
import pdb
import copy
from multiprocessing import Pool, Lock
import itertools
import datetime

import numpy as N

# from evac.plot.thumbnails import Thumbnails
from evac.plot.linegraph import LineGraph
from evac.stats.detscores import DetScores
from evac.stats.probscores import ProbScores
from evac.plot.scorecard import ScoreCard
from evac.datafiles.obs import Obs
import evac.utils as utils
from evac.utils.reproject_tools import reproject, WRF_native_grid, create_new_grid, VerifGrid


class Verif:
    """ A suite of verification scripts.

    This is a helper class to simplify the running of common
    verification patterns. Arguments to the class instantiation
    should be an ensemble object, followed by any number of arguments
    that are instances of verification classes. This class will
    search through these argumemts when verification is called
    on a relevant method. If the verification is ambiguous, it must
    be specified in the method's keywords.

    Common arguments and keyword arguments for most plotting scripts are:
        * `bblims`, is a dictionary of Nlim, Elim, Slim, and Wlim, setting
            bounding box limits for 2-D plots.
    
    Also note that anything listed in the dictionaries mplargs and 
    mplkwargs is passed to matplotlib when plotting.


    Todo:
        * Each method should all option to change bounding area of the plot.
        * Avoid doing `E.get()` numerous times for `Ensemble` instances, if
            looping over plotting times. Data should be loaded first as 
            3/4/5D array.
        * Compare domains of the ensemble, by reprojecting to a common
            domain that is contained inside all domains. Then a figure
            can be plotted that compares scores across given lead times
            ( like a scorecard?).

    
    Example:
        Set up the instance::
            from evac.datatypes.ensemble import Ensemble
            from evac.datatypes.da
            E = Ensemble(...)
            ST4 = StageIV(...)
            V = Verif(E,ST4)

    Args:
        ensemble: Class instance of `evac.datatypes.ensemble.Ensemble`.
        obs: instance, or list/tuple of instances, of Obs subclasses
        outdir (optional): path to output files
        datadir (optional): Path to which .npy data files will be saved
        reproject_dict (dict): if None, no reprojection. Otherwise,
            reprojection to a common domain is done with these options

    """


    def __init__(self,ensemble,obs,outdir=False,datadir=False,
                    reproject_dict=None):
        self.E = ensemble

        # A tuple of obs objects, or just one.
        # So a bunch of Radar objects for different times is no problem,
        # as we can search for the right one.

        if isinstance(obs,(list,tuple)):
            self.obs = obs
        else:
            self.obs = (obs,)
            # if isinstance(obs,Obs):
                # self.obs = (obs,)
        self.obdict = self.generate_obdict()

        # These are optional, for later use
        self.outdir = outdir
        self.datadir = datadir

        # Might attach self.E attributes such as times and member 
        # names to the class instance?

        # NOT IMPLEMENTED
        # self.fcst_times = {d:E.times(dom=d) for d in self.doms}

        # Mutual forecast times for each domain
        # NOT IMPLEMENTED
        # self.fcst_times['intersection'] = intersection(self.fcst_times)

        # Compute dom information
        self.doms, self.ndoms = self.check_ensemble_domains()

        # Dictionary of aribtrary WRFOut objects 
        self.arbdict = self.make_arbdict()

        # lats and lons for each domain
        self.lats = {}
        self.lons = {}
        for dom in self.doms:
            self.lats[dom],self.lons[dom] = self.E.get_latlons(dom=dom)

        # Reprojection - define a new domain to reproject to
        # Check all settings are in reproject_dict?
        self.RD = reproject_dict
        self.newgrid = self.create_newgrid()

        # Look up stats methods here
        self._STATS = self.generate_lookup_table()

    
    ##### INTERFACE METHODS - COMPUTE #####
   
    def compute_rmse(self,):
        pass

    def compute_stats(self,stats,vrbl,verif_times,dom=1,ncpus=1,
                        lvs=None,**kwargs):
        """ Megascript to control all stats computations.

        Args:
            stats: list/tuple of stats to compute. See 
                :func:`lookup_stats` for a list.
            verif_times: one, or list/tuple of datetime.datetimes.
            dom (int): domain(s) to loop over
            lvs: levels to loop over - default is None for surface etc
            ncpus (int): number of CPUs for parallel computing.
            kwargs: settings to pass to compute.

        Todo:
            Implement options, assign to self and then it can be seen
                by all compute_* methods.

        Returns:
            A numpy save file is created for each chunk looped over.

        """
        self.start_your_engines(ncpus)

        # This is written silly for debugging
        verif_times1 = self.ensure_list(verif_times)
        verif_times2 = self.ensure_integers(verif_times1)
        verif_times = verif_times2

        # add kwargs to class
        self.allowed_compute_kwargs = ['crps_thresholds',]
        for k,v in kwargs.items():
            if k in self.allowed_compute_kwargs:
                setattr(self,k,v)
            else:
                raise Exception("{} is not a valid keyword arg,".format(k))

        lvs = self.ensure_list(lvs)
        stats = self.ensure_list(stats)
        domlist = self.get_domlist(dom)

        # Create iterable of times and stats to iterate over
        itr = self.create_stats_iterable(vrbl=vrbl,verif_times=verif_times,
                                        doms=domlist,lvs=lvs)

        for stat in stats:
            print("Generating {} stats now.".format(stat))
            statfunc = self.lookup_statfunc(stat)
            if ncpus == 1:
                # Serial mode, good for debugging
                for chunk in itr:
                    #statfunc(next(itr))
                    statfunc(chunk)
            else:
                # result = self.pool.map(statfunc,itr,)#chunksize=1)
                # self.pool.close()
                # self.pool.join()
                for chunk in itr:
                    # self.pool.apply_async(statfunc<F3>,chunk)

                    fchr,dom,lv,vrbl = itr

                    # get_fcst_data
                    fcst = self.E.get(vrbl,fcsthr=fchr,Vdom=dom,level=lv,accum_hr=1)
                    fdata = self.reduce_data_dims(fcst)

                    # get_ob_data
                    utc = self.lookup_validtime(fchr)
                    obobj,obname = self.get_ob_instance(vrbl,return_name)
                    obdata = self.reduce_data_dims(obobj.get(utc,lv=lv))

                    # One parallel loop (async it)
                    # try to load reproject
                    rpj_fpath_f = check_for_reproj(vrbl=vrbl,model='fcst',lv=lv,utc=utc,
                                            dom=dom,ens='all',return_fpath=True)
                    if not rpj_fpath_f:
                        # reproject if not, PARALLEL, and save
                        flats, flons = self.E.get_latlons(dom=dom)
                        x
                        xfs = self.do_reprojection(fdata,flons,flats,save=rpj_fpath_f)
                    else:
                        xfs = N.load(rpj_fpath_f)

                    rpj_fpath_o = check_for_reproj(vrbl=vrbl,model=obname,lv=lv,utc=utc,
                                            dom=None,ens=None,return_fpath=True)
                    if not rpj_fpath_0:
                        oblats = obobj.lats
                        oblons = obobj.lons
                        xa = self.do_reprojection(obdata,oblons,oblats,save=rpk_fpath_o)
                    else:
                        xa = N.load(rpk_fpath_o)

                    # Next parallel loop
                    # Run statfunc and exit

                    statfunc = self.compute_crps_mp
                    kw = dict(xfs=xfs,xa=xa,vrbl=vrbl,lv=lv,fchr=fchr,dom=dom)
                    self.pool.apply_async(statfunc,args=chunk,kwargs=kw)
        return

    def compute_crps_mp(self,xfs,xa,vrbl,lv,fchr,dom,):
        P = ProbScores(xfs=xfs,xa=xa)
        crps = P.compute_crps(self.crps_thresholds)

        fpath = self.generate_npy_fname(vrbl,'CRPS',lv,fchr,
                            dom,'mean',fullpath=True)
        self.save_npy(crps,vrbl,'CRPS',lv,fchr,dom,'mean',fpath=fpath)
        return

    def compute_crps(self,itr):
        fchr,dom,lv,vrbl = itr

        # Check we have options 
        # This is something like 0 to 100 for 1-h QPF

        # TODO
        # Try again when not parallel
        # assert 'crps_thresholds' in self.__dict__

        # Get raw data

        # If this was a detscore, we'd need to loop over members.
        # Save for each.
        # TODO: implement fchr in E.get to look up accum_precip,
        # or instantaneous variables (validtime? utc?)
        # TODO: Not hard code QPF accum.
        fcst = self.E.get(vrbl,fcsthr=fchr,Vdom=dom,level=lv,accum_hr=1)
        fdata = self.reduce_data_dims(fcst)
        

        # Load verification
        utc = self.lookup_validtime(fchr)
        # Make lv ignored if it's nonsense
        # Need to implement get() in all obs type.
        # Do this in parent class and override.
        obobj = self.get_ob_instance(vrbl)
        obdata = self.reduce_data_dims(obobj.get(utc,lv=lv))

        # Reproject and compute
        flats, flons = self.E.get_latlons(dom=dom)
        xfs = self.do_reprojection(fdata,flons,flats)

        oblats = obobj.lats
        oblons = obobj.lons
        xa = self.do_reprojection(obdata,oblons,oblats)
        P = ProbScores(xfs=xfs,xa=xa)
        # for thresh in self.crps_thresholds:
        crps = P.compute_crps(self.crps_thresholds)

        # Save to disk
        # def save_npy(self,data.vrbl,score,lv,fchr,dom,ens,):
        fpath = self.generate_npy_fname(vrbl,'CRPS',lv,fchr,
                                dom,'mean',fullpath=True)
        self.save_npy(crps,vrbl,'CRPS',lv,fchr,dom,'mean',fpath=fpath)
        return

    ##### INTERFACE METHODS - PLOT #####

    def plot_domains(self,ensemble_domains='all',ob_domains='all'):
        """ Wrapper to plot domains of forecast or observation domains.
        """
        pass
        

    def plot_violin(self,vrbl,*args,**kwargs):
        """ Plot violin plot

        Args:
            vrbl (str): name of variable to evaluate

        This is an alternative to box-and-whisker.
        """
        self.__get_plot_options(*args,**kwargs)
        self.create_dict(models='all')
        VP = ViolinPlot()
        VP.plot()
        

    def plot_thumbnails(self,vrbl,ensmembers='all',ob='auto',
                        mplargs=None,mplkwargs=None,verif_first=True,
                        *args,**kwargs):
        """ Plot postage stamps or thumbnails of ensemble members listed
        (by default, this is all).

        Args:
            verif_first (bool): If True, the first axis (thumbnail/subplot)
                will be the observations. Method exits if verification does
                not exist in self.obdict.
        """
        self.__get_plot_options(*args,**kwargs)
        if verif_first:
            self.get_ob_instance(vrbl)
        TN = Thumbnails(verif_first=verif_first)
        B
        return

    def plot_scorecard(self,detscores='all',probscores='all'):
        """ To compare domains.
        """
        # Set up scorecard
        SC = ScoreCard()

        # Reproject domains to common grid

        # Evaluate detscores - save?

        # Evaluate probscores - save?

        # Plot to scorecard.

    def plot_all_detscores(self,scores='all',average='all',
                            *args,**kwargs):
        """ Plot all available deterministic scores.

        Currently available scores are:
            * CSI
            * POD...

        Args:
            scores (str,list): If `'all'` (default), this method will
                use all available scores hardcoded into the method.
            average (str): If 'all' (default), method will do mean across
                all cases and times with equal weighting. If None,
                scores will be plotted for every time.
        """
        args, kwargs = self.__get_plot_options(vrbl,*args,**kwargs)

        # if scores == 'all':
            # scores = ['CSI','POD',]

        scores = self.return_detscores(fcst=fcst,obs=obs,datafname=datafname)

        for score in scores:
            LG = LineGraph(self.outdir,fname=fname)
            LG.plot_score(xdata,ydata,hold=False)

    def return_detscores(self,fcst,obs,datafname=None):
        # Get data into right format.
        for utc in kwargs['list_of_times']:
            fcst = 0
            obs = 0
            DS = DetScores(fcst_arr=fcst,obs_arr=obs)

            scores = DS.compute_all(datadir=self.datadir,fname=datafname)
            return scores

    def plot_all_probscores(self,vrbl,scores='all',average='all',
                    *args,**kwargs):
        """ Plot all available probabilistic scores.

        Currently available scores are:
            * CRPS
            * ...

        Args:
            scores (str,list): If `'all'` (default), this method will
                use all available scores hardcoded into the method.
            average (str): If 'all' (default), method will do mean across
                all cases and times with equal weighting. If None,
                scores will be plotted for every time.
        """
        args, kwargs = self.__get_plot_options(vrbl,*args,**kwargs)

    ##### PRIVATE METHODS #####

    def __get_plot_options(self,vrbl,*args,**kwargs):
        """ Filter arguments and key-word arguments for plotting methods.

        Whatever is in dictionary will overwrite defaults in the plotting
        method.

        These may be
            * fhrs (forecast hour plotting times - or all)
            * ensmems (ensemble members, or all)


        """
        # Get plotting levels if not already given
        # TODO: rewrite this using hasattr() or something.
        S = Scales(vrbl)
        if not 'levels' in kwargs:
            kwargs['levels'] = S.clvs
        if not 'cmap' in kwargs:
            kwargs['cmap'] = S.cm

        # Specific things for certain variables
        if vrbl in ('REFL_10CM',"REFL_comp"):
            pass

        # Save all figures to a subdirectory
        if subdir in kwargs:
            utils.trycreate(subdir,is_folder=True)

        # What times are being plotted?
        # If all, return list of all times
        if 'utc' in kwargs:
            pass
        elif ('fchr' not in kwargs) or (kwargs['fchr'] == 'all'):
            kwargs['utc'] = E.list_of_times
        # Does this pick up on numpy arange?
        elif isinstance(kwargs['fchr'], (list,tuple,N.ndarray)):
            kwargs['utc'] = []
            for f in kwargs['fchr']:
                utc = self.inittime + datetime.timedelta(seconds=3600*f)
                kwargs['fchr'].append(utc)

        # Make domain smaller if requested

        # Save data before plotting
        clskwargs['save_data'] = kwargs.get('save_data',False)
        return clskwargs,plotkwargs,mplkwargs
            
    def generate_npy_fname(self,vrbl,score,lv,fchr,dom=None,ens=None,
                            fullpath=False,prefix=None,
                            *args,**kwargs):
        """ Save to disc with a naming scheme.

        A separate directory per init time (Ensemble instance) is
        needed to avoid overwriting.

        Note:
            ens is the ensemble member (name) or product (average, etc)


        Args:
            fullpath (bool): if True, give absolute path.
            args: a number of suffixes to add before the extension
                (in case A/B testing is needed on the same product)
        """
        def set_default(x,if_none,action=1):
            if x == None:
                return if_none
            else:
                if action == 1:
                    return str(x)
                elif action == 'dom':
                    return 'd{:02d}'.format(x)
                else:
                    return x

        # def lv_str(lv):
            # if lv == None:
                # return "sfc"
            # else:
                # return str(lv)

        vrblstr = vrbl
        scorestr = score
        lvstr = set_default(lv,'sfc',)
        # lvstr = lv_str(lv)
        fchrstr = '{}h'.format(fchr)
        domstr = set_default(dom,if_none='',action='dom')
        ensstr = set_default(ens,if_none='',)
        
        prefix = set_default(prefix,if_none='',)

        joinlist = [prefix,] + [vrblstr,scorestr,domstr,lvstr,fchrstr] + list(args) 
        fname = '_'.join(joinlist) + '.npy'

        if fullpath:
            return os.path.join(self.datadir,fname)
        else:
            return fname

    def lookup_validtime(self,fchr,intersect_only=True):
        """ Convert a forecast hour to a verification time.

        Args:
            intersect_only (bool): If True, raise error if 
                the time is not present for all domains.
                NOT IMPLEMENTED!
        """
        return self.E.initutc + datetime.timedelta(seconds=3600*fchr)


    def create_newgrid(self):
        """ Create neutral domain to reproject to.

        Done by creating generic Fields.

        Todo:
            * Refactor so that user selects domain to interp. to?
        """
        # DATA is the WRFOut or obs object.

        # Shouldn't all the obs/ensemble classes have a 
        # Field class attribute with the info of
        # x,y,lats,lons,bmapi,cen_lat,cen_lon,truelat1,
        # truelat2, Nlim/urlat, Elim etc

        # newgrid = Field(proj='lcc',...)

        # For now we'll use something else.
        # Make grid the smallest domain present.

        self.ng_nx = self.RD['nx']
        self.ng_ny = self.RD['ny']

        W = WRF_native_grid(self.arbdict[max(self.doms)])
        newgrid = VerifGrid(W,nx=self.RD['nx'],ny=self.RD['ny'])
        return newgrid


    def make_arbdict(self,):
        """

        Todo:
            * Make WRF_native_grid take objects or fpaths alike.
        """
        arbdict = {}
        for dom in self.doms:
            arbdict[dom] = self.E.arbitrary_pick(dom=dom,
                                    dataobj=False,give_path=True)
        return arbdict

    def check_ensemble_domains(self):
        # This should be a list of doms, not number.
        assert isinstance(self.E.doms,list)
        # doms should not be Python numbering~
        assert self.E.doms[0] == 1
        
        doms = self.E.doms
        ndoms = self.E.ndoms
        return doms,ndoms

    @staticmethod
    def ensure_list(a):
        """ Make sure the argument is a list

        Todo:
            * Move to utils.

        """
        if isinstance(a,list):
            return a
        elif isinstance(a,(tuple,N.ndarray)):
            return list(a)
        else:
            return (a,)


    def create_stats_iterable(self,vrbl,verif_times,doms,lvs):
        """ Generator for computing stats in parallel.

        Note:
            * All arguments must be lists/tuples.
            * vrbl is always the same. Limiting complexity...

        """
        for t,d,l in itertools.product(verif_times,
                                    doms,lvs):
            yield t,d,l,vrbl

    def start_your_engines(self,ncpus):
        """ Start multiprocessing pool.
        """
        self.pool = Pool(ncpus)
        return

    # @classmethod
    def lookup_statfunc(self,stat):
    # def lookup_statfunc(cls,stat):
        """ Dictionary lookup for stat-compute methods.

        Todo:
            * Why won't it work as classmethod?
        """
        lc_stat = stat.lower()
        return self._STATS[lc_stat]
        # return cls._STATS[stat]

    def generate_lookup_table(self):
        STATS = dict(
                contingency = self.compute_contingency,
                crps = self.compute_crps,
                rmse = self.compute_rmse,
                )
        return STATS

    def get_domlist(self,dom):
        """ Get list of domains.

        Could merge/use ensure_list?
        """
        if dom == 'all':
            domlist = self.doms
        elif isinstance(dom,int):
            domlist = (dom,)
        elif isinstance(dom,(list,tuple)):
            domlist = dom
        else:
            raise Exception

        return domlist

    def compute_contingency(self,dom='all'):
        """ Computed the 2x2 table over the whole ensemble
        """

        resultdict = self.analogue_dictionary(dom=dom)
        abcd = compute_contingency(fc,ob,th,'over')

    def analogue_dictionary(self,dom=1):
        """ Create a mirror-structure dictionary to
        the ensemble and domain requested.
        """


    def get_ob_instance(self,vrbl,return_name=False):
        """ Check original observation instances passed to init
        to see if it is available. Logic here could be easier using, e.g.
        `isinstance(k, StageIV)` but this involves importing a lot
        of unnecessary classes that should be instead invoked only
        by the user.

        Note:
            Must use lower case for all references. This might cause problems
                with saving plots/data though.
        """
        if vrbl in ('accum_precip','stageiv'):
            classname = 'stageiv'
        elif vrbl in ('REFL_comp','REFL_10CM','radar'):
            classname = 'radar'
        elif vrbl in ('stormreports','spclsr'):
            raise Exception("Not written yet.")
        else:
            raise Exception("Not implemented.")

        for k,v in self.obdict.items():
            # if isinstance(k, classname):
            # This needs to be something like
            if k == classname:
                if return_name:
                    return (v,classname)
                else:
                    return v
        return False
        

    def save_npy(self,data,vrbl,score,lv,fchr,dom,ens,fpath=False):
        # TODO
        # if not evac.number.int/float...
        # if not isinstance(data,N.ndarray):
            # data = N.array(data)
        if not fpath:
            fpath = self.generate_npy_fname(vrbl,score,lv,fchr,
                                dom,ens,fullpath=True)
        utils.trycreate(fpath)
        N.save(fpath,data)
        print("Saved data to {}".format(fpath))
        return


    def naming_for_reproj(self,vrbl,model,lv,utc,dom=None,
                            ens=None,fullpath=True):
        fpath = self.generate_npy_fname(vrbl=vrbl,score=model,lv=lv,fchr=utc,
                                        dom=dom,ens=ens,fullpath=fullpath,
                                        prefix='REPROJ')
        return fpath
    
    def check_for_reproj(self,vrbl,model,lv,utc,dom=None,ens=None):
        fpath = self.naming_for_reproj(vrbl=vrbl,model=model,lv=lv,utc=utc,
                                        dom=dom,ens=ens,fullpath=True,)
        if os.path.exists(fpath):
            print("Reprojection exists.")
            return True
        else:
            return False

    def save_reproj(self,data,vrbl,model,lv,utc,dom=None,ens=None)
        fpath = self.naming_for_reproj(vrbl=vrbl,model=model,lv=lv,utc=utc,
                                        dom=dom,ens=ens,fullpath=True,)
        utils.trycreate(fpath)
        N.save(fpath,data)
        print("Saved reprojection data to {}".format(fpath)
        return

    def load_reproj(self):
        arr = N.load(fpath)
        print("Loaded array from {}".format(fpath)
        return arr

    @staticmethod
    def _reproj_func(data):
        data_ng = reproject(data,xx_orig=data_ng_xx,
                yy_orig=data_ng_yy,xx_new=self.newgrid.xx,
                yy_new=self.newgrid.yy,)#method='linear')
        return data_ng

    def do_reprojection(self,data,lons,lats,save=None):
        """ Reproject data. If arguments are multiple, these
        are put onto common grid.

        Args:
            data: 2D array of data
            lons: longitudes of data
            lats: latitudes of data
        Returns:
            Data from new domain.
        """

        # New array for interpolated data.
        # x/y coordinates of reprojected data.
        data_ng_xx, data_ng_yy = self.newgrid.m(lons,lats)
        print("Reprojecting data.")
        if data.ndim == 3:
            data_ng = N.zeros((data.shape[0],self.ng_ny,self.ng_nx))
            for ens in range(data.shape[0]):
                data_ng[ens,:,:] = self._reproj_func(data[ens,:,:])
        elif data.ndim == 2:
            data_ng = self._reproj_func(data)
        else:
            raise Exception

        if save is not None:
            fpath = save
            self.save_reproj(fpath,data_ng)
        return data_ng

    def trim_radar_forever(self,bbdict):
        """ Make a copy of self.obdict[R], where R is the `Radar()`
        object, and return it with a smaller bounding box
        as specified.

        Todo:
            * Maybe we should put the copied, trimmed version back into
                the instance's obdict rather than returning it.
        """
        R2 = copy.copy(self.get_ob_instance(radar))
        R2.get_subdomain(overwrite=True,**bbdict)
        return R2

    @staticmethod
    def reduce_data_dims(data,lats=None,lons=None,level=None,utc=None):
        """ Get data from ensemble and deliver in 2-D chunks.

        Example:
            To return 
            reduce_data_dims(E) 
        """
        if data.ndim == 5:
            return data[:,0,0,:,:]
        elif data.ndim == 4:
            return data[0,0,:,:]
        else:
            raise Exception

    @staticmethod
    def ensure_integers(lst):
        # if isinstance(lst,N.ndarray):
            # return N.arra
        for n,x in enumerate(lst):
            lst[n] = int(x)
        return lst

    @staticmethod
    def ensure_floats(lst):
        for n,x in enumerate(lst):
            lst[n] = float(x)
        return lst
            
    def generate_obdict(self):
        obdict = {}
        for o in self.obs:
            n = self.get_object_name(o)
            obdict[n] = o
        return obdict

    @staticmethod
    def get_object_name(obj):
        """ This returns the name of the class instance, in lower case.

        There's surely a better way to do this...
        """
        a = str(obj.__class__).lower().split('.')[-1].split("'")[0]
        return a

""" 
DEPRECIATED - TOO BLOATED

Helper scripts to run stats over ensembles for a given
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

    User can run Verif in a mode that takes no ensemble forecast argument,
        but instead loads data saved to disc.

Todo:
    * Is this even essential? A lot of rewriting signatures to pass into functions
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
import time

import numpy as N

from evac.stats.objectbased import ObjectBased
from evac.plot.thumbnails import ThumbNails
from evac.plot.linegraph import LineGraph
from evac.plot.boxplot import BoxPlot
from evac.plot.violinplot import ViolinPlot
from evac.plot.histogram import Histogram
from evac.stats.detscores import DetScores
from evac.stats.probscores import ProbScores
from evac.plot.scorecard import ScoreCard
from evac.datafiles.radar import Radar
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
        obs: instance of ObsGroup, or an Obs subclass, or a list of either
        outdir (optional): path to output files
        datadir (optional): Path to which .npy data files will be saved
        reproject_dict (dict): if None, no reprojection. Otherwise,
            reprojection to a common domain is done with these options

    """


    def __init__(self,ensemble=None,obs=None,outdir=False,datadir=False,
                    newgrid=None):

        # FORECAST DATA
        self.E = ensemble
        if self.E:
            self.validtimes = self.E.validtimes
            self.initutc = self.E.initutc
            self.doms, self.ndoms = self.check_ensemble_domains()
            self.nmems = self.E.nmems
            self.arbdict = self.make_arbdict()
            self.lats, self.lons = self.get_latlons()
            
        # OBSERVATION DATA
        if isinstance(obs,(list,tuple)):
            self.obs = obs
        else:
            self.obs = (obs,)
        self.obdict = self.generate_obdict()

        # OPTIONAL SETTINGS FOR FIGURES AND SAVED DATA
        self.outdir = outdir
        self.datadir = datadir

        # OPTIONAL REPROJECTION DOMAIN
        self.newgrid = newgrid
    
    ##### INTERFACE METHODS - COMPUTE #####

    def compute_crps(self,xfs,xa,vrbl,utc=None,fchr=None,lv=None,dom=1):
        """ Compute CRPS.

        All arguments should be for a single time, variable, etc.
        Args:
            xfs: forecast data
            xa: observation data
            vrbl: variable to compute
            lv: level from forecast data. If None, use surface/lowest.
            utc (datetime.datetime): Datetime for plotting. If None,
                fchr must be specified.
            fchr: forecast hour. If None, utc must be specified.
            dom: domain. If None, use first.
        """
        fpath = self.generate_npy_fname(vrbl,'CRPS',lv=lv,fchr=fchr,
                            dom=dom,ens='mean',fullpath=True,)
        if not os.path.exists(fpath):
            P = ProbScores(xfs=xfs,xa=xa)
            crps = P.compute_crps(self.crps_thresholds)
            self.save_npy(crps,vrbl,'CRPS',lv=lv,fchr=fchr,dom=dom,ens='mean',fpath=fpath)
        else:
            print("Skipping - file already created.")
        return
   
    def compute_rmse(self,):
        pass

    def compute_stats(self,stats,vrbl,verif_times,dom=1,ncpus=1,
                        lvs=None,method=2,**kwargs):
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
            * Implement options, assign to self and then it can be seen
                by all compute_* methods.
            * Check if product has already been generated/saved. If not,
                also check that reprojection has been saved.

        Returns:
            A numpy save file is created for each chunk looped over.

        """
        self.start_your_engines(ncpus)

        # This is written silly for debugging
        verif_times1 = self.ensure_list(verif_times)
        verif_times2 = self.ensure_integers(verif_times1)
        verif_times = verif_times2

        # add kwargs to class
        self.allowed_compute_kwargs = ['crps_thresholds','det_thresholds']
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

        debug= 0
        for stat in stats:
            print("Generating {} stats now.".format(stat))
            statfunc = self.lookup_statfunc(stat)
            if method == 3:
                # Serial mode, good for debugging
                for chunk in itr:
                    statfunc(chunk)
            elif method == 1:
                for chunk in itr:
                    fchr,dom,lv,vrbl = chunk
    
                    # Check to see if stat has been computed/saved

                    # get_fcst_data
                    fcst = self.E.get(vrbl,fcsthr=fchr,dom=dom,level=lv,accum_hr=1)
                    # pdb.set_trace()
                    fdata = self.reduce_data_dims(fcst)

                    # get_ob_data
                    utc = self.lookup_validtime(fchr)
                    obobj, obname = self.get_ob_instance(vrbl,return_name=True)
                    obdata = self.reduce_data_dims(obobj.get(utc,lv=lv))

                    # One parallel loop (async it)
                    # try to load reproject
                    exists,rpj_fpath_f = self.check_for_reproj(vrbl=vrbl,model='fcst',lv=lv,fchr=fchr,
                                            dom=dom,ens='all',return_fpath=True,)
                    if not exists:
                        # reproject if not, PARALLEL, and save
                        flats, flons = self.E.get_latlons(dom=dom)
                        xfs = self.do_reprojection(fdata,flons,flats,save=rpj_fpath_f)
                    else:
                        if debug:
                            flats, flons = self.E.get_latlons(dom=dom)
                            xfs = self.do_reprojection(fdata,flons,flats,save=rpj_fpath_f)
                            pdb.set_trace()

                        xfs = N.load(rpj_fpath_f)

                    exists, rpj_fpath_o = self.check_for_reproj(vrbl=vrbl,model=obname,lv=lv,utc=utc,
                                            dom=None,ens=None,return_fpath=True,)
                    if not exists:
                        oblats = obobj.lats
                        oblons = obobj.lons
                        xa = self.do_reprojection(obdata,oblons,oblats,save=rpj_fpath_o)
                    else:
                        xa = N.load(rpj_fpath_o)

                    # Next parallel loop
                    # Run statfunc and exit

                    kw = dict(xfs=xfs,xa=xa,vrbl=vrbl,lv=lv,fchr=fchr,dom=dom,)
                    # self.pool.apply_async(statfunc,args=chunk,kwargs=kw)
                    statfunc(**kw)
            elif method == 2:
                codedict = {}
                for chunk in itr:
                    fchr,dom,lv,vrbl = chunk
                    code = self._get_loop_code(chunk)
                    utc = self.lookup_validtime(fchr)
                    obobj,obname = self.get_ob_instance(vrbl,return_name=True)

                    exists,rpj_fpath_f = self.check_for_reproj(vrbl=vrbl,model='fcst',lv=lv,fchr=fchr,
                                            dom=dom,ens='all',return_fpath=True)
                    if not exists:
                        codedict[code] = [rpj_fpath_f,]
                        

                    exists,rpj_fpath_o = self.check_for_reproj(vrbl=vrbl,model=obname,lv=lv,utc=utc,
                                            dom=None,ens=None,return_fpath=True)
                    if not exists:
                        codedict[code].append(rpj_fpath_o)

                # We now have a list of all chunks that need reprojection.
                # For fcst:
                for chunk in itr:
                    fchr,dom,lv,vrbl = chunk
                    code = self._get_loop_code(chunk)

                    if codedict[code][0]:
                        flats, flons = self.E.get_latlons(dom=dom)
                        fcst = self.E.get(vrbl,fcsthr=fchr,dom=dom,level=lv,accum_hr=1)
                        fdata = self.reduce_data_dims(fcst)
                        # self.do_reprojection(fdata,flons,flats,save=rpj_fpath_f)
                        ar = (fdata,flons,flats)
                        kar = dict(save=rpj_fpath_f)
                        result = self.pool.apply_async(self.do_reprojection,args=ar,kwargs=kar)
                    self.pool.join()


                # For obs:
                for chunk in itr:
                    fchr,dom,lv,vrbl = chunk
                    code = self._get_loop_code(chunk)

                    if codedict[code][1]:
                        utc = self.lookup_validtime(fchr)
                        oblats = obobj.lats
                        oblons = obobj.lons
                        obobj,obname = self.get_ob_instance(vrbl,return_name=True)
                        obdata = self.reduce_data_dims(obobj.get(utc,lv=lv))
                        # xa = self.do_reprojection(obdata,oblons,oblats,save=rpj_fpath_o)
                        ar = (obdata,oblons,oblats)
                        kar = dict(save=rpj_fpath_o)
                        result = self.pool.apply_async(self.do_reprojection,args=ar,kwargs=kar)
                    self.pool.join()

                    # self.pool.close()
                
                for chunk in itr:
                    xfs = N.load(rpj_fpath_f)
                    xa = N.load(rpj_fpath_o)
                    statfunc = self.compute_crps_mp
                    fchr,dom,lv,vrbl = chunk
                    kw = dict(xfs=xfs,xa=xa,vrbl=vrbl,lv=lv,fchr=fchr,dom=dom)
                    self.pool.apply_async(statfunc,args=chunk,kwargs=kw)
                self.pool.join()
                    

        return

    @staticmethod
    def _get_loop_code(itr):
        """ Merge iter items as one string for looking up later.
        """
        itr = list(itr)
        for n,i in enumerate(itr):
            if i is None:
                itr[n] = 'None'
            itr[n] = str(itr[n])
        s = '_'.join(itr)
        return s


    def compute_objectbased(self,vrbl1,vrbl2,fchrs='all',doms='all',lvs=None,
                                thresh='auto',footprint=500,suite=None,
                                obs_only=False,fcst_only=False,quickplot=False,
                                dx=None):
        """ Compute stats related to object-based scores.

        The procedure is:
            * Load arrays for all domains
            * Create ObjectBased instances
            * Pass these into SAL etc
            * Run other stats

        Args:
            suite (str): "SAL", etc
            vrbl (str): variable to identify objects 
            vrbl2 (str): variable for diagnostics
            fchrs: forecast hours to compute.
            footprint: in pixels if dx is None. If dx is used, 
                this is in km instead.
        """
        if doms == 'all':
            doms = self.doms
        ndoms = len(doms)

        if dx is None:
            fp = [footprint for d in doms]
        else:
            assert isinstance(dx,(tuple,list,N.ndarray))
            # dx = N.array(dx)
            fp_multi = [dx[0] * (dx[0]/dx[n])**2 for n in range(ndoms)]
            # fp = [footprint*((dx[0]/dx[n])**2)*dx[n] for n in range(ndoms)]
            fp = footprint * N.array(fp_multi)

        # if fchrs == 'all':
            # verif_times = 'all'
        # else:
            # verif_times = 
            # This needs computing - convert fchrs to datetimes.

        # pdb.set_trace()
        limdict = self.E.get_limits(dom=2)
        lats, lons = self.E.get_latlons(dom=1)

        def get_objects(all_objects,limdict,lats,lons):
            vrbl1_all = self.E.get(vrbl2,fcsthr=fchr,dom=dom,level=lv,members=member)[0,0,0,:,:]
            vrbl2_all = self.E.get(vrbl2,fcsthr=fchr,dom=dom,level=lv,members=member)[0,0,:,:,:]

            if dom == 1:
                vrbl1_arr, la_,lo_ = utils.return_subdomain(data=vrbl1_all,
                                        lats=lats,lons=lons,**limdict)
                vrbl2_arr, la_,lo_ = utils.return_subdomain(data=vrbl2_all,
                                        lats=lats,lons=lons,**limdict,)#axes=(1,2))
            else:
                vrbl1_arr = vrbl1_all
                vrbl2_arr = vrbl2_all
                la_ = None
                lo_ = None

            Obj = ObjectBased(vrbl1_arr,thresh=thresh,footprint=fp[ndom])
            uarr,_udict = Obj.get_updraughts(arr=vrbl2_arr)
            # pdb.set_trace()
            all_objects = N.hstack((all_objects,uarr))
            print("We found {} objects.".format(len(uarr)))
            return all_objects, Obj, la_, lo_

        verif_times = fchrs
        for ndom,dom in enumerate(self.doms):
            itr = self.create_stats_iterable(vrbl1,verif_times,doms=(dom,),lvs=lvs,first_time=False)
            all_objects = N.empty([0])
            for fchr,dom,lv,vrbl1 in itr:
                for member in self.E.member_names:
                    all_objects, Obj, la_, lo_ = get_objects(all_objects,limdict,lats,lons)

                    if quickplot and (member == self.E.member_names[0]):
                        fcmin = int(60*fchr)
                        W = self.E.return_datafile(member=member,dom=dom,utc=self.initutc)
                        fname = 'test_objs_d{:02d}_{:03d}min_{}.png'.format(dom,fcmin,member)
                        fpath = os.path.join(self.outdir,fname)
                        if dom == 1:
                            ld = limdict
                            W = None
                        else:
                            ld = dict()
                        lats = la_
                        lons = lo_
                        try:
                            Obj.plot(W=W,fpath=fpath,ld=ld,lats=lats,lons=lons)
                        except ValueError:
                            print("Skipping this time.")
                            continue
                        except IndexError:
                            print("No objects for this time.")
                            continue
                        # pdb.set_trace()
            fpath = self.generate_npy_fname(vrbl1,'object_{}_array'.format(vrbl2),
                                lv=None,fchr='all',
                                dom=dom,ens='all',fullpath=True,)
            print("Saving array to {}".format(fpath))
            utils.trycreate(fpath)
            N.save(arr=all_objects,file=fpath)

        return

    def plot_object_pdfs(self,vrbl1,vrbl2,fname=False,doms=(1,2)):
        """ Plot histograms of object diagnostics.
        """
        if not fname:
            fname = 'test_hist.png'

        # Load data
        arrs = []
        labels = list()
        for dom in doms:
            fpath = self.generate_npy_fname(vrbl1,'object_{}_array'.format(vrbl2),
                            lv=None,fchr='all',
                            dom=dom,ens='all',fullpath=True,)
            arr = self.load_data(fpath)
            arrs.append(arr)
            labels.append("Domain {}".format(dom))


        # pdb.set_trace()
        # Plot two domains cdf
        H = Histogram(fpath=os.path.join(self.outdir,fname),
                        data=arrs)
        H.plot(labels=labels)

        return

    def compute_pdf_stats(self,vrbl1,vrbl2,dom=1):
        """ Compute stats on pdf distributions of objects.
        """

        # Load data
        fpath = self.generate_npy_fname(vrbl1,'object_{}_array'.format(vrbl2),
                        lv=None,fchr='all',
                        dom=dom,ens='all',fullpath=True,)
        arr = self.load_data(fpath)

        # Run CRPS over pdf/cdfs.

        pass

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
        
        # If output exists, don't repeat

        fpath = self.generate_npy_fname(vrbl,'CRPS',lv,fchr,
                                dom,'mean',fullpath=True)
        if os.path.exists(fpath):
            print("Skipping - file already created.")
            return

        # Load verification
        utc = self.lookup_validtime(fchr)
        # Make lv ignored if it's nonsense
        # Need to implement get() in all obs type.
        # Do this in parent class and override.
        obobj,obname = self.get_ob_instance(vrbl,return_name=True)

        # Reproject and compute
        exists,rpj_fpath_f = self.check_for_reproj(vrbl=vrbl,model='fcst',lv=lv,fchr=fchr,
                                dom=dom,ens='all',return_fpath=True)
        if not exists:
            fcst = self.E.get(vrbl,fcsthr=fchr,dom=dom,level=lv,accum_hr=1)
            fdata = self.reduce_data_dims(fcst)
            flats, flons = self.E.get_latlons(dom=dom)
            xfs = self.do_reprojection(fdata,flons,flats,save=rpj_fpath_f)
        else:
            xfs = N.load(rpj_fpath_f)

        
        exists,rpj_fpath_o = self.check_for_reproj(vrbl=vrbl,model=obname,lv=lv,utc=utc,
                                dom=None,ens=None,return_fpath=True)
        if not exists:
            obdata = self.reduce_data_dims(obobj.get(utc,lv=lv))
            oblats = obobj.lats
            oblons = obobj.lons
            xa = self.do_reprojection(obdata,oblons,oblats,save=rpj_fpath_o)
        else:
            xa = N.load(rpj_fpath_o)

        P = ProbScores(xfs=xfs,xa=xa)
        # fr thresh in self.crps_thresholds:
        crps = P.compute_crps(self.crps_thresholds)

        # Save to disk
        self.save_npy(crps,vrbl,'CRPS',lv,fchr,dom,'mean',fpath=fpath)
        return

    def compute_detscores(self,xfs,xa,vrbl,lv,fchr,dom,debug=1):
        for thresh in self.det_thresholds:
            thstr = '{:02d}mmh'.format(thresh)
            for n,ens in enumerate(self.E.member_names):
                fpath = self.generate_npy_fname(vrbl,'contingency',lv=lv,fchr=fchr,
                                dom=dom,ens=ens,fullpath=True,suffix=thstr,npy_extension=False)
                # pdb.set_trace()
                if not os.path.exists(fpath+'.npz'):
                    DS = DetScores(fcst_arr=xfs[n,:,:],obs_arr=xa,
                                thresh=thresh,overunder='over')
                    scores = DS.compute_all(datadir=self.datadir,fname=fpath)
                else:
                    if debug:
                        DS = DetScores(fcst_arr=xfs[n,:,:],obs_arr=xa,
                                    thresh=thresh,overunder='over')
                        scores = DS.compute_all(datadir=self.datadir,fname=fpath)
                        print("initutc {}, dom {}, thresh {} = {}".format(self.initutc,
                                                    dom,thresh,scores['BIAS']))
                        # time.sleep(1)
                        # pdb.set_trace()
                    print("Skipping; scores exist for",fpath)
        return

    ##### INTERFACE METHODS - PLOT #####

    def plot_domains(self,ensemble_domains='all',ob_domains='all'):
        """ Wrapper to plot domains of forecast or observation domains.
        """
        pass
    
    def plot_trails(self,init_dict,score,vrbl,doms,interval=1,
                    use_hours=True,lv='sfc',ens='all',loop_prefix=None,
                    loop_suffix=None,mplargs=None,mplkwargs=None,
                    plotargs=None,plotkwargs=None,outdir='~/',fname='test.png',
                    by_leadtime=False):
        """ Plot stats for every x min for all domains.

        Args:
            init_dict (dict): dict(datetime.datetime:rootdir)
            score (str): 'crps', etc
            vrbl (str): 'accum_precip', etc
            doms: list, tuple, N.array of domain names/codes
            interval (int): interval in minutes (or hours) to plot
            use_hours (bool): if True, use hours for interval (for
                when loading data, when plotting, and setting interval)
            lv: By default, this is 'sfc'
            ens: By default, this is 'all'
            loop_suffix, loop_prefix (bool): when set (list of strings),
                script can loop over e.g. thresholds, etc
            by_leadtime (bool): if True, plot all trails by leadtime
                on the same axis.
            """ 
        # TODO: make this just class variable doms
        COLORS = (('red','tomato'),
                    ('blue','slateblue'),
                    ('darkgreen','forestgreen'),
                    ('darkorchid','mediumorchid'),)
        DATA = {}
        LG = LineGraph(outdir,fname=fname,vrbl=vrbl,
                        figsize=(6,6))
        # initutcs = init_dict.keys()
        for ninit,(initutc,rootdir) in enumerate(init_dict.items()):
            
            # (dom, score_at_valid_time)
            # number of valid times:
            nvts = 3 # needs to be dynamic
            validtimes = N.arange(1,nvts+1)
            actualtimes = [initutc + datetime.timedelta(seconds=3600*int(v)) for v in validtimes]
            #DATA[initutc] = N.zeros(ndoms,nvts,)
            for ndom,dom in enumerate(doms):
                data = N.zeros(nvts)
                for nt,t in enumerate(validtimes):
                    # Loops here for suffix, prefix
                    # TODO: npz extension and score extraction...
                    fname = self.generate_npy_fname(self,vrbl=vrbl,score=score,lv=lv,fchr=t,
                                dom=dom,ens=ens,fullpath=False,npy_extension=True)
                                # dom=dom,ens=ens,fullpath=False,suffix=thstr,npy_extension=False)
                    fpath = os.path.join(rootdir,fname)
                    data[nt] = self.load_data(fpath)
                plotkwargs['color'] = COLORS[ninit][ndom]
                # LG.plot_score(xdata=actualtimes,ydata=data,hold=True,mplkwargs=mplkwargs)
                LG.plot_score(use_plot_date=True,xdata=actualtimes,
                        ydata=data,hold=True,plotkwargs=plotkwargs)
                # LG.plot_score(xdata=validtimes,ydata=data,hold=True,mplkwargs=mplkwargs)
                # pdb.set_trace()


        LG.save()
        pdb.set_trace()

    def plot_stats(self,plot,
                    # init_dict,score,vrbl,
                    ndoms=1,outdir='./',
                    # interval=1,
                    # use_hours=True,lv='sfc',ens='all',loop_prefix=None,
                    # loop_suffix=None,
                    # outdir='~/',fname='test.png',
                    mplargs=None,mplkwargs=None,
                    plotargs=None,plotkwargs=None,
                    *args,**kwargs):
        """ Master script for plotting anything.
        """
        if not mplargs:
            mplargs = {}
        if not mplkwargs:
            mplkwargs = {}
        if not plotargs:
            plotargs = {}
        if not plotkwargs:
            plotkwargs = {}
        doms = N.arange(1,ndoms+1)
        plotfunc = self.lookup_plotfunc(plot)
        plotfunc(doms=doms,outdir=outdir,
                    mplargs=mplargs,mplkwargs=mplkwargs,
                    plotargs=plotargs,plotkwargs=plotkwargs,
                    *args,**kwargs)

    def lookup_plotfunc(self,plot):
        STATFUNCS = {}
        STATFUNCS['trails'] = self.plot_trails
        STATFUNCS['violin'] = self.plot_violin
        STATFUNCS['boxplot'] = self.plot_boxplot

        plot = plot.lower()
        return STATFUNCS[plot]

    def plot_boxplot(self,outdir):
        """ TODO: implement this.
        """
        raise Exception

    def plot_violin(self,outdir,
                    init_dict,score,vrbl,doms,
                    interval=1, use_hours=True,lv='sfc',ens='all',
                    loop_prefix=None,loop_suffix=None,
                    set_prefix=None,set_suffix=None,
                    nens=18,
                    *args,**kwargs):
        """ Plot violin plot

        Args:
            vrbl (str): name of variable to evaluate

        This is an alternative to box-and-whisker.
        """
        figsize = kwargs.get('figsize',(6,6))
        thresh = kwargs.get('thresh',5)

        nvts = 3 # needs to be dynamic
        validtimes = N.arange(1,nvts+1)
        # actualtimes = [initutc + datetime.timedelta(seconds=3600*int(v)) for v in validtimes]
        # Data is (inittime,validtime,dom,ensmember)
        data = {}
        for ninit,(initutc,rootdir) in enumerate(init_dict.items()):
            data[initutc] = {}
            for nt,t in enumerate(validtimes):
                data[initutc][t] = {}
                for nd,dom in enumerate(doms):
                    data[initutc][t][dom] = []
                    for ne,ens in enumerate(N.arange(1,nens+1)):
                        fname = self.generate_npy_fname(vrbl=vrbl,score='contingency',
                                            lv=lv,fchr=t,dom=dom,ens='m{:02d}'.format(ens),
                                            fullpath=False,npy_extension='npz',
                                            suffix='{:02d}mmh'.format(thresh))
                        fpath = os.path.join(rootdir,fname)
                        data[initutc][t][dom].append(self.load_data(fpath)[score])
                    data[initutc][t][dom] = N.array(data[initutc][t][dom])

        # Group ensemble members together?
        fname1 = 'all_violin_{}.png'.format(score)
        VP = ViolinPlot(outdir,fname=fname1,figsize=figsize)
        fname2 = 'all_boxplot_{}.png'.format(score)
        BX = BoxPlot(outdir,fname=fname2,figsize=figsize)
        # list of 1D arrays for each member
        plotlist = []
        pos_it = N.array([1,2,4,5,7,8])
        pos_lb = ['d01','d02'] * 12
        posns = pos_it.copy()
        for n,initutc in enumerate(init_dict):
            if n > 0:
                posns = N.hstack((posns,pos_it+(n*12)))
            for t in validtimes:
                for dom in doms:
                    plotlist.append(data[initutc][t][dom])
        # pdb.set_trace()
        pos = posns.flatten()
        VP.plot(dataset=plotlist,positions=pos,vert=True,)
        BX.plot(plotlist,positions=pos,vert=True,labels=pos_lb,
                autorange=True)

        # if score == 'CSI':
            # pdb.set_trace()

    def plot_thumbnails(self,vrbl,radardir=None,ensmembers='all',fchrs='all',
                        doms='all',limdict='auto',outdir=None,fname='auto',
                        mplargs=None,mplkwargs=None,verif_first=True,
                        lv=None,overwrite=True,
                        *args,**kwargs):
        """ Plot postage stamps or thumbnails of ensemble members listed
        (by default, this is all).

        Todo:
            * radardir should be more general to loading verification.
            * should be able to plot reprojections resulting from 
                score calculations.
            * Remove `radardir` to make this more general. The obs
                should be passed into the init script. How will this
                be done when Radar is only for one time?! We need a
                "obs ensemble" class.

        Args:
            verif_first (bool): If True, the first axis (thumbnail/subplot)
                will be the observations. Method exits if verification does
                not exist in self.obdict.
        """
        # self.__get_plot_options(*args,**kwargs)
        if outdir == None:
            outdir = self.outdir

        if fchrs == 'all':
            raise Exception("Need to implement this")
            # Use N.arange(1,4)

        if doms == 'all':
            doms = self.doms

    
        elif limdict is None:
            limdict = {}


        for fchr in fchrs:
            plotutc = self.lookup_validtime(fchr)
            for dom in doms:
                if limdict == 'auto':
                    ld = self.E.get_limits(dom=dom)
                else:
                    ld = limdict
                if fname == 'auto':
                    outfname = self.make_plot_fname('thumbnails',vrbl=vrbl,dom=dom,fchr=fchr,
                                                    lv=lv)
                if (overwrite is False) and (os.path.exists(outfname)):
                    print("Already created. Skipping.")
                    return
                TN = ThumbNails(outdir=outdir,fname=outfname)
                if verif_first:
                    OB = self.get_ob_instance(vrbl)
                    # R = Radar(plotutc,radardir)
                    # R.get_subdomain(**ld,overwrite=True)
                    # OB.get_subdomain(**ld,overwrite=True)
                    OB.set_subdomain(**ld)
                    TN.plot_verif(data=OB,utc=plotutc)

                arbW = self.get_arb(dom=dom,fpath_only=False)
                cen_lat,cen_lon = self.E.get_cenlatlons(dom=dom)
                fcstdata = self.reduce_data_dims(self.E.get(
                            vrbl,fcsthr=fchr,dom=dom,level=lv,accum_hr=1))
                TN.plot_fcst(data=fcstdata,vrbl=vrbl,W=arbW,
                        cen_lat=cen_lat,cen_lon=cen_lon,
                        titles=self.E.member_names)
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
            

    @staticmethod
    def _set_default(x,if_none,action=1):
        """ Default naming for fname generation methods.
        """
        if x == None:
            return if_none
        else:
            if action == 1:
                return str(x)
            elif action == 'dom':
                if isinstance(x,str):
                    return x
                else:
                    return 'd{:02d}'.format(x)
            else:
                return x

    def make_plot_fname(self,plottype,
                            # vrbl=None,dom=None,fcst=None,lv=None,
                            # fullpath=False,
                            *args,**kwargs):
        """ Wrapper for generate_fname.
        """
        f = self.generate_fname(score=plottype,*args,**kwargs)
        return f

    def generate_npy_fname(self,
                            # vrbl,score,lv,fchr=None,utc=None,dom=None,ens=None,
                            # fullpath=False,prefix=None,suffix=None,
                            *args,**kwargs):
        """ Wrapper for generate fname.
        """
        kwargs['npy_extension'] = kwargs.get('npy_extension',True)
        kwargs['datadir'] = self.datadir
        f = self.generate_fname(*args,**kwargs)
        return f

    @classmethod
    def generate_fname(cls,vrbl=None,score=None,lv=None,fchr=None,utc=None,
                            dom=None,ens=None,fullpath=False,prefix=None,
                            suffix=None,npy_extension=False,datadir=None):
        """ Save to disc with a naming scheme.

        A separate directory per init time (Ensemble instance) is
        needed to avoid overwriting.

        Note:
            ens is the ensemble member (name) or product (average, etc)


        Args:
            fullpath (bool): if True, give absolute path.
            prefix/suffix: a string to add at the start, or before the extension
            (in case A/B testing is needed on the same product)
        """

        # def lv_str(lv):
            # if lv == None:
                # return "sfc"
            # else:
                # return str(lv)

        vrblstr = vrbl
        scorestr = score
        lvstr = cls._set_default(lv,'sfc',)
        # lvstr = lv_str(lv)
        if utc:
            timestr = utils.string_from_time('output',utc)
        else:
            timestr = '{}h'.format(fchr)
        domstr = cls._set_default(dom,if_none='',action='dom')
        ensstr = cls._set_default(ens,if_none='',)
        
        prefixstr = cls._set_default(prefix,if_none='',)
        suffixstr = cls._set_default(suffix,if_none='',)

        joinlist = [vrblstr,scorestr,domstr,lvstr,timestr,ensstr] 
        if suffix:
            joinlist.append(suffix)
        if prefix:
            joinlist.insert(0,prefix)

        fname = '_'.join(joinlist) 

        if npy_extension:
            if npy_extension == True:
                ext = 'npy'
            elif isinstance(npy_extension,str):
                ext = npy_extension

            if not ext.startswith('.'):
                ext = '.' + ext
            fname = fname + ext

        if fullpath:
            return os.path.join(datadir,fname)
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


    def get_arb(self,dom,fpath_only=False):
        """ 
        Todo: 
            * Remove the arbdict and just use this to look
                up arbitrary WRFOuts or their paths.
        """
        if fpath_only:
            dataobj = False
            give_path = True
        else:
            dataobj = True
            give_path = False
        arb = self.E.arbitrary_pick(dom=dom,dataobj=dataobj,give_path=give_path)
        return arb

    def make_arbdict(self,):
        """

        Todo:
            * Make WRF_native_grid take objects or fpaths alike.
        """
        if self.E is None:
            return None
        arbdict = {}
        for dom in self.doms:
            arbdict[dom] = self.E.arbitrary_pick(dom=dom,
                                    dataobj=False,give_path=True)
        return arbdict

    def check_ensemble_domains(self):
        if self.E is None:
            return None, None
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


    def create_stats_iterable(self,vrbl,verif_times,doms,lvs='sfc',
                                first_time=True):
        """ Generator for computing stats in parallel.

        Note:
            * All arguments must be lists/tuples.
            * vrbl is always the same. Limiting complexity...

        """
        if doms == 'all':
            doms = self.doms
        if verif_times == 'all':
            verif_times = self.validtimes
            if not first_time:
                del verif_times[0]
        if not isinstance(lvs,(N.ndarray,list,tuple)):
            lvs = (lvs,)
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
                # crps = self.compute_crps,
                crps = self.compute_crps_mp,
                rmse = self.compute_rmse,
                detscores = self.compute_detscores,
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

    def naming_for_reproj(self,vrbl,model,lv,utc=None,fchr=None,dom=None,
                            ens=None,fullpath=True,):
        # utcstr = utils.string_from_time('output',utc,strlen='minute')
        fpath = self.generate_npy_fname(vrbl=vrbl,score=model,lv=lv,utc=utc,
                                        dom=dom,ens=ens,fullpath=fullpath,fchr=fchr,
                                        prefix='REPROJ')
        return fpath
    
    def check_for_reproj(self,vrbl,model,lv,utc=None,fchr=None,
                        dom=None,ens=None,return_fpath=False,
                        datadir=None):
        """ Check saved data exists.
        """
        fpath = self.naming_for_reproj(vrbl=vrbl,model=model,lv=lv,utc=utc,fchr=fchr,
                                        dom=dom,ens=ens,fullpath=True,)
        if os.path.exists(fpath):
            print("Reprojection exists.")
            ans = True
        else:
            ans = False

        if return_fpath:
            return ans,fpath
        else:
            return ans

    @classmethod
    def save_reproj(cls,data,vrbl=None,model=None,lv=None,utc=None,
                        dom=None,ens=None,fpath=None):
        if not fpath:
            fpath = cls.naming_for_reproj(vrbl=vrbl,model=model,lv=lv,utc=utc,
                                        dom=dom,ens=ens,fullpath=True,)
        utils.trycreate(fpath)
        N.save(fpath,data)
        print("Saved reprojection data to {}".format(fpath))
        return

    @staticmethod
    def load_reproj(fpath):
        arr = N.load(fpath)
        print("Loaded array from {}".format(fpath))
        return arr

    @staticmethod
    def load_data(fpath):
        """ Variation on load_reproj to load data.
        """
        arr = N.load(fpath)
        print("Loaded array from {}".format(fpath))
        return arr

    @staticmethod
    def _reproj_func(data,xx1,yy1,newgrid):
        # data_ng = reproject(data,xx_orig=data_ng_xx,
                # yy_orig=data_ng_yy,xx_new=self.newgrid.xx,
                # yy_new=self.newgrid.yy,)#method='linear')
        data_ng = reproject(data,xx_orig=xx1,yy_orig=yy1,
                    xx_new=newgrid.xx,yy_new=newgrid.yy)
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

        print("Reprojecting data.")

        if data.ndim == 3:
            data_ng = N.zeros((data.shape[0],self.newgrid.ny,self.newgrid.nx))
            for ens in range(data.shape[0]):
                data_ng[ens,:,:] = self._reproj_func(data[ens,:,:],
                                        data_ng_xx,data_ng_yy,
                                        self.newgrid)
        elif data.ndim == 2:
            data_ng = self._reproj_func(data,data_ng_xx,data_ng_yy,
                                        self.newgrid)
        else:
            raise Exception

        if save is not None:
            fpath = save
            self.save_reproj(data_ng,fpath=fpath)
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
            # In the event that the instance is ObsGroup
            if n == 'obsgroup':
                n = o.obs_type
            obdict[n] = o
        return obdict

    @staticmethod
    def get_object_name(obj):
        """ This returns the name of the class instance, in lower case.

        There's surely a better way to do this...
        """
        a = str(obj.__class__).lower().split('.')[-1].split("'")[0]
        return a

    def get_latlons():
        for dom in self.doms:
            lats[dom], lons[dom] = self.E.get_latlons(dom=dom)
        return lats, lons

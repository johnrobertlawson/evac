""" Helper scripts to run stats over ensembles for a given
initialisation time (looping over levels, ensemble members, valid times, and 
lat/lon if required)

The idea is that no actual plotting is done herein, but that this class
calls the relevant class methods outside of `Verif` to do this.

Todo:
    * Create a `evac.plot.thumbnails` class
    * Ensure all figure subtypes has a "print" call to notify of path to figure
    * Save numpy files of computed data for CPU-intense calculations.
    * Decorate all plotting/computation processes with timeme().
    * Be aware if saved data already exists - load if that's the case.
    * Wrap __get_plot_options into figure().

"""
import os
import pdb
import copy

from evac.plot.thumbnails import Thumbnails
from evac.plot.linegraph import LineGraph
from evac.stats.detscores import DetScores
from evac.stats.probscores import ProbScores
from evac.plot.scorecard import ScoreCard
import evac.utils as utils

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

    Args:
        ensemble: Class instance of `evac.datatypes.ensemble.Ensemble`.
        outdir: path to output files
        datadir (optional): Path to which .npy data files will be saved
        ST4 (optional): Class instance of `evac.datatypes.stageiv.StageIV`.

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

    """
    def __init__(self,ensembles,obs,outdir,datadir=False,):
        # A tuple of obs objects
        self.obs = obs
        # A tuple of forecast ensembles (for same times)
        self.E = ensembles    

        self.outdir = outdir
        self.obdict = *args
        self.datadir = datadir

        # Might attach self.E attributes such as times and member 
        # names to the class instance?

        # This should be a list of doms, not number.
        assert isinstance(self.E.doms,list)
        # doms should not be Python numbering~
        assert doms[0] == 1
        self.doms = self.E.doms
        self.ndoms = len(E.doms)

        # Get some useful things common to all plots

        # Dictionary of aribtrary WRFOut objects 
        self.arbW = {E.arbitrary_pick(dom=dom,dataobj=True) for dom in self.doms}

        # lats and lons for each domain
        self.lats = {}
        self.lons = {}
        for dom in self.doms:
            self.lats[dom],self.lons[dom] = self.E.get_latlons(dom=dom)

    def plot_domains(self,ensemble_domains='all',ob_domains='all'):
        """ Wrapper to plot domains of forecast or observation domains.
        """
        pass
        

    def get_ob_instance(self,vrbl):
        """ Check original observation instances passed to init
        to see if it is available. Logic here could be easier using, e.g.
        `isinstance(k, StageIV)` but this involves importing a lot
        of unnecessary classes that should be instead invoked only
        by the user.
        """
        if vrbl == ('accum_precip','stageiv'):
            classname = 'StageIV'
        elif vrbl in ('REFL_comp','REFL_10CM','radar'):
            classname = 'Radar':
        elif vrbl in ('stormreports','spclsr'):
            raise Exception("Not written yet.")

        for k in self.obdict.items():
            # if isinstance(k, classname):
            # This needs to be something like
            if k.__name__ == classname:
                return True
        return False

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
                        mplargs=None,mplkwargs=None,verif_first=True
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

    def generate_mutual_grid(self,):
        """ Create a new domain that all forecast and observation
        fields can be interpolated to.
        """

    def do_reprojection(self,):
        """ Reproject data. If arguments are multiple, these
        are put onto common grid.
        """
        # self.lats[dom] 
        pass

    def trim_radar_forever(self,bbdict):
        """ Make a copy of self.obdict[R], where R is the `Radar()`
        object, and return it with a smaller bounding box
        as specified.

        Todo:
            * Maybe we should put the copied, trimmed version back into
                the instance's obdict rather than returning it.
        """
        R2 = copy.copy(self.get_ob_instance(radar)
        R2.get_subdomain(**bbdict,overwrite=True)
        return R2

    @staticmethod
    def get_2D_data(lats=None,lons=None,level=None,utc=None):
        """ Get data from ensemble and deliver in 2-D chunks.

        Example:
            To return 
            get_2D_data(E) 
        """
        pass

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
        elif isinstance(kwargs['fchr'], (list,tuple)):
            kwargs['utc'] = []
            for f in kwargs['fchr']:
                utc = self.inittime + datetime.timedelta(seconds=3600*f)
                kwargs['fchr'].append(utc)

        # Make domain smaller if requested

        # Save data before plotting
        clskwargs['save_data'] = kwargs.get('save_data',False)
        return clskwargs,plotkwargs,mplkwargs
            

""" A catalogue of observation data files.

Only works on observation files that are subclasses of Obs, to ensure the correct
get() methods work.
"""

import pdb
import os
import glob

# import numpy as N

# import evac.utils as utils
from evac.datafiles.stageiv import StageIV
from evac.datafiles.radar import Radar
from evac.datafiles.obs import Obs

class ObsGroup:
    """ A group of observation data files.

    The instances or file types should all be the same.

    For instance, create a group of Stage IV data objects.

    Note:
        * obs_type strings are ['stageiv','radar',)...

    Args:
        fpath: Absolute path to directory or data file.
        obs_type: If "auto", the script looks for unambiguous
            directory contents (e.g., a bunch of Stage IV grib2 files).
            If there is ambiguity, obs_type needs to be a string
            that denotes which type is to be loaded. An exception
            is raised if this is left at "auto".
        load_objects (bool): If False, the catalogue
            dictionary is strings of filepaths and metadata only.
            Otherwise, objects are loaded into memory (how does
            the logic know what object to apply the class to?)
        st4_1h, st4_6h, st4_24h (bool): if True, and if the obs_type
            is set or determined to be Stage IV, these accumulation
            periods will be loaded into the catalogue.
        
    """
    def __init__(self,fpath,obs_type='auto',load_objects=False,
                    st4_1h=True,st4_6h=False,st4_24h=False):

        # By default, the full domain is returned
        self.return_subdomain = False

        self.load_objects = load_objects

        # Gather all options specific to certain datatypes.
        self.st4_opts = {k:v for k,v in locals().items() if k.startswith("st4")}

        # Check user options to see which accumulation times are wanted.
        self.hlist = self.return_st4_accumlist()
        assert self.hlist 
        
        self.INSTANCES = self.return_instances()
        self.PATTERNS = self.return_globpatterns()

        # A list of all valid times within. Unique values only
        self.utcs = set()

        # Determine if fpath is a directory or file
        if os.path.isdir(fpath):
            self.fdir = fpath
        else:
            self.fdir = os.path.dirname(fpath)

        # A dictionary catalogue of all obs files/instances
        self.catalogue = self.fill_catalogue(obs_type)

    @staticmethod
    def return_instances():
        INSTANCES = dict(stageiv = StageIV,
                        radar = Radar,
                        )
        return INSTANCES

    def lookup_instance(self,key=None):
        if not key:
            key = self.obs_type
        return self.INSTANCES[key]

    @staticmethod
    def return_globpatterns():
        PATTERNS = dict(stageiv = 'ST4*',
                        radar = "n0q_*.png",
                        )
        return PATTERNS

    def lookup_globpattern(self,key=None,reverse=False):
        """ 
        Args:
            key: if reverse is False, this is the obs_type string.
                if reverse is True, this is the glob pattern.
        """
        if (key is None) and (reverse is False):
            key = self.obs_type
        if reverse:
            PATTERNS = {v: k for k, v in self.PATTERNS.items()}
        else:
            PATTERNS = self.PATTERNS
        return PATTERNS[key]

    def fill_catalogue(self,obs_type):
        """ Gather observation files into dictionary.

        If self.load_objects is True, instances will be created.

        Args:
            obs_type: see Class docstring.
        """
        self.catalogue = {}
        if obs_type == 'auto':
            fs = glob.glob(os.path.join(self.fdir,'*'))
            for ot,pattern in self.PATTERNS.items():
                startstr = pattern.split("*")[0]
                print("Looking for the pattern {}".format(startstr))
                fs_check = [os.path.basename(f).startswith(startstr) for f in fs]
                if all(fs_check):
                    self.obs_type = ot
                    break
            else:
                print("No single obs_type found. Please specify.")
                raise Exception
        
        else:
            self.obs_type = obs_type
        globpattern = self.lookup_globpattern()
        fs = glob.glob(os.path.join(self.fdir,globpattern))

        self.instance = self.lookup_instance()
        if obs_type == 'stageiv':
            
            for h in self.hlist:
                fs_h = [f for f in fs if f.endswith(h)]

                for f in fs_h:
                    utc = self.instance.date_from_fname(f)
                    self.utcs.add(utc)
                    try:
                        self.catalogue[utc]
                    except KeyError:
                        self.catalogue[utc] = {h:{} for h in self.hlist}

                    # pdb.set_trace()
                    # if self.catalogue[utc]
                    self.catalogue[utc][h] = {'fpath':f}

                    if self.load_objects:
                        self.catalogue[utc][h]['loadobj'] = self.instance(f)

        elif obs_type == 'radar':
            # pdb.set_trace()
            for f in fs:
                utc = self.instance.date_from_fname(f)
                self.catalogue[utc] = {'fpath':f}
                if self.load_objects:
                    self.catalogue[utc]['loadobj'] = self.instance(f)

        return self.catalogue

    def catalogue_get(self,utc,auto_download=True,accum_hr=None):
                    #,loadobj=False):
        """ Load instance from catalogue.

        If it doesn't exist, download it.
        """
        if (accum_hr == None) and (self.obs_type == 'stageiv'):
            if len(self.hlist) == 1:
                accum_hr = hlist[0]
            else:
                raise Exception("Specify a single accum_hr that is needed.")

        try:
            self.catalogue[utc]
        except KeyError:
            if auto_download:
                self.catalogue[utc] = {'loadobj':self.instance(utc,self.fdir)}
                self.catalogue[utc]['fpath'] = self.catalogue[utc]['loadobj'].fpath
            else:
                assert 1==0
        finally:
            ob = self.catalogue[utc]

        if self.obs_type == 'stageiv':
            return ob[accum_hr]
        return self.instance(utc,ob['fpath'])
        

    def return_st4_accumlist(self,):
        hlist = []
        if self.st4_opts['st4_1h']:
            hlist.append("01h")
        if self.st4_opts['st4_6h']:
            hlist.append("06h")
        if self.st4_opts['st4_24h']:
            hlist.append("24h")
        return hlist

    def get(self,utc,*args,**kwargs):
        ob = self.catalogue_get(utc)
        data = ob.get(*args,**kwargs)
        if self.return_subdomain:
            data = ob.get_subdomain(data=data,**self.limdict)

        
    def arbitrary_pick(self,load_object=False):
        """ 
        Todo:
            * How to make this arbitrary lookup general to any data type?
        """
        raise Exception("Not implemented.")
        # arb_f = self.catalogue[self.utcs[0])
        # if load_object:
            # self.instance(arb_f)
            
    def set_subdomain(self,Nlim,Elim,Slim,Wlim,enable=True):
        """ Set limited domain. Whenever get is called, cut down domain.

        Call signature changed from Obs method, to
        enable a time to be chosen, and data to be
        returned (not overridden). These settings are
        kept for future get() commands.
        """
        self.Nlim = Nlim
        self.Elim = Elim
        self.Slim = Slim
        self.Wlim = Wlim

        self.limdict = dict(Nlim = Nlim,
                            Elim = Elim,
                            Slim = Slim,
                            Wlim = Wlim)

        if enable:
            self.return_subdomain = True
        return 

    def plot(self,utc,accum_hr=None,*args,**kwargs):
        """ Wrapper for members' plot methods.
        """
        ob = self.catalogue_get(utc,accum_hr=accum_hr)
        if self.return_subdomain:
            for k,v in self.limdict.items():
                kwargs[k] = v
        ob.plot(*args,**kwargs)
        return

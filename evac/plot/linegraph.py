""" Plotting basic line x-y graphs.

Todo:
    * __get_plot_options for all plots
    * Merge the types of plots
    * I seem to be typing the same stuff for `LineGraph`, `Verif`, etc...
        Can I not make it a method in `Figure`?
"""

import os
import pdb

from evac.plot.figure import Figure

class LineGraph(Figure):
    """ Basic graphs.

    To enable plotting of multiple lines, 

    For all plotting methods, there are three different args/kwargs:
        * plotkwargs, that are passed to ax.plot().
        * mplkwargs, that set various things (labels, titles, etc)
        * clskwargs, that sets general options relating to the class
    The private __get_plot_options finds defaults for all three kwargs.

    Args:
        x: ?
    
    Example:
        To plot forecast value::
            
            from evac.stats.forecastvalue import ForecastValue
            FV = ForecastValue()
            G = Graph()
            G.plot_fv(FV,hold=False)
            
    """
    def __init__(self,fpath=None,outdir=None,vrbl=None,fname=None,figsize=None,
                    mplkwargs=None,log=False,**kwargs):
        # From superclass
        #self.name = self.create_fname(fname=fname)
        # outfpath = os.path.join(self.outdir)
        if not fpath:
            fpath = os.path.join(outdir,fname)

        super().__init__(fpath=fpath,figsize=figsize,mplkwargs=mplkwargs,**kwargs)#ncol=1,nrow=1)
        self.vrbl = vrbl
        self.log = log
        
    def plot(self,xdata,ydata,use_plot_date=False, 
                    mplkwargs=None,plotkwargs=None,
                    # label=None,
                    *args,**kwargs):
        """ Plot e.g. deterministic scores.
        
        Examples:
            plot_score(xdata=range(len(threshs)),ydata=ydata,
                            labels=model,color=colors[model],linestyle=lines[model])

        Todo:
            * Can the labels argument be put into mplkwargs?
        """
        # self.mplkwargs.update(mplkwargs)
        # self.clskwargs.update(**kwargs)
        self.update_kwargs(mplkwargs,mpl=True)
        self.update_kwargs(kwargs,cls=True)
        if plotkwargs is None:
            plotkwargs = dict()

        if use_plot_date:
            plotfunc = self.ax.plot_date
            plotkwargs['xdate'] = True
            plotkwargs['marker'] = '.'
            plotkwargs['linestyle'] = 'solid'
        elif self.log == 'x':
            plotfunc = self.ax.semilogx
        elif self.log == 'y':
            plotfunc = self.ax.semilogy
        else:
            plotfunc = self.ax.plot
        plotfunc(xdata,ydata,**plotkwargs)

        # For legend, xlabel, title, xticks, etc
        self.ax.set(**self.mplkwargs)

        # if self.clskwargs['save']:
        if self.save_opt:
            self.save()

        return


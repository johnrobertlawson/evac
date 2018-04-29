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
    def __init__(self,outdir,vrbl,fname=None,figsize=None):
        # From superclass
        #self.name = self.create_fname(fname=fname)
        # outfpath = os.path.join(self.outdir)
        self.outdir = outdir
        self.fname = fname
        self.vrbl = vrbl

        super().__init__(self,figsize=figsize)#ncol=1,nrow=1)

        
    def plot_score(self,xdata,ydata,
                    use_plot_date=False,labels=None, 
                    mplargs=None,mplkwargs=None,
                    plotargs=None,plotkwargs=None,
                    *args,**kwargs):
        """ Plot e.g. deterministic scores.
        
        Examples:
            plot_score(xdata=range(len(threshs)),ydata=ydata,
                            labels=model,color=colors[model],linestyle=lines[model])

        Todo:
            * Can the labels argument be put into mplkwargs?
        """
        clskwargs, mplkwargs, plotkwargs = self._get_plot_options2(plotargs=plotargs,
                            plotkwargs=plotkwargs,mplargs=mplargs,mplkwargs=mplkwargs,
                                                # vrbl=self.vrbl,
                                                *args,**kwargs)

        if use_plot_date:
            plotfunc = self.ax.plot_date
            plotkwargs['xdate'] = True
            plotkwargs['marker'] = '.'
            plotkwargs['linestyle'] = 'solid'
        else:
            plotfunc = self.ax.plot
        plotfunc(xdata,ydata,**plotkwargs)

        # For legend, xlabel, title, xticks, etc
        self.ax.set(**mplkwargs)

        if clskwargs['save']:
            self.save()


""" Plot scores as a box plot.
"""

import pdb

from evac.plot.figure import Figure

class BoxPlot(Figure):
    """ Plot box-and-whisker in violin style.

    Examples:
        First, create a dictionary of scores::

            # Something like:
            # HRRR_scores = dict(scorename = score)

            VP = ViolinPlot(HRRR_scores,NEWSe_scores,
                    fname='vp.png',outdir='./')
            VP.plot()

    Args:
        *args: A number of dictionaries containing
            score data
    """
    def __init__(self,outdir,fname='test_violin',*args,**kwargs):
        self.fname = fname
        self.outdir = outdir
        super().__init__(self,*args,**kwargs)

    def plot(self,data,mplargs=None,mplkwargs=None,
                    plotargs=None,plotkwargs=None,
            *args,**kwargs):
        """ Plot the things.
        """
        clskwargs, mplkwargs, plotkwargs = self._get_plot_options2(plotargs=plotargs,
                            plotkwargs=plotkwargs,mplargs=mplargs,mplkwargs=mplkwargs,
                            # vrbl=self.vrbl,
                            *args,**kwargs)
        # self.get_options?
        self.ax.boxplot(data,*args,**kwargs)
        # TODO: implement a setter method.
        # self.pass_mpl_settings()

        if clskwargs['save']:
            self.save()


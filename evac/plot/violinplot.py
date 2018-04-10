""" Plot scores as a violin plot.
"""

import pdb

from evac.plot.figure import Figure

class ViolinPlot:
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
    def __init__(self,*args,**kwargs):
        super().__init__(self,*args,**kwargs)

    def plot(self,mplargs,mplkwargs,*args,**kwargs):
        """ Plot the things.
        """
        # self.get_options?
        self.ax.violinplot()
        self.pass_mpl_settings()

        if kwargs['save']:
            self.save()


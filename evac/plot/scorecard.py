""" Plot scores for numerous models side by side.

"""

import os
import pdb

from evac.plot.figure import Figure

class ScoreCard(Figure):
    """ Plot scores on a table-like scorecard.

    Note:
        * Args get passed straight to `evac.plot.figure.Figure`.

    Examples:
        First, create a dictionary of scores::
            
            # Something like:
            # HRRR_scores = dict(scorename = score)

            SC = ScoreCard(HRRR_scores,NEWSe_scores,
                    fname='sc.png',outdir='./')
            SC.plot()

    Args:
        *args: Any number of dictionaries with a given
            model's scores. See Examples.
        outdir: absolute path to output directory
        fname: file name of plot

    """
    def __init__(self,*args,**kwargs):
        pass

    def create_scorecard(self):
        pass

    def plot(self,):
        pass

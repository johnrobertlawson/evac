"""Briar and related scores/metrics.
"""

import os
import pdb

import numpy as N

from .detscores import DetScores
from .probscores import ProbScores
 
class Briar(ProbScores):
    def __init__():
        pass

    def compute_briar():
        """
        The better the forecast, the smaller the Briar score.
        Ideal limit is 0.

        Args:
        
        N   :   number of verifying observations
        gi  :   occurrence prob of event E for a verif date i
                    Must be between 0 and 1 inclusive.
        Oi  :   If E occurred, this is 1; otherwise 0.
                    Must be either 0 or 1.
        """
        BS = (1/self.N) * N.sum((gi - Oi)**2)

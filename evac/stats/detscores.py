import math
import pdb
import os

import numpy as N

from evac.stats.misc_stats import compute_contingency

class DetScores:
    """ Miscellaneous (and more simple) deterministic scores.

    Using Buizza 2001, MWR, as follows:
        * a = observed and forecast
        * b = forecast but not observed
        * c = observed but not forecast
        * d = neither observed nor forecast

    Users must supply the forecast and observation array (with threshold
    and over/under operation), or either arr or (a,b,c,d).

    Args:
        fcst_arr,obs_arr (N.ndarray)    :   Forecast and verification arrays,
                                                the same dimensions.
        thresh (float)  :   threshold
        overunder (str) :   over, under, overeq, undereq, equal
        arr2x2 (N.ndarray) :   Array representing 2x2 matrix.
                                Must be arranged so flatten array is
                                equal to a tuple of (a,b,c,d).
        a,b,c,d (int)   :   Number of elements in 2x2 matrix.
    """
    def __init__(self,arr2x2=None,a=None,b=None,c=None,d=None,
                    fcst_arr=None,obs_arr=None,thresh=None,overunder=None):
        if arr2x2 is not None:
            a,b,c,d = arr.flatten()
        elif fcst_arr is not None:
            a,b,c,d = compute_contingency(fcst_arr,obs_arr,thresh,overunder,
                                            fmt='tuple')
        else:
            for x in (a,b,c,d):
                assert x is not None

        # Number of total events
        self.n = sum((a,b,c,d))

        # Now add a tiny amount to zeros to avoid division error
        a,b,c,d = [N.finfo(N.float).eps if x == 0 else x for x in (a,b,c,d)]

        self.a = a
        self.b = b
        self.c = c
        self.d = d

        # Contingency table
        self.an = self.a/self.n
        self.bn = self.b/self.n
        self.cn = self.c/self.n
        self.dn = self.d/self.n

        # Lookup table for methods
        self.scores = self.assign_scores_lookup()

    def check_approx(self,x,y):
        if not math.isclose(x,y,rel_tol=0.0001):
            raise Exception("Values are not close.")

    def compute_pcli(self):
        """Observed frequency of event under consideration.
        Also known as base rate.
        """
        pcli = (self.a + self.c)/self.n
        return pcli

    def compute_hitrate(self,method=2):
        """ Accuracy.

        Eq. 6 in B01."""
        if method == 1:
            # This is Buizza version - same as POD
            HR = (self.a + self.d)/self.n
        else:
            # This is Jolliffe and Stephenson
            HR = self.a / (self.a + self.c)
        return HR

    def compute_falsealarmrate(self):
        # This was wrong was round with below!
        F = self.b / (self.b + self.d)
        return F

    def compute_falsealarmratio(self):
        """ Same as PFD"""
        FAR = self.b / (self.a + self.b)
        return FAR

    def compute_propcorrect(self):
        PC = (self.a + self.d)/self.n
        return PC

    def compute_E(self):
        # 3.17 in Jolliffe and Stephenson
        E = ( ( ((self.a+self.c)/self.n)*
                ((self.a+self.b)/self.n)) +
            ( ((self.b+self.d)/self.n)*
                ((self.c+self.d)/self.n)) )

        return E

    def compute_heidke(self):
        E = self.compute_E()
        HSS = (self.compute_propcorrect() - E)/(1-E)
        return HSS

    def compute_peirce(self):
        PSS = ( (self.a*self.d)-(self.b*self.c)/
                    ((self.b+self.d)*(self.a+self.c)) )
        return PSS

    def compute_csi(self):
        CSI = self.a / (self.a+self.b+self.c)
        return CSI

    def compute_ar(self):
        """Number of hits expected by pure chance
        """
        ar = ( (self.a+self.b)*(self.a + self.c))/self.n
        return ar

    def compute_gilbert(self):
        ar = self.compute_ar()
        GSS = (self.a - ar) /(
                (self.a - ar + self.b + self.c))
        return GSS

    def compute_yuleq(self):
        Q = (( (self.a*self.d) - (self.b*self.c) )/
             ( (self.a*self.d) + (self.b*self.c) ))
        return Q

    def compute_threat(self):
        """ Accuracy.

        Eq. 7 in B01."""
        TS = self.a/(self.a+self.b+self.c)
        return TS

    def compute_pod(self):
        """Accuracy - probability of detection.

        Eq. 8 in B01."""
        POD1 = self.a/(self.a+self.c)
        POD2 = (self.a/self.n) * (1/self.compute_pcli())
        self.check_approx(POD1,POD2)
        return POD1

    def compute_pfd(self):
        """Accuracy - probability of false detection.

        Eq. 9 in B01."""
        PFD1 = self.b/(self.b+self.d)
        PFD2 = (self.b/self.n) * (1/(1-self.compute_pcli()))
        self.check_approx(PFD1,PFD2)
        return PFD1

    def compute_bias(self):
        """ Eq. 10 in B01.
        """
        B = (self.a+self.b)/(self.a+self.c)
        return B

    def compute_kss(self):
        """ Eq. 13 in B01.
        """
        KSS1 = self.compute_pod() - self.compute_pfd()
        KSS2 = ( (self.a*self.d) - (self.b*self.c))/(
                 (self.a+self.c) * (self.b+self.d))
        self.check_approx(KSS1,KSS2)
        return KSS1

    def compute_successrate(self):
        """ 1-FAR.
        """
        SR = 1.0 - self.compute_falsealarmratio()
        return SR

    def assign_scores_lookup(self):
        scores = {'HR':self.compute_hitrate,
                'FAR':self.compute_falsealarmratio,
                'TS':self.compute_threat,
                'POD':self.compute_pod,
                'PFD':self.compute_pfd,
                'BIAS':self.compute_bias,
                'KSS':self.compute_kss,
                'CSI':self.compute_csi,
                'PCLI':self.compute_pcli,
                'F':self.compute_falsealarmrate,
                'PC':self.compute_propcorrect,
                'HSS':self.compute_heidke,
                'PSS':self.compute_peirce,
                'GSS':self.compute_gilbert,
                'Q':self.compute_yuleq,
                'SR':self.compute_successrate,
                }
        return scores

    def get(self,score):
        """ Getter method for obtaining a score via string.
        """
        if score in self.scores.keys():
            val = self.scores[score]()
            return val
        else:
            raise Exception("Specify from the following: \n{}".format(self.scores.keys()))

    def compute_all(self,datadir=None,fname=None,):
        """ If datadir/fname are not None, computed stats are saved.
        Else, they are just returned in the dictionary.
        """
        SCORES = dict(a=self.a, b=self.b, c=self.c, d=self.d)
        for score in self.scores.keys():
            SCORES[score] = self.get(score)
            
        if datadir:
            datapath = os.path.join(datadir,fname)
            print("Saving all data arrays to {}".format(datapath))
            N.savez(file=datapath,**SCORES)
            print("Saved.")
        return SCORES
        

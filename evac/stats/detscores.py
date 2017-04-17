""" Miscellaneous (and more simple) deterministic scores.

Using Buizza 2001, MWR.

a = observed and forecast
b = forecast but not observed
c = observed but not forecast
d = neither observed nor forecast
"""

from evac.stats.misc_stats import compute_contingency

class DetScores:
    def __init__(self,arr2x2=None,a=None,b=None,c=None,d=None,
                    fcst_arr=None,obs_arr=None,thresh=None,overunder=None):
        """
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

        self.a = a
        self.b = b
        self.c = c
        self.d = d

        # Contingency table
        self.an = self.a/self.n
        self.bn = self.b/self.n
        self.cn = self.c/self.n
        self.dn = self.d/self.n

    def compute_pcli(self):
        """Observed frequency of event under consideration.
        """
        pcli = (self.a + self.c)/self.n
        return pcli

    def compute_hitrate(self):
        """ Accuracy.

        Eq. 6 in B01."""
        HR = (self.a + self.d)/self.n
        return HR

    def compute_falsealarm(self):
        FAR = self.b / (self.a + self.b)
        return FAR

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
        assert POD1 == POD2
        return POD1

    def compute_pfd(self):
        """Accuracy - probability of false detection.

        Eq. 9 in B01."""
        PFD1 = self.b/(self.b+self.d)
        PFD2 = (self.b/self.n) * (1/(1-self.compute_pcli()))
        assert PFD1 == PFD2
        return PFD1

    def compute_bias(self):
        """ Eq. 10 in B01.
        """
        B = (self.a+self.b)/(self.a+self.c)
        return B

    def compute_kss(self):
        """ Eq. 13 in B01.
        """
        KSS = self.compute_pod() - self.compute_pfd()
        return KSS

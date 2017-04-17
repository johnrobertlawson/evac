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
        Also known as base rate.
        """
        pcli = (self.a + self.c)/self.n
        return pcli

    def compute_hitrate(self):
        """ Accuracy.

        Eq. 6 in B01."""
        # This is Buizza version?
        #HR = (self.a + self.d)/self.n
        # This is Jolliffe and Stephenson
        HR = self.a / (self.a + self.c)
        return HR

    def compute_falsealarmratio(self):
        F = self.b / (self.b + self.d)
        return F

    def compute_falsealarmrate(self):
        FAR = self.b / (self.a + self.b)
        return FAR

    def compute_propcorrect(self):
        PC = (self.a + self.d)/self.n
        return PC

    def compute_E(self):
        # 3.17 in Jolliffe and Stephenson
        E = ( ((self.a+self.c)/self.n)*
                ((self.a+self.b)/self.n)) +
            ( ((self.b+self.d)/self.n)*
                ((self.c+self.d)/self.n))

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
        GSS = (self.a - ar) /
                (self.a - ar + self.b + self.c)
        return GSS

    def compute_yuleq(self):
        Q = ( (self.a*self.d) - (self.b*self.c) )/
            ( (self.a*self.d) + (self.b*self.c) )
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

    def get(self,score):
        """ Getter method for obtaining a score via string.
        """
        scores = {'HR':self.compute_hitrate,
                    'FAR':self.compute_falsealarmratio,
                    'TS':self.compute_threat,
                    'POD':self.compute_pod,
                    'PFD':self.compute_pfd,
                    'BIAS':self.compute_bias,
                    'KSS':self.compute_bias,
                    'CSI',self.compute_csi,
                    'PCLI':self.compute_pcli,
                    'F':self.compute_falsealarmrate,
                    'PC':self.compute_propcorrect,
                    'HSS':self.compute_heidke,
                    'PSS':self.compute_peirce,
                    'GSS':self.compute_gilbert,
                    'Q':self.compute_yuleq,
                    }

        if score in scores.keys():
            return scores[score]()
        else:
            raise Exception("Specify from the following: \n{}".format(scores.keys()))

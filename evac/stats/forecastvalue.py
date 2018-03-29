import numpy as N
import pdb

from evac.utils.exceptions import FormatError, NonsenseError
from evac.stats.detscores import DetScores

class ForecastValue(DetScores):
    """ Compute FV, forecast value, as Buizza 2001.

    Example that returns FV for exceeding 0.1 units (e.g., QPF)
    using the default range of cost/loss ratios:
        FV = ForecastValue(fcst_arr=arr1,obs_arr=arr2,thresh=0.l,
                            overunder='over')

    Then, access the FV for a given cost/loss ratio with:
        FV.get(0.02)

    This will return None is the cost/lost ratio doesn't exist.

    Todo:
        * Shouldn't this be a mix-in of `ProbScores` and `DetScores`?

    Args:
        CLs (N.ndarray, tuple, list): Cost-loss ratios to compute
            FV as function of each ratio
        args,kwargs: Arguments sent to the superclass to compute 2x2 contigency.
            It's probably best to use keywords (`kwargs`).
    """
    def __init__(self,CLs=None,*args,**kwargs):
        super().__init__(*args,**kwargs)
        # super(DetScores,self).__init__(*args,**kwargs)
        self.pod = self.compute_pod()
        self.pfd = self.compute_pfd()
        self.pcli = self.compute_pcli()
        self.kss = self.compute_kss()

        if CLs is None:
            CLs = N.arange(0.01,1.01,0.01)
        elif isinstance(CLs,(int,float,N.ndarray,tuple,list)):
            CLs = N.array(CLs)
        else:
            raise FormatError("CLs argument must be one of the following: \n"
                              "int, float, N.ndarray, tuple, list.")
        self.CLs = CLs

        self.FVs = []
        for cl in self.CLs:
            self.FVs.append(self.compute_FV(cl))

    def compute_FV(self,cl):
        """ Compute FV.

        Note the variable "x" is the minimum of C/L ratio and P_cli.

        Args:
            cl (float)  :   Cost/loss ratio

        """
        # For ease of reading
        x = min(cl,self.pcli)

        # Eq. 26 in Buizza 2001
        FV = self.kss - (
            # Numerator
            ( ((1-self.pod) * (self.pcli-x)) + (self.pfd*(cl-x)) )/
            # Denominator
            (x - (self.pcli*cl))
        )

        # FV2 = self.kss - ( (1-self.pod)*(
                        # (self.pcli-cl)/
                        # (cl*(1-self.pcli))) )

        return FV

    def get_FV(self,cl):
        """ Return the FV for a given cost/loss ratio, otherwise None.
        """
        clidx = N.where(self.CLs == cl)
        if len(clidx) == 1:
            return self.FV[clidx]
        elif len(clidx) == 0:
            print("No C/L ratio found for that cl value")
            return None
        else:
            raise NonsenseError("C/L ratio is ambiguous.")

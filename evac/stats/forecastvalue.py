import pdb
import operator

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
    def __init__(self,fcst_arr,obs_arr,thresh,overunder,CLs=None):
        self.overunder = overunder
        self.fcst_arr = fcst_arr
        self.thresh = thresh
        self.obs_arr = obs_arr

        if fcst_arr.ndim == 3:
            if fcst_arr.shape[0] == 1:
                fcst_arr = fcst_arr[0,:,:]
                self.FVens = False
            else:
                self.arr2x2s, self.J = self.compute_FV_prob()
                self.FVens = True
        else:
            self.FVens = False

        if self.FVens:
            self.pod = dict()
            self.pfd = dict()
            self.pcli = dict()
            self.kss = dict()
            for j in self.J:
                DS = DetScores(arr2x2=self.arr2x2s[j])
                self.pod[j] = DS.compute_pod()
                self.pfd[j] = DS.compute_pfd()
                self.pcli[j] = DS.compute_pcli()
                self.kss[j] = DS.compute_kss()

        else:
            super().__init__(fcst_arr=fcst_arr,obs_arr=obs_arr,thresh=thresh,
                            overunder=overunder)
            
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
            if self.FVens:
                self.FVs.append(self.compute_FVens(cl))
            else:
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

    def compute_FVens(self,cl):
        fvs = []
        for j in self.J:
            x = min(cl, self.pcli[j])
            FV = self.kss[j] - (
                ( ((1-self.pod[j]) * (self.pcli[j]-x)) + (self.pfd[j]*(cl-x)) )/
                (x - (self.pcli[j]*cl)))
            fvs.append(FV)
        return max(fvs)

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


    def compute_FV_prob(self):
        self.nmems = self.fcst_arr.shape[0]
        # OU = dict(over = operator.gt,
                    # under = operator.lt,)
        OU = dict(over = N.greater,
                    under = N.less)
        prob_arr = N.sum(OU[self.overunder](self.fcst_arr,self.thresh),axis=0)/self.nmems
        yesno_arr = OU[self.overunder](self.obs_arr,self.thresh)
        X = dict()
        Y = dict()
        J = list()
        for idx in range(self.nmems):
            j = idx+1
            J.append(j)
            p0 = (j-1)/self.nmems
            p1 = (j/self.nmems)
            # yes = N.where(yesno_arr == True)
            # no = N.where(yesno_arr == False)
            # pr_idx = N.where((p1 > prob_arr) & (prob_arr >= p0))
            X[j] = N.where((p1 > prob_arr) & (prob_arr >= p0) & (yesno_arr == True))
            Y[j] = N.where((p1 > prob_arr) & (prob_arr >= p0) & (yesno_arr == False))
            # X[j] = N.intersect1d(pr_idx,yes)
            # Y[j] = N.intersect1d(pr_idx,no)
        A = dict()
        B = dict()
        C = dict()
        D = dict()
        arr2x2s = dict()
        for j in J:
            jidx = J.index(j)
            A[j] = sum([len(X[k][0]) for k in J[jidx+1:]])
            B[j] = sum([len(Y[k][0]) for k in J[jidx+1:]])
            C[j] = sum([len(X[k][0]) for k in J[1:jidx+1]])
            D[j] = sum([len(Y[k][0]) for k in J[1:jidx+1]])
            arr2x2s[j] = (A[j], B[j], C[j], D[j])
        # pdb.set_trace()
        return arr2x2s, J



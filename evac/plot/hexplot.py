import os
import pdb

import numpy as N
import seaborn as sns

#from evac.plot.figure import Figure

class HexPlot:
    def __init__(self,fpath):
        super().__init__(fpath=fpath)

    def plot(self,df,xname,yname,kind=hex,color='#4CB391'):
        hex = sns.jointplot(xname,yname,data=df,kind='hex')
        hex.savefig(self.fpath)
        print("Saved figure to",fpath)
        return

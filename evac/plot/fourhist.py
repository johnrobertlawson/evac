import os
import pdb

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='darkgrid')

class FourHist:
    def __init__(self,fpath):
        self.fpath = fpath

    def plot(self,df,xname,yname,dataname,bins=20):
        g = sns.FacetGrid(df,row=yname,col=xname,margin_titles=True)
        g.map(plt.hist, dataname, color='green',bins=bins)
        g.savefig(self.fpath)
        print("Saved figure to",self.fpath)
        return

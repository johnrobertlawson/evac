import os
import pdb

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

class PairPlot:
    def __init__(self,fpath):
        self.fpath = fpath

    def plot(self,df,color_name,vars=None):
        g = sns.pairplot(df,hue=color_name,palette='husl',vars=vars,
                            diag_kind='kde',kind='scatter',
                            plot_kws=dict(s=10,alpha=0.5,linewidth=0.1),
                            diag_kws=dict(shade=True,),
                            )
        g.savefig(self.fpath)
        print("Saved figure to",self.fpath)
        return

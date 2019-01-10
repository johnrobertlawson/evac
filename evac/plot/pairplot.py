import os
import pdb

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

class PairPlot:
    def __init__(self,fpath):
        self.fpath = fpath

    def plot(self,df,color_name,vars=None,palette='husl'):
        #current_palette = sns.color_palette(palette)
        #PAL = {}
        #for nc, c in enumerate(df[color_name].unique()):
        #    PAL[c] = current_palette[nc]

        g = sns.pairplot(df,hue=color_name,diag_kind='kde',kind='scatter',
                            palette=palette,vars=vars,
                            #vars=vars,palette=PAL,
                            plot_kws=dict(s=10,alpha=0.5,linewidth=0.1),
                            diag_kws=dict(shade=True,),
                            hue_order=sorted(df[color_name].unique()),
                            #legend=False,
                            )
        #g.add_legend(label_order = sorted(g._legend_data.keys(), key = str))
        g.savefig(self.fpath)
        print("Saved figure to",self.fpath)
        plt.close(g.fig)
        return

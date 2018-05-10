import os
import pdb

import matplotlib.pyplot as plt

class ColorBar(Figure):
    def __init__(self,fpath):
        fig = plt.figure()
        ax = plt.gca(fig)
        super().__init__(fpath=fpath,fig=fig,ax=ax)

    def create_colorbar(self,cf,label=None,ticks=False):
        """ Create colorbar.

        Todo:
            * Move to ColorBar(Figure)

        Inputs:
        fpath   :   path to file
        fname   :   filename
        cf      :   contour filling for legend
        label   :   colorbar label

        """
        CBax = self.fig.add_axes([0.15,0.05,0.7,0.02])
        CB = plt.colorbar(cf,cax=CBax,orientation='horizontal',ticks=ticks)
        CB.set_label(label)
        self.save(tight=False)
        return

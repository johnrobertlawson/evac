import pdb

# import

from evac.plot.figure import Figure

class Histogram(Figure):
    """ Plot a histogram.

    Args:
        data: can be 1D or 2D array. If 2D, then two
            datasets are plotted.
    """

    def __init__(self,data,fpath):
        self.data = data
        self.fpath = fpath

        super().__init__()

    def plot(self,labels=False,legend=True):
        self.ax.hist(self.data,label=labels,bins=12)
        if legend:
            self.ax.legend()
        self.save()
        return

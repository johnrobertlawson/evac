import pdb
import os
import datetime

import numpy as N
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from evac.datafiles.csvfile import CSVFile

class StormEvents(CSVFile):
    """ NWS storm events from (which website?).

    Todos:
        * Move or remove redundant plotting method

    Args:
        fpath: Absolute path to CSV file.
    """
    def __init__(self,fpath,):
        self.r = N.genfromtxt(fpath,dtype=None,
                names=True,delimiter=',',)#missing='',filling_values='none')
        # import pdb; pdb.set_trace()
        self.convert_times()

    def convert_times(self,):
        LOLtimes = [s for s in self.r['BEGIN_TIME']]
        LOLdates = [s.decode() for s in self.r['BEGIN_DATE']]
        # LOLtimes = self.r['BEGIN_TIME']
        padtimes = []
        for t in LOLtimes:
            intt = int(t)
            padtimes.append('{0:04d}'.format(intt))

        hours = [s[:-2] for s in padtimes]
        mins = [s[-2:] for s in padtimes]

        # pdb.set_trace()

        self.datetimes = N.array([datetime.datetime.strptime(s+h+m,'%m/%d/%Y%H%M')
                        for s,h,m in zip(LOLdates,hours,mins)])
        # import pdb; pdb.set_trace()
        # import numpy.lib.recfunctions
        # self.r = numpy.lib.recfunctions.append_fields(self.r,'datetimes',N.array(dates))

    def plot(self,reports,itime,ftime,fname,outdir,Nlim=False,
            Elim=False,Slim=False,Wlim=False,
            annotate=True,fig=None,ax=None,ss=80,color=None):
        if "+" in reports:
            reportlist = reports.split("+")
            hold = True
        else:
            reportlist = reports
            hold = False
        m = None
        for report in reportlist:
            LOLtypes = [s.decode() for s in self.r['EVENT_TYPE']]
            reportidx = N.array([n for n,t in zip(list(range(len(LOLtypes))),LOLtypes) if report in t])
            lateidx = N.where(self.datetimes > itime)
            earlyidx = N.where(self.datetimes < ftime)
            timeidx = N.intersect1d(earlyidx,lateidx,)#assume_unique=True)
            plotidx = N.intersect1d(reportidx,timeidx)

            # Need to edit if color is to be custom 
            assert report in ("Wind","Hail")
            if report == "Wind":
                color = "blue"
            elif report == "Hail":
                color = "lightgreen"
            else:
                assert True is False
            print("Plotting",report)

            if fig is None:
                # if this is the first time round a loop of reports, or custom fig arg
                fig,ax = plt.subplots(1,figsize=(6,8))

            if hold is True:
                if m is None:
                    m = Basemap(projection='merc',
                                llcrnrlat=Slim,
                                llcrnrlon=Wlim,
                                urcrnrlat=Nlim,
                                urcrnrlon=Elim,
                                lat_ts=(Nlim-Slim)/2.0,
                                resolution='h',
                                ax=ax)

                    m.drawcoastlines()
                    m.drawstates()
                    m.drawcountries()

            print("color =",color)
            m.scatter(self.r['BEGIN_LON'][plotidx],self.r['BEGIN_LAT'][plotidx],latlon=True,
                        c=color,s=ss,
                        marker='s',alpha=0.6,linewidths=1.5,
                        facecolors=color,edgecolors='black',
                        )
        fig.tight_layout()
        fpath = os.path.join(outdir,fname)
        plt.savefig(fpath)
        print("Saved to",fpath)

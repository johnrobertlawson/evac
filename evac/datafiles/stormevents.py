""" NWS storm events from <this website>...

"""

from .csvfile import CSVFile

class StormEvents(CSVFile):
    def __init__(self,fpath,):
        self.r = N.genfromtxt(fpath,dtype=None,
                names=True,delimiter=',',)#missing='',filling_values='none')
        # import pdb; pdb.set_trace()
        self.convert_times()

    def convert_times(self,):
        LOLtimes = self.r['BEGIN_TIME']
        padtimes = []
        for t in LOLtimes:
            intt = int(t)
            padtimes.append('{0:04d}'.format(intt))

        hours = [s[:-2] for s in padtimes]
        mins = [s[-2:] for s in padtimes]
        self.datetimes = N.array([datetime.datetime.strptime(s+h+m,'%m/%d/%Y%H%M')
                        for s,h,m in zip(self.r['BEGIN_DATE'],hours,mins)])
        # import pdb; pdb.set_trace()
        # import numpy.lib.recfunctions
        # self.r = numpy.lib.recfunctions.append_fields(self.r,'datetimes',N.array(dates))

    def plot(self,reports,itime,ftime,fname,outdir,Nlim=False,
            Elim=False,Slim=False,Wlim=False,
            annotate=True,fig=False,ax=False,ss=50,color='blue'):
        reportidx = N.array([n for n,t in zip(list(range(len(self.r['EVENT_TYPE']))),self.r['EVENT_TYPE']) if reports in t])
        lateidx = N.where(self.datetimes > itime)
        earlyidx = N.where(self.datetimes < ftime)
        timeidx = N.intersect1d(earlyidx,lateidx,)#assume_unique=True)
        plotidx = N.intersect1d(reportidx,timeidx)

        from mpl_toolkits.basemap import Basemap

        if fig==False:
            fig,ax = plt.subplots(1,figsize=(6,6))
        m = Basemap(projection='merc',
                    llcrnrlat=Slim,
                    llcrnrlon=Wlim,
                    urcrnrlat=Nlim,
                    urcrnrlon=Elim,
                    lat_ts=(Nlim-Slim)/2.0,
                    resolution='i',
                    ax=ax)

        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()

        m.scatter(self.r['BEGIN_LON'][plotidx],self.r['BEGIN_LAT'][plotidx],latlon=True,
                    marker='D',facecolors=color,edgecolors='black',s=ss)
        fig.tight_layout()
        plt.savefig(os.path.join(outdir,fname))

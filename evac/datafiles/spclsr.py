import pdb

from evac.datafiles.csvfile import CSVFile

class SPCLSR(CSVFile):
    """ Local storm reports from SPC site.

    Todo:
        * Work out the timing convention? Do we even need a `utc` argument?

    Args:
        utc (datetime.datetime): Time, but how does this work?!
        datapath: Absolute path to data CSV file.
        wind,hail,torn (bool): If `True`, use these data.
    """
    def __init__(self,utc,datapath,wind=True,hail=True,torn=True):
        self.utc = utc
        tt = utils.ensure_timetuple(utc)
        yr = str(tt[0])[-2:]
        mth = '{0:02d}'.format(tt[1])
        day = '{0:02d}'.format(tt[2])

        threats = []
        if wind:
            threats.append('wind')
        if hail:
            threats.append('hail')
        if torn:
            threats.append('torn')

        self.reports = {}
        for threat in threats:
            fname = '{0}{1}{2}_rpts_{3}.csv'.format(
                        yr,mth,day,threat)
            # Check to see if file exists
            fpath = os.path.join(datapath,fname)
            scan = glob.glob(fpath)
            if len(scan) == 0:
                url = 'http://www.spc.noaa.gov/climo/reports/'
                cmd = 'wget {0}{1} -P {2}'.format(
                                url,fname,datapath)
                os.system(cmd)

            if threat=='wind':
                names = ('time','speed','location','county',
                        'state','lat','lon')
                formats = ('S4','S4','S4','S4',
                            'S4','f4','f4')
            elif threat=='hail':
                names = ('time','size','location','county',
                        'state','lat','lon')
                formats = ('S4','S4','S4','S4',
                            'S4','f4','f4')
            elif threat=='torn':
                names = ('time','fscale','location','county',
                        'state','lat','lon')
                formats = ('S4','S4','S4','S4',
                            'S4','f4','f4')
            self.reports[threat] = N.loadtxt(fpath,dtype={'names':names,'formats':formats},
                                        skiprows=1,delimiter=',',usecols=list(range(8)))
            #times = reports['time']
            self.threats = threats

    def report_datenum(self,timestamp):
        """
        convert timestamp to datenum format.
        """
        tt = utils.ensure_timetuple(self.utc)
        itime_dn = calendar.timegm(tt[0],tt[1],tt[2],12,0,0)

        hr = int(timestamp[:2])
        mn = int(timestamp[2:])

        if hr<11:
            # After midnight UTC
            hr = hr+24

        tdelta = ((hr-12)*60*60) + (mn*60)
        return itime_dn + tdelta

    def plot_reports(self,fig=False,ax=False):
        plot_all = True
        for threat in threats:
            if plot_all:
                for t in self.reports[threat]['time']:
                    utc = self.report_datenum(t)

import pdb
import multiprocessing
from multiprocessing import Queue, Process, Pool

import numpy as N
import matplotlib as M
import matplotlib.pyplot as plt
import scipy.stats as SS

### 1D PLOT ###

### VARIABLES ###
R = 0 # Function, varying with scale t
dR = 0 # height of pulse - rainfatll intensity increment
rho = 3 # width of pulse - rainfall duration
ppcen = 0 # centre of pulse according to Poisson process 
mu = 0.025 # constant rate for Poisson process
rhotick = 0 # random duration
alpha = 5/3 # ?>>
H = 1/alpha # scaling exponent
lamb = 0 # ???
L0 = 0 # Interval min
L = 400 # Interval max
l = L/6 # Inner interval length
s = 8 # smoothing factor. S -> inf makes rectangular
u = 0 # distance from centre of echo
rho_i = 1 # minimum resolution or inner scale - in pixels
# rho_0 = mu*L # outer scale
rho_0 = 1200 # outer scale
n = 1.6e5 # number of pulses
# Enforce L/l between 3 and 6

def get_poisson(L,mu,size):
    return N.random.poisson(lam=L/mu,size=size)

def fractal_worker(size=False,high=False,rholist=False,xx=False,
                yy=False,delta=False,sigma=False,s=False,
                mu=False,pdf_full=False,R=False,q=False):
    centrex = N.random.random_integers(low=0,high=high,size=size)
    centrey = N.random.random_integers(low=0,high=high,size=size)
    count = 0
    for cenx, ceny in zip(centrex,centrey):
        rhotick = 0
        count += 1
        if not count%200:
            print('Doing plot #',count,"of",size)
        cenxP = get_poisson(cenx,mu,1)[0]
        cenyP = get_poisson(ceny,mu,1)[0]

        while rhotick < rho_i:
            rhotick = N.random.choice(rholist,p=pdf_full)
       
        dR = N.random.choice((-1,1)) * rhotick**H
        u = N.sqrt(((xx-cenx)**2)+((yy-ceny)**2))
        R += dR * N.exp(-1*((u**2/rhotick**2)-delta**2)/sigma**2)**(2*s)
    q.put(R)
    
if __name__ == '__main__':
    CPUs = multiprocessing.cpu_count()-1
    rho_i = 6
    rho_0 = 900
    l = 300
    L = 300
    # n = 160000
    # n = CPUs*9000
    n = CPUs*6000
    s = 2
    mu = 2.5
    # LAMB = 1.2 * N.pi
    LAMB = 1.2
    LAMBstar = N.sqrt(LAMB**2 - 1)
    delta = 0.5*(LAMBstar + LAMB)
    sigma = 0.5*(LAMB-LAMBstar)
    alpha = 5 # ?>>
    H = 1/alpha # scaling exponent


    rholist = N.arange(1,rho_0+1)
    cdf = 1/rholist
    pdf = N.diff(cdf)*-1
    pdf_full = N.append(pdf,1-pdf.sum())

    R = N.zeros([L,L],dtype=N.float64)
    xx, yy = N.mgrid[:L, :L]

    print("Starting computation now")
    # paralellise the processing
    npart = int(n/CPUs)
    # centrex = N.random.random_integers(low=0,high=L,size=n)
    # centrey = N.random.random_integers(low=0,high=L,size=n)
    Rs = []
    Rsum = N.zeros_like(R)


    # Pool = multiprocessing.Pool(CPUs)
    queues = [Queue() for c in range(CPUs)]
    kw = {}
    kw['size'] = npart
    kw['high'] = L
    kw['rholist'] = rholist
    kw['xx'] = xx
    kw['yy'] = yy
    kw['delta'] = delta
    kw['sigma'] = sigma
    kw['s'] = s
    kw['mu'] = mu
    kw['pdf_full'] = pdf_full
    kw['R'] = R
    jobs = []
    for q in queues:
        KW = kw.copy()
        KW['q'] = q
        jobs.append(Process(target=fractal_worker,kwargs=KW))

    for job in jobs:
        job.start()
    for qn,q in enumerate(queues):
        print("Getting Queue #",qn)
        Rq = q.get()
        Rsum += Rq
    for job in jobs:
        job.join()
    # for cpu in CPUs:
        # pool.apply(fractal_worker,args=(*listargs))
        # pool.close()
        # pool.join()
        # Rs.append(fractal_worker())
    # pdb.set_trace()
    R = Rsum
    fig,axes = plt.subplots(2,figsize=(12,12))

    ax = axes.flat[0]
    cc = ax.pcolormesh(xx,yy,R,cmap=M.cm.Greys_r,)
    ax.set_aspect('equal')
    ax = axes.flat[1]
    # vals, edges = N.histogram(R,bins=20)
    ax.hist(R.flatten(),100)
    fig.colorbar(cc)
    fig.show()
    fig.savefig('lovejoy_example.png')


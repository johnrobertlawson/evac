import pdb
import copy
import datetime
import os
import glob
import pickle
import itertools

import numpy as N
import matplotlib as M
M.use('Agg')
import matplotlib.pyplot as plt

from evac.stats.detscores import DetScores
from evac.plot.performance import Performance
from evac.stats.forecastvalue import ForecastValue

datadir = '/work/john.lawson/HRRRvNEWSe/stats/'
outdir = '/home/john.lawson/HRRRvNEWSe/pyoutput'

intersects = True

closestcomp = False
FSS = False; save_sigtest = False
Casati = False
rankhist = False
SAL_allmembers = False
salth = 15 #(15,22,30)
compute_abcd = False
plot_abcd = False; allscores = False
compute_fvens = False
plot_FV_matrix = False
plot_performance = True

if intersects == True:
    outdir = '/home/john.lawson/HRRRvNEWSe/pyoutput/intersects/'

models = ("NEWSe","HRRRo","HRRRx")
ensmems = ["ENS_MEM_{0}".format(n) for n in range(1,19)]
fchrs = (1,2,3)

#### FUNCTIONS ####
def proceed_intersect(cf,ix_dts,model=None,countin=None,intersects=True,return_model=False,
                        return_membername=False,return_utc=False):
    """
    Returns:
        (whether to proceed, count)
    """

    returns = []
    proceed = True
    utc = None
    member = None

    model_of_file = [m for m in models if os.path.basename(cf).startswith(m)][0]
    # pdb.set_trace()
    if (model != None) and (model != model_of_file):
        proceed = False
    else:
        if (model_of_file == 'NEWSe'):
            member = return_member(cf)
        else:
            member = None

        if intersects:
        # checks the case/time is common to all three models

            utc = time_from_name(cf,model=model_of_file,prod=member)
            if utc not in ix_dts:
                proceed = False

    returns.append(proceed)
    if countin is not None:
        if proceed:
            if isinstance(countin,dict):
                # countout = countin[model_of_file]
                countin[model_of_file] += 1
            else:
                # countout = countin
                countin += 1

            # countout += 1

        # countin/countout are connected
        returns.append(countin)

    # if model_of_file != 'NEWSe':
        # pdb.set_trace()
    if return_model:
        returns.append(model_of_file)

    if return_utc:
        returns.append(utc)

    if return_membername:
        returns.append(member)

    return returns

def find_intersects(ls,fmt='dict'):
    # Find times, init times, and fcst times common to all three
    lsf = [os.path.basename(l) for l in ls]

    # Store date/time in lists
    TIMES = {m:{} for m in models}
    for l in lsf:
        model,case,itime,*_ = l.split('_')
        case = case[:8]
        itime = itime[:4]
        if case not in TIMES[model].keys():
            TIMES[model][case] = set()
        # if itime not in TIMES[model][case]:
        TIMES[model][case].add(itime)

    cases = TIMES[models[0]].keys()
    for model in models:
        assert TIMES[model].keys() == cases

    # DTS[case] = list of mutual init times
    DTS = {}
    for case in cases:
        DTS[case] = set.intersection(*[TIMES[m][case] for m in models])

    if fmt == 'dict':
        dts = DTS
    elif fmt == 'datetime':
        dtfmt = "%Y%m%d%H%M"
        dts = []
        for case in cases:
            for i in DTS[case]:
                st = ''.join((case,i))
                dt = datetime.datetime.strptime(st,dtfmt)
                dts.append(dt)

    return dts

def savepickle(ob,fname,fdir='/work/john.lawson/HRRRvNEWSe/stats/'):
    f = os.path.join(fdir,fname)
    with open(f,'wb') as fo:
        pickle.dump(ob,fo)
    return

def loadpickle(f):
    with open(f,'rb') as fo:
        data = pickle.load(fo)
        return data

def generate_fileloop(fs):
    for f in fs:
        yield f

def generate_fssloop(threshs,neighs,fchrs,ndicts):
    for th in threshs:
        for ne in neighs:
            for fc in fchrs:
                for nd in range(ndicts):
                    yield th,ne,fc,nd

def time_from_name(fp,model=None,prod='mean'):
    f, *_ = os.path.basename(fp).split('.')
    model, *_ = f.split('_')
    if model == 'NEWSe':
        if 'saldata' in f:
            fmt = "{}_%Y%m%d_%H%M_saldata_{}".format(model,prod)
        else:
            fmt = "{}_%Y%m%d_%H%M_{}".format(model,prod)
    else:
        if 'saldata' in f:
            fmt = "{}_%Y%m%d_%H%M_saldata".format(model)
        else:
            fmt = "{}_%Y%m%d_%H%M".format(model)


    t = datetime.datetime.strptime(f,fmt)
    return t


def return_closest(plot_hist=False,return_dict=True):
    cc_dir = os.path.join(datadir,'closest')
    cc_ls = glob.glob(os.path.join(cc_dir,'*'))
    cc_outdir = outdir

    cc1 = plot_hist
    cc2 = return_dict

    if cc1:
        stats = N.zeros([4,18]) # four times, 18 members
        for c in cc_ls:
            lst = loadpickle(c)
            if len(lst) != 4:
                continue
            for tidx in range(4):
                if lst[tidx] == None:
                    continue
                mem = lst[tidx]
                ensn =  ensmems.index(mem)
                stats[tidx,ensn] += 1

        for tidx in range(4):
            fig,ax = plt.subplots(1)
            fname = "closest_members_{}h".format(tidx)
            ax.bar(range(18),stats[tidx,:])
            ax.set_xticks(range(18))
            ax.set_xticklabels(["M{:02d}".format(e) for e in range(1,19)])
            ax.set_title("Hour {}".format(tidx))
            outfpath = os.path.join(outdir,fname)
            fig.savefig(outfpath)
            print("Saved figure {}".format(outfpath))

    if cc2: 
        CM = {}
        for c in cc_ls:
            t = time_from_name(c)
            CM[t] = {}
            lst = loadpickle(c)
            for tidx in range(len(lst)):
                if lst[tidx] == None:
                    CM[t][tidx] = None
                else:
                    CM[t][tidx] = lst[tidx]

    return CM

def return_member(fname):
    fname = os.path.basename(fname)
    found_idx = False
    for substr in ("ENS","mean","probs"):
        try:
            idx = fname.index(substr)
        except:
            pass
        else:
            found_idx = True
            break
    if not found_idx:
        raise Exception
    return fname.split('.')[0][idx:]

#### PROCEDURE ####
if closestcomp:
    return_closest()

if compute_abcd:
    C_dir = os.path.join(datadir,'ABCD')
    C_ls = glob.glob(os.path.join(C_dir,'*'))
    C_outdir = os.path.join(outdir,"ABCD")
    cs = ('a','b','c','d')
    threshs = sorted(loadpickle(C_ls[0])[0].keys())
    hrs = range(3)
    models = ['NEWSe','HRRRx','HRRRo']
    counts = {m:0 for m in models}
    modelsplus = models + ['NEWSe_closest']
    scores = ['POD','TS','HR','FAR','PFD','BIAS','KSS','CSI',
                'PCLI','F','PC','HSS','PSS','GSS','Q','SR']
    allscores = scores + ['FV',]
    scoredata = {sc: {m: {h: {th: [] for th in threshs} for h in hrs} 
                for m in modelsplus} for sc in allscores}

    CM = return_closest()
    # pdb.set_trace()
    print("Starting loading and computation")
    ix_dts = find_intersects(C_ls,fmt='datetime')
    for cf in C_ls:
        proceed,counts,model,utc,member = proceed_intersect(cf,ix_dts,countin=counts,intersects=intersects,return_model=True,return_utc=True,return_membername=True)
        if proceed:
            if intersects:
                print("Computing for {}".format(os.path.basename(cf)))
            lst = loadpickle(cf)
            # List of times, each with dictionary of 2x2
            if (len(lst) != 3) or (all(lst) == False):
                # Skip all but the full runs of NEWS-e
                continue
            for tidx in hrs:
                closest_switch = False
                if (model == 'NEWSe'):
                    if (CM[utc][tidx] == member):
                        closest_switch = True
                for th in threshs:
                    KW = {c:lst[tidx][th][cix] for cix,c in enumerate(cs)}

                    # FV is forecast value, for 2x2 cont. matrix
                    # Default C/L ratios are 0.01 to 1.00 every 0.01
                    FV = ForecastValue(**KW)
                    scoredata['FV'][model][tidx][th].append(FV.FVs)
                    if closest_switch:
                        scoredata['FV']['NEWSe_closest'][tidx][th].append(FV.FVs)

                    # DS is deterministic scores, like POD
                    DS = DetScores(**KW)
                    for s in scores:
                        scoredata[s][model][tidx][th].append(DS.get(s))
                        if closest_switch:
                            scoredata[s]['NEWSe_closest'][tidx][th].append(DS.get(s))

    
    for model in models:
        print("{} has {} files.".format(model,counts[model]))
    if intersects:
        picklename = 'detscore_intersects.pickle'
    else:
        picklename = 'detscore_dict.pickle'

    # pdb.set_trace()

    print("About to save pickle.")
    savepickle(scoredata,picklename)

if compute_fvens:
    # TODO: Logic to only pick FVens files for intersecting times
    raise Exception
    R_dir = os.path.join(datadir,'ROC')
    R_ls = glob.glob(os.path.join(R_dir,'*'))
    R_outdir = os.path.join(outdir,"ABCD") # Save with other FV stuff
    cs = ('a','b','c','d')
    ix_dts = find_intersects(R_ls,fmt='datetime')
    count = 0
    R_ls2 = []
    for r in R_ls:
        assert os.path.basename(r).startswith("NEWSe")
        proceed,count = proceed_intersect(r,ix_dts,countin=count,intersects=intersects,model='NEWSe')
        if proceed:
            R_ls2.append(r)
            rr = loadpickle(r)
            if rr[0] != None:
                threshs = sorted(rr[0].keys())
                ws = N.linspace(0,1,len(rr[0][1][0]))
                break
    print("Found {} NEWSe times.".format(count))
    hrs = range(3)
    model = 'NEWSe'

    FVens_data = {h: {th: {w: [] for w in ws} for th in threshs} for h in hrs} 

    print("Starting loading and computation")
    for rf in R_ls2:
        lst = loadpickle(rf)
        # List of times, each with dictionary of 2x2
        if (len(lst) != 3) or (all(lst) == False):
            # Skip all but the full runs of NEWS-e
            continue
        for tidx in hrs:
            for th in threshs:
                for widx,w in enumerate(ws):
                    if lst[tidx] == None:
                        continue
                    # print('== {} h ======= {} th ======== {} w =='.format(tidx,th,w))
                    KW = {c:lst[tidx][th][cix][widx] for cix,c in enumerate(cs)}

                    # FV is forecast value, for 2x2 cont. matrix
                    # Default C/L ratios are 0.01 to 1.00 every 0.01
                    FV = ForecastValue(**KW)
                    fvarr = N.array(FV.FVs)
                    fvarr[fvarr < 0.0] = 0.0
                    FVens_data[tidx][th][w].append(fvarr)
                    pdb.set_trace()

    print("About to save pickle.")
    if intersects:
        savepickle(FVens_data,'FVens_intersects.pickle')
    else:
        savepickle(FVens_data,'FVens_dict.pickle')

if plot_FV_matrix:
    if intersects:
        scoredata = loadpickle(os.path.join(datadir,'detscore_intersects.pickle'))
    else:
        scoredata = loadpickle(os.path.join(datadir,'detscore_dict.pickle'))

    threshs = sorted(scoredata['HR']['NEWSe'][1].keys())
    hrs = (1,2,3)
    nthresh = len(threshs)

    # X-axis of cost-loss ratio
    cls = N.arange(0.01,1.01,0.01)
    # cls_lab = N.arange(0.0,1.05,0.05)
    cls_lab = [0.01,0.1,0.5,1.0]

    models = ('NEWSe','HRRRx','HRRRo')
    # Data for indiv. members
    ncases = {}
    for model in models:
        ncases[model] = {}
        for hidx,hr in enumerate(hrs):
            ncases[model][hidx] = len(scoredata['FV'][model][hidx][threshs[0]])
    ncls = len(cls)
    FVdata = N.zeros([3,3,ncls,nthresh])
    clvs = N.arange(-0.15,0.16,0.01)

    for hidx, hr in enumerate(hrs):
        for thidx, th in enumerate(threshs):
            for nm, m in enumerate(models):
                ncs = ncases[m][hidx]
                FVall = N.zeros([ncs,ncls])
                #for c in range(ncs):
                FVall[:,:] = scoredata['FV'][m][hidx][th]
                FVall[FVall<0] = 0
                FVdata[hidx,nm,:,thidx] = N.mean(FVall,axis=0)

    for i,m in enumerate(('NEWSe','HRRRx')):
        fig,axes = plt.subplots(nrows=3,ncols=1,figsize=(8,8))
        for hidx,hr in enumerate(hrs):
            ax = axes.flat[hidx]
            ax.set_xscale('log')
            if m == 'NEWSe':
                cm = M.cm.seismic
            else:
                cm = M.cm.PRGn_r
            fname = 'FV_matrix_{0}.png'.format(m)
            #fname = 'FV_matrix_{0}_{1}h.png'.format(m,hr)
            FVplot = N.transpose(FVdata[hidx,i,:,:] - FVdata[hidx,i+1,:,:])
            cb = ax.contourf(cls,N.arange(nthresh),FVplot,levels=clvs,cmap=cm)
            #pdb.set_trace()
            #cb = ax.pcolormesh(FVplot,vmin=clvs.min(),vmax=clvs.max(),cmap=cm)
            # ax.set_xticks([0.01,0.1,0.5,1.0])
            ax.set_xticks(cls_lab)
            ax.set_xticklabels(cls_lab)
            #ax.get_xaxis().set_major_formatter(M.ticker.ScalarFormatter())
            #ax.set_xticklabels(cls_lab)
            #ax.set_ylim([0,0.6])
            ax.set_yticks(N.arange(nthresh)[::2])
            ax.set_yticklabels(threshs[::2])
            ax.set_ylabel("Threshold (mm/hr)")
        ax.set_xlabel("C/L ratio")
        outpath = os.path.join(outdir,fname)
        fig.tight_layout()

        fig.subplots_adjust(bottom=0.17)
        cbar_ax = fig.add_axes([0.18,0.06,0.7,0.025])
        cbx = plt.colorbar(cb,cax=cbar_ax,orientation='horizontal')
        cbx.set_label("Difference in FV")

        fig.savefig(outpath)
        plt.close(fig)


if plot_abcd:
    hrs = range(3)
    models = ['NEWSe','HRRRx','HRRRo']
    modelsplus = models + ['NEWSe_closest']
    C_outdir = os.path.join(outdir,"ABCD")

    if allscores:
        scores = ['POD','TS','HR','FAR','PFD','BIAS','KSS','CSI',
                    'PCLI','F','PC','HSS','PSS','GSS','Q']
    else:
        scores = ['CSI',]

    if intersects:
        scoredata = loadpickle(os.path.join(datadir,'detscore_intersects.pickle'))
    else:
        scoredata = loadpickle(os.path.join(datadir,'detscore_dict.pickle'))
    threshs = sorted(scoredata['HR']['NEWSe'][1].keys())
    colors = {'NEWSe':'red','HRRRo':'green','HRRRx':'blue','NEWSe_closest':'red'}
    lines = {'NEWSe':'-','HRRRo':'-','HRRRx':'-','NEWSe_closest':'--'}

    # POD, FAR, etc
    for score in scores:
        print("Plotting all figures for {}".format(score))
        for hidx in hrs:
            fname = "{}_{}h.png".format(score,hidx+1)
            outfpath = os.path.join(C_outdir,fname)
            if os.path.isfile(outfpath):
                print("Already plotted.")
                continue
            # fig,axes = plt.subplots(3,figsize=(10,12))
            fig,ax = plt.subplots()
            # for model in models:
            for model in modelsplus:
                ydata = []
                for th in threshs:
                    ydata.append(N.mean(scoredata[score][model][hidx][th][:]))
                ax.plot(range(len(threshs)),ydata,label=model,color=colors[model],linestyle=lines[model])
                #ax.plot(threshs,ydata,label=model,color=colors[model],linestyle=lines[model])
            ax.set_xlabel("Threshold (mm/hr)")
            ax.set_ylabel("Score")
            # ax.set_xticks(threshs)
            ax.set_xticks(range(len(threshs))[::2])
            ax.set_xticklabels(threshs[::2])
            #ax.set_title("{}_{}h.png".format(score,hidx+1))
            ax.legend()
            fig.tight_layout()
            fig.savefig(outfpath)
            plt.close(fig)


    # FV for individual members and ensemble members

    # FV for ensemble probs
    if intersects:
        FVens = loadpickle(os.path.join(datadir,'FVens_intersects.pickle'))
    else:
        FVens = loadpickle(os.path.join(datadir,'FVens_dict.pickle'))
    ws = sorted(FVens[0][1].keys())
    nws = len(ws)
    #pdb.set_trace()

    # X-axis of cost-loss ratio
    cls = N.arange(0.01,1.01,0.01)
    # cls_lab = N.arange(0.0,1.05,0.05)
    cls_lab = [0.01,0.1,0.5,1.0]

    # Data for indiv. members
    ncases = {}
    for model in modelsplus:
        ncases[model] = {}
        for hidx in hrs:
            ncases[model][hidx] = len(scoredata['FV'][model][hidx][threshs[0]])
    ncls = len(cls)
    assert ncls == len(scoredata['FV'][models[0]][0][threshs[0]][0])

    method = 3
    print("Now plotting FV.")
    for hidx in hrs:
        for th in threshs:
            fname = "FV_{}h_{}mm.png".format(hidx+1,th,method)
            # fname = "FV_{}h_{}mm_method{}.png".format(hidx+1,th,method)
            outfpath = os.path.join(C_outdir,fname)
            if os.path.isfile(outfpath):
                print("Already plotted.")
                continue
            # fig,axes = plt.subplots(3,figsize=(10,12))
            print("Plotting FV to {}".format(outfpath))
            fig,ax = plt.subplots()
            ax.set_xscale('log')

            modelsplot = ['NEWSe determ.' if m.startswith('NEWSe') else m for m 
                                in models] + ['NEWSe_closest',]
                                #in models] + ['NEWSe probab','NEWSe_closest']
            for model in modelsplot:
                if model.endswith('probab'):
                    continue
                    # FVens is function of C/L
                    # for each w, find max FV along C/L axis. or for each c/l, find max FV along w axis?
                    # FVens_data[0][1][1.0][0]   ---> 100, all C/L vals
                    # FVens_data[0][1][1.0]      ---> 104, all cases
                    # (104,100)
                    
                    # FVens_data[tidx][th][w].append(FV.FVs)
                    FVarr = N.zeros([100,len(FVens[hidx][th][0.0]),nws])
                    for widx,w in enumerate(ws):
                        FVarr[:,:,widx] = N.array(FVens[hidx][th][w]).T
                        # FVarr[:,:,w][FVarr[:,:,w] < 0.0] = 0.0
                    # (100,104)
                    # if th == 1:
                        # pdb.set_trace()
                    if method in (1,2):
                        # Max along the probability axis 
                        FVarr2 = N.max(FVarr,axis=2)
                        # Average max along all cases
                        FVarr2 = N.mean(FVarr2,axis=1)
                        ax.plot(cls,FVarr2,label=model,color='red')
                        for widx,w in enumerate(ws):
                            ax.plot(cls,N.mean(FVarr[:,:,widx],axis=1),lw=0.2,color='red')
                    elif method == 8:
                        FVarr11 = N.max(N.mean(FVarr,axis=0),axis=1)
                        ax.plot(cls,FVarr11,label=model,color='red')
                    elif method == 3:
                        FVarr5 = N.mean(FVarr,axis=1)
                        FVarr6 = N.max(FVarr5,axis=1)
                        ax.plot(cls,FVarr6,label=model,color='black')
                        for widx,w in enumerate(ws):
                            ax.plot(cls,N.mean(FVarr[:,:,widx],axis=1),lw=0.2,color='black')
                    elif method == 4:
                        FVarr8 = N.mean(N.max(FVarr,axis=2),axis=1)
                        ax.plot(cls,FVarr8,label=model,color='red')
                        for widx,w in enumerate(ws):
                            ax.plot(cls,N.mean(FVarr[:,:,widx],axis=1),lw=0.2,color='red')

                    if method == 2:
                        # FVarr3 = N.mean(FVarr,axis=1)
                        # FVarr = N.max(FVarr,axis=1)
                        for w in ws:
                            arr = N.zeros([100,len(FVens[hidx][th][w])])
                            for n in range(len(FVens[hidx][th][w])):
                                arr[:,n] = FVens[hidx][th][w][n]
                            # arr[arr < 0.0] = 0.0
                            # pdb.set_trace()
                            ax.plot(cls,N.mean(arr,axis=1),lw=0.4,color=colors[model],
                                                linestyle=lines[model])

                # elif model.endswith('closest'):
                    # pass
                else:
                    if model.endswith('closest'):
                        pass
                    elif model.startswith('NEWSe'):
                        model = 'NEWSe'
                    ydata_all = N.zeros([ncases[model][hidx],ncls])
                    for cl in cls:
                        ydata_all[:,:] = scoredata['FV'][model][hidx][th]
                    ydata_all[ydata_all < 0.0] = 0.0
                    ydata = N.mean(ydata_all,axis=0)
                    ax.plot(cls,ydata,label=model,linestyle=lines[model],color=colors[model])
            ax.set_xlabel("C/L ratio")
            ax.set_xticks(cls_lab)
            # ax.get_xaxis().set_major_formatter(M.ticker.ScalarFormatter())
            ax.set_xticklabels(cls_lab)
            ax.set_ylim([0,0.4])
            ax.set_ylabel("Forecast value")
            # ax.set_title("{}h_{}mm.png".format(hidx+1,th))
            ax.legend()
            fig.tight_layout()
            fig.savefig(outfpath)
            plt.close(fig)
            

    


if Casati:
    #assert True == False
    # Need to fix intersects logic.

    C_dir = os.path.join(datadir,'CASATI')
    C_ls = glob.glob(os.path.join(C_dir,'*'))
    C_outdir = os.path.join(outdir,"casati")
    scores = ('MSE','SS')
    models = ('NEWSe','HRRRx','HRRRo')
    hrs = range(1,4)

    thresholds = sorted(loadpickle(C_ls[0])[0][0].keys())
    Ls = sorted(loadpickle(C_ls[0])[0][0][1].keys()) # assuming we have a 1 mm/hr thresh!
    # DATA = {score:{th:{hr:[] for hr in hrs} 
                # for th in thresholds} for score in scores}

    # Loop over each model
    FM = {model:[] for model in models}
    [FM[m].append(loadpickle(f)) for f,m in itertools.product(C_ls, models)
                            if os.path.basename(f).startswith(m)]
    avedata = {model:{score:{} for score in scores} for model in models}
    for model in models:
        DATA = {f:{score:{hr:N.zeros([len(thresholds),len(Ls)]) for hr in hrs} 
                            for score in scores} for f in range(len(FM[model]))}

        print("Starting to collate data for {}.".format(model))
        # Count number of pickle files with 3 times
        f3 = []
        for f,data in enumerate(FM[model]):
#def proceed_intersect(cf,ix_dts,model=None,countin):
            # lst[number of times][(MSE,SS)][thresh][L]
            # data = loadpickle(f)
            ntimes = len(data)
            if ntimes == 3:
                f3.append(f)
                for (thidx,th), hr,(sidx,score),(Lidx,L) in itertools.product(
                        enumerate(thresholds),hrs,enumerate(scores),enumerate(Ls)):
                    hridx = hr-1
                    DATA[f][score][hr][thidx,Lidx] = data[hridx][sidx][th][L]
            else:
                del DATA[f]

        print("Done - moving to plotting.")
        # Plot
        xx = N.arange(len(thresholds))
        yy = N.arange(len(Ls))
        xxtx = xx
        yytx = yy
        xxtl = thresholds
        yytl = Ls
        CMAPS = {'MSE':M.cm.plasma,'SS':M.cm.inferno}
        # CMAPS = {'MSE':M.cm.inferno,'SS':M.cm.inferno_r}
        LEVELS = {'MSE':N.arange(0,0.13,0.01),'SS':N.arange(-5,1.5,0.5)}
        # levels = ( N.arange(0.01,MSEarr.max()+0.01,0.01), N.arange(N.floor(SSarr.min()),N.ceil(SSarr.max()),1),) 

        for score,hr in itertools.product(scores,hrs):
            # Mean data
            alldata = N.zeros([len(f3),len(thresholds),len(Ls)])
            for fidx,f in enumerate(f3):
                if len(DATA[f][score].keys())==3:
                    alldata[fidx,:,:] = DATA[f][score][hr]
            alldata[alldata == N.inf] = N.nan
            alldata[alldata == -N.inf] = N.nan
            alldata[alldata < -5] = -5
            # if score == 'SS' and model == 'HRRRo' and hr == 3:
                # pdb.set_trace()
            avedata[model][score][hr] = N.nanmean(alldata,axis=0)

    skip_mse = 1
    for score,hr in itertools.product(scores,hrs):
        if skip_mse:
            continue
        # Figure time
        fig,axes = plt.subplots(ncols=1,nrows=3,figsize=(10,10))
        fname = '{}_{}.png'.format(score,hr)
        print("Plotting {}".format(fname))
        cmap = CMAPS[score]

        for ax,model in zip(axes.flat,models):
            
            title = '{} {} score at {} h (hourly QPF)'.format(model,score,hr)
            # title = '{} {} score at {} h (hourly QPF)'.format(model,score,hr)
            cf = ax.contourf(xx,yy,avedata[model][score][hr].T,cmap=cmap,levels=LEVELS[score])
            ax.set_title(title)
            ax.set_xticks(xxtx)
            ax.set_xticklabels([str(x) for x in xxtl])
            ax.set_yticks(yytx)
            ax.set_yticklabels([str(2**y) for y in yytl])
            ax.set_xlabel("Thresholds (mm/hr)")
            ax.set_ylabel("Spatial scale (km)")
            # plt.colorbar(cf)
            fig.tight_layout()
            fig.savefig(os.path.join(C_outdir,fname))
        
    # Plot differences
    # NEWSe - HRRRx
    lablist = ('NEWSe','HRRRx')
    levels = N.arange(-0.016,0.017,0.001)
    for hr in hrs:
        diff1 = avedata['NEWSe']['MSE'][hr].T - avedata['HRRRx']['MSE'][hr].T
        diff2 = avedata['HRRRx']['MSE'][hr].T - avedata['HRRRo']['MSE'][hr].T

        for n,diff in enumerate((diff1, diff2)):
            fig,ax = plt.subplots(1,figsize=(8,5))
            fname = 'MSE_valueadded_{}_{}h.png'.format(lablist[n],hr)
            print("Plotting {}".format(fname))

            # Color bars opposite to FV_matrix because lower error is better
            if n == 0:
                cm = M.cm.seismic_r
            else:
                cm = M.cm.seismic_r
                #cm = M.cm.PRGn
            
            cf = ax.contourf(xx,yy,diff,cmap=cm,levels=levels)
            ax.set_xticks(xxtx)
            ax.set_xticklabels([str(x) for x in xxtl])
            ax.set_yticks(yytx)
            ax.set_yticklabels([str(2**y) for y in yytl])
            ax.set_xlabel("Thresholds (mm/hr)")
            ax.set_ylabel("Spatial scale (km)")
            plt.colorbar(cf,orientation='horizontal')
            fig.tight_layout()
            fig.savefig(os.path.join(C_outdir,fname))
            print("Saved",os.path.join(C_outdir,fname))

    

if FSS:
    FSSdir = os.path.join(datadir,'FSS')
    FSSls = glob.glob(os.path.join(FSSdir,'*'))
    ix_dts = find_intersects(FSSls,fmt='datetime')
    # 0-1h, 1-2h, 2-3h, 0.5-1.5h
    # Maybe re-run half-hour stats so it's 0h-1h?
    DICTS = {}
    MEANS = {}
    missed = {} 
    # [times, models (and members), threshs,neighs, cases]
    SIGTEST = N.ones([3,20,7,62,51]) + 98
    for _T, fidx in enumerate(range(3)):
        MEANS[fidx] = {}
        missed[fidx] = 0
        fig,axes = plt.subplots(3,figsize=(10,12))
        counts = {m:0 for m in models}
        for m,ax in zip(models,axes.flat):
            MEANS[fidx][m] = {}
            ax.set_title(m)
            ax.set_ylim([0,1])

            
            # f_iter = generate_fileloop(FSSls)
            # dict_all = [loadpickle(f) for f in f_iter if os.path.basename(f).startswith(m)]
#def proceed_intersect(cf,ix_dts,model=None,countin):
            dict_all = []
            for f in FSSls:
                proceed, counts = proceed_intersect(f,ix_dts,model=m,countin=counts,intersects=intersects)
                if proceed:
                    this_pickle = loadpickle(f)
                    dict_all.append(this_pickle)

                    if save_sigtest:
                        if m == 'NEWSe':
                            mem_name = return_member(f)
                            NMODEL = int(mem_name.split('_')[-1]) - 1
                        else:
                            if m == 'HRRRo':
                                NMODEL = 18
                            elif m == 'HRRRx':
                                NMODEL = 19
                            mem_name = m
                        # [times, models (and members), threshs, cases]
                        this_case = time_from_name(f,prod=mem_name)
                        NCASE = ix_dts.index(this_case)
                        for NTHRESH, th in enumerate(sorted(this_pickle[0].keys())):
                            for NNEIGH,neigh in enumerate(sorted(this_pickle[0][1].keys())):
                                for NTIME in range(3):
                                    if this_pickle[NTIME] is None:
                                        val = N.nan
                                    else:
                                        val = this_pickle[NTIME][th][neigh]
                                    SIGTEST[NTIME,NMODEL,NTHRESH,NNEIGH,NCASE] = val

            #print(dict_all.__len__())
            #pdb.set_trace()
            DICTS[m] = [d for d in dict_all if len(d) == 3]
            print("{0}files found for {1}".format(len(DICTS[m]),m))
            thresholds = sorted(DICTS[m][0][0].keys())
            neighs = sorted(DICTS[m][0][0][1].keys()) # assuming a 1 mm threshold exists!
            for th in thresholds:
                data = N.zeros([len(neighs),len(DICTS[m])])
                for n,dicn in itertools.product(range(len(neighs)),range(len(DICTS[m]))):
                    try:
                        data[n,dicn] = DICTS[m][dicn][fidx][th][neighs[n]]
                    except TypeError:
                        missed[fidx] += 1
                        data[n,dicn] = N.nan
                MEANS[fidx][m][th] = N.nanmean(data,axis=1)
                # std_th = N.std(data,axis=1)
                ax.plot(neighs,MEANS[fidx][m][th],label="Thresh = {}".format(th))
                # if m == 'NEWSe': pdb.set_trace()
                #pdb.set_trace()
            ax.grid('on')
        axes.flat[0].legend(ncol=3)
        fname = 'FSS_score_all_{}h.png'.format(fidx+1)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir,fname))
        print("{} missed times.".format(missed[fidx]))
        print("{} total items loaded.".format(counts[m]))

        if save_sigtest:
            f_out = os.path.join(datadir,'FSS_Sigtest.npy')
            N.save(f_out,SIGTEST)
            assert 1==0
            save_sigtest = False
            # This has been done, don't do it again.

    # Loop over neighbourhoods and plot change with time for three models, for a few thresholds.
    # [models, times, neighbourhoods, thresholds]
    print("Creating time-evolution array.")

    PLOT = {'NEWSe':'red','HRRRo':'green','HRRRx':'blue',
            0.5:'-',4:'--',16:'-.'}

    neighpick = neighs[:10:2]
    npidxs = N.searchsorted(neighs,neighpick)
    nnp = len(neighpick)
    threshpick = (0.5,4,16)
    ntp = len(threshpick)
    comparearr = N.zeros([3,3,nnp,ntp])
    for fidx, (midx,m), (tpidx,tp) in itertools.product(range(3),enumerate(models),enumerate(threshpick)):
        comparearr[fidx,midx,:,tpidx] = MEANS[fidx][m][tp][npidxs]
    fig,axes = plt.subplots(1,nnp,figsize=(8,7))
    fname = 'FSSmean_time_ev.png'
    for (midx,m),(nidx,n),(tidx,t) in itertools.product(enumerate(models),enumerate(neighpick),enumerate(threshpick)):
        ax = axes.flat[nidx]
        if nidx == 0:
            newlabs = {'HRRRo':"HRRRv1",'HRRRx':"HRRRv2",'NEWSe':"NEWSe"}
            lab = '{}\n{} mm/h\n({:.2f} in/h)'.format(newlabs[m],t,(t*0.03937))
        else:
            lab = ''
        ax.plot(comparearr[:,midx,nidx,tidx],color=PLOT[m],linestyle=PLOT[t],label=lab)
        # Place annotation at the middle model's line
        # mmid = N.where(comparearr[-1,:,nidx,tidx] == sorted(comparearr[-1,:,nidx,tidx])[1])[0][0]
        if midx == 0:
        # if midx == mmid:
            # ax.annotate(str(n),xy=(2.03,comparearr[2,midx,nidx,tidx]))
            ax.set_xlim([0,2])
            ax.set_xticks(range(3))
            ax.set_xticklabels([str(n+1) for n in range(3)])
            ax.set_xlabel("Time")
            #ax.set_title("Neighborhood\n(sq. length) = {}".format(n))
    axes[0].set_ylabel("FSS score")
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.31)
    # -0.225    
    axes[0].legend(bbox_to_anchor=(3.35,-0.43),loc=8,borderaxespad=0.05,ncol=3)
    fpathout = os.path.join(outdir,fname)
    #pdb.set_trace()
    fig.savefig(fpathout)
    print("Saved FSS figure to ",fpathout)
    
if rankhist:
    # 20 samples per fcst time x 3 fcst times x 12 cases x ~10 init times ~ 6500
    # Maybe increase samples if assuming independence
    ls = glob.glob(os.path.join(datadir,'RANKHIST','*'))
    arr_total = N.zeros([3,19],dtype=N.int64)
    arr_halfhour = N.zeros([1,19],dtype=N.int64)
    count = 0
    ix_dts = find_intersects(ls,fmt='datetime')
    for f in ls:
        if os.path.basename(f).startswith('NEWSe'):
            proceed, count = proceed_intersect(f,ix_dts,countin=count,intersects=intersects)
            if proceed:
                print("Loading ",f)
                arr_f = N.load(f)
                assert arr_f.shape[1] == 19
                if arr_f.shape[0] == 3:
                    arr_total += arr_f
                    # print(arr_total)
                    # pdb.set_trace()
                elif arr_f.shape[0] == 1:
                    arr_halfhour += arr_f
                else:
                    raise Exception
                assert N.all(arr_f > -1)
    print("Counted {} files".format(count))
    for hridx in range(5):
        fig,ax = plt.subplots()
        if hridx == 4:
            ax.set_title("NEWS-e rank histogram for all hourly QPF periods")
            fname = 'NEWSe_hist_all.png'
            arr_all = N.sum(arr_total,axis=0) + arr_halfhour
            ax.bar(range(19),arr_all[0,:])
        elif hridx == 3:
            ax.set_title("NEWS-e rank histogram for 0.5-1.5h QPF (half-hour init)")
            fname = 'NEWSe_hist_halfhr.png'
            ax.bar(range(19),arr_halfhour[0,:])
        else:
            ax.set_title("NEWS-e rank histogram for {}-{}h QPF".format(hridx,hridx+1))
            fname = 'NEWSe_hist_{0}h.png'.format(hridx+1)
            ax.bar(range(19),arr_total[hridx,:])
        ax.set_xticks(N.arange(19))
        ax.set_xticklabels(N.arange(19))
        fpath = os.path.join(outdir,fname)
        fig.savefig(fpath)

if SAL_allmembers:
    """ TODO -  Highlight the closest-to-mean NEWSe members.
             -  Plot the mean, prob-matched mean.
    """
    salstr = 'SAL_{}'.format(salth)
    ls = glob.glob(os.path.join(datadir,salstr,'*'))
    models = ('NEWSe','HRRRx','HRRRo')
    comps = ('SS','AA','LL')
    colors = ('red','blue','green')
    times = list(range(1,4))
    
    # Dictionary for all data
    SALDICT = {}
    # Dictionary for ignored file/time stats
    IGNORED = {}
    # Dictionary for investigating outliers    
    WEIRD = {}

    counts = {m:0 for m in models}
    print("Initialising dictionaries.")
    for model,color in zip(models,colors):
        WEIRD[model] = {}
        IGNORED[model] = {k:0 for k in comps}
        IGNORED[model]['files'] = 0
        SALDICT[model] = {k:{} for k in comps}
        SALDICT[model]['col'] = color
        SALDICT[model]['filecount'] = 0
        for comp in comps:
            SALDICT[model][comp] = {'all':[],'plotall':[],
                                    1:[],2:[],3:[],}
            IGNORED[model][comp] = {1:0,2:0,3:0}
            WEIRD[model][comp] = {1:[],2:[],3:[]}

    print("Loading data.")
    ix_dts = find_intersects(ls,fmt='datetime')
    for nf, f in enumerate(ls): 
#def proceed_intersect(cf,ix_dts,model=None,countin):
        proceed,counts = proceed_intersect(f,ix_dts,countin=counts,intersects=intersects)
        if proceed:
            datafname = os.path.basename(f)
            for k in models:
                if datafname.startswith(k):
                    model = k
                    SALDICT[model]['filecount'] += 1
            salarr = N.load(f)

            # Sanity check for SAL components
            assert salarr.shape[1] == 3
            # Only include 3-h forecasts
            if salarr.shape[0] != 3:
                IGNORED[model]['files'] += 1
                continue

            for tidx,t in enumerate(times):
                plotthis = 1
                for compidx,comp in enumerate(comps):
                    s = salarr[tidx,compidx]
                    SALDICT[model][comp]['all'].append(s)
                    if (s<-1.8) or (s>1.8):
                        WEIRD[model][comp][t].append(datafname)
                    if (s<-1.99) or (s>1.99) or (s==0.0) or (s==N.nan):
                        plotthis = 0
                        IGNORED[model][comp][t] += 1
                    if plotthis:
                        SALDICT[model][comp]['plotall'].append(s) 
                        SALDICT[model][comp][t].append(s)
          
    # Stats
    for model in models:
        print("We ignored {} files for {}".format(IGNORED[model]['files'],model))
        # print("We plotted data from {} files for {}".format(SALDICT[model]['filecount']
                                                                # -IGNORED[model]['files'],model))
        for t in times:
            allignored = sum([IGNORED[model][c][t] for c in comps])
            print("We ignored {} times for {} {}h-{}h forecasts".format(allignored,model,t-1,t))
            # print("We plotted data from {} files for {}".format(len(SALDICT[model][t])
                                                                    # -IGNORED[model][t],model))
            

    # Plot scatter points
    for tidx, t in enumerate(times+['plotall']):
        fig,ax = plt.subplots(figsize=(8,8))
        for model in models:
            SALDICT[model]['label'] = False
            print("Plotting data for {} at time {}.".format(model,t))
            for ss, aa in zip(SALDICT[model]['SS'][t],SALDICT[model]['AA'][t]):
                lab = model if not SALDICT[model]['label'] else ""
                ax.scatter(ss,aa,c=SALDICT[model]['col'],
                            marker='s',s=30,edgecolor='k',linewidth=0.25,
                            zorder=500,alpha=0.8,label=lab)
                SALDICT[model]['label'] = True

        # Plot median/IQR box
        for model in SALDICT.keys():
            # SALDICT[model][comp] = {}
            for comp in comps:
                median = N.median(N.array(SALDICT[model][comp][t]))
                col = SALDICT[model]['col']
                medkwargs = {'color':col,'zorder':700,'linewidth':0.8,'linestyle':':'}
                if comp is 'SS':
                    ax.axvline(median,**medkwargs)
                    lbS = N.percentile(SALDICT[model][comp][t],25)
                    ubS = N.percentile(SALDICT[model][comp][t],75)
                elif comp is 'AA':
                    ax.axhline(median,**medkwargs)
                    lbA = N.percentile(SALDICT[model][comp][t],25)
                    ubA = N.percentile(SALDICT[model][comp][t],75)
            width = ubS-lbS
            height = ubA-lbA
            ax.add_patch(M.patches.Rectangle((lbS,lbA),width,height,
                            facecolor=col,alpha=0.2,linewidth=0.5,
                            zorder=100,edgecolor='k'))

        ax.set_facecolor('lightgrey')
        ax.legend()
        plt.axhline(0, color='k')
        plt.axvline(0, color='k')
        ax.set_xlim([-2,2])
        ax.set_ylim([-2,2])
        ax.set_xlabel('Structural component')
        ax.set_ylabel('Amplitude component')
        ax.set_aspect('equal',adjustable='box')

        textkwargs = {'va':'center','ha':'center'}
        ax.text(1.4,1.4,'Too wet \n Too stratiform',**textkwargs)
        ax.text(-1.4,-1.4,'Too dry \n Too convective',**textkwargs)
        ax.text(-1.4,1.4,'Too wet \n Too convective',**textkwargs)
        ax.text(1.4,-1.4,'Too dry \n Too stratiform',**textkwargs)

        fig.tight_layout()
        tstr = '{}h'.format(t)
        fname = '{}_allmembers_{}.png'.format(salstr,tstr)
        fpath = os.path.join(outdir,fname)
        fig.savefig(fpath)
        plt.close(fig)

    # Histogram of data
    
    for tidx, t in enumerate(times+['all']):
        fig,axes = plt.subplots(3,3,figsize=(10,10))
        for ax,(model,comp) in zip(axes.flatten(),itertools.product(models,comps)):
            ax.hist(SALDICT[model][comp][t],20)
            ax.set_title('{} {}'.format(model,comp))
            if comp == 'LL':
                ax.set_xlim([0,1])
            else:
                ax.set_xlim([-2,2])
                ax.axvline(0,color='k')
            # ax.axvline(medians[comp][model],**medkwargs)
        fname = '{}_hists_{}h.png'.format(salstr,t)
        fig.savefig(os.path.join(outdir,fname))
        fig.tight_layout()
        plt.close(fig)
    for m in models:
        print('{} files found for {}'.format(counts[m],m))

if plot_performance:
    hrs = range(3)
    models = ['NEWSe','HRRRx','HRRRo']
    modelsplus = models + ['NEWSe_closest']
    C_outdir = os.path.join(outdir,"Performance")

    scores = ['POD','FAR']

    if intersects:
        scoredata = loadpickle(os.path.join(datadir,'detscore_intersects.pickle'))
    else:
        scoredata = loadpickle(os.path.join(datadir,'detscore_dict.pickle'))
    threshs = sorted(scoredata['HR']['NEWSe'][1].keys())
    colors = {'NEWSe':'red','HRRRo':'green','HRRRx':'blue',}#'NEWSe_closest':'darkred'}
    shapes = {0.5:'o',4:'D',16:'s'}
    #shapes = {0.25:'D',1:'o',4:'*'}
    # alphas = {1:0.9,2:0.6,3:0.3} 
    # alphas = {1:0.6,3:0.99} 

    # fname = "{}_{}h.png".format(hidx+1)

    lk = {'bbox_to_anchor':(-0.15,1.0),'ncol':1,}#'fontsize':11}

    for hidx,hr in enumerate(N.arange(1,4)):
        fname = "all_performance_{}h.png".format(hr)
        PD = Performance(C_outdir,fname,legendkwargs=lk)
        for model in modelsplus:
            for th in sorted(shapes.keys()):
                if model.endswith('closest'):
                    continue
                pod = N.mean(scoredata['POD'][model][hidx][th][:])
                # pod = scoredata['POD'][model][hidx][th][:]
                far = N.mean(scoredata['FAR'][model][hidx][th][:])
                # far = scoredata['FAR'][model][hidx][th][:]
                # sr = scoredata['SR'][model][hidx][th][:]
                # print(pod,far)
                count = 0
                # for p,f in zip(pod,far):
                #for p,s in zip(pod,sr):
                #    if not count:
                newlabs = {'HRRRo':"HRRRv1",'HRRRx':"HRRRv2",'NEWSe':"NEWSe"}
                label = "{} {}h {}mm/hr".format(newlabs[model],hr,th)
                #        count = 1
                #    else:
                #        label = None

                p = pod
                f = far
                # pldic = {'marker':shapes[th],'c':colors[model],'alpha':alphas[hr],
                pldic = {'marker':shapes[th],'c':colors[model],
                            'edgecolor':'black','s':50}
                PD.plot_data(pod=pod,far=far,plotkwargs=pldic,label=label)
                # PD.plot_data(pod=p,sr=s,plotkwargs=pldic,label=label)
        PD.save()

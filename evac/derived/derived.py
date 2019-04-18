""" For computing derived variables from existing fields.

Todo:
    * Steal a lot of these from WRFOut.
    * There might be some overlap from stats.
    * Resolve the API. Do all functions start with a parent to inherit
        from? What if the user wants to call without the parent (e.g.,
        using arrays of data that already exist)?

For every class, parent is the class from which to "inherit" from.
"""
import os
import pdb
import itertools

import numpy as N
from scipy.interpolate import interp1d, interpn, RegularGridInterpolator

import evac.utils.met_constants as mc
import evac.utils.utils as utils

def cold_pool_strength(parent,X,time,swath_width=100,env=0,dz=0):
    """
    Returns array the same shape as WRF domain.

    Args:
        X   :   cross-section object with given path
                    This path goes front-to-back through a bow
        km  :   width in the line-normal direction
        env :   (x,y) for location to sample environmental dpt
    """

    # Set up slices
    tidx = parent.get_time_idx(time)
    lvidx = False
    lonidx = False
    latidx = False

    # Get wind data
    wind10 = parent.get('wind10',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    T2 = parent.get('T2',tidx,lvidx,lonidx,latidx)[0,0,:,:]

    # This is the 2D plane for calculation data
    coldpooldata = N.zeros(wind10.shape)

    # Compute required C2 fields to save time
    dpt = parent.get('dpt',tidx,lvidx,lonidx,latidx)[0,:,:,:]
    Z = parent.get('Z',tidx,lvidx,lonidx,latidx)[0,:,:,:]
    HGT = parent.get('HGT',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    heights = Z-HGT
    # pdb.set_trace()

    if isinstance(env,collections.Sequence):
        xx_env = N.arange(env[0]-2,env[0]+3)
        yy_env = N.arange(env[1]-2,env[1]+3)
        dpt_env = N.mean(dpt[:,yy_env,xx_env],axis=1)

    # All cross-sections (parallel)

    # xxx = [X.xx,]
    # yyy = [X.yy,]
    X.translate_xs(-swath_width/2)
    for n in range(swath_width):
        print(("Cross section #{0}".format(n)))
    # for xx, yy in zip(xxx, yyy):
        X.translate_xs(1)
        xx = X.xx
        yy = X.yy
        xx = xx.astype(int)
        yy = yy.astype(int)
        wind_slice = wind10[yy,xx]
        T2_slice = T2[yy,xx]
        slice_loc = parent.find_gust_front(wind_slice,T2_slice,X.angle)

        gfx = xx[slice_loc]
        gfy = yy[slice_loc]
        gf_pts = N.intersect1d(N.where(xx==gfx)[0],N.where(yy==gfy)[0])
        gf_pt = gf_pts[int(len(gf_pts)/2.0)]
        xx_cp = xx[:gf_pt]
        yy_cp = yy[:gf_pt]
        # pdb.set_trace()

        # Compute enviromental dpt at each height
        # Average all levels from the location of gust front
        # forwards to the end of the cross-section.
        if not env:
            xx_env = xx[gf_pt+1:]
            yy_env = yy[gf_pt+1:]
            dpt_env = N.mean(dpt[:,yy_env,xx_env],axis=1)
        # pdb.set_trace()

        for x,y in zip(xx_cp,yy_cp):
        #for x,y in zip(xx,yy):
            if dz:
                coldpooldata[y,x] = parent.compute_cpdz(x,y,dpt[:,y,x],heights[:,y,x],dpt_env)
                # pdb.set_trace()
            else:
                coldpooldata[y,x] = N.sqrt(parent.compute_C2(x,y,dpt[:,y,x],heights[:,y,x],dpt_env))

    return coldpooldata

def compute_cpdz(parent,x,y,dpt,heights,dpt_env):
    """
    Cold pool depth

    Args:
        x       :   x location in domain
        y       :   y location in domain
        dpt     :   density potential temperature slice
        heights :   height AGL slice
        dpt_env :   environmental dpt, column
    """

    dz, zidx = parent.cold_pool_depth(dpt,heights,dpt_env)
    # import pdb; pdb.set_trace()
    return dz

def compute_C2(parent,x,y,dpt,heights,dpt_env):
    """
    C^2 as found in James et al. 2006 MWR

    Args:
        x       :   x location in domain
        y       :   y location in domain
        dpt     :   density potential temperature slice
        heights :   height AGL slice
        dpt_env :   environmental dpt, column
    """

    dz, zidx = parent.cold_pool_depth(dpt,heights,dpt_env)
    C2 = -2*mc.g*((dpt[zidx]-dpt_env[zidx])/dpt_env[zidx])*dz
    # print("dpt = {0} ... dpt_env = {1} ... C2 = {2}".format(dpt[0],dpt_env[0],C2))
    return C2

def cold_pool_depth(parent,dpt,heights,dpt_env):
    dz = 0
    # thresh = -2.0
    thresh = -1.0
    for d,z, de in zip(dpt[1:],heights[1:],dpt_env[1:]):
        dptp = d - de
        if dptp > thresh:
            break
        dz = z

    if isinstance(dz,float):
        zidx = N.where(heights==dz)[0]
    else:
        zidx = 0
    # import pdb; pdb.set_trace()
    return dz, zidx

def find_gust_front(parent,wind_slice,T2_slice,angle,method=3):
    """
    Find location of maximum shear in the horizontal wind along a
    1D slice.

    Args:
        wind_slice      :   1D numpy array
        T2_slice        :   temp 2m slice
        angle           :   angle of slice cross-section
        method          :   way to locate gust front
    """

    shp = wind_slice.shape

    # Compute gradient quantities
    shear = N.zeros(shp)
    T2grad = N.zeros(shp)

    for n in range(shp[0]):
        if n == 0 or n == shp[0]-1:
            shear[n] = 0
        else:
            len1 = abs(parent.dx / N.sin(angle))
            len2 = abs(parent.dx / N.cos(angle))
            hyp = min((len1,len2))
            # In kilometres:
            shear[n] = ((wind_slice[n+1]-wind_slice[n-1])/(2*hyp))*1000.0
            T2grad[n] = ((T2_slice[n+1]-T2_slice[n-1])/(2*hyp))*1000.0
            # N.W.DX

    # pdb.set_trace()

    # Go from B to A
    # Find first location where T2 drops and shear is ?

    if method==1:
        ### METHOD 1: USING THRESHOLDS
        # By default
        gfidx = shp/2
        for n, s, t in zip(list(range(shp))[::-1],shear[::-1],T2grad[::-1]):
            if (abs(s)>2.0) and (t<2.0):
                gfidx = n
                break

    elif method==2 or method==3:

        ### METHOD 2: FINDING MAX GRADIENTS AND AVERAGING
        shear = abs(shear)
        T2grad = abs(T2grad)

        xsh_idx = N.where(shear == shear.max())[0][0]
        xtg_idx = N.where(T2grad == T2grad.max())[0][0]
        print(("Max shear loc: {0} ... max tempgrad loc: {1}".format(xsh_idx,xtg_idx)))

        if method==2:
            gfidx = int((xsh_idx + xtg_idx)/2.0)
        else:
            gfidx = max([xsh_idx,xtg_idx])

    return gfidx
    # maxshearloc[0][0] returns the integer

def compute_frontogenesis(parent,time,level):
    """
    Note that all variables fetched with parent.get have been
    destaggered and are at the same location.

    Returns:
        Front       :   Frontgenesis in Kelvin per second.
    """
    # import pdb; pdb.set_trace()
    #ds = 1 # Number of grid point to compute horizontal gradient over
    #dx = ds
    #dy = ds
    #dz = 1 # Normal for vertical
    tidx = parent.get_time_idx(time)
    tidxs = (tidx-1,tidx,tidx+1)

    if (tidx == 0) or (tidx == parent.wrf_times.shape[0]-1):
        Front = None
    else:
        nt,nl,ny,nx = parent.get('U',utc=tidx,level=1).shape
        U = N.zeros([3,3,ny,nx])
        V = N.zeros_like(U)
        W = N.zeros_like(U)
        T = N.zeros_like(U)
        if level == 2000:
            # Use the bottom three model levels
            P = N.zeros_like(U)
            for n, t in enumerate(tidxs):
                U[n,...] = parent.get('U',utc=t,level=1)
                V[n,...] = parent.get('V',utc=t,level=1)
                W[n,...] = parent.get('W',utc=t,level=1)
                T[n,...] = parent.get('T',utc=t,level=N.arange(3))
                P[n,...] = parent.get('pressure',utc=t,level=N.arange(3))
                # Average different in pressure between model levels
                # This field is passed into the gradient
                # THIS IS NOT USED RIGHT NOW
                dp = N.average(abs(0.5*(0.5*(P[2,2,:,:]-P[2,0,:,:]) +
                            0.5*(P[0,2,:,:]-P[0,0,:,:]))))

        elif isinstance(level,int):
            dp = 15 # hPa to compute vertical gradients

            for n, t in enumerate(tidxs):
                U[n,...] = parent.get_p('U',t,level)
                V[n,...] = parent.get_p('V',t,level)
                W[n,...] = parent.get_p('W',t,level)

                # 3D array has dimensions (vertical, horz, horz)
                T[n,...] = parent.get_p('T',t,(level-dp,level,level+dp))

                # Compute omega
                # P = rho* R* drybulb
                # drybulb = T/((P0/P)^(R/cp)

        drybulb = 273.15 + (T/((100000.0/(level*100.0))**(mc.R/mc.cp)))
        rho = (level*100.0)/(mc.R*drybulb)
        omega = -rho * mc.g * W

        # Time difference in sec
        dt = parent.utc[tidx+1]-parent.utc[tidx]
        dTdt, dTdz, dTdy, dTdx = N.gradient(T,dt,dp*100.0,parent.dy, parent.dx)
        # Gradient part
        grad = (dTdx**2 + dTdy**2)**0.5
        # Full derivative - value wrong for dgraddz here
        dgraddt, dgraddz, dgraddy, dgraddx = N.gradient(grad,dt,dp*100.0,parent.dy, parent.dx)
        # Full equation
        Front = dgraddt[1,1,:,:] + U[1,1,:,:]*dgraddx[1,1,:,:] + V[1,1,:,:]*dgraddy[1,1,:,:] # + omega[1,1,:,:]*dgraddz[1,1,:,:]
    return Front

def compute_accum_rain(parent,utc,accum_hr):
    """
    Needs to be expanded to include other precip
    """
    dn = utils.ensure_datenum(utc)
    idx0 = parent.get_time_idx(dn-(3600*accum_hr))
    idx1 = parent.get_time_idx(dn)
    # range_idx = range(start_idx,end_idx+1)
    total0 = (parent.get('RAINNC',utc=idx0) +
                parent.get('RAINC',utc=idx0))
    total1 = (parent.get('RAINNC',utc=idx1) +
                parent.get('RAINC',utc=idx1))
    accum = total1 - total0
    return accum

def compute_satvappres(parent,tidx,lvidx,lonidx,latidx,other):
    t = parent.get('drybulb',utc=tidx,level=lvidx,lons=lonidx,lats=latidx,other='C')
    es = 6.1*N.exp(0.073*t)
    return es

def compute_vappres(parent,tidx,lvidx,lonidx,latidx,other):
    RH = parent.get('RH',utc=tidx,level=lvidx,lons=lonidx,lats=latidx)
    es = parent.get('es',utc=tidx,level=lvidx,lons=lonidx,lats=latidx)
    e = RH*es
    return e

def compute_spechum(parent,tidx,lvidx,lonidx,latidx,other):
    es = parent.get('es',utc=tidx,level=lvidx,lons=lonidx,lats=latidx)
    p = parent.get('pressure',utc=tidx,level=lvidx,lons=lonidx,lats=latidx)
    q = 0.622*(es/p)
    return q

def compute_Td_2(parent,tidx,lvidx,lonidx,latidx,other='C'):
    """
    Another version of Td
    From p70 Djuric Weather Analysis
    """
    e = parent.get('e',utc=tidx,level=lvidx,lons=lonidx,lats=latidx)
    Td = 273*((N.log(e) - N.log(6.1))/(19.8 - (N.log(e) - N.log(6.1))))
    if other == 'K':
        Td =+ 273.15
    return Td

def compute_derivatives(parent,U,V,axis=None):
    dargs = (parent.dx,parent.dx)
    dkwargs = {'axis':None}
    if len(U.shape) == 2:
        dudx, dudy = N.gradient(U,parent.dx,parent.dx)
        dvdx, dvdy = N.gradient(V,parent.dx,parent.dx)
    elif len(U.shape) == 3:
        nt = U.shape[0]
        dudx = N.zeros_like(U)
        dvdx = N.zeros_like(U)
        dudy = N.zeros_like(U)
        dvdy = N.zeros_like(U)
        for t in range(nt):
            dudx[t,...], dudy[t,...] = N.gradient(U[t,...],parent.dx,parent.dx)
            dvdx[t,...], dvdy[t,...] = N.gradient(V[t,...],parent.dx,parent.dx)
    elif len(U.shape) == 4:
        nt,nz,nlat,nlon = U.shape
        dudx = N.zeros_like(U)
        dvdx = N.zeros_like(U)
        dudy = N.zeros_like(U)
        dvdy = N.zeros_like(U)
        for t in range(nt):
            for z in range(nz):
                dudx[t,z,:,:], dudy[t,z,:,:] = N.gradient(U[t,z,:,:],parent.dx,parent.dx)
                dvdx[t,z,:,:], dvdy[t,z,:,:] = N.gradient(V[t,z,:,:],parent.dx,parent.dx)
    return dudx, dudy, dvdx, dvdy

def compute_stretch_deformation(parent,U,V):
    dudx, dudy, dvdx, dvdy = compute_derivatives(parent,U,V)
    Est = dudx - dvdy
    return Est

def compute_shear_deformation(parent,U,V):
    dudx, dudy, dvdx, dvdy = compute_derivatives(parent,U,V)
    Esh = dudy + dvdx
    return Esh

def compute_total_deformation(parent,U,V):
    Esh = parent.compute_shear_deformation(U,V)
    Est = parent.compute_stretch_deformation(U,V)
    E = (Est**2 + Esh**2)**0.5
    return E

def compute_vorticity(parent,U,V):
    """
    Returns vertical vorticity. In future, change to allow all thre
    components to be be returned.
    """
    dudx, dudy, dvdx, dvdy = compute_derivatives(parent,U,V,)
    zeta = dvdx - dudy
    return zeta

def return_vorticity(parent,tidx,lvidx,lonidx,latidx,other):
    # pdb.set_trace()
    U = parent.get('U',tidx,lvidx,lonidx,latidx)[:,0,:,:]
    V = parent.get('V',tidx,lvidx,lonidx,latidx)[:,0,:,:]
    zeta = parent.compute_vorticity(parent,U,V)
    return zeta

def compute_fluid_trapping_diagnostic(parent,tidx,lvidx,lonidx,latidx,other):
    U = parent.get('U10',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    V = parent.get('V10',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    E = parent.compute_total_deformation(U,V)
    zeta = parent.compute_vorticity(parent,U,V)
    omega2 = 0.25*(E**2 - zeta**2)
    return omega2

def compute_divergence(parent,U,V):
    dudx, dudy, dvdx, dvdy = compute_derivatives(parent,U,V)
    div = dudx + dvdy
    return div

def compute_instantaneous_local_Lyapunov(parent,tidx,lvidx,lonidx,latidx,other):
    # import pdb; pdb.set_trace()
    U = parent.get('U',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    V = parent.get('V',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    E = parent.compute_total_deformation(U,V)
    zeta = parent.compute_vorticity(U,V)
    div = parent.compute_divergence(U,V)
    EzArr = (E**2  - zeta**2)**0.5
    Ez_nonan = N.nan_to_num(EzArr)
    # D =  0.5*(div + (E**2  - zeta**2)**0.5)
    D =  0.5*(div + Ez_nonan)
    # import pdb; pdb.set_trace()
    return D

def return_axis_of_dilatation_components(parent,tidx,lvidx=False,lonidx=False,
                                            latidx=False,other=False):
    U = parent.get('U10',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    V = parent.get('V10',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    Esh = parent.compute_shear_deformation(U,V)
    Est = parent.compute_stretch_deformation(U,V)
    E = parent.compute_total_deformation(U,V)
    zeta = parent.compute_vorticity(U,V)

    psi1 = 0.5 * N.arctan2(Esh,Est)
    # chi1 = psi1 + 0.5*N.nan_to_num(N.arcsin(zeta/E))
    chi1 = psi1 + 0.5*(N.arcsin(zeta/E))

    return N.cos(chi1), N.sin(chi1)

def compute_omega(parent,tidx,lvidx,lonidx,latidx,other):
    # Rising motion in Pa/s
    # dp/dt of air parcel
    W = parent.get('W',tidx,lvidx,lonidx,latidx)[0,:,:,:]
    rho = parent.get('density',tidx,lvidx,lonidx,latidx)[0,:,:,:]
    omega = -rho * -mc.g * W # I think it's meant to be minus g?
    # import pdb; pdb.set_trace()
    return omega

def compute_density(parent,tidx,lvidx,lonidx,latidx,other):
    drybulb = parent.get('drybulb',tidx,lvidx,lonidx,latidx,other='K')
    P = parent.get('pressure',tidx,lvidx,lonidx,latidx)
    rho = P/(mc.R*drybulb)
    # drybulb = 273.15 + (T/((100000.0/(level*100.0))**(mc.R/mc.cp)))
    return rho

def compute_CAPE(parent,tidx,lvidx,lonidx,latidx,other):
    """CAPE, taken from Pat Skinner's cookbook.
    """
    kwargs = dict(utc=tidx,level=lvidx,lons=lonidx,lats=latidx,other=other)
    # Environmental temperature
    t_env = parent.get('drybulb',**kwargs)
    # Lifted parcel temperature
    t_parc = parent.get('LPT',**kwargs)
    # Environmental pressure
    p = parent.get('pressure',**kwargs)
    # LCL pressure
    lcl_p = parent.get('LCL_P',**kwargs)
    # Vertical grid spacing at each grid point
    dz = parent.get('dz',**kwargs)

    t_diff = t_parc-t_env
    t_diff[t_diff > 0] = 0

    for lv in range(1,t_diff.shape[0]):
        t_diff[lv,:,:] = N.where(p[lv,:,:] > lcl_p, 0.0, t_diff[lv,:,:])

    CAPE = mc.g * N.trapz((t_diff/t_env),dx=dz[:-1,:,:],axis=0)
    return CAPE

def compute_lifted_parcel_temp(parent,tidx,lvidx,lonidx,latidx,other):
    assert True is False
    # TODO....

    kwargs = dict(utc=tidx,level=lvidx,lons=lonidx,lats=latidx,other=other)
    # t = parent.get('drybulb',**kwargs)
    # th = parent.get('theta',**kwargs)
    p = parent.get('pressure',**kwargs)

    # 2D slices for the layer to be lifted.
    the_base = 0
    p_base = 0
    t_base = 0

    # t_parcel = N.zeros_like()
    # th_parcel = N.zeros_like()

    return

def __compute_CAPE(parent,tidx,lvidx,lonidx,latidx,other):
    """
    INCOMPLETE!

    CAPE method based on GEMPAK's pdsuml.f function

    Inputs:

    tidx,lvidx,lonidx,latidx  :   dictionary of level/time/lat/lon



    Outputs:
    CAPE    :   convective available potential energy
    CIN     :   convective inhibition
    """
    # Make sure all levels are obtained
    #tidx,lvidx,lonidx,latidx['lv'] = slice(None,None)

    totalCAPE = 0
    totalCIN = 0

    theta = parent.get('theta',tidx,lvidx,lonidx,latidx)
    Z = parent.get('Z',tidx,lvidx,lonidx,latidx)

    for lvidx in range(theta.shape[1]-1):
        if lvidx < 20:
            continue
        # This should loop over the levels?
        """
        z1      :   bottom of layer (index)
        z2      :   top of layer (index)
        th1     :   theta (environ) at z1
        th2     :   theta (environ) at z2
        thp1    :   theta (parcel) at z1
        thp2    :   theta (parcel) at z2
        """

        z1 = Z[0,lvidx,...]
        z2 = Z[0,lvidx+1,...]

        th1 = theta[0,lvidx,...]
        th2 = theta[0,lvidx+1,...]

        thp1 = 0
        thp2 = 0

        capeT = 0.0
        cinT = 0.0

        dt2 = thp2 - th2
        dt1 = thp1 - th1

        dz = z2 - z1

        dt1_pos_ma = N.ma.masked_greater(dt1,0)
        dt1_neg_ma = N.ma.masked_less(dt1,0)

        dt2_pos_ma = N.ma.masked_greater(dt2,0)
        dt2_neg_ma = N.ma.masked_less(dt2,0)

        dt1_pos = N.select([dt1>0],[dt1])
        dt1_neg = N.select([dt1<0],[dt1])
        dt2_pos = N.select([dt2>0],[dt1])
        dt2_neg = N.select([dt2<0],[dt1])

        pdb.set_trace()

        if (dt1 > 0) and (dt2 > 0):
            capeT = ((dt2 + dt1)*dz)/(th2+th1)
        elif dt1 > 0:
            ratio = dt1/(dt1-dt2)
            zeq = z1 + (dz*ratio)
            teq = th1 + ((th2-th1)*ratio)
            capeT = (dt1*(zeq-z1))/(th1+teq)
            cinT = (dt2*(z2-zeq))/(th2+teq)
        elif dt2 > 0:
            ratio = dt2/(dt2-dt1)
            zfc = z2-(dz*ratio)
            tfc = th2-((th2-th1)*ratio)
            capeT = (dt2*(z2-zfc)/(tfc+th2))
            cinT = (dt1*(zfc-z1)/(tfc+th1))
        else:
            cinT = ((dt2+dt1)*dz)/(th2+th1)

        if capeT > 0:
            CAPE = capeT * cc.g
        else:
            CAPE = 0

        if cinT < 0:
            CIN = cinT * cc.g
        else:
            CIN = 0

            totalCAPE += CAPE
            totalCIN += CIN

    return totalCAPE,totalCIN

def compute_RH(parent,tidx,lvidx,lonidx,latidx,other):

    T = parent.get('drybulb',tidx,lvidx,lonidx,latidx,other='C')
    Td = parent.get('Td',tidx,lvidx,lonidx,latidx)
    RH = N.exp(0.073*(Td-T))
    # pdb.set_trace()
    return RH*100.0

def compute_temp_advection(parent,tidx,lvidx,lonidx,latidx,other):
    U = parent.get('U',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    V = parent.get('V',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    T = parent.get('drybulb',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    dTdx, dTdy = N.gradient(T,parent.DX,parent.DY)
    field = -U*dTdx - V*dTdy
    # pdb.set_trace()
    return field

def compute_PMSL_gradient(parent,tidx,lvidx,lonidx,latidx,other):
    P = parent.get('PMSL',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    dPdx, dPdy = N.gradient(P,parent.dx,parent.dy)
    field = N.sqrt(dPdx**2 + dPdy**2)
    # import pdb; pdb.set_trace()
    return field

def compute_T2_gradient(parent,tidx,lvidx,lonidx,latidx,other):
    T2 = parent.get('T2',tidx,lvidx,lonidx,latidx)[0,0,:,:]
    dTdx, dTdy = N.gradient(T2,parent.dx,parent.dy)
    field = N.sqrt(dTdx**2 + dTdy**2)
    # import pdb; pdb.set_trace()
    return field

def compute_dryairmass(parent,tidx,lvidx,lonidx,latidx,other):
    MU = parent.get('MU',tidx,lvidx,lonidx,latidx)
    MUB = parent.get('MUB',tidx,lvidx,lonidx,latidx)
    return MU + MUB

def compute_pmsl(parent,tidx,lvidx,lonidx,latidx,other):
    P = parent.get('PSFC',tidx,lvidx,lonidx,latidx)
    T2 = parent.get('T2',tidx,lvidx,lonidx,latidx)
    HGT = parent.get('HGT',tidx,lvidx,lonidx,latidx)

    temp = T2 + (6.5*HGT)/1000.0
    pmsl = P*N.exp(9.81/(287.0*temp)*HGT)

    #sm = kwargs.get('smooth',1)
    #data = pmsl[0,::sm,::sm]
    #return data
    return pmsl

def compute_buoyancy(parent,tidx,lvidx,lonidx,latidx,other=False):
    """
    Method from Adams-Selin et al., 2013, WAF
    """
    theta = parent.get('theta',tidx,lvidx,lonidx,latidx)
    thetabar = N.mean(theta)
    qv = parent.get('QVAPOR',tidx,lvidx,lonidx,latidx)
    qvbar = N.mean(qv)

    B = cc.g * ((theta-thetabar)/thetabar + 0.61*(qv - qvbar))
    return B

def compute_mixing_ratios(parent,tidx,lvidx,lonidx,latidx,other=False):
    qv = parent.get('QVAPOR',tidx,lvidx,lonidx,latidx)
    qc = parent.get('QCLOUD',tidx,lvidx,lonidx,latidx)
    qr = parent.get('QRAIN',tidx,lvidx,lonidx,latidx)

    try:
        qi = parent.get('QICE',tidx,lvidx,lonidx,latidx)
    except KeyError:
        print("MP scheme has no ice data.")
        qi = 0

    try:
        qs = parent.get('QSNOW',tidx,lvidx,lonidx,latidx)
    except KeyError:
        print("MP scheme has no snow data.")
        qs = 0

    try:
        qg = parent.get('QGRAUP',tidx,lvidx,lonidx,latidx)
    except KeyError:
        print("MP scheme has no graupel data.")
        qg = 0

    rh = qc + qr + qi + qs + qg
    rv = qv

    return rh, rv

def compute_qtotal(parent,tidx,lvidx,lonidx,latidx,other):
    qtotal, _ = parent.compute_mixing_ratios(tidx,lvidx,lonidx,latidx)
    return qtotal

def compute_dptp(parent,tidx,lvidx,lonidx,latidx,other):
    dpt = parent.get('dpt',tidx,lvidx,lonidx,latidx)
    dpt_mean = N.mean(dpt)
    dptp = dpt - dpt_mean
    return dptp

def compute_T2_pertub(parent,tidx,lvidx,lonidx,latidx,other):
    T2 = parent.get('T2',tidx,lvidx,lonidx,latidx)
    T2_mean = N.mean(T2)
    T2p = T2-T2_mean
    return T2p

def compute_Q_pert(parent,tidx,lvidx,lonidx,latidx,other):
    Q = parent.get('QVAPOR',tidx,lvidx,lonidx,latidx)
    Q_mean = N.mean(Q)
    Qp = Q-Q_mean
    return Qp

def compute_dpt(parent,tidx,lvidx,lonidx,latidx,other):
    """
    Potential: if surface level is requested, choose sigma level just
    about the surface. I don't think this affects any
    other dictionaries around...
    """
    # if tidx,lvidx,lonidx,latidx['lv'] == 0:
        # tidx,lvidx,lonidx,latidx['lv'] = 0
    theta = parent.get('theta',tidx,lvidx,lonidx,latidx)
    rh, rv = parent.compute_mixing_ratios(tidx,lvidx,lonidx,latidx)

    dpt = theta * (1 + 0.61*rv - rh)
    return dpt

def compute_geopotential_height(parent,tidx,lvidx,lonidx,latidx,other):
    geopotential = parent.get('PH',tidx,lvidx,lonidx,latidx) + parent.get('PHB',tidx,lvidx,lonidx,latidx)
    Z = geopotential/9.81
    return Z

def compute_geopotential(parent,tidx,lvidx,lonidx,latidx,other):
    geopotential = parent.get('PH',tidx,lvidx,lonidx,latidx) + parent.get('PHB',tidx,lvidx,lonidx,latidx)
    return geopotential

def compute_wind10(parent,tidx,lvidx,lonidx,latidx,other):
    u = parent.get('U10',tidx,lvidx,lonidx,latidx)
    v = parent.get('V10',tidx,lvidx,lonidx,latidx)
    data = N.sqrt(u**2 + v**2)
    return data

def compute_pressure(parent,tidx,lvidx,lonidx,latidx,other):
    PP = parent.get('P',tidx,lvidx,lonidx,latidx)
    PB = parent.get('PB',tidx,lvidx,lonidx,latidx)
    pressure = PP + PB
    return pressure

def compute_drybulb(parent,tidx,lvidx,lonidx,latidx,other='K'):
    theta = parent.get('theta',tidx,lvidx,lonidx,latidx)
    P = parent.get('pressure',tidx,lvidx,lonidx,latidx)

    # Theta-e at level 2
    drybulb = theta*((P/100000.0)**(287.04/1004.0))
    if other=='K':
        return drybulb
    elif other=='C':
        return drybulb-273.15

def compute_theta(parent,tidx,lvidx,lonidx,latidx,other):
    theta = parent.get('T',tidx,lvidx,lonidx,latidx)
    Tbase = 300.0
    theta = Tbase + theta
    return theta

def compute_wind(parent,tidx,lvidx,lonidx,latidx,other):
    # pdb.set_trace()
    u = parent.get('U',tidx,lvidx,lonidx,latidx)
    v = parent.get('V',tidx,lvidx,lonidx,latidx)
    data = N.sqrt(u**2 + v**2)
    return data

def compute_shear(parent,tidx,lvidx,lonidx,latidx,other=False):
    """
    Args:
        other   :      dictionary of 'top' and 'bottom' levels, km

    Todos:
        * Could make this faster with numpy.digitize()?
    """
    if not other:
        print("No shear heights specified. Using 0-6 km by default.")
        topm = 6000.0
        botm = 0.0
        # print("Choose top and bottom for shear calc.")
        # raise Exception
    else:
        topm = other['top']*1000
        botm = other['bottom']*1000

    u = parent.get('U',tidx,lvidx,lonidx,latidx)
    v = parent.get('V',tidx,lvidx,lonidx,latidx)
    Z = parent.get('Z',tidx,lvidx,lonidx,latidx)

    topidx = N.zeros((parent.y_dim,parent.x_dim))
    botidx = N.zeros((parent.y_dim,parent.x_dim))
    ushear = N.zeros((parent.y_dim,parent.x_dim))
    vshear = N.zeros((parent.y_dim,parent.x_dim))

    for j in range(parent.x_dim):
        for i in range(parent.y_dim):
            # import pdb; pdb.set_trace()
            topidx[i,j] = round(N.interp(
                            topm,Z[0,:,i,j],list(range(parent.z_dim))))
            botidx[i,j] = round(N.interp(
                            botm,Z[0,:,i,j],list(range(parent.z_dim))))
            ushear[i,j] = u[0,topidx[i,j],i,j] - u[0,botidx[i,j],i,j]
            vshear[i,j] = v[0,topidx[i,j],i,j] - v[0,botidx[i,j],i,j]

    # Find indices of bottom and top levels
    # topidx = N.where(abs(Z-topm) == abs(Z-topm).min(axis=1))
    # botidx = N.where(abs(Z-botm) == abs(Z-botm).min(axis=1))

    # ushear = u[0,:,topidx] - u[0,:,botidx]
    # vshear = v[0,topidx,:,:] - v[0,botidx,:,:]

    shear = N.sqrt(ushear**2 + vshear**2)
    # pdb.set_trace()
    return shear

def compute_thetae(parent,tidx,lvidx,lonidx,latidx,other):
    P = parent.get('pressure',tidx,lvidx,lonidx,latidx) # Computed
    Drybulb = parent.get('temp',tidx,lvidx,lonidx,latidx)
    Q = parent.get('Q',tidx,lvidx,lonidx,latidx)

    thetae = (Drybulb + (Q * cc.Lv/cc.cp)) * (cc.P0/P) ** cc.kappa
    return thetae

def compute_olr(parent,tidx,lvidx,lonidx,latidx,other):
    OLR = parent.get('OLR',tidx,lvidx,lonidx,latidx)
    sbc = 0.000000056704
    ir = ((OLR/sbc)**0.25) - 273.15
    return ir

def compute_REFL_comp(parent,tidx,lvidx,lonidx,latidx,other=False):
    # lvidx = None
    # pdb.set_trace()
    refl = parent.get('REFL_10CM',tidx,lvidx,lonidx,latidx,other)[:,:,:,:]
    refl_comp = N.max(refl,axis=1)
    return refl_comp

def compute_comp_ref(parent,tidx,lvidx,lonidx,latidx,other):
    """Amend this so variables obtain at start fetch only correct date, lats, lons
    All levels need to be fetched as this is composite reflectivity
    """
    T2 = parent.get('T2',tidx,False,lonidx,latidx)
    # QR = parent.nc.variables['QRAIN'][PS['t'],:,PS['la'],PS['lo']]
    QR = parent.get('QRAIN',tidx,False,lonidx,latidx) # This should get all levels
    PSFC = parent.get('PSFC',tidx,False,lonidx,latidx)

    try:
        QS = parent.get('QSNOW',tidx,False,lonidx,latidx)
    except:
        QS = N.zeros(N.shape(QR))
    rhor = 1000.0
    rhos = 100.0
    rhog = 400.0
    rhoi = 917.0

    no_rain = 8.0E6
    # How do I access this time?
    no_snow = 2.0E6 * N.exp(-0.12*(T2-273.15))
    no_grau = 4.0E6

    density = N.divide(PSFC,(287.0 * T2))
    Qra_all = QR[0,...]
    Qsn_all = QS[0,...]

    for j in range(len(Qra_all[1,:,1])):
        curcol_r = []
        curcol_s = []
        for i in range(len(Qra_all[1,1,:])):
                maxrval = N.max(Qra_all[:,j,i])
                maxsval = N.max(Qsn_all[:,j,i])
                curcol_r.append(maxrval)
                curcol_s.append(maxsval)
        N_curcol_r = N.array(curcol_r)
        N_curcol_s = N.array(curcol_s)
        if j == 0:
            Qra = N_curcol_r
            Qsn = N_curcol_s
        else:
            Qra = N.row_stack((Qra, N_curcol_r))
            Qsn = N.row_stack((Qsn, N_curcol_s))

    # Calculate slope factor lambda
    lambr = (N.divide((3.14159 * no_rain * rhor), N.multiply(density, Qra)+N.nextafter(0,1))) ** 0.25
    lambs = N.exp(-0.0536 * (T2 - 273.15))

    # Calculate equivalent reflectivity factor
    Zer = (720.0 * no_rain * (lambr ** -7.0)) * 1E18
    Zes = (0.224 * 720.0 * no_snow * (lambr ** -7.0) * (rhos/rhoi) ** 2) * 1E18
    Zes_int = N.divide((lambs * Qsn * density), no_snow)
    Zes = ((0.224 * 720 * 1E18) / (3.14159 * rhor) ** 2) * Zes_int ** 2

    Ze = N.add(Zer, Zes)
    dBZ = N.nan_to_num(10*N.log10(Ze))
    return dBZ

def compute_simref_atlevel(parent,level=1):
    pass
    return data

def compute_DCP(parent):
    """
    Derecho Composite Parameter (Evans and Doswell, 2001, WAF)
    And info from SPC Mesoanalyses
    """
    DCAPE = parent.get('DCAPE')
    MUCAPE = parent.get('MUCAPE')
    shear_0_6 = parent.get('shear',0,6)
    meanwind_0_6 = parent.get('meanwind',0,6)
    DCP = (DCAPE/980.0)*(MUCAPE/2000.0)*(shear_0_6/20.0)*(meanwind_0_6/16.0)
    return DCP

def compute_DCAPE(parent):

    pass

def compute_thetae(parent,tidx,lvidx,lonidx,latidx,other):
    P = parent.get('pressure',tidx,lvidx,lonidx,latidx)
    T = parent.get('drybulb',tidx,lvidx,lonidx,latidx,units='K')
    Td = parent.get('Td',tidx,lvidx,lonidx,latidx)
    p2, t2 = thermo.drylift(P,T,Td)
    x = thermo.wetlift(p2,t2,100.0)
    thetae = thermo.theta(100.0, x, 1000.0)
    return thetae


def compute_Td(parent,tidx,lvidx,lonidx,latidx,other):
    """
    Using HootPy equation
    """
    Q = parent.get('QVAPOR',tidx,lvidx,lonidx,latidx)
    P = parent.get('pressure',tidx,lvidx,lonidx,latidx)
    w = N.divide(Q, N.subtract(1,Q))
    e = N.divide(N.multiply(w,P), N.add(0.622,w))/100.0
    a = N.multiply(243.5,N.log(N.divide(e,6.112)))
    b = N.subtract(17.67,N.log(N.divide(e,6.112)))
    Td = N.divide(a,b)
    # pdb.set_trace()
    return Td



def compute_strongest_wind(parent,tidx,lvidx,lonidx,latidx,other):
    """
    Pass the array of time indices and it will find the max
    along that axis.
    """
    if 'WSPD10MAX' in parent.fields:
        ww = parent.get('WSPD10MAX',tidx,lvidx,lonidx,latidx)
        if ww.max() > 0.1:
            print("Using WSPD10MAX data")
            wind = ww
        else:
            print("Using wind10 data")
            wind = parent.get('wind10',tidx,lvidx,lonidx,latidx)
    else:
        print("Using wind10 data")
        wind = parent.get('wind10',tidx,lvidx,lonidx,latidx)
    wind_max = N.amax(wind,axis=0)
    # wind_max_smooth = parent.test_smooth(wind_max)
    # return wind_max_smooth

    return wind_max

def return_updraught_helicity_02(parent,tidx,lvidx,lonidx,latidx,other):
    return return_updraught_helicity(parent,tidx,lvidx,lonidx,latidx,other,z0=0,z1=2000)

def return_updraught_helicity_25(parent,tidx,lvidx,lonidx,latidx,other):
    return return_updraught_helicity(parent,tidx,lvidx,lonidx,latidx,other,z0=2000,z1=5000)

def return_updraught_helicity(parent,tidx,lvidx,lonidx,latidx,other,z0=2000,z1=5000):
    def method1():
        # Optimisation no. 1:
        # Simplest!
        __z = parent.get("Z",tidx,None,lonidx,latidx)[:,:,:,:]
        nt, nlv, nlat, nlon = __z.shape

        # Closest model levels for each lat/lon point to desired AGL-km levels
        idx0 = N.argmin(N.abs(__z-z0),axis=1)[0,:,:]
        idx1 = N.argmin(N.abs(__z-z1),axis=1)[0,:,:]

        # Now, subset xi and w to be only levels of interest
        bot = int(N.median(idx0))
        top = int(N.median(idx1))
        zs = slice(bot,top+1,1)

        # Vorticity and updraft
        _u = parent.get("U",tidx,None,lonidx,latidx)[:,zs,:,:].data
        _v = parent.get("V",tidx,None,lonidx,latidx)[:,zs,:,:].data
        _w = parent.get("W",tidx,None,lonidx,latidx)[:,zs,:,:].data
        utils.enforce_same_dimensions(_u,_v,_w)
        _xi = compute_vorticity(parent=parent,U=_u,V=_v)

        # Interpolate _w, _xi in the vertical
        nt, nlv, nlat, nlon = _w.shape
        _tidx = N.arange(nt)
        _lvidx = N.arange(nlv)
        i_lvidx = N.arange(0.5,nlv,1.0)
        _latidx = N.arange(nlat)
        _lonidx = N.arange(nlon)
        # oldidx = (_tidx,_lvidx,_latidx,_lonidx)
        # RGI = RegularGridInterpolator((_tidx,_lvidx,_latidx,_lonidx),_w)


        #newidx = N.array([_tidx,i_lvidx,_latidx,_lonidx])
        #w = RGI(newidx)
        #xi = RGI(newidx)

        interp_func = interp1d(x=_lvidx,y=_w,axis=1)
        w = interp_func(i_lvidx)

        pdb.set_trace()
        #xi = interp1d(points=oldidx,values=_xi,xi=newidx)



        dz = N.diff(__z[:,zs,:,:],axis=1)
        # Final UH computation:
        assert xi.ndim == 4
        UH = N.sum(xi*w*dz,axis=1)[0,:,:]
        return UH

    def method2():
        pass
        # no. 2:
        # Need to optimise with a new scheme that splits domain into
        # regions (use scikit?) and computes UH in areas separately

    def method3():
        # no. 3:
        # For all permutations of bot_uniq:top_uniq,
        # find the indices that match this and compute UH

        # Tip: find all points that have e.g. 10 levels, those with 11, 12, etc
        # Then extract only those points, compute, and place the values in UH array
        bot_unique = N.unique(idx0)
        top_unique = N.unique(idx1)

        for bot, top in itertools.product(bot_unique,top_unique):
            pass


    # Meta-data
    # nt, nlv, nlat, nlon = z.shape

    method = 1
    METHODS = {1:method1,2:method2,3:method3}
    UH = METHODS[method]()


    return UH

def return_maxcol_updraught(parent,tidx,lvidx,lonidx,latidx,other):
    w = parent.get("W",tidx,lvidx,lonidx,latidx)
    wmax = N.max(w,axis=1)
    return wmax[:,N.newaxis,:,:]

def compute_updraught_helicity(parent,u,v,w,dz=None,z=None,z0=2000,z1=5000):
    """Eq. 11 & 12 from Kain et al 2008, WAF.

    Requires either z or dz to be specified.

    Args:
        u,v,w (N.ndarray)           :   3-D fields of wind, all same size.
        dz (N.ndarray)              :   1-D vector of difference in z
        z (N.ndarray)               :   1-D vector height above ground level (m)
        z0 (int,float,optional)     :   bottom level in m
        z1 (int,float,optional)     :   top level in m

    Todo:
        * Not sure if this works.
    """

    # Need to do mid-point approximation here
    # If dz, don't need to do anything.
    # If z, need to do differences.
    UH = (vort_z * w) * dz

    # Logic for smoothing?
    return UH

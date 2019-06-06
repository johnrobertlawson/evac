import pdb
import os

import numpy as N

import evac.derived.pat as pat

def compute_CAPE(nc=None,P=None,PB=None,PH=None,PHB=None,T=None,QVAPOR=None,QCLOUD=None,
                    QRAIN=None,QICE=None,QSNOW=None,QGRAUP=None,PSFC=None,tidx=None):
    """
    You need either:
        * nc
        * tidx
    Or:
        * P
        * PB
        * PH
        * PHB
        * T
        * QVAPOR
        * QCLOUD
        * QRAIN
        * QICE
        * QSNOW
        * QGRAUP
        * PSFC

    The second group are the WRF "keys" for various fields as a mnemonic. These are 3D numpy
    fields, apart from PSFC, which is 2D.

    Find the WRF "key" meanings here: 
    http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap5.htm#fields

    ============
    THE FOLLOWING IS JUST FOR JRL - NOTES
    ============
    We need:
        * t_env - 3d numpy array [z,y,x] of environmental temperature (K)
        * t_parc - 3d numpy array [z,y,x] of lifted parcel temperature (K)
        * p - 3d numpy array [z,y,x] of environmental pressure (hPa)
        * lcl_p - 2d numpy array [y,x] of LCL pressure (hPa)
        * dz - 3d numpy array [y,x] of vertical grid spacing for each grid point (m)


    Pat does:
        calc_cape(t_v, t_100mb_parcel, p, lcl_p_100mb, dz)    
    Where:
        * t_v = calc_thv(temp, qv, qt)
            * temp is calc_t(th,p)
                * p = (p + pb) / 100
                    * p and pb (P, PB) from WRF - NEED
                * th is T from WRF, + 300 - NEED
            * qv is QVAPOR from WRF - NEED
            * qt is sum of QCLOUD, QRAIN, QICE, QSNOW, QGRAUP - NEED ALL
        * t_100mb_parcel = calc_parcel_dj(p, th_e_100mb[f,:,:], t_v_100mb, p_100mb)
            * (p)
            * th_e_100mb[f,:,:], lcl_t_100mb = calc_the_bolt(p_100mb, t_100mb, qv_100mb)
                * f is what? forecast time I think...
                * p_100mb = N.ma.average(masked_p, axis=0, weights=dz)
                    * masked_p = N.ma.masked_where((p_sfc_temp - p) > 100., (p))
                        * p_sfc_temp =  PSFC from WRF [0,:,:] / 100 - NEED
                        * (p)
                    * z, dz = calc_height(ph, phb)
                        * PH, PHB from WRF
                * t_100mb = N.ma.average(masked_temp, axis=0, weights=dz)
                    * masked_temp = N.ma.masked_where((p_sfc_temp - p) > 100., (temp))
                * qv_100mb = N.ma.average(masked_qv, axis=0, weights=dz)
                    * masked_qv = N.ma.masked_where((p_sfc_temp - p) > 100., (qv))
            * t_v_100mb = N.ma.average(masked_t_v, axis=0, weights=dz)
                * masked_t_v = N.ma.masked_where((p_sfc_temp - p) > 100., (t_v))
            * (p_100mb)
        * (p)
        * lcl_p_100mb =  1000. / (th_v_100mb[f,:,:] / lcl_t_100mb)**(1004. / 287.)
            * th_v_100mb[f,:,:] = N.ma.average(masked_th_v, axis=0, weights=dz)
                * masked_th_v = N.ma.masked_where((p_sfc_temp - p) > 100., (th_v))
                    * th_v = calc_thv(th, qv, qt)
                        * (th) - don't forget the +300
                        * (qv)
                        * (qt)
            * (lcl_t_100mb)
        * (dz)
            
    """
    if nc is not None:
        # only one time
        assert isinstance(tidx,int)
        P = nc.variables['P'][tidx,...]
        PB = nc.variables['PB'][tidx,...]
        PH = nc.variables['PH'][tidx,...]
        PHB = nc.variables['PHB'][tidx,...]
        T = nc.variables['T'][tidx,...]
        QVAPOR = nc.variables['QVAPOR'][tidx,...]
        QCLOUD = nc.variables['QCLOUD'][tidx,...]
        QRAIN = nc.variables['QRAIN'][tidx,...]
        QICE = nc.variables['QICE'][tidx,...]
        QSNOW = nc.variables['QSNOW'][tidx,...]
        QGRAUP = nc.variables['QGRAUP'][tidx,...]
        PSFC = nc.variables['PSFC'][tidx,...]

    qt = QCLOUD + QRAIN + QICE + QSNOW + QGRAUP
    th = T+300.0
    p = (P+PB)/100.0
    z, dz = pat.calc_height(PH, PHB)
    # pdb.set_trace()
    p_sfc_temp = PSFC/ 100.0
    temp = pat.calc_t(th,p)

    t_v = pat.calc_thv(temp, QVAPOR, qt)
    th_v = pat.calc_thv(th, QVAPOR, qt)

    masked_p = N.ma.masked_where((p_sfc_temp - p) > 100., (p))
    masked_temp = N.ma.masked_where((p_sfc_temp - p) > 100., (temp))
    masked_qv = N.ma.masked_where((p_sfc_temp - p) > 100., (QVAPOR))
    masked_t_v = N.ma.masked_where((p_sfc_temp - p) > 100., (t_v))
    masked_th_v = N.ma.masked_where((p_sfc_temp - p) > 100., (th_v))

    p_100mb = N.ma.average(masked_p, axis=0, weights=dz)
    t_100mb = N.ma.average(masked_temp, axis=0, weights=dz)
    qv_100mb = N.ma.average(masked_qv, axis=0, weights=dz)
    t_v_100mb = N.ma.average(masked_t_v, axis=0, weights=dz)

    th_e_100mb, lcl_t_100mb = pat.calc_the_bolt(p_100mb, t_100mb, qv_100mb)
    th_v_100mb = N.ma.average(masked_th_v, axis=0, weights=dz)

    t_100mb_parcel = pat.calc_parcel_dj(p, th_e_100mb, t_v_100mb, p_100mb)

    lcl_p_100mb = 1000. / (th_v_100mb / lcl_t_100mb)**(1004. / 287.)

    CAPE = pat.calc_cape(t_v,t_100mb_parcel,p,lcl_p_100mb,dz)
    # pdb.set_trace()
    # Not masked
    return CAPE.data

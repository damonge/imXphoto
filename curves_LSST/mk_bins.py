import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,InterpolatedUnivariateSpline
from scipy.special import erf
from scipy.optimize import brentq
import py_cosmo_mad as csm

plot_stuff=False

def get_zarr(s_photoz,n_width,z_max) :
    dz_half=s_photoz*n_width*0.5
    zc_arr=[dz_half/(1-dz_half)]
    zw_arr=[dz_half*(1+zc_arr[0])]
    z_last=zc_arr[0]
    while z_last<=z_max :
        zc_new=(z_last+dz_half*(2+z_last))/(1-dz_half)
        zw_new=dz_half*(1+zc_new)
        zc_arr.append(zc_new)
        zw_arr.append(zw_new)
        z_last=zc_new
    zc_arr=np.array(zc_arr)
    zw_arr=np.array(zw_arr)
    indices=np.where(zc_arr<z_max)
    zc_arr_out=zc_arr[indices]
    zw_arr_out=zw_arr[indices]
    sz_arr_out=s_photoz*(1+zc_arr_out)

    return zc_arr_out,zw_arr_out,sz_arr_out

def plot_bins(fname,z_max,nz_func) :
    nz_plot=128
    data=np.loadtxt(fname,unpack=True)
    zc_arr=(data[1]+data[0])/2
    zw_arr=(data[1]-data[0])/2
    sz_arr=data[2]
    for i in np.arange(len(zc_arr)) :
        zc=zc_arr[i]
        zw=zw_arr[i]
        sz=sz_arr[i]
        z0=zc-zw
        zf=zc+zw
        zmax=zf+3*zw
        if z0<3*zw :
            zmin=0
        else :
            zmin=z0-3*zw
        denom=1./(np.sqrt(2)*sz)
        zarr=zmin+(zmax-zmin)*(np.arange(nz_plot)+0.5)/nz_plot
        nzarr=np.array([nz_func(z) for z in zarr])
        wzarr=np.array([erf((zf-z)*denom)-erf((z0-z)*denom) for z in zarr])
        wzarr/=np.amax(wzarr)
        nwarr=nzarr*wzarr
        if plot_stuff :
            plt.plot(zarr,nwarr)
    z_global_arr=z_max*1.2*(np.arange(nz_plot)+0.5)/nz_plot
    nz_global_arr=np.array([nz_func(z) for z in z_global_arr])
    if plot_stuff :
        plt.plot(z_global_arr,nz_global_arr,'k--')
        plt.xlim([0,z_max*1.2])
        plt.show()

def get_kmax(pcs,dz,st) :
    def sig_minus_target(lk) :
        k=10.**lk
        return np.sqrt(pcs.sig0_L(1./k,1./k,"SharpK","SharpK"))*dz-st
    return 10.**brentq(sig_minus_target,-3.,1.)

def get_bins(s_photoz,z_max,n_photoz_width,nzfile=None,fname_out=None,sigma_thr=None,lmax=2000) :
    zc,zw,sz=get_zarr(s_photoz,n_photoz_width,z_max)
    print "%d bins"%(len(zc))

    if fname_out!=None :
        data_write=np.zeros([6,len(zc)])
        data_write[0]=zc-zw
        data_write[1]=zc+zw
        data_write[2]=sz
        data_write[3,:]=0
        data_write[4,:]=0
        data_write[5,:]=lmax
        
        pcs=csm.PcsPar()
        pcs.set_verbosity(0)
        pcs.background_set(0.3,0.7,0.05,-1.0,0.0,0.7,2.7255)
        pcs.set_linear_pk("EH",-4.,3.,0.01,0.96,0.8)
        
        if sigma_thr!=None :
            for i in np.arange(len(zc)) :
                z=zc[i]
                gf=pcs.growth_factor(1/(1+z))/pcs.growth_factor(1)
                chi=pcs.radial_comoving_distance(1./(1+z))
                kmax=get_kmax(pcs,gf,sigma_thr)
                ell=min(int(chi*kmax),lmax)
                print z,ell
                data_write[5,i]=ell

        np.savetxt(fname_out,np.transpose(data_write),fmt="%lf %lf %lf %d %d %d",header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

    if nzfile!=None :
        z,nz=np.loadtxt(nzfile,unpack=True)
        nz=InterpolatedUnivariateSpline(z,nz,k=3)
        plot_bins(fname_out,z_max,nz)

#Bins for red galaxies
n_photoz=3
#s_photoz_red =0.02; z_max_red =1.4; nzfile_red =None;#"nz_red.txt";
#s_photoz_blue=0.05; z_max_blue=3.0; nzfile_blue=None;#"nz_blue.txt";
#s_photoz_gold=0.05; z_max_gold=3.0; nzfile_gold=None;#"nz_shear_fiducial.txt";
s_photoz_gold=0.03; z_max_gold=3.0; nzfile_gold="nz_shear_fiducial.txt";

pcs=csm.PcsPar()
pcs.set_verbosity(0)
pcs.background_set(0.3,0.7,0.05,-1.0,0.0,0.7,2.7255)
pcs.set_linear_pk("EH",-4.,3.,0.01,0.96,0.8)
gf=pcs.growth_factor(1./(1+0.))/pcs.growth_factor(1)
sthr_arr=(np.arange(50)+1.)*0.021
kmax_arr=np.array([get_kmax(pcs,gf,st) for st in sthr_arr])
if plot_stuff :
    plt.plot(sthr_arr,kmax_arr);
    plt.xlabel("$\\sigma_{\\rm thr}$",fontsize=16)
    plt.ylabel("$k_{\\rm max}(z=0)\\,[h\\,{\\rm Mpc}^{-1}]$",fontsize=16)
    plt.show()
np.savetxt('kvals_st.txt',np.transpose([sthr_arr,kmax_arr]))

get_bins(s_photoz_gold,z_max_gold,n_photoz,nzfile=nzfile_gold,
         fname_out="bins_gold_sthr0p50.txt",sigma_thr=0.50)
get_bins(s_photoz_gold,z_max_gold,n_photoz,nzfile=nzfile_gold,
         fname_out="bins_gold_sthr0p75.txt",sigma_thr=0.75)
get_bins(s_photoz_gold,z_max_gold,n_photoz,nzfile=nzfile_gold,
         fname_out="bins_gold_lmax2000.txt",sigma_thr=None)

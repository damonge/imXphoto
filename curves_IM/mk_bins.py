import numpy as np
from scipy.optimize import brentq
import py_cosmo_mad as csm

plot_stuff=False
KMAX=0.2
nuHI=1420.405

pcs=csm.PcsPar()
pcs.set_verbosity(1)
pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.725)
pcs.set_linear_pk("EH",-4.,3.,0.01,0.96,0.8)

def get_kmax(z,st) :
    dz=pcs.growth_factor(1./(1+z))/pcs.growth_factor(1.)
    def sig_minus_target(lk) :
        k=10.**lk
        return np.sqrt(pcs.sig0_L(1./k,1./k,"SharpK","SharpK"))*dz-st
    return 10.**brentq(sig_minus_target,-3.,1.)

def mk_binfiles(prefix_out,nu0,nuf,nnu,sigma_thr=None,lmax=2000) :
    zf=nuHI/nu0-1
    z0=nuHI/nuf-1
    dz=(zf-z0)/nnu
    z0arr=z0+dz*np.arange(nnu)
    zfarr=z0arr+dz
    zmarr=z0arr+0.5*dz
    szarr=np.zeros_like(z0arr)
    larr=np.ones_like(z0arr)*lmax

    if sigma_thr!=None :
        for i in np.arange(nnu) :
            z=zmarr[i]
            chi=pcs.radial_comoving_distance(1./(1+z))
            kmax=get_kmax(z,sigma_thr)
            ell=min(int(chi*kmax),lmax)
            print z,ell
            larr[i]=ell

    np.savetxt(prefix_out+".txt",np.transpose([z0arr,zfarr,szarr,larr]),
               fmt="%.3lE %.3lE %.3lE 0 0 %d",
               header='[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax')

nbins=[50,100,200]
for nb in nbins :
    mk_binfiles("bins_HIRAX_sthr0p75_n%d"%nb,400.,800.,nb,sigma_thr=0.75)
    mk_binfiles("bins_HIRAX_sthr0p50_n%d"%nb,400.,800.,nb,sigma_thr=0.50)
    mk_binfiles("bins_HIRAX_lmax2000_n%d"%nb,400.,800.,nb)

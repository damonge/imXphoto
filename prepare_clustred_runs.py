import numpy as np
import py_cosmo_mad as csm
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import experiments as xp
import os

nuHI=1420.405

#Save single-bin photo-z file
#Save multi-bin IM file

pcs=csm.PcsPar()
pcs.background_set(0.3,0.7,0.05,-1,0,0.7,2.7255)
pcs.set_linear_pk("EH",-4.,3.,0.01,0.96,0.8)

def get_dnu(z,dchi) :
    return nuHI*pcs.hubble(1/(1+z))*dchi/(1+z)**2

def get_lmax(z,st) :
    dz=pcs.growth_factor(1/(1+z))/pcs.growth_factor(1)
    chi=pcs.radial_comoving_distance(1/(1+z))
    def sig_minus_target(lk) :
        k=10.**lk
        return np.sqrt(pcs.sig0_L(1./k,1./k,"SharpK","SharpK"))*dz-st
    kmax=10.**brentq(sig_minus_target,-3.,1.)
    lmax=int(kmax*chi)

    return lmax

def write_param_file(xp_photo,xp_im,bins_photo,bz_photo,bins_im,bz_im,output_dir,parfile,fsky,
                     inc_align='no',inc_rsd='no',inc_mag='no',inc_gr='no') :
    stout= "#Cosmological parameters\n"
    stout+="[och2]\n"
    stout+="x= 0.1197\n"
    stout+="dx= 0.001\n"
    stout+="is_free= yes\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[obh2]\n"
    stout+="x= 0.02222\n"
    stout+="dx= 0.0001\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[hh]\n"
    stout+="x= 0.69\n"
    stout+="dx= 0.01\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[w0]\n"
    stout+="x= -1.0\n"
    stout+="dx= 0.01\n"
    stout+="is_free= no\n"
    stout+="onesided= 1\n"
    stout+="\n"
    stout+="[wa]\n"
    stout+="x= 0.0\n"
    stout+="dx= 0.05\n"
    stout+="is_free= no\n"
    stout+="onesided= 1\n"
    stout+="\n"
    stout+="[A_s]\n"
    stout+="x= 2.1955\n"
    stout+="dx= 0.01\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[ns]\n"
    stout+="x= 0.9655\n"
    stout+="dx= 0.005\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[tau]\n"
    stout+="x= 0.06\n"
    stout+="dx= 0.01\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[rt]\n"
    stout+="x= 0.00\n"
    stout+="dx= 0.005\n"
    stout+="is_free= no\n"
    stout+="onesided= 1\n"
    stout+="\n"
    stout+="[mnu]\n"
    stout+="x= 60.\n"
    stout+="dx= 10.\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[nnu]\n"
    stout+="x= 2.046\n"
    stout+="dx= 0.1\n"
    stout+="is_free= no\n"
    stout+="onesided= 0\n"
    stout+="\n"
    stout+="[pan]\n"
    stout+="x= 0.\n"
    stout+="dx= 0.02\n"
    stout+="is_free= no\n"
    stout+="onesided= 1\n"
    stout+="\n"
    stout+="[Tracer 1]\n"
    stout+="tracer_name= "+xp_photo['name']+"\n"
    stout+="tracer_type= gal_clustering\n"
    stout+="bins_file= "+bins_photo+"\n"
    stout+="nz_file= "+xp_photo['nzfi']+"\n"
    stout+="bias_file= "+bz_photo+"\n"
    stout+="sbias_file= "+xp_photo['szfi']+"\n"
    stout+="ebias_file= "+xp_photo['ezfi']+"\n"
    stout+="use_tracer= yes\n"
    stout+="\n"
    stout+="[Tracer 2]\n"
    stout+="tracer_name= "+xp_im['name']+"\n"
    stout+="tracer_type= intensity_mapping\n"
    stout+="bins_file= "+bins_im+"\n"
    stout+="nz_file= "+xp_im['nzfi']+"\n"
    stout+="bias_file= "+bz_im+"\n"
    stout+="sbias_file= "+xp_im['szfi']+"\n"
    stout+="ebias_file= "+xp_im['ezfi']+"\n"
    stout+="tz_file= "+xp_im['tzfi']+"\n"
    stout+="dish_size= %.3lf\n"%xp_im['dish_size']
    stout+="t_inst= %.3lf\n"%xp_im['t_inst']
    stout+="t_total= %.3lf\n"%xp_im['t_total']
    stout+="n_dish= %d\n"%xp_im['n_dish']
    stout+="area_efficiency= %.3lf\n"%xp_im['area_eff']
    stout+="fsky_im= %.3lf\n"%fsky
    stout+="is_single_dish= "+xp_im['is_sd']+"\n"
    stout+="base_file= "+xp_im['base_file']+"\n"
    stout+="use_tracer= yes\n"
    stout+="\n"
    stout+="[CLASS parameters]\n"
    stout+="lmax_cmb= 5000\n"
    stout+="lmax_lss= 2000\n"
    stout+="lmin_limber= 1000.\n"
    stout+="include_alignment= "+inc_align+"\n"
    stout+="include_rsd= "+inc_rsd+"\n"
    stout+="include_magnification= "+inc_mag+"\n"
    stout+="include_gr_vel= "+inc_gr+"\n"
    stout+="include_gr_pot= "+inc_gr+"\n"
    stout+="exec_path= ./class_mod\n"
    stout+="use_nonlinear= yes\n"
    stout+="f_sky= %.3lf\n"%fsky
    stout+="\n"
    stout+="[Output parameters]\n"
    stout+="output_dir= "+output_dir+"\n"
    stout+="output_spectra= run0\n"
    stout+="output_fisher= Fisher\n"
    stout+="\n"
    stout+="[Behaviour parameters]\n"
    stout+="model= LCDM\n"
    stout+="save_cl_files= yes\n"
    stout+="save_param_files= yes\n"

    f=open(parfile,"w")
    f.write(stout)
    f.close()

def prepare_bin(z0,zf,sz,lmx,predir,ibin,lmax_im=2000,sthr=None,fac_sigma=3) :
    fname_bin_single=predir+"/bins_photoz_b%d.txt"%ibin
    data_line=np.array([z0,zf,sz,lmx])
    np.savetxt(fname_bin_single,data_line[None,:],
               fmt="%lf %lf %lf 1 1 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

    z0_im=np.amax([z0-fac_sigma*sz,0])
    zf_im=zf+fac_sigma*sz
    nu0=nuHI/(1+zf_im); nuf=nuHI/(1+z0_im); dnu=get_dnu(0.5*(z0+zf),20.);
    nnu=2*3*fac_sigma; #nnu=int((nuf-nu0)/dnu);
    dnu=(nuf-nu0)/nnu;
    nuf_arr=nuf-(np.arange(nnu)+0)*dnu
    nu0_arr=nuf-(np.arange(nnu)+1)*dnu
    z0_arr=nuHI/nuf_arr-1; zf_arr=nuHI/nu0_arr-1; zm_arr=0.5*(z0_arr+zf_arr)
    lmax_arr=np.ones_like(z0_arr)*lmax_im;
    if sthr!=None :
        for i in np.arange(nnu) :
            lmax_arr[i]=np.amin([lmax_arr[i],get_lmax(zm_arr[i],sthr)])

    fname_im=predir+"/bins_im_b%d.txt"%ibin
    np.savetxt(fname_im,np.transpose([z0_arr,zf_arr,lmax_arr]),
               fmt="%lf %lf 0.0 0 0 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

def prepare_files(xp_photo,xp_im,bins_photo,sth_im,predir,fsky) :
    os.system('mkdir -p '+predir)
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        prepare_bin(z0_ph_arr[i],zf_ph_arr[i],sz_ph_arr[i],lmx_arr[i],predir,i,sthr=sth_im)
        zm=0.5*(z0_ph_arr[i]+zf_ph_arr[i])
        z,bz,ms=np.loadtxt(xp_photo['bzfi'],unpack=True); bzf=interp1d(z,bz);
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_photoz_b%d.txt"%i,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        z,bz,ms=np.loadtxt(xp_im['bzfi'],unpack=True); bzf=interp1d(z,bz);
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_im_b%d.txt"%i    ,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        write_param_file(xp_photo,xp_im,
                         predir+"/bins_photoz_b%d.txt"%i,predir+"/bz_photoz_b%d.txt"%i,
                         predir+"/bins_im_b%d.txt"%i,predir+"/bz_im_b%d.txt"%i,
                         predir+"/output_b%d"%i,predir+"/params_b%d.ini"%i,fsky)

prepare_files(xp.phoz_LSSTgold,xp.im_HIRAX_32_6,"curves_LSST/bins_gold_lmax2000.txt",None,"test_bins",0.4)

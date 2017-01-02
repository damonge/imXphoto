import numpy as np
import py_cosmo_mad as csm
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import experiments as xp
import os

JUST_CLS="no"
nuHI=1420.405
run_online=False
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

def write_param_file(xp_photo,xp_spec,bins_photo,bz_photo,bins_spec,bz_spec,output_dir,parfile,fsky,
                     inc_align='no',inc_rsd='no',inc_mag='no',inc_gr='no',fishname=None) :
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
    stout+="tracer_name= "+xp_spec['name']+"\n"
    stout+="tracer_type= gal_clustering\n"
    stout+="bins_file= "+bins_spec+"\n"
    stout+="nz_file= "+xp_spec['nzfi']+"\n"
    stout+="bias_file= "+bz_spec+"\n"
    stout+="sbias_file= "+xp_spec['szfi']+"\n"
    stout+="ebias_file= "+xp_spec['ezfi']+"\n"
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
    if run_online :
        stout+="exec_path= ./class_mod\n"
    else :
        stout+="exec_path= mpisubonepernodecmb \"GoFish\" 1x12 ./class_mod\n"
    stout+="use_nonlinear= yes\n"
    stout+="f_sky= %.3lf\n"%fsky
    stout+="\n"
    stout+="[Output parameters]\n"
    stout+="output_dir= "+output_dir+"\n"
    stout+="output_spectra= run0\n"
    if fishname==None :
        stout+="output_fisher= Fisher\n"
    else :
        stout+="output_fisher= "+fishname+"\n"
    stout+="\n"
    stout+="[Behaviour parameters]\n"
    stout+="model= LCDM\n"
    stout+="save_cl_files= yes\n"
    stout+="save_param_files= yes\n"
    stout+="just_run_cls= "+JUST_CLS+"\n"

    f=open(parfile,"w")
    f.write(stout)
    f.close()

def prepare_bin(z0,zf,sz,lmx,predir,ibin,lmax_spec=2000,sthr=None,fac_sigma=3) :
    fname_bin_single=predir+"/bins_photoz_b%d.txt"%ibin
    data_line=np.array([z0,zf,sz,lmx])
    np.savetxt(fname_bin_single,data_line[None,:],
               fmt="%lf %lf %lf 1 1 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

    z0_spec=np.amax([z0-fac_sigma*sz,0])
    zf_spec=zf+fac_sigma*sz
    nz=2*3*fac_sigma; #nnu=int((nuf-nu0)/dnu);
    dz=(zf_spec-z0_spec)/nz
    z0_arr=z0_spec+(np.arange(nz)+0)*dz
    zf_arr=z0_spec+(np.arange(nz)+1)*dz
    zm_arr=0.5*(z0_arr+zf_arr)
    lmax_arr=np.ones_like(z0_arr)*lmax_spec;
    if sthr!=None :
        for i in np.arange(nz) :
            lmax_arr[i]=np.amin([lmax_arr[i],get_lmax(zm_arr[i],sthr)])

    fname_spec=predir+"/bins_spec_b%d.txt"%ibin
    np.savetxt(fname_spec,np.transpose([z0_arr,zf_arr,lmax_arr]),
               fmt="%lf %lf 1.0E-5 0 0 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

def prepare_files(xp_photo,xp_spec,bins_photo,sth_spec,predir,fsky) :
    os.system('mkdir -p '+predir)
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        prepare_bin(z0_ph_arr[i],zf_ph_arr[i],sz_ph_arr[i],lmx_arr[i],predir,i,sthr=sth_spec)
        zm=0.5*(z0_ph_arr[i]+zf_ph_arr[i])
        z,bz,ms=np.loadtxt(xp_photo['bzfi'],unpack=True); bzf=interp1d(z,bz)
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_photoz_b%d.txt"%i,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        z,bz=np.loadtxt(xp_spec['bzfi'],unpack=True); bzf=interp1d(z,bz)
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_spec_b%d.txt"%i    ,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        write_param_file(xp_photo,xp_spec,
                         predir+"/bins_photoz_b%d.txt"%i,predir+"/bz_photoz_b%d.txt"%i,
                         predir+"/bins_spec_b%d.txt"%i,predir+"/bz_spec_b%d.txt"%i,
                         predir+"/output_b%d"%i,predir+"/params_b%d.ini"%i,fsky)
        os.system("python main.py "+predir+"/params_b%d.ini"%i)

def run_fsky(xp_photo,xp_spec,bins_photo,sth_spec,predir,fsky) :
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        parname=predir+"/params_fs%.3lf"%fsky+"_b%d.ini"%i
        write_param_file(xp_photo,xp_spec,
                         predir+"/bins_photoz_b%d.txt"%i,predir+"/bz_photoz_b%d.txt"%i,
                         predir+"/bins_spec_b%d.txt"%i,predir+"/bz_spec_b%d.txt"%i,
                         predir+"/output_b%d"%i,parname,fsky,fishname="Fisher_fs%.3lf"%fsky)
        os.system("python main.py "+parname)

prepare_files(xp.phoz_LSSTgold,xp.spec_DESI  ,"curves_LSST/bins_gold_lmax2000.txt",None,"DESI_fs0p20"  ,0.20)
prepare_files(xp.phoz_LSSTgold,xp.spec_Euclid,"curves_LSST/bins_gold_lmax2000.txt",None,"Euclid_fs0p40",0.40)
prepare_files(xp.phoz_LSSTgold,xp.spec_WFIRST,"curves_LSST/bins_gold_lmax2000.txt",None,"WFIRST_fs0p05",0.05)

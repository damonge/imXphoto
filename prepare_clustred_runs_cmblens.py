import numpy as np
import py_cosmo_mad as csm
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import experiments as xp
import os

nuHI=1420.405
run_cls=False
marg_A=False

def write_param_file(xp_photo,xp_cmbl,bins_photo,bz_photo,output_dir,parfile,fsky,fishname,
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
    if run_cls :
        stout+="is_free= yes\n"
    else :
        if marg_A :
            stout+="is_free= yes\n"
        else :
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
    stout+="is_photometric= yes\n"
    stout+="use_tracer= yes\n"
    stout+="\n"
    stout+="[Tracer 2]\n"
    stout+="tracer_name= CMBLens\n"
    stout+="tracer_type= cmb_lensing\n"
    stout+="sigma_t= %.3lf\n"%xp_cmbl['sigma_t']
    stout+="beam_amin= %.3lf\n"%xp_cmbl['beam_amin']
    stout+="lmin= %d\n"%xp_cmbl['lmin']
    stout+="lmax= 3000\n"
    if xp_cmbl['name']=='none' :
        stout+="use_tracer= no\n"
    else :
        stout+="use_tracer= yes\n"
    stout+="\n"
    stout+="[CLASS parameters]\n"
    stout+="lmax_cmb= 5000\n"
    stout+="lmax_lss= 3000\n"
    stout+="lmin_limber= 100.\n"
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
    stout+="output_fisher= "+fishname+"\n"
    stout+="\n"
    stout+="[Behaviour parameters]\n"
    stout+="model= LCDM\n"
    stout+="save_cl_files= yes\n"
    stout+="save_param_files= yes\n"
    if run_cls :
        stout+="just_run_cls= yes\n"
    else :
        stout+="just_run_cls= no\n"

    f=open(parfile,"w")
    f.write(stout)
    f.close()

def get_bin_fnames(predir,expname,ibin,sthr) :
    fisher_prefix="Fisher_"
    if marg_A :
        fisher_prefix+="wA_"
    else :
        fisher_prefix+="woA_"

    if sthr==None :
        fname_bin_single=predir+"/bins_photoz_b%d_lmax3000_l3000.txt"%ibin
        fname_params=predir+"/params_b%d_"%ibin+expname+"_lmax3000"
        fname_fisher=fisher_prefix+expname+"_lmax3000"
    else :
        fname_bin_single=predir+"/bins_photoz_b%d_"%ibin+"sthr%.3lf.txt"%sthr
        fname_params=predir+"/params_b%d_"%ibin+expname+"_sthr%.3lf"%sthr
        fname_fisher=fisher_prefix+expname+"_sthr%.3lf"%sthr

    return fname_bin_single,fname_params,fname_fisher

def prepare_bin(z0,zf,sz,lmx,predir,expname,ibin,sthr=None,fac_sigma=3,fac_sample_sigma=2) :
    fname_bin_single,fname_params,fname_fisher=get_bin_fnames(predir,expname,ibin,sthr)
    data_line=np.array([z0,zf,sz,lmx])
    np.savetxt(fname_bin_single,data_line[None,:],
               fmt="%lf %lf %lf 1 1 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

def prepare_files(xp_photo,xp_cmbl,bins_photo,sth_im,predir,fsky) :
    os.system('mkdir -p '+predir)
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        fname_bins_photo,fname_params,fname_fisher=get_bin_fnames(predir,xp_cmbl['name'],i,sth_im)
        prepare_bin(z0_ph_arr[i],zf_ph_arr[i],sz_ph_arr[i],lmx_arr[i],predir,xp_cmbl['name'],i,sthr=sth_im)
        zm=0.5*(z0_ph_arr[i]+zf_ph_arr[i])
        z,bz,ms=np.loadtxt(xp_photo['bzfi'],unpack=True); bzf=interp1d(z,bz);
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_photoz_b%d.txt"%i,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        write_param_file(xp_photo,xp_cmbl,
                         fname_bins_photo,predir+"/bz_photoz_b%d.txt"%i,
                         predir+"/output_b%d"%i,fname_params+".ini",fsky,fname_fisher)
        if run_cls :
#            os.system("addqueuecmb \"1h bin %d\" 1x12 /usr/local/shared/python/2.7.6-gcc/bin/python main.py "%i+fname_params+".ini")
            os.system("addqueue -q cmb -s -n 1x12 -m 1 -c \"1h bin %d\" /usr/local/shared/python/2.7.6-gcc/bin/python main.py "%i+fname_params+".ini")
#            os.system("python main.py "+fname_params+".ini")
        else :
            os.system("python main.py "+fname_params+".ini")

exper_array=[xp.cmb_S3_opt,xp.cmb_S4_opt,xp.cmb_S3,xp.cmb_S4,xp.cmb_none]
if run_cls :
    prepare_files(xp.phoz_LSSTgold,xp.cmb_S3_opt,"curves_LSST/bins_gold_lmax3000_l3000.txt",None,
                  "runs/CMBL",xp.cmb_S3_opt['fsky'])
else :
#    for exper in exper_array :
    for exper in exper_array :
        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_lmax3000_l3000.txt",None,
                      "runs/CMBL",exper['fsky'])
        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p50.txt",0.50,
                      "runs/CMBL",exper['fsky'])
        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p75.txt",0.75,
                      "runs/CMBL",exper['fsky'])
        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr1p00.txt",1.00,
                      "runs/CMBL",exper['fsky'])

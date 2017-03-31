import numpy as np
import py_cosmo_mad as csm
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import experiments as xp
import os

nuHI=1420.405
run_cls=False
include_fg=False
run_wedge=True
a_fg=1.0
xi_fg=10.0
marg_A=True

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

def mk_xp_generic(sigmaT,dmin,dmax,fsky,name) :
    xpr = {
        "name" : name,
        "nzfi" : "curves_IM/nz_HI.txt",
        "bzfi" : "curves_IM/bz_HI.txt",
        "szfi" : "curves_IM/sz_HI.txt",
        "ezfi" : "curves_IM/ez_HI.txt",
        "tzfi" : "curves_IM/tz_HI.txt",
        "dish_size" : 13.5,
        "t_inst" : sigmaT,
        "t_total" : 4000.,
        "n_dish" : 64,
        "area_eff" : 1.0,
        "im_type" : "generic",
        "base_file" : "none",
        "base_min" : dmin,
        "base_max" : dmax,
        "fsky" : fsky
        }

    return xpr;

def write_param_file(xp_photo,xp_im,bins_photo,bz_photo,bins_im,bz_im,output_dir,parfile,fsky,fishname,
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
    stout+="tracer_name= IM\n"
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
    stout+="instrument_type= "+xp_im['im_type']+"\n"
    stout+="base_file= "+xp_im['base_file']+"\n"
    stout+="baseline_min= %.3lf\n"%xp_im['base_min']
    stout+="baseline_max= %.3lf\n"%xp_im['base_max']
    if include_fg :
        stout+="include_foregrounds= yes\n"
        stout+="fit_foregrounds= no\n"
        stout+="A_fg= %.3lf\n"%a_fg
        stout+="alpha_fg=-2.7\n"
        stout+="beta_fg=-2.0\n"
        stout+="xi_fg=%.3lf\n"%xi_fg
        stout+="nux_fg=130.\n"
        stout+="lx_fg=1000.\n"
    elif run_wedge :
        stout+="include_foregrounds= yes\n"
        stout+="include_wedge= yes\n"
        stout+="fit_foregrounds= no\n"
        stout+="A_fg= %.3lf\n"%a_fg
        stout+="alpha_fg=-2.7\n"
        stout+="beta_fg=-2.0\n"
        stout+="xi_fg=10.\n"
        stout+="nux_fg=130.\n"
        stout+="lx_fg=1000.\n"
    else :
        stout+="include_foregrounds= no\n"
        stout+="fit_foregrounds= no\n"
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
    if include_fg :
        fisher_prefix+="wFG_a%.3lf_"%a_fg+"xi%.3lf_"%xi_fg
    elif run_wedge :
        fisher_prefix+="wFG_a%.3lf_"%a_fg+"wedge_"
    else :
        fisher_prefix+="woFG_"

    if sthr==None :
        fname_bin_single=predir+"/bins_photoz_b%d_lmax2000.txt"%ibin
        fname_im=predir+"/bins_im_b%d_lmax2000.txt"%ibin
        fname_params=predir+"/params_b%d_"%ibin+expname+"_lmax2000"
        fname_fisher=fisher_prefix+expname+"_lmax2000"
    else :
        fname_bin_single=predir+"/bins_photoz_b%d_"%ibin+"sthr%.3lf.txt"%sthr
        fname_im=predir+"/bins_im_b%d_"%ibin+"sthr%.3lf.txt"%sthr
        fname_params=predir+"/params_b%d_"%ibin+expname+"_sthr%.3lf"%sthr
        fname_fisher=fisher_prefix+expname+"_sthr%.3lf"%sthr

    return fname_bin_single,fname_im,fname_params,fname_fisher

def prepare_bin(z0,zf,sz,lmx,predir,expname,ibin,lmax_im=2000,sthr=None,fac_sigma=3,fac_sample_sigma=2) :
    fname_bin_single,fname_im,fname_params,fname_fisher=get_bin_fnames(predir,expname,ibin,sthr)
    data_line=np.array([z0,zf,sz,lmx])
    np.savetxt(fname_bin_single,data_line[None,:],
               fmt="%lf %lf %lf 1 1 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

    z0_im=np.amax([z0-fac_sigma*sz,0])
    zf_im=zf+fac_sigma*sz
    nz=int(fac_sample_sigma*(zf_im-z0_im)/(sz-1E-5)); dz=(zf_im-z0_im)/nz
    z0_arr=z0_im+(np.arange(nz)+0)*dz; zf_arr=z0_im+(np.arange(nz)+1)*dz; zm_arr=0.5*(z0_arr+zf_arr)
    lmax_arr=np.ones_like(z0_arr)*lmax_im;
    if sthr!=None :
        for i in np.arange(nz) :
            lmax_arr[i]=np.amin([lmax_arr[i],get_lmax(zm_arr[i],sthr)])

    np.savetxt(fname_im,np.transpose([z0_arr,zf_arr,lmax_arr]),
               fmt="%lf %lf 0.0 0 0 %d",
               header="[1]-z0 [2]-zf [3]-sz [4]-marg_sz [5]-marg_bz [6]-lmax")

def prepare_files(xp_photo,xp_im,bins_photo,sth_im,predir,fsky) :
    os.system('mkdir -p '+predir)
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        fname_bins_photo,fname_bins_im,fname_params,fname_fisher=get_bin_fnames(predir,xp_im['name'],i,sth_im)
        prepare_bin(z0_ph_arr[i],zf_ph_arr[i],sz_ph_arr[i],lmx_arr[i],predir,xp_im['name'],i,sthr=sth_im)
        zm=0.5*(z0_ph_arr[i]+zf_ph_arr[i])
        z,bz,ms=np.loadtxt(xp_photo['bzfi'],unpack=True); bzf=interp1d(z,bz);
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_photoz_b%d.txt"%i,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        z,bz,ms=np.loadtxt(xp_im['bzfi'],unpack=True); bzf=interp1d(z,bz);
        zarr=np.array([z[0],zm,z[-1]]); bzarr=bzf(zarr); mask=np.array([0,1,0])
        np.savetxt(predir+"/bz_im_b%d.txt"%i    ,np.transpose([zarr,bzarr,mask]),fmt='%lf %lf %d')
        write_param_file(xp_photo,xp_im,
                         fname_bins_photo,predir+"/bz_photoz_b%d.txt"%i,
                         fname_bins_im,predir+"/bz_im_b%d.txt"%i,
                         predir+"/output_b%d"%i,fname_params+".ini",fsky,fname_fisher)
        if run_cls :
            os.system("addqueuecmb \"1h bin %d\" 1x12 /usr/local/shared/python/2.7.6-gcc/bin/python main.py "%i+fname_params+".ini")
        else :
            os.system("python main.py "+fname_params+".ini")

def run_fsky(xp_photo,xp_im,bins_photo,sth_im,predir,fsky) :
    z0_ph_arr,zf_ph_arr,sz_ph_arr,dum1,dum2,lmx_arr=np.loadtxt(bins_photo,unpack=True)
    for i in np.arange(len(z0_ph_arr)) :
        fname_bins_photo,fname_bins_im,fname_params,fname_fisher=get_bin_fnames(predir,xp_im['name'],i,sth_im)
        parname=fname_params+"_fs%.3lf.ini"%fsky
        fishname=fname_fisher+"_fs%.3lf"%fsky
        write_param_file(xp_photo,xp_im,
                         fname_bins_photo,predir+"/bz_photoz_b%d.txt"%i,
                         fname_bins_im,predir+"/bz_im_b%d.txt"%i,
                         predir+"/output_b%d"%i,parname,fsky,fishname)
        os.system("python main.py "+parname)

exper_array=[xp.im_HIRAX_32_6,
             xp.im_SKA_SD,    xp.im_SKA_IF,    xp.im_SKA,
             xp.im_MeerKAT_SD,xp.im_MeerKAT_IF,xp.im_MeerKAT]
#exper_array=[xp.im_SKA_SD,xp.im_SKA_IF,xp.im_SKA,xp.im_MeerKAT_SD,xp.im_MeerKAT_IF,xp.im_MeerKAT_IF]
#exper_array=[xp.im_HIRAX_32_6]
#fsky_arr=[0.01,0.02,0.05,0.1,0.2,0.4]
if run_cls :
    prepare_files(xp.phoz_LSSTgold,xp.im_HIRAX_32_6,"curves_LSST/bins_gold_lmax2000.txt",None,
                  "runs/IMAP",xp.im_HIRAX_32_6['fsky'])
else :
#    for exper in exper_array :
    for exper in [xp.im_HIRAX_32_6,xp.im_SKA,xp.im_MeerKAT] :
#        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_lmax2000.txt",None,
#                      "runs/IMAP",exper['fsky'])
#        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p50.txt",0.50,
#                      "runs/IMAP",exper['fsky'])
#        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p75.txt",0.75,
#                      "runs/IMAP",exper['fsky'])
        prepare_files(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr1p00.txt",1.00,
                      "runs/IMAP",exper['fsky'])
#        for fs in fsky_arr :
#            run_fsky(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_lmax2000.txt",None,"runs/IMAP",fs)
#            run_fsky(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p50.txt",0.50,"runs/IMAP",fs)
#            run_fsky(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr0p75.txt",0.75,"runs/IMAP",fs)
#            run_fsky(xp.phoz_LSSTgold,exper,"curves_LSST/bins_gold_sthr1p00.txt",1.00,"runs/IMAP",fs)
#
#    #Study noise level
#    for texp in [1.0,1.5,2.0,2.5,3.0,3.5] :
#        tinst=10.**(-texp); name="gen_sT%.3lf"%texp
#        xpr=mk_xp_generic(tinst,6.,180.,0.4,name)
#        prepare_files(xp.phoz_LSSTgold,xpr,"curves_LSST/bins_gold_sthr1p00.txt",1.00,
#                      "runs/IMAP",xpr['fsky'])
#
#    #Study minimum baseline for interferometers
#    for dmin in [0.,1.5,3.,6.,12.,24.,48.] :
#        name="gen_dmn%.3lf"%dmin
#        xpr=mk_xp_generic(1E-3,dmin,1000.,0.4,name)
#        prepare_files(xp.phoz_LSSTgold,xpr,"curves_LSST/bins_gold_sthr1p00.txt",1.00,
#                      "runs/IMAP",xpr['fsky'])
#
#    #Study maximum baseline for single-dish
#    for dmax in [7.5,15.,30.,60.,120.,240.] :
#        name="gen_dmx%.3lf"%dmax
#        xpr=mk_xp_generic(1E-3,0,dmax,0.4,name)
#        prepare_files(xp.phoz_LSSTgold,xpr,"curves_LSST/bins_gold_sthr1p00.txt",1.00,
#                      "runs/IMAP",xpr['fsky'])

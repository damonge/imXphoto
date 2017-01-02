import numpy as np
import matplotlib.pyplot as plt
import fishplot as fshr

names=np.array(['och2','bias_gal','sphz','bphz','bias_im'])
nptot=len(names)

def find_param(name) :
    for i in np.arange(nptot) :
        if names[i]==name :
            return i

def read_fish(prefix,lmin=-1,lmax=-1) :
    fish_full_l=np.load(prefix+"/fisher_raw_l.npy")

    if lmin<0 :
        lmn=2
    else :
        lmn=lmin
    if lmax<0 :
        lmx=-1
    else :
        lmx=lmax

    fish_full=np.sum(fish_full_l[:,:,lmn:lmx],axis=2)
    
#    errors=np.sqrt(np.diag(np.linalg.inv(fish_full)))
#    print prefix
#    for i in np.arange(nptot) :
#        print names[i],errors[i]
#    print ""
    return fish_full

vpar={'och2'    :0.1197,
      'bias_gal':1.0,
      'sphz'    :0.03,
      'bphz'    :0.00,
      'bias_im' :1.0
      }
rpar={'och2'    :None,
      'bias_gal':None,
      'sphz'    :None,
      'bphz'    :None,
      'bias_im' :None
      }
ppar={'och2'    :-1,
      'bias_gal':-1,
      'sphz'    :-1,
      'bphz'    :-1,
      'bias_im' :-1
      }
lpar={'och2'    :"$\\omega_c$",
      'bias_gal':"$b_{\\rm gal}$",
      'sphz'    :"$\\sigma_z$",
      'bphz'    :"$\\Delta z$",
      'bias_im' :"$b_{\\rm IM}$"
      }
plpar={'och2'    :True,
       'bias_gal':True,
       'sphz'    :True,
       'bphz'    :True,
       'bias_im' :True
       }
mpar={'och2'    :True,
      'bias_gal':True,
      'sphz'    :True,
      'bphz'    :True,
      'bias_im' :True
      }

#for fs in [0.05,0.1,0.2,0.4] :
#    fish=read_fish("HIRAX_32_6/output_b4/Fisher_fs%.3lf"%fs)

def compute_errors(dirname,col,z_min,z_max,label,fishname="Fisher") :
    i_sthz=find_param('sphz')
    i_bthz=find_param('bphz')
    z0,zf,dum,dum,dum,dum=np.loadtxt("curves_LSST/bins_gold_lmax2000.txt",unpack=True);
    zarr=0.5*(z0+zf)
    ind=np.where((zarr>=z_min) & (zarr<=z_max))[0]
    err_sthz=np.zeros_like(zarr)
    err_bthz=np.zeros_like(zarr)
    for b in np.arange(len(z0)) :
        fish=read_fish(dirname+"/output_b%d/"%b+fishname)
        error=np.sqrt(np.diag(np.linalg.inv(fish)))
        err_sthz[b]=error[i_sthz]
        err_bthz[b]=error[i_bthz]
    plt.plot(zarr[ind],err_sthz[ind],col+'-',lw=2,label=label)
    plt.plot(zarr[ind],err_bthz[ind],col+'--',lw=2)
    return zarr,err_sthz,err_bthz

plt.figure()
zarr_DESI  ,err_sthz_DESI_arr  ,err_bthz_DESI_arr  =compute_errors("DESI_fs0p20"  ,'r',
                                                                   0.00,1.85,label="DESI")
zarr_Euclid,err_sthz_Euclid_arr,err_bthz_Euclid_arr=compute_errors("Euclid_fs0p40",'g',
                                                                   0.65,2.15,label="Euclid")
zarr_WFIRST,err_sthz_WFIRST_arr,err_bthz_WFIRST_arr=compute_errors("WFIRST_fs0p05",'b',
                                                                   0.95,2.85,"WFIRST")
zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("HIRAX_32_6"   ,'y',
                                                                   0.00,3.00,"HIRAX-6m")
zarr_SKASD ,err_sthz_SKASD_arr ,err_bthz_SKADS_arr =compute_errors("SKA_SD"       ,'m',
                                                                   0.00,3.00,"SKA single dish")
plt.xlim([0,2.7])
plt.ylim([5E-5,1E-2])#6E-3])
plt.yscale('log')
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$\\sigma(\\sigma_z,\\Delta z)$',fontsize=16)
plt.plot(zarr_SKASD,1E-3*(1+zarr_SKASD),'k-',label='LSST req.')
plt.plot([-1,-1],[-1,-2],'k-',lw=2,label="$\\sigma_z$")
plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\Delta z$")
plt.legend(loc='upper left',frameon=False,labelspacing=0.1,ncol=2)
plt.savefig("photoz_calib_compare.pdf",bbox_inches='tight')

plt.figure()
zarr_HIRAXa,err_sthz_HIRAXa_arr,err_bthz_HIRAXa_arr=compute_errors("HIRAX_32_6",'r',0.00,3.00,
                                                                   "HIRAX-6m, $\\ell<2000$",
                                                                   fishname="Fisher")
zarr_SKASDa,err_sthz_SKASDa_arr,err_bthz_SKADSa_arr=compute_errors("SKA_SD"    ,'g',0.00,3.00,
                                                                   "SKA single dish, $\\ell<2000$",
                                                                   fishname="Fisher")
zarr_HIRAXb,err_sthz_HIRAXb_arr,err_bthz_HIRAXb_arr=compute_errors("HIRAX_32_6",'b',0.00,3.00,
                                                                   "HIRAX-6m, $\\sigma_{\\rm thr}=0.75$",
                                                                   fishname="Fisher_sthr0p75")
zarr_SKASDb,err_sthz_SKASDb_arr,err_bthz_SKADSb_arr=compute_errors("SKA_SD"    ,'y',0.00,3.00,
                                                                   "SKA single dish, $\\sigma_{\\rm thr}=0.75$",
                                                                   fishname="Fisher_sthr0p75")
plt.xlim([0,2.7])
plt.ylim([5E-5,1E-2])#6E-3])
plt.yscale('log')
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$\\sigma(\\sigma_z,\\Delta z)$',fontsize=16)
plt.plot(zarr_SKASDa,1E-3*(1+zarr_SKASDa),'k-',label='LSST req.')
plt.plot([-1,-1],[-1,-2],'k-',lw=2,label="$\\sigma_z$")
plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\Delta z$")
plt.legend(loc='upper left',frameon=False,labelspacing=0.1,ncol=2)
plt.show()

exit(1)
fshp=[]
for ip in np.arange(nptot) :
    pname=names[ip]
    fshp.append(fshr.ParamFisher(vpar[pname],ppar[pname],pname,lpar[pname],True,plpar[pname],ppar[pname]))
fshp=np.array(fshp)

fshr.plot_fisher_all(fshp,[fish],['none'],[2],['solid'],['red'],['Default'],2,rpar,'fishtest.pdf',do_1D=False,fsize=14)
plt.show()

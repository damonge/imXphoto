import numpy as np
import matplotlib.pyplot as plt
import fishplot as fshr
import os
import sys
import py_cosmo_mad as csm
import experiments as xpr
from scipy.interpolate import interp1d
from scipy.integrate import quad

NU_21=1420.405751786
CLIGHT=299.792458
FWHM2G=0.42466090014
col_HIRAX="#004C99"
col_SKA="#CC0000"
col_MKT="#D8BB00"

if len(sys.argv)!=2 :
    print "Usage plot_utils.py figname"
    exit(1)
whichfig=sys.argv[1]

def read_fisher(prefix,integrate_l=False,lmin=-1,lmax=-1) :
    data=np.load(prefix+"/fisher_raw.npz")

    if integrate_l :
        if lmin>0 : lmn=lmin
        else : lmn=0
        if lmax>0 : lmx=lmax
        else : lmx=-1

        fish_cls=np.sum(data['fisher_l'][:,:,lmn:lmx])
        fish_tot=fish_cls+data['fisher_bao']
    else :
        fish_tot=data['fisher_tot']
    
    fish_package={'names' :data['names'],
                  'values':data['values'],
                  'labels':data['labels'],
                  'fisher':fish_tot}
    return fish_package

def get_sigma(fishp,names) :
    sigmarr=np.sqrt(np.diag(np.linalg.inv(fishp['fisher'])))
    result=[]
    for name in names :
        i=np.where(fishp['names']==name)[0]
        result.append(sigmarr[i])

    return np.array(result)

def compute_errors(dirname,fishname,col,z_min,z_max,label,plot_bphz=True,plot_sphz=True,lt='-') :
    z0,zf,dum,dum,dum,dum=np.loadtxt("curves_LSST/bins_gold_lmax2000.txt",unpack=True);
    zarr=0.5*(z0+zf)
    ind=np.where((zarr>=z_min) & (zarr<=z_max))[0]
    err_sthz=np.zeros_like(zarr)
    err_bthz=np.zeros_like(zarr)
    for b in np.arange(len(z0)) :
        fish=read_fisher(dirname+"/output_b%d/"%b+fishname)
        errors=get_sigma(fish,['bphz_LSST_gold_node0','sphz_LSST_gold_node0'])
        err_bthz[b]=errors[0]
        err_sthz[b]=errors[1]
    if plot_bphz :
        if plot_sphz :
            plt.plot(zarr[ind],err_bthz[ind],'--',color=col,lw=2)
        else :
            plt.plot(zarr[ind],err_bthz[ind],lt,color=col,lw=2,label=label)
    if plot_sphz :
        plt.plot(zarr[ind],err_sthz[ind],lt,color=col,lw=2,label=label)
    return zarr,err_sthz,err_bthz

if whichfig=='fig1' or whichfig=='all':
    plt.figure()
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,"HIRAX-6m",plot_bphz=False)
    zarr_SKAFL ,err_sthz_SKAFL_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA",plot_bphz=False)
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT",plot_bphz=False)
    zarr_DESI  ,err_sthz_DESI_arr  ,err_bthz_DESI_arr  =compute_errors("runs/DESI"  ,"Fisher_woA_sthr1.000"                ,'#009900',0.00,1.85,label="DESI",plot_bphz=False,lt='--')
    zarr_Euclid,err_sthz_Euclid_arr,err_bthz_Euclid_arr=compute_errors("runs/Euclid","Fisher_woA_sthr1.000"                ,'#CC00CC',0.65,2.15,label="Euclid",plot_bphz=False,lt='--')
    zarr_WFIRST,err_sthz_WFIRST_arr,err_bthz_WFIRST_arr=compute_errors("runs/WFIRST","Fisher_woA_sthr1.000"                ,'#00CCCC',0.95,2.85,"WFIRST",plot_bphz=False,lt='--')
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=16)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=16)
    plt.plot(zarr_HIRAX,1E-3*(1+zarr_HIRAX),'k-',label='LSST req.')
#plt.plot([-1,-1],[-1,-2],'k-',lw=2,label="$\\sigma_z$")
#plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\Delta z$")
    plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=3)
    plt.savefig("compare_spec.pdf",bbox_inches='tight')

if whichfig=='fig2' or whichfig=='all':
    plt.figure()
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,"HIRAX-6m")
    zarr_SKAFL ,err_sthz_SKAFL_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA full")
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT full")
    plt.plot([-1,-1],[-1,-2],'k-',lw=2,label="$\\sigma_z$")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\Delta z$")
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=16)
    plt.ylabel('$\\sigma(\\sigma_z,\\Delta z)$',fontsize=16)
    plt.plot(zarr_HIRAX,1E-3*(1+zarr_HIRAX),'k-',label='LSST req.')
    plt.legend(loc='lower right',frameon=False,labelspacing=0.1,ncol=2)
    plt.savefig("compare_wbias.pdf",bbox_inches='tight')

if whichfig=='fig3' or whichfig=='all':
    plt.figure()
    zarr_SKAFL ,err_sthz_SKAFL_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA full",plot_bphz=False)
    zarr_SKAFL ,err_sthz_SKASD_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_SD_sthr1.000"    ,col_SKA,0.00,3.00,"SKA single dish",lt='--',plot_bphz=False)
    zarr_SKAFL ,err_sthz_SKAIF_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_IF_sthr1.000"    ,col_SKA,0.00,3.00,"SKA interferometer",lt='-.',plot_bphz=False)
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT full",plot_bphz=False)
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_SD_sthr1.000",col_MKT,0.00,3.00,"MeerKAT single dish",lt='--',plot_bphz=False)
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_IF_sthr1.000",col_MKT,0.00,3.00,"MeerKAT interferometer",lt='-.',plot_bphz=False)
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=16)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=16)
    plt.plot(zarr_SKAFL,1E-3*(1+zarr_SKAFL),'k-')#,label='LSST req.')
    plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=2)
    plt.savefig("compare_if_sd.pdf",bbox_inches='tight')

if whichfig=='fig4' or whichfig=='all':
    plt.figure()
    plt.plot([-1,-1],[-1,-2],'k:' ,lw=2,label="$\\sigma_{\\rm thr}=0.5$")
    plt.plot([-1,-1],[-1,-2],'k-.',lw=2,label="$\\sigma_{\\rm thr}=0.75$")
    plt.plot([-1,-1],[-1,-2],'k-' ,lw=2,label="$\\sigma_{\\rm thr}=1$")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\ell\\leq 2000$")
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_sthr0.500",col_HIRAX,0.00,3.00,None      ,lt=':',plot_bphz=False)
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_sthr0.750",col_HIRAX,0.00,3.00,None      ,lt='-.',plot_bphz=False)
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,"HIRAX-6m",lt='-',plot_bphz=False)
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_HIRAX_32_6_lmax2000" ,col_HIRAX,0.00,3.00,None      ,lt='--',plot_bphz=False)
    zarr_SKA ,err_sthz_SKA_arr ,err_bthz_SKA_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr0.500",col_SKA,0.00,3.00,None      ,lt=':',plot_bphz=False)
    zarr_SKA ,err_sthz_SKA_arr ,err_bthz_SKA_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr0.750",col_SKA,0.00,3.00,None      ,lt='-.',plot_bphz=False)
    zarr_SKA ,err_sthz_SKA_arr ,err_bthz_SKA_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_sthr1.000",col_SKA,0.00,3.00,"SKA"     ,lt='-',plot_bphz=False)
    zarr_SKA ,err_sthz_SKA_arr ,err_bthz_SKA_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_SKA_lmax2000" ,col_SKA,0.00,3.00,None      ,lt='--',plot_bphz=False)
    zarr_MeerKAT ,err_sthz_MeerKAT_arr ,err_bthz_MeerKAT_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr0.500",col_MKT,0.00,3.00,None      ,lt=':',plot_bphz=False)
    zarr_MeerKAT ,err_sthz_MeerKAT_arr ,err_bthz_MeerKAT_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr0.750",col_MKT,0.00,3.00,None      ,lt='-.',plot_bphz=False)
    zarr_MeerKAT ,err_sthz_MeerKAT_arr ,err_bthz_MeerKAT_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_sthr1.000",col_MKT,0.00,3.00,"MeerKAT"     ,lt='-',plot_bphz=False)
    zarr_MeerKAT ,err_sthz_MeerKAT_arr ,err_bthz_MeerKAT_arr =compute_errors("runs/IMAP"  ,"Fisher_woA_woFG_MeerKAT_lmax2000" ,col_MKT,0.00,3.00,None      ,lt='--',plot_bphz=False)
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=16)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=16)
    plt.plot(zarr_HIRAX,1E-3*(1+zarr_HIRAX),'k-')#,label='LSST req.')
    plt.legend(loc='lower right',frameon=False,labelspacing=0.1,ncol=2)
    plt.savefig("compare_nlin.pdf",bbox_inches='tight')

if whichfig=='fig5' or whichfig=='all' :
    nu0=650.
    z0=NU_21/nu0-1
    pcs=csm.PcsPar()
    pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.7255)
    lmax=10000

    def get_noisepower_imap(fname_tz,xp) :
        z,tz=np.loadtxt(fname_tz,unpack=True)
        tofz=interp1d(z,tz); tbg=tofz(z0);
        tsys=(xp['t_inst']+60.*(nu0/300.)**(-2.5))*1000
        sigma2_noise=(tsys/tbg)**2
        sigma2_noise*=4*np.pi*xp['fsky']/(3.6E9*xp['t_total']*NU_21)

        beam_fwhm=CLIGHT/(xp['dish_size']*nu0)
        l_arr=np.arange(lmax)
        if xp['im_type']=='interferometer' :
            lam0=CLIGHT/nu0
            dist,nbase=np.loadtxt(xp['base_file'],unpack=True)
            ndistint=interp1d(dist,nbase*dist*2*np.pi,bounds_error=False,fill_value=0.)
            norm=0.5*xp['n_dish']*(xp['n_dish']-1.)/quad(ndistint,dist[0],dist[-1])[0]
            nbase*=norm; ndist=interp1d(dist,nbase,bounds_error=False,fill_value=0.)
            n_baselines=ndist(l_arr*lam0/(2*np.pi))
            factor_beam_if=n_baselines*(lam0/beam_fwhm)**2
        else :
            factor_beam_if=1E-16*np.ones_like(l_arr)
        if xp['im_type']=='single_dish' :
            beam_rad=beam_fwhm*FWHM2G
            factor_beam_sd=xp['n_dish']*np.exp(-l_arr*(l_arr+1.)*beam_rad**2)
        else :
            factor_beam_sd=1E-16*np.ones_like(l_arr)
        factor_beam=factor_beam_sd+factor_beam_if
        
        cl_noise=sigma2_noise/factor_beam
        chi=pcs.radial_comoving_distance(1./(1+z0))
        pk_noise=cl_noise*(chi*(1+z0))**2/pcs.hubble(1./(1+z0))
        k_arr=l_arr/chi
        return k_arr,pk_noise

    def get_noisepower_spec(xp,zx) :
        z,nz=np.loadtxt(xp['nzfi'],unpack=True)
        nofz=interp1d(z,nz,bounds_error=False,fill_value=0)
        nzval=nofz(zx)*(180.*60/np.pi)**2;
        chi=pcs.radial_comoving_distance(1./(1+zx))
        k_arr=np.arange(lmax)/chi
        pk_noise=chi**2/(nzval*pcs.hubble(1./(1+zx)))*np.ones_like(k_arr)
        return k_arr,pk_noise
    
    k_HIRAX ,pk_HIRAX =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_HIRAX_32_6)
    k_SKASD ,pk_SKASD =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_SKA_SD)
    k_SKAIF ,pk_SKAIF =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_SKA_IF)
    k_MKTSD ,pk_MKTSD =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_MeerKAT_SD)
    k_MKTIF ,pk_MKTIF =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_MeerKAT_IF)
    k_DESI  ,pk_DESI  =get_noisepower_spec(xpr.spec_DESI  ,z0)
    k_Euclid,pk_Euclid=get_noisepower_spec(xpr.spec_Euclid,z0)
    k_WFIRST,pk_WFIRST=get_noisepower_spec(xpr.spec_WFIRST,z0)
    plt.figure()
    plt.plot(k_HIRAX ,pk_HIRAX ,'-' ,color=col_HIRAX,lw=2,label='HIRAX')
    plt.plot(k_SKASD ,pk_SKASD ,'-' ,color=col_SKA  ,lw=2,label='SKA (S.D.)')
    plt.plot(k_SKAIF ,pk_SKAIF ,'--',color=col_SKA  ,lw=2,label='SKA (interf.)')
    plt.plot(k_MKTSD ,pk_MKTSD ,'-' ,color=col_MKT  ,lw=2,label='MeerKAT (S.D.)')
    plt.plot(k_MKTIF ,pk_MKTIF ,'--',color=col_MKT  ,lw=2,label='MeerKAT (interf.)')
    plt.plot(k_DESI  ,pk_DESI  ,'k-' ,lw=2,label='DESI')
    plt.plot(k_Euclid,pk_Euclid,'k--',lw=2,label='Euclid')
    plt.plot(k_WFIRST,pk_WFIRST,'k-.',lw=2,label='WFIRST')
    plt.legend(loc='upper left',frameon=False,labelspacing=0.1)
    plt.loglog()
    plt.xlim([1E-3,3])
    plt.ylim([1E2,1E6])
plt.show()
exit(1)

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



for i in np.arange(15) :
    f1=read_fisher("runs/IMAP/output_b%d/Fisher_woA_woFG_HIRAX_32_6_lmax2000"%i)
    f2=read_fisher("runs/IMAP/output_b%d/Fisher_wA_woFG_HIRAX_32_6_lmax2000"%i)
    print i, get_sigma(f1,['sphz_LSST_gold_node0'])/get_sigma(f2,['sphz_LSST_gold_node0'])

exit(1)

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

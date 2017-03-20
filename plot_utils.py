import numpy as np
import matplotlib.pyplot as plt
import fishplot as fshr
import os
import sys
import py_cosmo_mad as csm
import experiments as xpr
from scipy.interpolate import interp1d
from scipy.integrate import quad
from mpl_toolkits.axes_grid1.inset_locator import mark_inset,zoomed_inset_axes
import copy

NU_21=1420.405751786
CLIGHT=299.792458
FWHM2G=0.42466090014
col_HIRAX="#004C99"
col_SKA="#CC0000"
col_MKT="#D8BB00"

zreq=3.0*np.arange(256)/255.

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

def get_sigma(fishp,names,fix_params=[]) :
    fish=(fishp['fisher']).copy()
    for name in fix_params :
        i=np.where(fishp['names']==name)[0]
        fish[i,i]+=1E16**2
    sigmarr=np.sqrt(np.diag(np.linalg.inv(fish)))
    result=[]
    for name in names :
        i=np.where(fishp['names']==name)[0]
        result.append(sigmarr[i])

    return np.array(result)

def compute_errors(dirname,fishname,col,z_min,z_max,label,plot_bphz=True,plot_sphz=True,lt='-',
                   l5000=False,fix_sphz=False,do_plot=True) :
    z0,zf,dum,dum,dum,dum=np.loadtxt("curves_LSST/bins_gold_lmax2000.txt",unpack=True);
    zarr=0.5*(z0+zf)
    ind=np.where((zarr>=z_min) & (zarr<=z_max))[0]
    err_sthz=np.zeros_like(zarr)
    err_bthz=np.zeros_like(zarr)
    if fix_sphz :
        plot_sphz=False
        namefix=['sphz_LSST_gold_node0','och2']
    else :
        namefix=[]
    for b in np.arange(len(z0)) :
        if l5000 :
            fish=read_fisher(dirname+"/output_l5000_b%d/"%b+fishname)
        else :
            fish=read_fisher(dirname+"/output_b%d/"%b+fishname)
        errors=get_sigma(fish,['bphz_LSST_gold_node0','sphz_LSST_gold_node0'],fix_params=namefix)
        err_bthz[b]=errors[0]
        err_sthz[b]=errors[1]
    if do_plot :
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
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,"HIRAX",plot_bphz=False)
    zarr_SKAFL ,err_sthz_SKAFL_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA",plot_bphz=False)
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT",plot_bphz=False)
    zarr_DESI  ,err_sthz_DESI_arr  ,err_bthz_DESI_arr  =compute_errors("runs/DESI"  ,"Fisher_wA_sthr1.000"                ,'#000000',0.00,1.85,label="DESI",plot_bphz=False,lt='-')
    zarr_Euclid,err_sthz_Euclid_arr,err_bthz_Euclid_arr=compute_errors("runs/Euclid","Fisher_wA_sthr1.000"                ,'#000000',0.65,2.15,label="Euclid",plot_bphz=False,lt='--')
    zarr_WFIRST,err_sthz_WFIRST_arr,err_bthz_WFIRST_arr=compute_errors("runs/WFIRST","Fisher_wA_sthr1.000"                ,'#000000',0.95,2.85,"WFIRST",plot_bphz=False,lt='-.')
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.plot(zreq,1E-3*(1+zreq),'k-',label='LSST req.')
    plt.plot(zreq,2E-3*(1+zreq),'k--',label='DES req.')
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=3,fontsize=16)
    plt.savefig("bak/compare_spec.pdf",bbox_inches='tight')

if whichfig=='fig2' or whichfig=='all':
    plt.figure()
    zarr_HIRAX ,err_sthz_HIRAX_arr ,err_bthz_HIRAX_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,"HIRAX")
    zarr_SKAFL ,err_sthz_SKAFL_arr ,err_bthz_SKAFL_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA")
    zarr_MKTFL ,err_sthz_MKTFL_arr ,err_bthz_MKTFL_arr =compute_errors("runs/IMAP"  ,"Fisher_wA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT")
    plt.plot([-1,-1],[-1,-2],'k-',lw=2,label="$\\sigma_z$")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\Delta z$")
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z,\\Delta z)$',fontsize=18)
    plt.plot(zreq,1E-3*(1+zreq),'k-',label='LSST req.')
    plt.legend(loc='lower right',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_wbias.pdf",bbox_inches='tight')

if whichfig=='fig3' or whichfig=='all':
    plt.figure()
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_sthr1.000"       ,col_SKA,0.00,3.00,"SKA full",plot_bphz=False)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_SD_sthr1.000"    ,col_SKA,0.00,3.00,"SKA sing. dish",lt='--',plot_bphz=False)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_IF_sthr1.000"    ,col_SKA,0.00,3.00,"SKA interf.",lt='-.',plot_bphz=False)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr1.000"   ,col_MKT,0.00,3.00,"MeerKAT full",plot_bphz=False)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_SD_sthr1.000",col_MKT,0.00,3.00,"MeerKAT sing. dish",lt='--',plot_bphz=False)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_IF_sthr1.000",col_MKT,0.00,3.00,"MeerKAT interf.",lt='-.',plot_bphz=False)
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_if_sd.pdf",bbox_inches='tight')

if whichfig=='fig4' or whichfig=='all':
    plt.figure()
    plt.plot([-1,-1],[-1,-2],'k:' ,lw=2,label="$\\sigma_{\\rm thr}=0.5$")
    plt.plot([-1,-1],[-1,-2],'k-.',lw=2,label="$\\sigma_{\\rm thr}=0.75$")
    plt.plot([-1,-1],[-1,-2],'k-' ,lw=2,label="$\\sigma_{\\rm thr}=1$")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\ell\\leq 5000$")
    compute_errors("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_sthr0.500_l5000",col_HIRAX,0.00,3.00,None,lt=':',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_sthr0.750_l5000",col_HIRAX,0.00,3.00,None,lt='-.',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_sthr1.000_l5000",col_HIRAX,0.00,3.00,"HIRAX",lt='-',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_lmax5000_l5000" ,col_HIRAX,0.00,3.00,None,lt='--',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_sthr0.500_l5000",col_SKA,0.00,3.00,None,lt=':',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_sthr0.750_l5000",col_SKA,0.00,3.00,None,lt='-.',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_sthr1.000_l5000",col_SKA,0.00,3.00,"SKA",lt='-',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_lmax5000_l5000" ,col_SKA,0.00,3.00,None,lt='--',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr0.500_l5000",col_MKT,0.00,3.00,None,lt=':',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr0.750_l5000",col_MKT,0.00,3.00,None,lt='-.',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr1.000_l5000",col_MKT,0.00,3.00,"MeerKAT",lt='-',plot_bphz=False,l5000=True)
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_lmax5000_l5000" ,col_MKT,0.00,3.00,None,lt='--',plot_bphz=False,l5000=True)
    plt.xlim([0,2.7])
    plt.ylim([5E-5,1E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='lower right',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_nlin.pdf",bbox_inches='tight')

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
    plt.xlabel('$k_\\perp\\,[h\\,{\\rm Mpc}^{-1}]$',fontsize=18)
    plt.ylabel('$N(k_\\perp)\\,\\,[{\\rm Mpc}\\,h^{-1}]^3$',fontsize=18)
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.xlim([1E-3,3])
    plt.ylim([1E2,1E6])
    plt.savefig("bak/nk.pdf",bbox_inches='tight')

if whichfig=='test_sigmaT' :
    nu0=650.
    z0=NU_21/nu0-1
    pcs=csm.PcsPar()
    pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.7255)
    lmax=10000

    def get_noisepower_imap(fname_tz,xp) :
        tinst=xp['t_inst']*1000
        sigma2_noise=tinst**2*4*np.pi*xp['fsky']/(3.6E9*xp['t_total'])
        lam0=CLIGHT/nu0

        beam_fwhm=CLIGHT/(xp['dish_size']*nu0)
        l_arr=np.arange(lmax)
        if xp['im_type']=='interferometer' :
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
        d_arr=l_arr*lam0/(2*np.pi)
        return d_arr,np.sqrt(cl_noise)
    
    k_HIRAX ,pk_HIRAX =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_HIRAX_32_6)
    k_SKASD ,pk_SKASD =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_SKA_SD)
    k_SKAIF ,pk_SKAIF =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_SKA_IF)
    k_MKTSD ,pk_MKTSD =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_MeerKAT_SD)
    k_MKTIF ,pk_MKTIF =get_noisepower_imap("curves_IM/tz_HI.txt",xpr.im_MeerKAT_IF)
    plt.figure()
    plt.plot(k_HIRAX ,pk_HIRAX ,'-' ,color=col_HIRAX,lw=2,label='HIRAX')
    plt.plot(k_SKASD ,pk_SKASD ,'-' ,color=col_SKA  ,lw=2,label='SKA (S.D.)')
    plt.plot(k_SKAIF ,pk_SKAIF ,'--',color=col_SKA  ,lw=2,label='SKA (interf.)')
    plt.plot(k_MKTSD ,pk_MKTSD ,'-' ,color=col_MKT  ,lw=2,label='MeerKAT (S.D.)')
    plt.plot(k_MKTIF ,pk_MKTIF ,'--',color=col_MKT  ,lw=2,label='MeerKAT (interf.)')
    plt.legend(loc='upper left',frameon=False,labelspacing=0.1)
    plt.loglog()

if whichfig=='fig6' or whichfig=='all':
    plt.figure()
    plt.plot([-1,-1],[-1,-2],'k-' ,lw=2,label="${\\rm No\\,foregrounds}$")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$\\xi=10$")
    plt.plot([-1,-1],[-1,-2],'k-.',lw=2,label="$\\xi=1$")
    plt.plot([-1,-1],[-1,-2],'k:' ,lw=2,label="$\\xi=0.1$")
    compute_errors("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_sthr1.000"               ,col_HIRAX,0.00,3.00,'HIRAX',plot_bphz=False,lt='-')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi10.000_HIRAX_32_6_sthr1.000",col_HIRAX,0.00,3.00,None   ,plot_bphz=False,lt='--')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi1.000_HIRAX_32_6_sthr1.000" ,col_HIRAX,0.00,3.00,None   ,plot_bphz=False,lt='-.')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi0.100_HIRAX_32_6_sthr1.000" ,col_HIRAX,0.00,3.00,None   ,plot_bphz=False,lt=':')
    compute_errors("runs/IMAP","Fisher_wA_woFG_SKA_sthr1.000"               ,col_SKA,0.00,3.00,'SKA',plot_bphz=False,lt='-')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi10.000_SKA_sthr1.000",col_SKA,0.00,3.00,None ,plot_bphz=False,lt='--')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi1.000_SKA_sthr1.000" ,col_SKA,0.00,3.00,None ,plot_bphz=False,lt='-.')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi0.100_SKA_sthr1.000" ,col_SKA,0.00,3.00,None ,plot_bphz=False,lt=':')
    compute_errors("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr1.000"               ,col_MKT,0.00,3.00,'MeerKAT',plot_bphz=False,lt='-')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi10.000_MeerKAT_sthr1.000",col_MKT,0.00,3.00,None     ,plot_bphz=False,lt='--')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi1.000_MeerKAT_sthr1.000" ,col_MKT,0.00,3.00,None     ,plot_bphz=False,lt='-.')
    compute_errors("runs/IMAP","Fisher_wA_wFG_a1.000_xi0.100_MeerKAT_sthr1.000" ,col_MKT,0.00,3.00,None     ,plot_bphz=False,lt=':')
    plt.xlim([0,2.7])
    plt.ylim([1E-4,3E-2])
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(bbox_to_anchor=[0.7,1.0],frameon=False,labelspacing=0.1,ncol=2)
    plt.savefig("bak/compare_FG.pdf",bbox_inches='tight')
    
if whichfig=='fig7' or whichfig=='all' :
    cols=['#003300','#006600','#00CC00','#33FF33','#99FF99']

    plt.figure()
    icol=0
    for texp in [1,2,3,4] :
        z,es,eb=compute_errors("runs/IMAP","Fisher_wA_woFG_gen_sT%.3lf_sthr1.000"%(float(texp)),
                               cols[icol],0.00,3.00,
                               "$10^{-%d}\\,[{\\rm mK\\,rad\\,MHz}^{1/2}]$"%texp,plot_bphz=False)
        icol+=1
    texp=6
    z,es,eb=compute_errors("runs/IMAP","Fisher_wA_woFG_gen_sT%.3lf_sthr1.000"%(float(texp)),
                           cols[icol],0.00,3.00,
                           "Noiseless case",plot_bphz=False)
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.xlim([0,2.7])
    plt.ylim([3E-5,3E-1])
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='upper right',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_noise.pdf",bbox_inches='tight');

if whichfig=='fig8' or whichfig=='all' :
    cols=['#003300','#006600','#009900','#00CC00','#00FF00','#33FF33','#66FF66','#99FF99','#CCFFCC']

    plt.figure()
    icol=0
    for dmin in [200.,100.,50.,25.,12.,6.,3.,0.] :
        if dmin==7.5 :
            num="%.1lf"%dmin
        else :
            num="%d"%(int(dmin))
        z,es,eb=compute_errors("runs/IMAP","Fisher_wA_woFG_gen_dmn%.3lf_sthr1.000"%dmin,
                               cols[icol],0.00,3.00,'$d_{\\rm min}='+num+'\\,{\\rm m}$',plot_bphz=False)
        icol+=1
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.xlim([0,2.7])
    plt.ylim([4E-5,1E-2])
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_dmin.pdf",bbox_inches='tight');

if whichfig=='fig9' or whichfig=='all' :
    cols=['#003300','#006600','#009900','#00CC00','#00FF00','#33FF33','#66FF66','#99FF99','#CCFFCC']

    plt.figure()
    icol=0
    for dmax in [7.5,15.,30.,60.,120.,240.] :
        if dmax==7.5 :
            num="%.1lf"%dmax
        else :
            num="%d"%(int(dmax))
        z,es,eb=compute_errors("runs/IMAP","Fisher_wA_woFG_gen_dmx%.3lf_sthr1.000"%dmax,
                               cols[icol],0.00,3.00,'$d_{\\rm max}='+num+'\\,{\\rm m}$',plot_bphz=False)
        icol+=1
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\sigma_z)$',fontsize=18)
    plt.xlim([0,2.7])
    plt.ylim([1E-4,1E-2])
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='upper left',frameon=False,labelspacing=0.1,ncol=2,fontsize=16)
    plt.savefig("bak/compare_dmax.pdf",bbox_inches='tight');

if whichfig=='fig10' or whichfig=='all' :
    def plot_nz(fname,area,lt,label) :
        z,nz=np.loadtxt(fname,unpack=True)
        plt.plot(z,nz*area*60**2/area,'k'+lt,lw=2,label=label)

    plt.figure()
    plot_nz("curves_spec/nz_DESI.txt"  ,14000.,'-' ,'DESI')
    plot_nz("curves_spec/nz_Euclid.txt",15000.,'--','Euclid')
    plot_nz("curves_spec/nz_WFIRST.txt", 2000.,'-.','WFIRST')
    plt.yscale('log')
    plt.xlim([0,3]); plt.ylim([1E2,3E4]);
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$dN/dz\\,d\\Omega\\,[{\\rm deg}^{-2}]$',fontsize=18)
    plt.legend(loc='upper left',frameon=False,labelspacing=0.1,fontsize=16)
    plt.savefig("bak/nz_spec.pdf",bbox_inches='tight')

if whichfig=='fig11' or whichfig=='all':
    plt.figure()
    zarr_nofg ,err_sthz_nocmb_arr ,err_bthz_nofg_arr =compute_errors("runs/CMBL/","Fisher_woA_none_sthr1.000",'#009900',0.00,3.00,'Self-calibration',fix_sphz=True)
    zarr_nofg ,err_sthz_cmbS3_arr ,err_bthz_cmbS3_arr=compute_errors("runs/CMBL/","Fisher_woA_S3_opt_sthr1.000"  ,'#6600CC',0.00,3.00,'$8.5\\,\\mu{\\rm K}\\,{\\rm arcmin}$ CMB lensing',fix_sphz=True)
    zarr_nofg ,err_sthz_cmbS4_arr ,err_bthz_cmbS4_arr=compute_errors("runs/CMBL/","Fisher_woA_S4_opt_sthr1.000"  ,'#009999',0.00,3.00,'$1\\,\\mu{\\rm K}\\,{\\rm arcmin}$ CMB lensing',fix_sphz=True)
    plt.xlim([0,2.7])
    plt.ylim([1E-3,1E-1])
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.yscale('log')
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\Delta z)$',fontsize=18)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(loc='upper right',frameon=False,labelspacing=0.1,ncol=1,fontsize=16)
    plt.savefig("bak/compare_CMBlens.pdf",bbox_inches='tight')

if whichfig=='fig12' or whichfig=='all':
    plt.figure()
    plt.plot([-1,-1],[-1,-2],'k-' ,lw=2,label="Nominal")
    plt.plot([-1,-1],[-1,-2],'k--',lw=2,label="$f_{\\rm sky}/2$")
    plt.plot([-1,-1],[-1,-2],'k-.',lw=2,label="$f_{\\rm sky}/4$")
    plt.plot([-1,-1],[-1,-2],'k:' ,lw=2,label="$f_{\\rm sky}/8$")
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_HIRAX_32_6_sthr1.000_fs0.400",col_HIRAX,0.00,3.00,"HIRAX",plot_bphz=False,lt='-')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_HIRAX_32_6_sthr1.000_fs0.200",col_HIRAX,0.00,3.00,None,plot_bphz=False,lt='--')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_HIRAX_32_6_sthr1.000_fs0.100",col_HIRAX,0.00,3.00,None,plot_bphz=False,lt='-.')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_HIRAX_32_6_sthr1.000_fs0.050",col_HIRAX,0.00,3.00,None,plot_bphz=False,lt=':')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_SKA_sthr1.000_fs0.400",col_SKA,0.00,3.00,"SKA",plot_bphz=False,lt='-')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_SKA_sthr1.000_fs0.200",col_SKA,0.00,3.00,None,plot_bphz=False,lt='--')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_SKA_sthr1.000_fs0.100",col_SKA,0.00,3.00,None,plot_bphz=False,lt='-.')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_SKA_sthr1.000_fs0.050",col_SKA,0.00,3.00,None,plot_bphz=False,lt=':')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_MeerKAT_sthr1.000_fs0.400",col_MKT,0.00,3.00,"MeerKAT",plot_bphz=False,lt='-')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_MeerKAT_sthr1.000_fs0.200",col_MKT,0.00,3.00,None,plot_bphz=False,lt='--')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_MeerKAT_sthr1.000_fs0.100",col_MKT,0.00,3.00,None,plot_bphz=False,lt='-.')
    z,es,eb=compute_errors("runs/IMAP","Fisher_woA_woFG_MeerKAT_sthr1.000_fs0.050",col_MKT,0.00,3.00,None,plot_bphz=False,lt=':')
    plt.xlim([0,2.7])
    plt.ylim([1E-4,1E-2])
    plt.yscale('log')
    plt.gca().tick_params(axis='x', labelsize=14)
    plt.gca().tick_params(axis='y', labelsize=14)
    plt.xlabel('$z$',fontsize=18)
    plt.ylabel('$\\sigma(\\Delta z)$',fontsize=18)
    plt.plot(zreq,1E-3*(1+zreq),'k-')
    plt.legend(frameon=False,labelspacing=0.1,ncol=2,bbox_to_anchor=[0.65,1.0])
    plt.savefig("bak/compare_fsky.pdf",bbox_inches='tight')

def plot_nz_err(fname,col,label,ax=None,plot_smooth=False) :
    data=np.load(fname)
    zm=data['z']
    phi=data['phi']
    ephi=data['err_phi']

    dz=zm[-1]-zm[0]
    if plot_smooth :
        zarr=zm[0]-dz/2+2*dz*np.arange(256)/255.
        phif=interp1d(zm,phi,kind='cubic',bounds_error=False,fill_value=0)
        phiarr=phif(zarr)
        ax.plot(zarr,phiarr,'k-',lw=1)
    ibad=np.where(ephi>0.1)[0]
    if ibad!=[] :
        ephi[ibad]=ephi[ibad[0]-1]*np.exp(5*(zm[ibad]-zm[ibad[0]-1]))
    ax.errorbar(zm,phi,yerr=ephi,fmt='o',ecolor=col,color=col,label=label,ms=0.1,elinewidth=2)
    return zm[0]-dz/10,zm[-1]+dz/5

if whichfig=='fig13' or whichfig=='all':
    ibin=4
    plt.figure()
    ax=plt.gca()
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_MeerKAT.npz",col_MKT,'MeerKAT',ax=ax,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_DESI.npz",'#00FF00','DESI',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,None,ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=ax)
    plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize=16)
    axins=zoomed_inset_axes(ax,2.5,loc=1)
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_MeerKAT.npz",col_MKT,'MeerKAT',ax=axins,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_DESI.npz",'#00FF00','DESI',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,None,ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=axins)
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    axins.set_xlim(0.48,0.52)
    axins.set_ylim(0.12,0.17)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    ax.set_xlim(zr)
    ax.set_ylim([-0.02,0.2])
    ax.get_yaxis().set_ticks([])
    ax.set_xlabel('$z$',fontsize=18)
    ax.set_ylabel('$\\phi(z)$',fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.savefig("bak/compare_nonpar_bin%d.pdf"%ibin,bbox_inches='tight')

    ibin=8
    plt.figure()
    ax=plt.gca()
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_MeerKAT.npz",col_MKT,'MeerKAT',ax=ax,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_DESI.npz",'#00FF00','DESI',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=ax)
    plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize=16)
    axins=zoomed_inset_axes(ax,2.5,loc=1)
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_MeerKAT.npz",col_MKT,'MeerKAT',ax=axins,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_DESI.npz",'#00FF00','DESI',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=axins)
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    axins.set_xlim(1.09,1.18)
    axins.set_ylim(0.12,0.17)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    ax.set_xlim(zr)
    ax.set_ylim([-0.02,0.2])
    ax.get_yaxis().set_ticks([])
    ax.set_xlabel('$z$',fontsize=18)
    ax.set_ylabel('$\\phi(z)$',fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.savefig("bak/compare_nonpar_bin%d.pdf"%ibin,bbox_inches='tight')

    ibin=11
    plt.figure()
    ax=plt.gca()
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=ax,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_WFIRST.npz",'#00FF00','WFIRST',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=ax)
    plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize=16)
    axins=zoomed_inset_axes(ax,2.5,loc=1)
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=axins,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_WFIRST.npz",'#00FF00','WFIRST',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=axins)
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    axins.set_xlim(1.75,1.85)
    axins.set_ylim(0.13,0.155)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    ax.set_xlim(zr)
    ax.set_ylim([-0.02,0.17])
    ax.get_yaxis().set_ticks([])
    ax.set_xlabel('$z$',fontsize=18)
    ax.set_ylabel('$\\phi(z)$',fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.savefig("bak/compare_nonpar_bin%d.pdf"%ibin,bbox_inches='tight')

    ibin=14
    plt.figure()
    ax=plt.gca()
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=ax,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_WFIRST.npz",'#00FF00','WFIRST',ax=ax)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=ax)
    plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize=16)
    axins=zoomed_inset_axes(ax,2.5,loc=1)
    zr=plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_SKA.npz",col_SKA,'SKA',ax=axins,
                   plot_smooth=True)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,'HIRAX',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_WFIRST.npz",'#00FF00','WFIRST',ax=axins)
    plot_nz_err("runs_non_parametric/result_b%d_"%ibin+"lmax2000_HIRAX_32_6.npz",col_HIRAX,None,ax=axins)
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    axins.set_xlim(2.6,2.735)
    axins.set_ylim(0.124,0.155)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    ax.set_xlim(zr)
    ax.set_ylim([-0.02,0.17])
    ax.get_yaxis().set_ticks([])
    ax.set_xlabel('$z$',fontsize=18)
    ax.set_ylabel('$\\phi(z)$',fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.savefig("bak/compare_nonpar_bin%d.pdf"%ibin,bbox_inches='tight')

'''
pars_all=np.array([{'name':'tau' ,'val': 6.0000E-02,'label':'$\\tau$'},
                   {'name':'mnu' ,'val': 6.0000E+01,'label':'$\\Sigma m_\\nu$'},
                   {'name':'och2','val': 1.1970E-01,'label':'$\\omega_c$'},
                   {'name':'hh'  ,'val': 6.9000E-01,'label':'$h$'},
                   {'name':'obh2','val': 2.2220E-02,'label':'$\\omega_b$'},
                   {'name':'ns'  ,'val': 9.6550E-01,'label':'$n_s$'},
                   {'name':'A_s' ,'val': 2.1955E+00,'label':'$A_s\\times10^9$'},
                   {'name':'wa'  ,'val': 0.0000E+00,'label':'$w_a$'},
                   {'name':'w0'  ,'val':-1.0000E+00,'label':'$w_0$'},
                   {'name':'bias_LSST_gold_node1','val':1.0000E+00,'label':'$b_1$'},
                   {'name':'bias_LSST_gold_node2','val':1.0000E+00,'label':'$b_2$'},
                   {'name':'bias_LSST_gold_node3','val':1.0000E+00,'label':'$b_3$'},
                   {'name':'bias_LSST_gold_node4','val':1.0000E+00,'label':'$b_4$'},
                   {'name':'bias_LSST_gold_node5','val':1.0000E+00,'label':'$b_5$'},
                   {'name':'bias_LSST_gold_node6','val':1.0000E+00,'label':'$b_6$'},
                   {'name':'bias_LSST_gold_node7','val':1.0000E+00,'label':'$b_7$'},
                   {'name':'bias_LSST_gold_node8','val':1.0000E+00,'label':'$b_8$'},
                   {'name':'bias_LSST_gold_node9','val':1.0000E+00,'label':'$b_9$'},
                   {'name':'sphz_LSST_gold_node1','val':1.0000E+00,'label':'$\\sigma_z^1$'},
                   {'name':'sphz_LSST_gold_node2','val':1.0000E+00,'label':'$\\sigma_z^2$'},
                   {'name':'sphz_LSST_gold_node3','val':1.0000E+00,'label':'$\\sigma_z^3$'},
                   {'name':'sphz_LSST_gold_node4','val':1.0000E+00,'label':'$\\sigma_z^4$'},
                   {'name':'sphz_LSST_gold_node5','val':1.0000E+00,'label':'$\\sigma_z^5$'},
                   {'name':'sphz_LSST_gold_node6','val':1.0000E+00,'label':'$\\sigma_z^6$'},
                   {'name':'sphz_LSST_gold_node7','val':1.0000E+00,'label':'$\\sigma_z^7$'},
                   {'name':'sphz_LSST_gold_node8','val':1.0000E+00,'label':'$\\sigma_z^8$'},
                   {'name':'sphz_LSST_gold_node9','val':1.0000E+00,'label':'$\\sigma_z^9$'},
                   {'name':'bphz_LSST_gold_node1','val':1.0000E+00,'label':'$\\Delta z^1$'},
                   {'name':'bphz_LSST_gold_node2','val':1.0000E+00,'label':'$\\Delta z^2$'},
                   {'name':'bphz_LSST_gold_node3','val':1.0000E+00,'label':'$\\Delta z^3$'},
                   {'name':'bphz_LSST_gold_node4','val':1.0000E+00,'label':'$\\Delta z^4$'},
                   {'name':'bphz_LSST_gold_node5','val':1.0000E+00,'label':'$\\Delta z^5$'},
                   {'name':'bphz_LSST_gold_node6','val':1.0000E+00,'label':'$\\Delta z^6$'},
                   {'name':'bphz_LSST_gold_node7','val':1.0000E+00,'label':'$\\Delta z^7$'},
                   {'name':'bphz_LSST_gold_node8','val':1.0000E+00,'label':'$\\Delta z^8$'},
                   {'name':'bphz_LSST_gold_node9','val':1.0000E+00,'label':'$\\Delta z^9$'}])
p_all=np.array([fshr.ParamFisher(p['val'],0,p['name'],p['label'],True,True,1E4) for p in pars_all[:]])
'''

def get_fisher(prefix,priors_sz,priors_bz) :
    data=np.load(prefix+"/fisher_raw.npz")
    p_all=np.array([fshr.ParamFisher(data['values'][i],0,data['names'][i],data['labels'][i],True,False,1E4)
                    for i in np.arange(len(data['names']))])
    fisher=data['fisher_tot']

    for i in np.arange(len(data['names'])) :
        for inode in np.arange(9) :
            if priors_sz[inode]>0 :
                if p_all[i].name=='sphz_LSST_gold_node%d'%inode :
                    fisher[i,i]+=1./priors_sz[inode]**2
            if priors_bz[inode]>0 :
                if p_all[i].name=='bphz_LSST_gold_node%d'%inode :
                    fisher[i,i]+=1./priors_bz[inode]**2

    return p_all,fisher

def plot_subset(fish_arr,pars_arr,dict_plot,dict_prior,properties,names_arr,fac_sigma=2,fs=18,nticks=4,fname='none') :
    nfish=len(fish_arr)
    pars_arr_here=copy.deepcopy(pars_arr)
    fish_arr_here=[]
    
    for ifish in np.arange(nfish) :
        fish_arr_here.append(fish_arr[ifish])
        for key in dict_plot.keys() :
            ipar=0
            for par in pars_arr_here[ifish] :
                if par.name==key :
                    pars_arr_here[ifish][ipar].do_plot=dict_plot[key]
                    break
                ipar+=1
        for key in dict_prior.keys() :
            ipar=0
            for par in pars_arr_here[ifish] :
                if par.name==key :
                    fish_arr_here[ifish][ipar,ipar]+=1./dict_prior[key]**2
                    break
                ipar+=1
                
        cov=np.linalg.inv(fish_arr_here[ifish])

        print "Case "+names_arr[ifish]+":"
        for key in dict_plot.keys() :
            ipar=0
            for par in pars_arr_here[ifish] :
                if par.name==key :
                    print key+" %lE"%(np.sqrt(cov[ipar,ipar]))
                    break
                ipar+=1
                
        npar=len(pars_arr[ifish])
        covde=np.zeros([2,2])
        ipw0=-1; ipwa=-1;
        for ipar in np.arange(npar) :
            if pars_arr_here[ifish][ipar].name=="w0" :
                ipw0=ipar
            if pars_arr_here[ifish][ipar].name=="wa" :
                ipwa=ipar
        covde[0,0]=cov[ipw0,ipw0]
        covde[0,1]=cov[ipw0,ipwa]
        covde[1,0]=cov[ipwa,ipw0]
        covde[1,1]=cov[ipwa,ipwa]
        print "FoM : %lE"%(1./np.sqrt(np.linalg.det(covde)))
        print ""

    fshr.plot_fisher_all(pars_arr_here[0],fish_arr_here,properties,
                         names_arr,fac_sigma,fname,do_1D=False,fs=fs,nticks=nticks)

if whichfig=='fig14' or whichfig=='all':
    data=np.loadtxt('curves_LSST/bins_gold_conservative_lmax2000.txt',unpack=True);
    zm=0.5*(data[0]+data[1])
    
    def get_errz(prefix,fisher) :
        z,es,eb=compute_errors(prefix,fisher,'#FFFFFF',0.00,3.00,'none',do_plot=False)
        esf=interp1d(z,es)
        ebf=interp1d(z,eb)
        return esf,ebf
    
    es_HIRAX_f,eb_HIRAX_f=get_errz("runs/IMAP","Fisher_wA_woFG_HIRAX_32_6_sthr1.000")
    es_HIRAX=es_HIRAX_f(zm)
    eb_HIRAX=eb_HIRAX_f(zm)

    es_SKA_f,eb_SKA_f=get_errz("runs/IMAP","Fisher_wA_woFG_SKA_sthr1.000")
    es_SKA=es_SKA_f(zm)
    eb_SKA=eb_SKA_f(zm)

    es_MeerKAT_f,eb_MeerKAT_f=get_errz("runs/IMAP","Fisher_wA_woFG_MeerKAT_sthr1.000")
    es_MeerKAT=es_MeerKAT_f(zm)
    eb_MeerKAT=eb_MeerKAT_f(zm)

    es_DESI_f,eb_DESI_f=get_errz("runs/DESI","Fisher_wA_sthr1.000")
    es_DESI=es_DESI_f(zm)
    eb_DESI=eb_DESI_f(zm)

    p0,f0=get_fisher("outputs_FisherPerBin/Fisher_clust_shear_cmb",
                     -1*np.ones_like(zm),-1*np.ones_like(zm))
    p1,f1=get_fisher("outputs_FisherPerBin/Fisher_clust_shear_cmb",
                     es_HIRAX,eb_HIRAX)
    p2,f2=get_fisher("outputs_FisherPerBin/Fisher_clust_shear_cmb",
                     1E-16*(1+zm),1E-16*(1+zm))
    plot_subset([f0,f1,f2],[p0,p1,p2],{'w0':True,'wa':True,'mnu':True},
                {},
                [{'ls':'solid','col':'blue','alpha':0.5,'lw':2},
                 {'ls':'solid','col':'red','alpha':1.0,'lw':2},
                 {'ls':'dashed','col':'black','alpha':-1 ,'lw':2}],
                ['Self-calibrated','Calibrated with 21cm','No photo-$z$ uncert.'],
                fname='bak/compare_constraints.pdf')

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import erf
import experiments as xpr
import emcee as mc
import corner as cn

NU_21=1420.405751786
CLIGHT=299.792458
FWHM2G=0.42466090014
prefix_ell="lmax2000"
exp_arr=[xpr.im_HIRAX_32_6,xpr.im_SKA,xpr.im_MeerKAT,xpr.spec_DESI,xpr.spec_Euclid,xpr.spec_WFIRST]
n_bins_photo=15

def run_lj(prefix):
    stout ="#Params for LimberJack"
    stout+="#Cosmological parameters\n"
    stout+="omega_m= 0.3\n"
    stout+="omega_l= 0.7\n"
    stout+="omega_b= 0.05\n"
    stout+="w0= -1.\n"
    stout+="wa= 0.\n"
    stout+="h= 0.7\n"
    stout+="ns= 0.96\n"
    stout+="s8= 0.8\n"
    stout+="\n"
    stout+="#Radial resolution\n"
    stout+="d_chi= 2.\n"
    stout+="r_smooth= 0.1\n"
    stout+="z_kappa= 0.1\n"
    stout+="z_isw= 0.1\n"
    stout+="\n"
    stout+="#Maximum multipole\n"
    stout+="l_max= 2000\n"
    stout+="dl= 1\n"
    stout+="l_limber_min= 0\n"
    stout+="\n"
    stout+="#Behavioural flags (include number counts? lensing shear? CMB lensing? ISW?)\n"
    stout+="do_nc= 1\n"
    stout+="has_nc_dens= 1\n"
    stout+="has_nc_rsd= 0\n"
    stout+="has_nc_lensing= 0\n"
    stout+="has_nc_lognorm= 0\n"
    stout+="do_shear= 0\n"
    stout+="has_sh_intrinsic= 0\n"
    stout+="do_cmblens= 0\n"
    stout+="do_isw= 0\n"
    stout+="\n"
    stout+="#Angular correlation function\n"
    stout+="do_w_theta= 0\n"
    stout+="use_logbin= 0\n"
    stout+="theta_min= 0\n"
    stout+="theta_max= 10.\n"
    stout+="n_bins_theta= 15\n"
    stout+="n_bins_decade= 5\n"
    stout+="\n"
    stout+="#File names (window function, bias, magnification bias, power spectrum)\n"
    stout+="window_1_fname= runs_non_parametric/nz1.txt\n"
    stout+="window_2_fname= runs_non_parametric/nz2.txt\n"
    stout+="bias_fname= runs_non_parametric/bias.txt\n"
    stout+="sbias_fname= bullshit\n"
    stout+="abias_fname= bullshit\n"
    stout+="pk_fname= EH\n"
    stout+="\n"
    stout+="#Output prefix\n"
    stout+="prefix_out= "+prefix+"\n"

    param_fname=prefix+"_par.ini"
    f=open(param_fname,"w")
    f.write(stout)
    f.close()

    command="./LimberJack/LimberJack "+param_fname+" > "+prefix+"_log"
    os.system(command)
    os.system("rm "+param_fname+" "+prefix+"_log")

def compute_cij_single(zedg1,zedg2,prefix) :
    nz=32
    def get_nzarr(z0,zf) :
        dz=zf-z0
        zarr=z0-2*dz+5*dz*np.arange(5*nz)/(5*nz-1)
        narr=np.zeros_like(zarr);
        narr[np.where((zarr<=zf) & (zarr>z0))]=1.
        return zarr,narr
    z1,n1=get_nzarr(zedg1[0],zedg1[1])
    z2,n2=get_nzarr(zedg2[0],zedg2[1])
    zbarr=3.*np.arange(1024)/1023.
    bbarr=np.ones_like(zbarr)

    np.savetxt("runs_non_parametric/nz1.txt",np.transpose([z1,n1]))
    np.savetxt("runs_non_parametric/nz2.txt",np.transpose([z2,n2]))
    np.savetxt("runs_non_parametric/bias.txt",np.transpose([zbarr,bbarr]))
    
    run_lj(prefix);

def compute_cij_full(nbins) :
    for ib in np.arange(nbins) :
        fname="runs/IMAP/bins_im_b%d_lmax2000.txt"%ib
        data=np.loadtxt(fname,unpack=True)
        z0_arr=data[0];
        zf_arr=data[1];
        nb=len(z0_arr)
        for i in np.arange(nb) :
            print ib, i
            compute_cij_single([z0_arr[i],zf_arr[i]],[z0_arr[i],zf_arr[i]],'runs_non_parametric/thin_b%d_'%ib+'sb%d'%i)

def pdf_photo(z,z0,zf,sz) :
    denom=1./np.sqrt(2*sz*sz)
    return 0.5*(erf((zf-z)*denom)-erf((z0-z)*denom))

def compute_noise_cl(exper,z0ar,zfar,ls) :
    nbin=len(z0ar)
    nij=np.zeros([nbin,len(ls)])
    if exper['type']=='gal_clustering' :
        znz,nnz=np.loadtxt(exper['nzfi'],unpack=True)
        nnz*=(180.*60/np.pi)**2
        nzf=interp1d(znz,nnz,bounds_error=False,fill_value=0)
        for i in np.arange(nbin) :
            ndens=max(1E-16,quad(nzf,z0ar[i],zfar[i])[0])
            nij[i,:]=1./ndens
    elif exper['type']=='intensity_mapping' :
        nu0_arr=NU_21/(1+zfar)
        nuf_arr=NU_21/(1+z0ar)
        nu_arr=0.5*(nu0_arr+nuf_arr)
        dnu_arr=nuf_arr-nu0_arr

        z,tz=np.loadtxt(exper['tzfi'],unpack=True)
        tofz=interp1d(z,tz)
        tbg_arr=np.array([tofz(z) for z in 0.5*(z0ar+zfar)])

        tsys_arr=(exper['t_inst']+60.*(nu_arr/300.)**(-2.5))*1000
        sigma2_noise=(tsys_arr/tbg_arr/exper['area_eff'])**2
        sigma2_noise*=4*np.pi*exper['fsky']/(3.6E9*exper['t_total']*dnu_arr)
        beam_fwhm=CLIGHT/(exper['dish_size']*nu_arr)
        
        if ((exper['im_type']=="single_dish") or (exper['im_type']=="hybrid")) :
            beam_rad=beam_fwhm*FWHM2G
            factor_beam_sd=exper['n_dish']*np.exp(-(ls*(ls+1))[None,:]*(beam_rad**2)[:,None])
        else :
            factor_beam_sd=np.zeros([len(nu_arr),len(ls)])
        if ((exper['im_type']=="interferometer") or (exper['im_type']=="hybrid")) :
            lambda_arr=CLIGHT/nu_arr
            dist,nbase=np.loadtxt(exper['base_file'],unpack=True)
            ndistint=interp1d(dist,nbase*dist*2*np.pi,bounds_error=False,fill_value=0.)
            norm=0.5*exper['n_dish']*(exper['n_dish']-1.)/quad(ndistint,dist[0],dist[-1])[0]
            nbase*=norm; ndist=interp1d(dist,nbase,bounds_error=False,fill_value=0.)
            n_baselines=ndist(ls[None,:]*lambda_arr[:,None]/(2*np.pi))
            factor_beam_if=n_baselines[:,:]*((lambda_arr/beam_fwhm)**2)[:,None]
        elif exper['im_type']=="generic" :
            lambda_arr=CLIGHT/nu_arr
            dist_arr=(ls[None,:]*lambda_arr[:,None]).flatten()
            f=np.exp(-(dist_arr*FWHM2G/np.fmax(exper['base_max'],1E-1))**2)
            f[np.where(dist_arr<2*np.pi*exper['base_min'])]=0.
            factor_beam_if=np.reshape(f,[len(lambda_arr),len(ls)])
            sigma2_noise=(exper['t_inst']/tbg_arr)**2/dnu_arr
        else :
            factor_beam_if=np.zeros([len(nu_arr),len(ls)])
        factor_beam=np.fmax(factor_beam_sd,factor_beam_if)
        for i in np.arange(nbin) :
            nij[i,:]=sigma2_noise[i]/np.fmax(factor_beam[i,:],1E-16)

    return nij
            
def compute_errors_single(ibin,nzf,bzphf,exper) :
    fsky=exper['fsky']

    data0=np.loadtxt(exper['bzfi'],unpack=True)
    bzspf=interp1d(data0[0],data0[1],bounds_error=False,fill_value=1.)

    data1=np.loadtxt("runs/IMAP/bins_photoz_b%d_"%ibin+prefix_ell+".txt",unpack=True)
    z0_ph=data1[0]; zf_ph=data1[1]; sz_ph=data1[2]; zm_ph=0.5*(z0_ph+zf_ph)

    data2=np.loadtxt("runs/IMAP/bins_im_b%d_"%ibin+prefix_ell+".txt",unpack=True)
    z0_sp=data2[0]; zf_sp=data2[1]; sz_sp=data2[2]; zm_sp=0.5*(z0_sp+zf_sp); lm_sp=data2[5]
    nb=len(z0_sp)

    def nzph(z) :
        return nzf(z)*pdf_photo(z,z0_ph,zf_ph,sz_ph)
    ndens=quad(nzph,z0_ph-5*sz_ph,zf_ph+5*sz_ph)[0]
    phiph=nzph(zm_sp); phiph/=np.sum(phiph)
    bph=bzphf(zm_sp)
    bsp=bzspf(zm_sp)

    cij=[]
    for i in np.arange(nb) :
        l,cl=np.loadtxt('runs_non_parametric/thin_b%d_'%ibin+'sb%d_cl_dd.txt'%i,unpack=True)
        cij.append(cl)
    cij=np.array(cij)

    noi=compute_noise_cl(exper,z0_sp,zf_sp,l)
    for i in np.arange(nb) :
        noi[i,np.where(l>lm_sp[i])[0]]=1E16

    a00=np.sum((phiph*bph)[:,None]*cij[:,:],axis=0)+1./ndens
    a0i=(phiph*bph*bsp)[:,None]*cij[:,:]
    aii=(bsp**2)[:,None]*cij[:,:]+noi[:,:]
    ader=(bsp*bph)[:,None]*cij[:,:]

    rcorr=a0i/np.sqrt(a00[None,:]*aii)
    s=1./(1-np.sum(rcorr**2,axis=0))

    fish_l=np.zeros([nb,nb,len(l)])
    prefac=fsky*(l+0.5)*s/a00
    for i1 in np.arange(nb) :
        cont_diag=1./aii[i1,:]
        for i2 in np.arange(nb) :
            v=2*s*np.sqrt(rcorr[i1,:]*rcorr[i2,:]/(aii[i1,:]*aii[i2,:]))
            if i1==i2 :
                v+=1./aii[i1,:]
            fish_l[i1,i2,:]=prefac*v*ader[i1,:]*ader[i2,:]

    fish=np.sum(fish_l,axis=2)
    covar=np.linalg.inv(fish)

#    plt.errorbar(zm_sp,phiph,yerr=np.sqrt(np.diag(covar)),label=exper['name'])
    np.savez("runs_non_parametric/result_b%d_"%ibin+prefix_ell+"_"+exper['name'],
             z=zm_sp,phi=phiph,err_phi=np.sqrt(np.diag(covar)),covar=covar)

def get_sigma_errors(ibin,nzf,exper) :
    data1=np.loadtxt("runs/IMAP/bins_photoz_b%d_"%ibin+prefix_ell+".txt",unpack=True)
    z0_ph=data1[0]; zf_ph=data1[1]; sz_ph=data1[2]; zm_ph=0.5*(z0_ph+zf_ph)

    data=np.load("runs_non_parametric/result_b%d_"%ibin+prefix_ell+"_"+exper['name']+".npz")
    zarr=data['z']
    phiarr=data['phi']
    errarr=data['err_phi']

    def theo(dz,sz) :
        phi=nzf(zarr)*pdf_photo(zarr,z0_ph+dz,zf_ph+dz,sz)
        norm=np.sum(phi)
        if norm<=0 :
            return 1E16
        else :
            return phi/norm

    def logprob(pars,dum) :
        phith=theo(pars[0],pars[1])
        dphi=phith-phiarr
        return -0.5*np.sum((dphi/errarr)**2)

    p0=np.array([0.0,sz_ph])
    ng=512
    res=0.01
    p0arr=p0[0]-res+2*res*(np.arange(ng)+0.5)/ng
    p1arr=p0[1]-res+2*res*(np.arange(ng)+0.5)/ng
    gr=np.array([[logprob([pp0,pp1],0) for pp0 in p0arr] for pp1 in p1arr])
    prob=np.exp(gr); prob/=np.sum(prob)
    sig_dz=np.sqrt(np.sum(prob[:,:]*((p0arr-p0[0])**2)[None,:]))
    sig_sz=np.sqrt(np.sum(prob[:,:]*((p1arr-p0[1])**2)[:,None]))
    print sig_dz,sig_sz

#    plt.imshow(prob,interpolation='nearest',origin='lower',extent=[-res,res,-res,res])
#    plt.colorbar()
#    plt.show()
#    nwalkers=100
#    nsteps_burn=100
#    nsteps_per_chain=1000
#    npar=2
#    par0=(p0)[None,:]+(0.001*np.random.randn(nwalkers))[:,None]
#    sampler=mc.EnsembleSampler(nwalkers,npar,logprob,args=[0.0])
#    print "burn"
#    pos,prob,stat=sampler.run_mcmc(par0,nsteps_burn)
#    sampler.reset()
#    print "run"
#    sampler.run_mcmc(pos,nsteps_per_chain)
#    print np.shape(sampler.chain)
#    samples=(sampler.chain)[:,nsteps_per_chain/2:,:].reshape((-1,npar))
#    print np.mean(samples,axis=0)
#    print np.std(samples,axis=0)
#    print("Mean acceptance fraction: {0:.3f}"
#          .format(np.mean(sampler.acceptance_fraction)))
#    
#    fig=cn.corner(samples,labels=['$\\Delta z$','$\\sigma_z$'],truths=[0.0,sz_ph])
#    plt.show()
    
    np.savetxt("runs_non_parametric/sigmas_b%d_"%ibin+prefix_ell+"_"+exper['name']+".txt",
               [sig_dz,sig_sz])

#compute_cij_full(n_bins_photo)
#zph,nzph=np.loadtxt("curves_LSST/nz_shear_fiducial.txt",unpack=True)
#zphb,bphb,dum=np.loadtxt("curves_LSST/bz_gold.txt",unpack=True)
#nzphf=interp1d(zph,nzph*(180*60/np.pi)**2,bounds_error=False,fill_value=0)
#bzphf=interp1d(zphb,bphb,bounds_error=False,fill_value=1.)

#for ib in np.arange(n_bins_photo) :
#    plt.figure()
#    plt.title("%d"%ib)
#    for xper in exp_arr :
#        compute_errors_single(ib,nzphf,bzphf,xper)
#        get_sigma_errors(ib,nzphf,xper)
#    plt.legend(loc='upper left')
#    plt.show()

cols=['k','r','b','g','y','c','m']
icol=0
for xper in exp_arr :
    zarr=np.zeros(n_bins_photo); szarr=np.zeros(n_bins_photo); eszarr=np.zeros(n_bins_photo); edzarr=np.zeros(n_bins_photo)
    for ib in np.arange(n_bins_photo) :
        data=np.loadtxt("runs/IMAP/bins_photoz_b%d_"%ib+prefix_ell+".txt",unpack=True)
        zarr[ib]=(data[0]+data[1])/2
        szarr[ib]=data[2]
        data=np.loadtxt("runs_non_parametric/sigmas_b%d_"%ib+prefix_ell+"_"+xper['name']+".txt",unpack=True)
        edzarr[ib]=data[0]
        eszarr[ib]=data[1]
    np.savetxt("runs_non_parametric/sigmas_all_"+prefix_ell+"_"+xper['name']+".txt",
               np.transpose([zarr,szarr,edzarr,eszarr]),header='[0]-z [1]-sz [2]-e_dz [3]-e_sz')
    plt.plot(zarr,eszarr,cols[icol]+'-',label=xper['name'])
    plt.plot(zarr,edzarr,cols[icol]+'--')
    icol+=1
plt.legend(loc='upper right')
plt.show()

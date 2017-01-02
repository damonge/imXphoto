import numpy as np
import matplotlib.pyplot as plt
import py_cosmo_mad as csm
from scipy.interpolate import interp1d

plot_stuff=True
nz=256
zmax=3.5
bDESI=0.84
bEuclid=0.76
bWFIRST=0.76

pcs=csm.PcsPar()
pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.7255)

zDESI_in,nzDESI_in=np.loadtxt("nz_DESI_deg.txt",unpack=True); 
nzf=interp1d(zDESI_in,nzDESI_in,bounds_error=False,fill_value=0)
zDESI=zmax*np.arange(nz)/(nz-1.)
nzDESI=nzf(zDESI)/60.**2
bzDESI=bDESI*pcs.growth_factor(1.)/np.array([pcs.growth_factor(1./(1+z)) for z in zDESI])
np.savetxt("nz_DESI.txt",np.transpose([zDESI,nzDESI]))
np.savetxt("bz_DESI.txt",np.transpose([zDESI,bzDESI]))
np.savetxt("sz_DESI.txt",np.transpose([zDESI,0.4*np.ones(nz)]))
np.savetxt("ez_DESI.txt",np.transpose([zDESI,np.zeros(nz)]))

zEuclid_in,nzEuclid_in=np.loadtxt("nz_Euclid_deg.txt",unpack=True); 
nzf=interp1d(zEuclid_in,nzEuclid_in,bounds_error=False,fill_value=0)
zEuclid=zmax*np.arange(nz)/(nz-1.)
nzEuclid=nzf(zEuclid)/60.**2
bzEuclid=bEuclid*pcs.growth_factor(1.)/np.array([pcs.growth_factor(1./(1+z)) for z in zEuclid])
np.savetxt("nz_Euclid.txt",np.transpose([zEuclid,nzEuclid]))
np.savetxt("bz_Euclid.txt",np.transpose([zEuclid,bzEuclid]))
np.savetxt("sz_Euclid.txt",np.transpose([zEuclid,0.4*np.ones(nz)]))
np.savetxt("ez_Euclid.txt",np.transpose([zEuclid,np.zeros(nz)]))

zWFIRST_in,nzWFIRST_in=np.loadtxt("nz_WFIRST_deg.txt",unpack=True); 
nzf=interp1d(zWFIRST_in,nzWFIRST_in,bounds_error=False,fill_value=0)
zWFIRST=zmax*np.arange(nz)/(nz-1.)
nzWFIRST=nzf(zWFIRST)/60.**2
bzWFIRST=bWFIRST*pcs.growth_factor(1.)/np.array([pcs.growth_factor(1./(1+z)) for z in zWFIRST])
np.savetxt("nz_WFIRST.txt",np.transpose([zWFIRST,nzWFIRST]))
np.savetxt("bz_WFIRST.txt",np.transpose([zWFIRST,bzWFIRST]))
np.savetxt("sz_WFIRST.txt",np.transpose([zWFIRST,0.4*np.ones(nz)]))
np.savetxt("ez_WFIRST.txt",np.transpose([zWFIRST,np.zeros(nz)]))


if plot_stuff :
    plt.figure()
    plt.plot(zDESI  ,nzDESI  ,'r-')
    plt.plot(zEuclid,nzEuclid,'g-')
    plt.plot(zWFIRST,nzWFIRST,'b-')
    
    plt.figure()
    plt.plot(zDESI  ,bzDESI  ,'r-')
    plt.plot(zEuclid,bzEuclid,'g-')
    plt.plot(zWFIRST,bzWFIRST,'b-')
    
    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def get_grid_square(ns,dl) :
    pos_1=(np.arange(ns)+0.5)*dl
    pos=np.array([(pos_1[:,None]*np.ones([ns,ns])).flatten(),
                  (pos_1[None,:]*np.ones([ns,ns])).flatten()])

    return pos

def get_basedist(pos,nbins,dmin,dmax,fname=None) :
    dpos=pos[:,:,None]-pos[:,None,:];
    dpos_b=np.array([dpos[0].flatten(),dpos[1].flatten()])
    dposr=np.sqrt(dpos_b[0]**2+dpos_b[1]**2); ind0=np.where(dposr!=0)[0]; dposr=dposr[ind0]
    counts,bins=np.histogram(dposr,bins=nbins,range=[dmin,dmax])#2*ns,range=[0,2*lside])
    rm=(bins[1:]+bins[:-1])*0.5
    vol=np.pi*(bins[1:]**2-bins[:-1]**2)
    n_d=counts/vol/2.;
    print rm

    if fname!=None :
        np.savetxt(fname,np.transpose([rm,n_d]))
    return rm,n_d

def get_basedist_grid(nside,d_min,fac=1,fname=None) :
    pos=get_grid_square(nside,d_min)
    return get_basedist(pos,int(fac*nside*0.75),0,1.5*nside*d_min,fname=fname)

def plotarr(d,n) :
    nf=interp1d(d,n,bounds_error=False,fill_value=0,kind='cubic')
    darr=d[0]+(d[-1]-d[0])*(np.arange(256)+0.5)/256
    narr=nf(darr)
    plt.plot(darr,narr)

darr1,nbarr1=get_basedist_grid(32,7.,fname="baseline_file_HIRAX_7m.txt")
darr2,nbarr2=get_basedist_grid(32,6.,fname="baseline_file_HIRAX_6m.txt")
darr3,nbarr3=get_basedist_grid(64,7.,fname="baseline_file_HIRAX64_7m.txt")
darr4,nbarr4=get_basedist_grid(64,6.,fname="baseline_file_HIRAX64_6m.txt")
fac_l=2*np.pi*400./300.
plotarr(fac_l*darr1,nbarr1)
plotarr(fac_l*darr2,nbarr2)
plotarr(fac_l*darr3,nbarr3)
plotarr(fac_l*darr4,nbarr4)
plt.plot(fac_l*darr1,nbarr1)
plt.plot(fac_l*darr2,nbarr2)
plt.plot(fac_l*darr3,nbarr3)
plt.plot(fac_l*darr4,nbarr4)
plt.show()

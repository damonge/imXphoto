import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plot_stuff=True

def get_grid_square(ns,dl) :
    pos_1=(np.arange(ns)+0.5)*dl
    pos=np.array([(pos_1[:,None]*np.ones([ns,ns])).flatten(),
                  (pos_1[None,:]*np.ones([ns,ns])).flatten()])

    return pos

def get_basedist(pos,nbins,dmin=None,dmax=None,fname=None,dec=90.,weigh_dist=False) :
    fac_dec=np.sin(np.pi*dec/180)
    dpos=pos[:,:,None]-pos[:,None,:];
    dpos_b=np.array([dpos[0].flatten(),dpos[1].flatten()])
    dposr=np.sqrt(dpos_b[0]**2+dpos_b[1]**2)*fac_dec;
    ind0=np.where(dposr!=0)[0]; dposr=dposr[ind0]
    if dmin==None :
        dmin=np.amin(dposr)
    if dmax==None :
        dmax=np.amax(dposr)
    counts,bins=np.histogram(dposr,bins=nbins,range=[dmin,dmax])
    rm=(bins[1:]+bins[:-1])*0.5
    if weigh_dist :
        hr,bins=np.histogram(dposr,bins=nbins,range=[dmin,dmax],weights=dposr)
        rm[np.where(counts>0)]=(hr/counts)[np.where(counts>0)]
    vol=np.pi*(bins[1:]**2-bins[:-1]**2)
    n_d=counts/vol/2.;
    print np.amin(rm), np.amax(rm), np.amin(dposr), np.amax(dposr), np.sum(n_d*vol*2)-len(pos[0])*(len(pos[0])-1.)

    if fname!=None :
        np.savetxt(fname,np.transpose([rm,n_d]))
    return rm,n_d,counts

def get_basedist_grid(nside,d_min,fac=1,fname=None,dec=90.,weigh_dist=False) :
    pos=get_grid_square(nside,d_min)
    return get_basedist(pos,int(fac*nside*0.75),dmin=0,dmax=1.5*nside*d_min,fname=fname,dec=dec,weigh_dist=weigh_dist)

def get_basedist_file(fname_in,dmin,nbins,fname=None,dec=90.,dmax=None,weigh_dist=False) :
    x,y=np.loadtxt(fname_in,unpack=True)
    pos=np.array([x,y])
    return get_basedist(pos,nbins,dmin=dmin,dmax=dmax,fname=fname,dec=dec,weigh_dist=weigh_dist)

def plotarr(d,n) :
    nf=interp1d(d,n,bounds_error=False,fill_value=0,kind='cubic')
    darr=d[0]+(d[-1]-d[0])*(np.arange(256)+0.5)/256
    narr=nf(darr)
    plt.plot(darr,narr)

darr1,nbarr1,c=get_basedist_grid(32,7.,fname="baseline_file_HIRAX_7m.txt")
darr2,nbarr2,c=get_basedist_grid(32,6.,fname="baseline_file_HIRAX_6m.txt")
darr3,nbarr3,c=get_basedist_grid(64,7.,fname="baseline_file_HIRAX64_7m.txt")
darr4,nbarr4,c=get_basedist_grid(64,6.,fname="baseline_file_HIRAX64_6m.txt")

dmax_ska=158000.; dmax_mkt=7700.; db=25.; nd_ska=int(dmax_ska/db)+1; nd_mkt=int(dmax_mkt/db)+1;
decmax=80.; decstep=0.5; ndec=int(decmax/decstep)+1
nbarr5=np.zeros(nd_ska); darr5=np.zeros(nd_ska); dcumul5=np.zeros(nd_ska); cumul5=np.zeros(nd_ska);
nbarr6=np.zeros(nd_mkt); darr6=np.zeros(nd_mkt); dcumul6=np.zeros(nd_mkt); cumul6=np.zeros(nd_mkt);
decarr=90.-decmax*np.arange(ndec)/(ndec-1.)
for dec in decarr :
    d,n,c=get_basedist_file("positions_ska.txt",0.,nd_ska,dmax=dmax_ska,dec=dec,weigh_dist=True);
    nbarr5+=n; darr5+=d; cumul5+=c; dcumul5+=d*c; 
    d,n,c=get_basedist_file("positions_mkt.txt",0.,nd_mkt,dmax=dmax_mkt,dec=dec,weigh_dist=True);
    nbarr6+=n; darr6+=d; cumul6+=c; dcumul6+=d*c;
darr5/=ndec; darr5[np.where(cumul5>0)]=(dcumul5/cumul5)[np.where(cumul5>0)]; nbarr5/=ndec;
darr6/=ndec; darr6[np.where(cumul6>0)]=(dcumul6/cumul6)[np.where(cumul6>0)]; nbarr6/=ndec;

np.savetxt("baseline_file_SKA.txt"    ,np.transpose([darr5,nbarr5]))
np.savetxt("baseline_file_MeerKAT.txt",np.transpose([darr6,nbarr6]))

fac_l=2*np.pi*400./300.
if plot_stuff :
    plt.plot(fac_l*darr1,nbarr1)
    plt.plot(fac_l*darr2,nbarr2)
    plt.plot(fac_l*darr3,nbarr3)
    plt.plot(fac_l*darr4,nbarr4)
    plt.plot(fac_l*darr5,nbarr5)
    plt.plot(fac_l*darr6,nbarr6)
    plt.show()

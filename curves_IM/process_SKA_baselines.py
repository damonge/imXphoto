import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

radius_karoo=6371831.
DTOR=np.pi/180

def get_physical_coords(lat,lon,rad) :
    #Get spherical coordinates
    cth_arr=np.cos(DTOR*(90-lat))
    sth_arr=np.sin(DTOR*(90-lat))
    phi_arr=DTOR*lon

    #Compute antena positions and subtract center
    vec_arr=np.array([sth_arr*np.cos(phi_arr),sth_arr*np.sin(phi_arr),cth_arr])
    vec_mean=np.mean(vec_arr,axis=1)
    vec_arr-=vec_mean[:,None]

    #Compute normal vectors in the theta and phi directions
    the_mean,phi_mean=hp.vec2ang(vec_mean)
    n_the=np.array([np.cos(the_mean)*np.cos(phi_mean),np.cos(the_mean)*np.sin(phi_mean),-np.sin(the_mean)])
    n_phi=np.array([-np.sin(phi_mean)                ,np.cos(phi_mean)                 ,0                ])
    
    #Project on theta and phi directions
    x_arr=rad*np.dot(np.transpose(n_phi),vec_arr)
    y_arr=rad*np.dot(np.transpose(n_the),vec_arr)[0]

    return x_arr,y_arr

coords_ang=np.genfromtxt("positions_SKA.txt",dtype=[('name','S6'),('lon','<f8'),('lat','<f8')])

x_ska,y_ska=get_physical_coords(coords_ang['lat'],coords_ang['lon'],radius_karoo)
x_sko=x_ska[64:]; y_sko=y_ska[64:]
x_mkt=x_ska[:64]; y_mkt=y_ska[:64]

np.savetxt("positions_ska.txt",np.transpose([x_ska,y_ska]))
np.savetxt("positions_mkt.txt",np.transpose([x_mkt,y_mkt]))

plt.figure()
plt.scatter(x_mkt,-y_mkt,c='r')
plt.scatter(x_sko,-y_sko,c='b')
plt.xlim([-5700,5300])
plt.ylim([-4200,4400])
plt.show()

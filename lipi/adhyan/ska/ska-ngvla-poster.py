import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pyproj
import scipy.spatial.transform as te
from astropy.io import fits

def convert_ecef2latlong(x,y,z):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
    return lon, lat, alt

def convert_ecef2enu(x,y,z,lat0, lon0, alt0):
    transformer = pyproj.Transformer.from_crs({"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
                                {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},)
    x_org, y_org, z_org = transformer.transform(lon0,lat0,alt0,radians=False)
    vec=np.array([[ x-x_org, y-y_org, z-z_org]]).T
    rot1 =  te.Rotation.from_euler('x', -(90-lat0), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  te.Rotation.from_euler('z', -(90+lon0), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rotMatrix = rot1.dot(rot3)
    enu = rotMatrix.dot(vec).T.ravel()
    return enu.T.reshape(x.shape[0],3)

ska_mid=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/ska1mid.tm/layout.txt')
ska_low=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/ska1low.tm/layout.txt')
ngvla=np.loadtxt('/home/rohit/simulations/ska-ngvla/ngvla-revD.main_ants.cfg')
ngvla=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/ngvla-revC.tm/layout.txt')
mwa=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/mwa.phase1.tm/layout.txt')
lofar=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/lofar.tm/layout.txt')
vlac=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/vla.c.tm/layout.txt')
meerkat=np.loadtxt('/home/rohit/karabo/karabo-pipeline/karabo/data/meerkat.tm/layout.txt')

tel_param={'ngvla':[34,-107.6,1000]}
x, y, z=ngvla[:,0],ngvla[:,1],ngvla[:,2]
lon, lat, alt = convert_ecef2latlong(x,y,z)
lat0,lon0,alt0=tel_param['ngvla']
enu_coord=convert_ecef2enu(x,y,z,lat0, lon0, alt0)
enu_x,enu_y,enu_z=enu_coord[:,0],enu_coord[:,1],enu_coord[:,2]
np.savetxt(np.array([enu_x,enu_y]).T,'/home/rohit/karabo/karabo-pipeline/karabo/data/ngvla-revD.tm/layout.txt')

#---------------------------------------

f,ax=plt.subplots(1,1)
ax.scatter(ska_mid[:,0]/1.e3,ska_mid[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()

f,ax=plt.subplots(1,1)
ax.scatter(ska_low[:,0]/1.e3,ska_low[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()

f,ax=plt.subplots(1,1)
ax.scatter(ngvla[:,0]/1.e3+1600,ngvla[:,1]/1.e3+5000,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()

f,ax=plt.subplots(1,1)
ax.scatter(lofar[:,0]/1.e3,lofar[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
inset_ax = inset_axes(ax, width="40%", # width = 30% of parent_bbox
                    height=1., # height : 1 inch
                    loc=2)
inset_ax.scatter(lofar[:,0]/1.e3,lofar[:,1]/1.e3,color='k')
inset_ax.set_ylim(-3,3);inset_ax.set_xlim(-3,3)
plt.show()

f,ax=plt.subplots(1,1)
ax.scatter(vlac[:,0]/1.e3,vlac[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()

f,ax=plt.subplots(1,1)
ax.scatter(mwa[:,0]/1.e3,mwa[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()


f,ax=plt.subplots(1,1)
ax.scatter(meerkat[:,0]/1.e3,meerkat[:,1]/1.e3,color='k')
ax.set_ylabel('Coordinate-Y (km)')
ax.set_xlabel('Coordinate-X (km)')
plt.show()


aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_mwa0.fits')
mwa_snap=aa[0].data[0][0]
print('MWA snap RMS:',np.nanstd(mwa_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(mwa_snap)/1.e4, '/ DR:',np.nanmax(mwa_snap)/np.nanstd(mwa_snap[0:500,0:500]))

aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_lofar0.fits')
lofar_snap=aa[0].data[0][0]
print('LOFAR snap RMS:',np.nanstd(lofar_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(lofar_snap)/1.e4, '/ DR:',np.nanmax(lofar_snap)/np.nanstd(lofar_snap[0:500,0:500]))

aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_ska-low0.fits')
skalow_snap=aa[0].data[0][0]
print('SKA low snap RMS:',np.nanstd(skalow_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(skalow_snap)/1.e4, '/ DR:',np.nanmax(skalow_snap)/np.nanstd(skalow_snap[0:500,0:500]))

print('-----------------------')

aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_meerkat0.fits')
meerkat_snap=aa[0].data[0][0]
print('Meerkat snap RMS:',np.nanstd(meerkat_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(meerkat_snap)/1.e4, '/ DR:',np.nanmax(meerkat_snap)/np.nanstd(meerkat_snap[0:500,0:500]))


aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_ska_mid0.fits')
skamid_snap=aa[0].data[0][0]
print('SKA mid snap RMS:',np.nanstd(skalow_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(skamid_snap)/1.e4, '/ DR:',np.nanmax(skamid_snap)/np.nanstd(skamid_snap[0:500,0:500]))

aa=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/data/solar/final_maps/solar_ngvlad0.fits')
ngvlad_snap=aa[0].data[0][0]
print('ngVLAd snap RMS:',np.nanstd(ngvlad_snap[0:500,0:500])/1.e4,'/ Peak:',np.nanmax(ngvlad_snap)/1.e4, '/ DR:',np.nanmax(ngvlad_snap)/np.nanstd(ngvlad_snap[0:500,0:500]))


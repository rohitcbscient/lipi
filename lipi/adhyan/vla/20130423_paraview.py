import numpy as np
from pyevtk.hl import imageToVTK
import sys
#nx, ny, nz = 6, 6, 2
#ncells = nx * ny * nz
#npoints = (nx + 1) * (ny + 1) * (nz + 1)
#pressure = np.random.rand(ncells).reshape( (nx, ny, nz), order = 'C')
#temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1))
#imageToVTK("./image", cellData = {"pressure" : pressure}, pointData = {"temp" : temp} )


from tvtk.api import tvtk, write_data
from scipy.io import readsav
from astropy.io import fits
import matplotlib.pyplot as plt
from sunpy.map import Map
import astropy.units as u
# Unpack velocity information

#aa=readsav('/media/rohit/VLA/20130423/hmi.M_720s.20130423_202220.E34N12CR.CEA.NAS.sav')
aa=readsav('/media/rohit/VLA/20130423/hmi.M_720s.20130423_202220.E31N13CR.CEA.NAS.sav')

vx=aa['box']['bx'][0]#[0:50,120:-120,120:-120]#flow['vx']
vy=aa['box']['by'][0]#[0:50,120:-120,120:-120]#flow['vy']
vz=aa['box']['bz'][0]#[0:50,120:-120,120:-120]#flow['vz']
babs=np.sqrt(vx*vx+vy*vy+vz*vz)
bspk=babs[:,172,77];r=np.arange(400)*1.4
plt.plot(r,bspk,'o-',label='Magnetic Field')
plt.show()

h=np.arange(400)*1400/1.e3
omega_spk=bspk/2.8
ne=1.e9+1.16e17*np.exp(-1*(h*1000/500.))
fp=9000*np.sqrt(ne)/1.e6
R=np.linspace(2000,280.e3,1000)
nkne=4.2e4 *10**(4.32/((695700.+h*1.e6)/695700.))
th=10*np.pi/180.
saitone=(3.09e8*(1-0.5*np.sin(th))/((695700.+h*1.e6)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+h*1.e6)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+h*1.e6)/695700.)**2.5)
nkfp=9000*np.sqrt(nkne)/1.e6
saitofp=9000*np.sqrt(saitone)/1.e6


f,ax=plt.subplots(2,2)
ax[0,0].plot(h,bspk,'o-',label='Magnetic Field')
ax[0,0].set_ylabel('|B| (Gauss)')
ax[0,1].plot(h,ne,'o-',label='Electron density')
ax[0,1].set_ylabel(' n$_e$ ($cm^{-3}$)')
ax[1,0].plot(h,omega_spk,'o-',label='Gyroresonce Frequency')
ax[1,0].set_ylabel('Frequency (MHz)')
ax[1,1].plot(h,fp,'-',label='Plasma Freqeuncy')
ax[1,1].plot(h,omega_spk,'o-',label='Gyroresonce Frequency')
ax[1,1].set_ylabel('Frequency (MHz)')
ax[1,0].set_xlabel('Coronal Height (Mm)');ax[1,1].set_xlabel('Coronal Height (Mm)')
ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()
plt.show()



dim=vx.shape

xc=680;yc=250
hmi=fits.open('/media/rohit/VLA/20130423/hmi.m_45s.2013.04.23_20_35_15_TAI.magnetogram.fits')
hmidata=hmi[0].data
hmihead=hmi[0].header
hmimap=Map('/media/rohit/VLA/20130423/hmi.m_45s.2013.04.23_20_35_15_TAI.magnetogram.fits')
hmimap.pixel_to_world(694*u.pix,1550*u.pix)

sys.exit()
#########################

nhmi=150
zz,yy,xx=np.mgrid[0:nhmi,0:nhmi,0:10]
zz1,yy1,xx1=np.mgrid[-1:1:150j, -1:1:150j, 0:1:10j]
pts = np.empty((150,150,10) + (3,), dtype=np.float64)
pts[..., 0] = xx1
pts[..., 1] = yy1
pts[..., 2] = zz1
pts.shape = pts.size // 3, 3

vectors = np.empty((150,150,10) + (1,), dtype=np.float64)
vectors[:,:,0, 0] = hmidata[619:769,1475:1625]
#vectors[:,:,5,0]=hmidata[-1200:-1000,-1200:-1000]
vectors.shape = vectors.size // 1, 1
sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)

sg.point_data.scalars = vectors[:,0]
sg.point_data.scalars.name = 'HMI'

write_data(sg, '20130323_hmi_1d.vtk')

sys.exit()

############ Magnetic Field ###########################
# Generate the grid
zz,yy,xx=np.mgrid[0:dim[0],0:dim[1],0:dim[2]]
pts = np.empty(vx.shape + (3,), dtype=np.int32)
pts[..., 0] = xx
pts[..., 1] = yy
pts[..., 2] = zz

vectors = np.empty(vx.shape + (3,), dtype=np.float64)
vectors[..., 0] = vx
vectors[..., 1] = vy
vectors[..., 2] = vz

# We reorder the points and vectors so this is as per VTK's
# requirement of x first, y next and z last.
pts = pts.transpose(2, 1, 0, 3).copy()
pts.shape = pts.size // 3, 3

vectors = vectors.transpose(2, 1, 0, 3).copy()
vectors.shape = vectors.size // 3, 3

sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)

sg.point_data.vectors = vectors
sg.point_data.vectors.name = 'Magnetic Field'

absv=np.sqrt(vectors[:,0]**2+vectors[:,1]**2+vectors[:,2]**2)
sg.point_data.scalars = absv
sg.point_data.scalars.name = 'Magnitude'

write_data(sg, '20130423_mapfull_high_resolution.vtk')




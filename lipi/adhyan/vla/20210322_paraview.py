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
# Unpack velocity information

aa=readsav('/sdata/20210322/magnetic_fields_exp_20210322.sav')

vx=aa['box'][0][0][0]#[0:50,120:-120,120:-120]#flow['vx']
vy=aa['box'][0][0][1]#[0:50,120:-120,120:-120]#flow['vy']
vz=aa['box'][0][0][2]#[0:50,120:-120,120:-120]#flow['vz']

dim=vx.shape

hmi=fits.open('/media/rohit/VLA/20160409_EUV/hmi/hmi.m_45s.2016.04.09_18_45_45_TAI.magnetogram.fits')
hmidata=hmi[0].data
hmimap=Map(hmidata,hmi[0].header)
nhmi=400
zz,yy,xx=np.mgrid[0:nhmi,0:nhmi,0:10]
zz1,yy1,xx1=np.mgrid[-1:1:400j, -1:1:400j, 0:1:200j]
pts = np.empty((400,400,200) + (3,), dtype=np.float64)
pts[..., 0] = xx1
pts[..., 1] = yy1
pts[..., 2] = zz1
pts.shape = pts.size // 3, 3

vectors = np.empty((400,400,200) + (1,), dtype=np.float64)
vectors[:,:,0, 0] = hmidata[-1200:-1000,-1200:-1000]
#vectors[:,:,5,0]=hmidata[-1200:-1000,-1200:-1000]
vectors.shape = vectors.size // 1, 1
sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)

sg.point_data.scalars = vectors[:,0]
sg.point_data.scalars.name = 'HMI'

write_data(sg, 'hmi_1d.vtk')

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

write_data(sg, 'mapfull.vtk')




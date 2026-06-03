import numpy as np
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
import os

#on py39 on bhargav

# path for the carrington rotation here: /media/rohit/ray/predsci/data/runs/
# cr2286-high | hmi_masp_mas_std_0201 hmi_mast_mas_std_0101 hmi_mast_mas_std_0201 | corona.zip helio.zip
'''HMI1:  
hmi_mast_std_0101 (Heating model 1): The 'standard' heating model from Lionello et al. 2009. 
- This model typically gives hotter QS temperatures, typically between 1.4-1.8 MK.
- Coronal holes are wider/more open.
- Active regions aren't as visible due to the properties of the magnetic maps and heating model at this resolution.
    ;
HMI2:
;hmi_mast_std_0201 (Heating model 2): used for most recent eclipse predictions.
- This model typically gives cooler QS temperatures, typically between 1.2-1.5 MK.
- The coronal base density is slightly higher everywhere.
- Coronal holes are smaller and slightly more dense than heating model 1.
- Active regions are heated more and are hotter/more visible.
HMI3:
hmi_masp_std_0201
GONG2:
also available as gong_mast_std_0201'''

# corona
# bp002.hdf  br_r0.hdf  cs002.hdf  jr002.hdf  omas      rho002.hdf  va002.hdf  vr002.hdf
# br002.hdf  bt002.hdf  jp002.hdf  jt002.hdf  p002.hdf  t002.hdf    vp002.hdf  vt002.hdf
# helio
# bp002.hdf  br_r0.hdf  cs002.hdf  jr002.hdf  omasip    rho002.hdf  t002.hdf  va002.hdf  vr002.hdf  vt002.hdf
# br002.hdf  bt002.hdf  jp002.hdf  jt002.hdf  p002.hdf  rho_r0.hdf  t_r0.hdf  vp002.hdf  vr_r0.hdf
#-----------------------
#par_cor=['bp002','br_r0','cs002','jr002','omas','rho002','t002','va002','vr002','vt002','br002','bt002','jp002.hdf  jt002','p002','vp002']
#par_hel=['bp002','br_r0','cs002','jr002','omasip','rho002','t002','va002','vr002','vt002','br002','bt002','jp002.hdf  jt002','p002','rho_r0','t_r0','vp002','vr_r0']

path = '/media/rohit/ray/predsci/data/runs/cr2286-high/hmi_masp_mas_std_0201/'
path_cor = path + 'corona/'
path_hel = path + 'helio/'
hdflist_cor = [f.split('.hd')[0] for f in os.listdir(path_cor) if f.endswith('.hdf')]
hdflist_hel = [f.split('.hd')[0] for f in os.listdir(path_hel) if f.endswith('.hdf')]
print('Reading Corona...',path)
hdfdict_cor = {f: [s.select(k)[:] for k in s.datasets()]
                       for f in hdflist_cor
                                      for s in [SD(path_cor + f + '.hdf', SDC.READ)]}
print('Reading Helio...',path)
hdfdict_hel = {f: [s.select(k)[:] for k in s.datasets()]
                       for f in hdflist_hel
                                      for s in [SD(path_hel + f + '.hdf', SDC.READ)]}
#'Dictionary index 0,1,2,3 | 3 is the quantity & others are axes
print('Parameters available:')
print('Corona:',hdfdict_cor.keys())
print('Helio:',hdfdict_hel.keys())
bp_cor=hdfdict_cor['bp002'][-1] # mas.bp
br_cor=hdfdict_cor['br002'][-1] # mas.br
bt_cor=hdfdict_cor['bt002'][-1] # mas.bt
jp_cor=hdfdict_cor['jp002'][-1]
jr_cor=hdfdict_cor['jr002'][-1]
jt_cor=hdfdict_cor['jt002'][-1]
va_cor=hdfdict_cor['va002'][-1]
vr_cor=hdfdict_cor['vr002'][-1] # mas.vr
vp_cor=hdfdict_cor['vp002'][-1] # mas.vp
vt_cor=hdfdict_cor['vt002'][-1] # mas.vt
print('Reading rho')
rho_cor=hdfdict_cor['rho002'][-1] # mas.dens
t_cor=hdfdict_cor['t002'][-1] # mas.theta
p_cor=hdfdict_cor['p002'][-1] # mas.phase
cs_cor=hdfdict_cor['cs002'][-1] 
br_r0_cor=hdfdict_cor['br_r0'][-1]

dx=0.0;dy=0.0;xcenter=0;ycenter=0;xxmin=0.0;xxmax=0.0
yymin=0.0;yymax=0.0;rheight=0.0

cmer = -90 # central meridian longitude of observer
phio = 0 # shift to impose on model structure -- note this
#               allows us to make a Carrington map in fixed
#               Carrington coordinates and move the structure around, or
#               fit its location to data
ngrid = 256
dang=180.0/ngrid
nx=int(2.0*ngrid)
ny=int(ngrid)
nxny = np.meshgrid(np.arange(nx),np.arange(ny))
rheight=1
r=nxny[0]*0 + rheight
th=nxny[0]*0
ph=nxny[0]*0
mdtor = np.pi/180
for i in range(ny):
    ph[i]=i*dang + (cmer-phio)

ph = ph/np.max(ph)*np.pi
if(dang<1):
    dstop=np.round(1.0/dang)
if(dang>1):
    dstop=1.0

west=0;east=1;cmer=1
if(west):
    for j in range(int(ny-dstop)):
        th[:,j] = 180.0 + j*dang

if(east):
    for j in range(int(ny-dstop)):
        th[:,j] = 180.0 - j*dang

if(cmer==1):
    for j in range(int(ny-dstop)):
        #th[:,j] = 0
        r[:,j] = rheight*(np.cos((180.0-j*dang)))

#th=th*mdtor
#ph=ph*mdtor
xcenter=180.0
ycenter=90.0

xxmin=0.0
xxmax=360.0
yymin=0.0
yymax=180.0

dx = (xxmax-xxmin)/(nx-1)
dy = (yymax-yymin)/(ny-1)
dx = 360.0/(nx-1) # pixel size in degrees, longitude
dy = 180.0/(ny-1) # pixel size in degrees, latitude

xcenter=(xxmax+xxmin)/2.0
ycenter=(yymax+yymin)/2.0





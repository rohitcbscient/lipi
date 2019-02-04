import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
import pickle
from scipy import signal
from surya.utils import main as ut
import glob

def get_forward(list_):
    freq=[0]*len(list_)
    Tb=[0]*len(list_)
    flux=[0]*len(list_)
    ll=0
    for l in list_:
        data=readsav(l)
        freq[ll]=float(l.split('_')[1].split('MHz')[0])
        dx=data['quantmap'][0][4]*16*60
        Tb_=data['stokesstruct'][0][0]
        size=dx*Tb_.shape[0]
        Tb[ll]=np.mean(Tb_)
        flux[ll]=np.mean(ut.Tb2flux(Tb[ll], size, size, freq[ll]/1000.))
        ll=ll+1
    return freq,flux,Tb

dir_='/home/i4ds1807205/20151203/'
data_100MHz=pickle.load(open(dir_+'Tb_1133149192-084-085.p','r'))
data_160MHz=pickle.load(open(dir_+'Tb_1133149192-125-126.p','r'))
data_240MHz=pickle.load(open(dir_+'Tb_1133149192-187-188.p','r'))

mwa_100MHz=np.mean(np.array(data_100MHz[0]),axis=0)
mwa_160MHz=np.mean(np.array(data_160MHz[0]),axis=0)
mwa_240MHz=np.mean(np.array(data_240MHz[0])[:-2,:,:],axis=0)

bmin_100MHz=data_100MHz[3]*60 # In arcmin
bmax_100MHz=data_100MHz[4]*60 # In arcmin
bmin_160MHz=data_160MHz[3]*60 # In arcmin
bmax_160MHz=data_160MHz[4]*60 # In arcmin
bmin_240MHz=data_240MHz[3]*60 # In arcmin
bmax_240MHz=data_240MHz[4]*60 # In arcmin

flist=[108,120,132,145,161,179,196,217,240]
mwa_flux=[2.52,3.45,4.50,5.49,6.59,8.32,11.52,16.0,18.64]
emwa_flux=[0.34,0.38,1.30,0.57,0.68,0.86,1.19,1.64,1.93]

list_=sorted(glob.glob(dir_+'*psimas.sav'))
fall,fluxall,Tball=get_forward(list_)

forward_100MHz=readsav(dir_+'20151203_100MHz_psimas.sav')
forward_160MHz=readsav(dir_+'20151203_160MHz_psimas.sav')
forward_240MHz=readsav(dir_+'20151203_240MHz_psimas.sav')

pres=forward_100MHz['modsolstruct'][0][0]
dens=forward_100MHz['modsolstruct'][0][1]
temp=forward_100MHz['modsolstruct'][0][2]
br=forward_100MHz['modsolstruct'][0][3]
bth=forward_100MHz['modsolstruct'][0][4]
bph=forward_100MHz['modsolstruct'][0][5]
vr=forward_100MHz['modsolstruct'][0][6]
vth=forward_100MHz['modsolstruct'][0][7]
vph=forward_100MHz['modsolstruct'][0][8]

forward_time=forward_100MHz['quantmap'][0][5]
forward_dx=forward_100MHz['quantmap'][0][3]*16*60 # In arcsec
forward_dy=forward_100MHz['quantmap'][0][4]*16*60 # In arcsec; originally in Rsun

I_100MHz=forward_100MHz['stokesstruct'][0][0]
I_160MHz=forward_160MHz['stokesstruct'][0][0]
I_240MHz=forward_240MHz['stokesstruct'][0][0]

fwhmx_100MHz=int(bmin_100MHz*60/forward_dx)
fwhmy_100MHz=int(bmax_100MHz*60/forward_dy)
fwhmx_160MHz=int(bmin_160MHz*60/forward_dx)
fwhmy_160MHz=int(bmax_160MHz*60/forward_dy)
fwhmx_240MHz=int(bmin_240MHz*60/forward_dx)
fwhmy_240MHz=int(bmax_240MHz*60/forward_dy)

beam_101MHz=ut.makeGaussian(60, fwhmx_100MHz, fwhmy_100MHz , center=None)
beam_160MHz=ut.makeGaussian(60, fwhmx_160MHz, fwhmy_160MHz , center=None)
beam_240MHz=ut.makeGaussian(60, fwhmx_240MHz, fwhmy_240MHz , center=None)

I_100MHz_convolved = signal.convolve(I_100MHz,beam_101MHz, mode='same')/np.sum(beam_101MHz)
I_160MHz_convolved = signal.convolve(I_160MHz,beam_160MHz, mode='same')/np.sum(beam_160MHz)
I_240MHz_convolved = signal.convolve(I_240MHz,beam_240MHz, mode='same')/np.sum(beam_240MHz)

mwa_x=np.linspace(-2500,2500,100)
mwa_y=np.linspace(-2500,2500,100)

forward_x=np.linspace(-1440,1400,256)
forward_y=np.linspace(-1440,1400,256)

id_l=ut.find_nearest(mwa_x,-1440)[0]
id_r=ut.find_nearest(mwa_x,1440)[0]
del_id=id_r-id_l
nforward_100MHz=np.zeros((del_id,del_id))
for i in range(del_id):
    for j in range(del_id):
        nforward_100MHz[i,j]=I_100MHz_convolved[ut.find_nearest(mwa_x,id_l+i)[0],ut.find_nearest(mwa_y,id_l+j)[0]]

############### Optical depth        

id240=np.where(mwa_240MHz==np.max(mwa_240MHz))
id160=np.where(mwa_160MHz==np.max(mwa_160MHz))
id100=np.where(mwa_100MHz==np.max(mwa_100MHz))
maxTb_240MHz=np.mean(mwa_240MHz[id240[0][0]-2:id240[0][0]+2,id240[1][0]-2:id240[1][0]+2])
maxTb_160MHz=np.mean(mwa_160MHz[id160[0][0]-2:id160[0][0]+2,id160[1][0]-2:id160[1][0]+2])
maxTb_100MHz=np.mean(mwa_100MHz[id100[0][0]-2:id100[0][0]+2,id100[1][0]-2:id100[1][0]+2])

fw240=np.where(I_240MHz_convolved==np.max(I_240MHz_convolved))
fw160=np.where(I_160MHz_convolved==np.max(I_160MHz_convolved))
fw100=np.where(I_100MHz_convolved==np.max(I_100MHz_convolved))
maxfw_240MHz=np.mean(I_240MHz_convolved[fw240[0][0]-4:fw240[0][0]+4,fw240[1][0]-4:fw240[1][0]+4])
maxfw_160MHz=np.mean(I_160MHz_convolved[fw160[0][0]-4:fw160[0][0]+4,fw160[1][0]-4:fw160[1][0]+4])
maxfw_100MHz=np.mean(I_100MHz_convolved[fw100[0][0]-4:fw100[0][0]+4,fw100[1][0]-4:fw100[1][0]+4])

tau240MHz=np.log((maxfw_240MHz/maxTb_240MHz)-1)
tau160MHz=np.log((maxfw_160MHz/maxTb_160MHz)-1)
tau100MHz=np.log((maxfw_100MHz/maxTb_100MHz)-1)

mwa=np.array([mwa_240MHz,mwa_160MHz,mwa_100MHz])

fwdx=forward_100MHz['quantmap'][0][3]*16*60 # Rsun to arcsec
fwdy=forward_100MHz['quantmap'][0][4]*16*60 

size=fwdx*I_100MHz.shape[0]
Tb_100MHz=np.mean(I_100MHz_convolved)
Tb_160MHz=np.mean(I_160MHz_convolved)
Tb_240MHz=np.mean(I_240MHz_convolved)
fw_flux_100MHz=ut.Tb2flux(Tb_100MHz, size, size, 0.1)
fw_flux_160MHz=ut.Tb2flux(Tb_160MHz, size, size, 0.16)
fw_flux_240MHz=ut.Tb2flux(Tb_240MHz, size, size, 0.24)


print '##### FORWARD ###'
print 100,np.mean(fw_flux_100MHz),np.sqrt(size*size/np.pi)/60,np.mean(Tb_100MHz)/1.e6
print 160,np.mean(fw_flux_160MHz),np.sqrt(size*size/np.pi)/60,np.mean(Tb_160MHz)/1.e6
print 240,np.mean(fw_flux_240MHz),np.sqrt(size*size/np.pi)/60,np.mean(Tb_240MHz)/1.e6
mwa_size=len(np.where(mwa_240MHz!=0)[0])
print '##### MWA ####'
print 100,mwa_flux[0],np.sqrt(len(np.where(mwa_100MHz!=0)[0]))*50/60
print 160,mwa_flux[4],np.sqrt(len(np.where(mwa_160MHz!=0)[0]))*50/60
print 240,mwa_flux[-1],np.sqrt(len(np.where(mwa_240MHz!=0)[0]))*50/60


plot_mlab=0
if(plot_mlab==1):
    z,x,y=np.mgrid[0:2:3j,-2500:2500:100j,-2500:2500:100j]
    iso=mlab.contour3d(z,x,y,mwa,vmin=mwa.min(),vmax=mwa.max(),opacity=0.3)
    iso.contour.number_of_contours = 15
    mlab.show()

############### PLOT
plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')

plot_fluxes=0
if(plot_fluxes==1):
    plt.plot(fall,np.array(fluxall),'o-',color='red',label='FORWARD')
    plt.errorbar(flist,mwa_flux,emwa_flux,color='blue')
    plt.plot(flist,mwa_flux,'o',color='blue',label='MWA')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend()
    plt.xlim([100,250])
    plt.show()

plot_composite=0
if(plot_composite==1):
    levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    f, ax = plt.subplots(2, 3, figsize=(20,10),sharex=True)
    im00=ax[0,0].imshow(mwa_100MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im00,ax=ax[0,0],label='T$_{B}$(MK)')
    ax[0,0].contour(mwa_100MHz/np.max(mwa_100MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,0].set_title('100 MHz')
    im10=ax[1,0].imshow(I_100MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440,1440,-1440,1440],origin=True)
    f.colorbar(im10,ax=ax[1,0],label='T$_{B}$(MK)')
    ax[1,0].contour(I_100MHz_convolved/np.max(I_100MHz_convolved),levels=levels_,extent=[-1440,1440,-1440,1440],colors='red')
    ax[1,0].set_xlim([-2500,2500])
    ax[1,0].set_ylim([-2500,2500])
    ax[0,0].grid(True)
    ax[1,0].grid(True)
    ax[1,0].set_xlabel('arcsec')
    ax[1,0].set_ylabel('arcsec')
    ax[0,0].set_ylabel('arcsec')
    ##
    im01=ax[0,1].imshow(mwa_160MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im01,ax=ax[0,1],label='T$_{B}$(MK)')
    ax[0,1].contour(mwa_160MHz/np.max(mwa_160MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,1].set_title('160 MHz')
    im10=ax[1,1].imshow(I_160MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440,1440,-1440,1440],origin=True)
    f.colorbar(im10,ax=ax[1,1],label='T$_{B}$(MK)')
    ax[1,1].contour(I_160MHz_convolved/np.max(I_160MHz_convolved),levels=levels_,extent=[-1440,1440,-1440,1440],colors='red')
    ax[1,1].set_xlim([-2500,2500])
    ax[1,1].set_ylim([-2500,2500])
    ax[0,1].grid(True)
    ax[1,1].grid(True)
    ax[1,1].set_xlabel('arcsec')
    ##
    im02=ax[0,2].imshow(mwa_240MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im02,ax=ax[0,2],label='T$_{B}$(MK)')
    ax[0,2].contour(mwa_240MHz/np.max(mwa_240MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,2].set_title('240 MHz')
    im12=ax[1,2].imshow(I_240MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440,1440,-1440,1440],origin=True)
    f.colorbar(im12,ax=ax[1,2],label='T$_{B}$(MK)')
    ax[1,2].contour(I_240MHz_convolved/np.max(I_240MHz_convolved),levels=levels_,extent=[-1440,1440,-1440,1440],colors='red')
    ax[1,2].set_xlim([-2500,2500])
    ax[1,2].set_ylim([-2500,2500])
    ax[0,2].grid(True)
    ax[1,2].grid(True)
    ax[1,2].set_xlabel('arcsec')
    f.tight_layout()
    f.show()

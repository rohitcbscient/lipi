import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import pickle
import numpy as np
import glob
from scipy.io import readsav
from mpl_toolkits.mplot3d import proj3d
from surya.utils import main as ut
from surya.plot import main as pt
from scipy import signal
from astropy.io import fits
import matplotlib as mpl
from scipy import interpolate
from matplotlib.colors import LogNorm
from matplotlib import patches
from surya.radio import get_maps as tb


def fit_cubic(y):
    x=np.arange(y.shape[0])
    ynew=np.zeros((y.shape[0],y.shape[1],y.shape[2]))
    snew=np.zeros((y.shape[0],y.shape[1],y.shape[2]))
    for i in range(y.shape[1]):
        for j in range(y.shape[2]):
            f = np.polyfit(x, y[:,i,j], 3)
            p = np.poly1d(f)
            snew[:,i,j]=p(x)
            ynew[:,i,j]=y[:,i,j]-p(x)
    return ynew,snew

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')
#freq=[109.0,121.0,134.0,147.0,162.0,180.0,198.0,218.0,241.0]
freq=[109.0,134.0,147.0,162.0,180.0,198.0,218.0,241.0]
#fband=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
fband=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
data=[0]*len(freq)
udata_dirty=[0]*len(freq)
udata=[0]*len(freq)
Tb=[0]*len(freq)
Tb_all=[0]*len(freq)
sTb_all=[0]*len(freq)
diff_Tb_all=[0]*len(freq)
Tb_ch=[0]*len(freq)
Tb_ar=[0]*len(freq)
Tb_qs=[0]*len(freq)
udata_all=[0]*len(freq)
diff_udata_all=[0]*len(freq)
sudata_all=[0]*len(freq)
udata_dirty_all=[0]*len(freq)
udata_ch=[0]*len(freq)
udata_ar=[0]*len(freq)
udata_qs=[0]*len(freq)
udata_ch_dirty=[0]*len(freq)
udata_ar_dirty=[0]*len(freq)
udata_qs_dirty=[0]*len(freq)
udata_ch_dirty_rms=[0]*len(freq)
udata_ar_dirty_rms=[0]*len(freq)
udata_qs_dirty_rms=[0]*len(freq)
Tb_rms=[0]*len(freq)
Tbmax=[0]*len(freq)
Tbrmsmax=[0]*len(freq)
bmin=[0]*len(freq)
bmax=[0]*len(freq)
flux=[0]*len(freq)
flux_std=[0]*len(freq)
flux_ar=[0]*len(freq)
flux_ch=[0]*len(freq)
flux_qs=[0]*len(freq)
flux_ar_res=[0]*len(freq)
flux_ch_res=[0]*len(freq)
flux_qs_res=[0]*len(freq)
nonsun_rms=[0]*len(freq)
sun_rms=[0]*len(freq)
sun_rms_time=[0]*len(freq)
sun_res_rms=[0]*len(freq)
SEFD=[0]*len(freq)
NEFD=[0]*len(freq)
snoise_F=[0]*len(freq)
snoise_T=[0]*len(freq)
nnoise=[0]*len(freq)
fact=[0]*len(freq)
print 'Reading files...'

Tsky_=np.array([665,483,371,290,225,184,179,181,122,94])
Trec_=np.array([30,28,26,24,21,20,21,23,27,32])
f=np.array([103,117,131,148,167,189,213,240,272,299])
finters=interpolate.interp1d(f, Tsky_,kind='cubic')
finterr=interpolate.interp1d(f, Trec_,kind='cubic')
Tsky=finters(freq)
Trec=finterr(freq)
baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
flux_path='/media/rohit/VLA/20151203_MWA/pickle/flux_V1_1133149192-%b'
img_path='/media/rohit/VLA/20151203_MWA/images_all/1133149192-%b'
outdir='/media/rohit/VLA/20151203_MWA/Tb_new/'
centre_file='/home/i4ds1807205/20151203/20151203_nasa_horizon.dat'
ids='1133149192'

get_flux=1
if(get_flux):
    ff=0
    Ssun_all,time_str,timesec=tb.mean_flux(flux_path,fband[ff],baseline_filelist,0.5)
    Ssun_mean=np.mean(Ssun_all,axis=(0,1))
    Ssun_std=np.std(Ssun_all,axis=(0,1))
    

get_Tb=1
if(get_Tb==1):
    print 'Starting Tb computation..'
    del_=50
    angle=-15.4 # Negative angle imples clockwise rotation and vice versa
    res=50 # arcsec
    for ff in range(len(fband)):
        print str(freq[ff])+' MHz'
        imglist=sorted(glob.glob(img_path+fband[ff]+'*.image.FITS'))
        reslist=sorted(glob.glob(img_path+fband[ff]+'*.residual.FITS'))
        img_time=[0]*len(imglist)
        Tb=[0]*len(imglist)
        resi=[0]*len(imglist)
        resi_sun=[0]*len(imglist)
        Tb_resi=[0]*len(imglist)
        flux_resi=[0]*len(imglist)
        Tb_sun=[0]*len(imglist)
        flux_sun=[0]*len(imglist)
        Tb_sun_resi=[0]*len(imglist)
        Tb_fac=[0]*len(imglist)
        flux_sun_resi=[0]*len(imglist)
        flux_fac=[0]*len(imglist)
        polTb=[0]*len(imglist)
        flux=[0]*len(imglist)
        bmaj=[0]*len(imglist)
        bmin=[0]*len(imglist)
        bpa=[0]*len(imglist)
        ndata=[0]*len(imglist)
        offsun_mean=[0]*len(imglist)
        offsun_std=[0]*len(imglist)
        centre=[0]*len(imglist)
        j=0
        for fitsfile in imglist:
            print j
            xc,yc,img_time[j]=tb.solar_center_pixel(fitsfile,centre_file)
            centre[j]=[xc,yc]
            Tb[j],flux[j],offsun_mean[j],offsun_std[j],bmaj[j],bmin[j],bpa[j],ndata[j]=tb.compute_Tb(fitsfile,xc,yc,del_,angle,res,freq[ff]*1.e6,5,Ssun_mean[j])
            resi_sun[j],resi[j]=tb.get_residuals(fitsfile,xc,yc,del_,angle)
            polTb[j]=ut.cart2polar(Tb[j])
            Tb_sun_resi[j],flux_sun_resi[j],Tb_fac[j],flux_fac[j]=tb.scale_residuals(resi_sun[j],flux[j],ndata[j],res)
            Tb_sun[j]=resi[j]*Tb_fac[j]
            flux_sun[j]=resi[j]*flux_fac[j]
            j=j+1
        Tb=np.array(Tb)
        Tb_fac=np.array(Tb_fac)
        flux_fac=np.array(flux_fac)
        Tb_sun=np.array(Tb_sun)
        Tb_sun_resi=np.array(Tb_sun_resi)
        Tb_resi=np.array(Tb_resi)
        flux=np.array(flux)
        flux_fac=np.array(flux_fac)
        flux_sun=np.array(flux_sun)
        flux_resi=np.array(flux_resi)
        flux_sun_resi=np.array(flux_sun_resi)
        resi=np.array(resi)
        resi_sun=np.array(resi_sun)
        polTb=np.array(polTb)
        flux=np.array(flux)
        ndata=np.array(ndata)
        print 'Mean Tb: ',np.max(Tb[0])
        pickle.dump([Tb,offsun_mean,offsun_std,bmaj,bmin,bpa,centre,flux,ndata,polTb],open(outdir+'Tb_'+str(ids)+'-'+str(fband[ff])+'.p','w'))
        pickle.dump([resi,resi_sun,Tb_sun_resi,flux_sun_resi,Tb_fac,flux_fac],open(outdir+'res_'+str(ids)+'-'+str(fband[ff])+'.p','w'))



Tb=[0]*len(fband)
flux=[0]*len(fband)
bmaj=[0]*len(fband)
bmin=[0]*len(fband)
bpa=[0]*len(fband)
ndata=[0]*len(fband)
offsun_mean=[0]*len(fband)
offsun_std=[0]*len(fband)
centre=[0]*len(fband)
for ff in range(len(fband)):
    Tb[ff],offsun_mean[ff],offsun_std[ff],bmaj[ff],bmin[ff],bpa[ff],centre[ff],flux[ff]=pickle.load(open(outdir+'Tb_'+str(ids)+'-'+str(fband[ff])+'.p'),'r')

sys.exit()
for i in range(len(freq)):
    print fband[i]
    #data[i]=pickle.load(open('/media/rohit/VLA/20151203_MWA/Tb/images_all/Tb_1133149192-'+fband[i]+'.p','r'))
    data[i]=pickle.load(open('/home/i4ds1807205/20151203/images_all/Tb_1133149192-'+fband[i]+'.p','r'))
    #udata[i],udata_dirty[i]=pickle.load(open('/media/rohit/VLA/20151203_MWA/Tb/images_all/res_1133149192-'+fband[i]+'.p','r'))
    udata[i],udata_dirty[i]=pickle.load(open('/home/i4ds1807205/20151203/images_all/res_1133149192-'+fband[i]+'.p','r'))
    # dirty means residual, udata is the image
    Tb_all[i]=np.array(data[i][0])[0:577] 
    udata_all[i]=np.array(udata[i][0:577])
    diff_Tb_all[i],sTb_all[i]=fit_cubic(Tb_all[i])
    diff_udata_all[i],sudata_all[i]=fit_cubic(udata_all[i])
    udata_mean=np.mean(udata_all[i],axis=0)
    udata_dirty_all[i]=np.array(udata_dirty[i][0:577]) 
    nonsun_rms[i]=np.std(np.array(data[i][2][0:577])) # Time non-solar RMS in Jy/beam
    sun_rms[i]=np.std(udata_all[i],axis=(1,2)) # Solar Spatial RMS in Jy/beam 
    sun_rms_time[i]=np.std(udata_all[i],axis=0) # Solar time RMS in Jy/beam 
    sun_res_rms[i]=np.std(udata_dirty_all[i],axis=(1,2)) # Spatial RMS in residual Jy/beam
    Tb[i]=np.mean(Tb_all[i],axis=0)
    fact[i]=np.sum(Tb_all[i][0])/np.sum(udata_mean)
    Tb_rms[i]=np.std(Tb_all[i],axis=0)
    Tb[i][np.isnan(Tb[i])]=0
    Tbmax[i]=Tb[i].max()
    Tbrmsmax[i]=Tb_rms[i].max()
    Tb[i][np.where(Tb[i]==0)]=np.nan
    Tb_ch[i]=np.mean(Tb_all[i][:,55:61,48:54],axis=(1,2))
    Tb_ar[i]=np.mean(Tb_all[i][:,41:46,58:63],axis=(1,2))
    Tb_qs[i]=np.mean(Tb_all[i][:,40:45,42:47],axis=(1,2))
    udata_ch[i]=np.mean(udata_all[i][:,55:61,48:54],axis=(1,2))*fact[i]
    udata_ar[i]=np.mean(udata_all[i][:,41:46,58:63],axis=(1,2))*fact[i]
    udata_qs[i]=np.mean(udata_all[i][:,40:45,42:47],axis=(1,2))*fact[i]
    udata_ch_dirty[i]=np.mean(udata_dirty_all[i][:,55:61,48:54],axis=(1,2))*fact[i]
    udata_ar_dirty[i]=np.mean(udata_dirty_all[i][:,41:46,58:63],axis=(1,2))*fact[i]
    udata_qs_dirty[i]=np.mean(udata_dirty_all[i][:,40:45,42:47],axis=(1,2))*fact[i]
    udata_ch_dirty_rms[i]=np.std(udata_dirty_all[i][:,55:61,48:54],axis=(1,2))*fact[i]
    udata_ar_dirty_rms[i]=np.std(udata_dirty_all[i][:,41:46,58:63],axis=(1,2))*fact[i]
    udata_qs_dirty_rms[i]=np.std(udata_dirty_all[i][:,40:45,42:47],axis=(1,2))*fact[i]
    bmin[i]=data[i][3]*3600
    bmax[i]=data[i][4]*3600
    pick=pickle.load(open('/home/i4ds1807205/20151203/pickle/flux_V1_1133149192-%b'+fband[i]+'_T000-009.p','rb'))
    pick=pick[17][3][0]
    pick[np.isnan(pick)]=0
    flux[i]=np.mean(pick)
    flux_std[i]=np.std(pick)
    flux_ar[i]=ut.Tb2flux(Tb_ar[i], bmin[i], bmax[i], freq[i]/1000.)
    flux_qs[i]=ut.Tb2flux(Tb_qs[i], bmin[i], bmax[i], freq[i]/1000.)
    flux_ch[i]=ut.Tb2flux(Tb_ch[i], bmin[i], bmax[i], freq[i]/1000.)
    flux_ar_res[i]=ut.Tb2flux(udata_ar_dirty_rms[i], bmin[i], bmax[i], freq[i]/1000.)
    flux_qs_res[i]=ut.Tb2flux(udata_qs_dirty_rms[i], bmin[i], bmax[i], freq[i]/1000.)
    flux_ch_res[i]=ut.Tb2flux(udata_ch_dirty_rms[i], bmin[i], bmax[i], freq[i]/1000.)
    #### Self Noise ####
    N=128
    dnu=2.e6 # in Hz
    dtau=0.5 # in Sec
    M=np.sqrt(N*(N-1)*dnu*dtau)
    Aeff=21.5 # in m**2
    SEFD[i]=0.5*1.38*1.e-23*1.e22*Tsky[i]/Aeff # in SFU
    NEFD[i]=0.5*1.38*1.e-23*1.e22*Trec[i]/Aeff # in SFU
    SEFD[i]=SEFD[i]/M
    SEFD[i]=NEFD[i]/M
    snoise_T[i]=np.mean(Tb[i][np.isfinite(Tb[i])])/np.sqrt(2*64*40000*0.5)
    snoise_F[i]=flux[i]/((Tb[i][np.isfinite(Tb[i])].shape[0])*np.sqrt(2*64*40000*0.5))
    


Tb=np.array(Tb)
nonsun_rms=np.array(nonsun_rms)
#nonsun_rms_mean=np.mean(nonsun_rms,axis=1)
sun_rms=np.array(sun_rms)
sun_rms_time=np.array(sun_rms_time)
sun_res_rms=np.array(sun_res_rms)
#sun_res_rms_mean=np.mean(sun_res_rms,axis=1)
#sun_rms_mean=np.mean(sun_rms,axis=1)
Tb_all=np.array(Tb_all)
Tb_ch=np.array(Tb_ch)
flux_ch=np.array(flux_ch)
flux_ar=np.array(flux_ar)
flux_qs=np.array(flux_qs)

diff_udata_all=np.array(diff_udata_all)
diff_Tb_all=np.array(diff_Tb_all)
sTb_all=np.array(sTb_all)
sudata_all=np.array(sudata_all)
udata_ch_dirty=np.array(udata_ch_dirty)
udata_ar_dirty=np.array(udata_ar_dirty)
udata_qs_dirty=np.array(udata_qs_dirty)
udata_ch_dirty_rms=np.array(udata_ch_dirty_rms)
udata_ar_dirty_rms=np.array(udata_ar_dirty_rms)
udata_qs_dirty_rms=np.array(udata_qs_dirty_rms)
udata_ch=np.array(udata_ch)
udata_all=np.array(udata_all)
udata_dirty_all=np.array(udata_dirty_all)
udata_ar=np.array(udata_ar)
udata_qs=np.array(udata_qs)
udata_ch_mean=np.mean(udata_ch,axis=1)
udata_qs_mean=np.mean(udata_qs,axis=1)
udata_ar_mean=np.mean(udata_ar,axis=1)
udata_ch_rms=np.std(udata_ch,axis=1)
udata_qs_rms=np.std(udata_qs,axis=1)
udata_ar_rms=np.std(udata_ar,axis=1)
udata_ch_res_rms_t=np.std(udata_ch_dirty,axis=1)
udata_qs_res_rms_t=np.std(udata_qs_dirty,axis=1)
udata_ar_res_rms_t=np.std(udata_ar_dirty,axis=1)
udata_ch_res_mean=np.std(udata_ch_dirty_rms,axis=1)
udata_qs_res_mean=np.std(udata_qs_dirty_rms,axis=1)
udata_ar_res_mean=np.std(udata_ar_dirty_rms,axis=1)
udata_ch_res_rms_st=np.std(udata_ch_dirty_rms,axis=1)
udata_qs_res_rms_st=np.std(udata_qs_dirty_rms,axis=1)
udata_ar_res_rms_st=np.std(udata_ar_dirty_rms,axis=1)

Tb_ch_mean=np.mean(Tb_ch,axis=1)
Tb_qs_mean=np.mean(Tb_qs,axis=1)
Tb_ar_mean=np.mean(Tb_ar,axis=1)
Tb_ch_rms=np.std(Tb_ch,axis=1)
Tb_qs_rms=np.std(Tb_qs,axis=1)
Tb_ar_rms=np.std(Tb_ar,axis=1)
flux_ch_mean=np.mean(flux_ch,axis=1)
flux_qs_mean=np.mean(flux_qs,axis=1)
flux_ar_mean=np.mean(flux_ar,axis=1)
flux_ch_rms=np.std(flux_ch,axis=1)
flux_qs_rms=np.std(flux_qs,axis=1)
flux_ar_rms=np.std(flux_ar,axis=1)
#Tb_ch_rms_dirty=np.std(Tb_ch_dirty,axis=1)
#Tb_qs_rms_dirty=np.std(Tb_qs_dirty,axis=1)
#Tb_ar_rms_dirty=np.std(Tb_ar_dirty,axis=1)


#### EUV #####
aiafile='/home/i4ds1807205/20151203/saia_00193_fd_20151203_192129.fts'
h=fits.open(aiafile)
d=h[0].data
x1=np.arange(d.shape[0])
y1=np.arange(d.shape[1])
X1,Y1=np.meshgrid(x1,y1)
X11=(X1-d.shape[0]*0.5)*h[0].header['CDELT1']
Y12=(Y1-d.shape[1]*0.5)*h[0].header['CDELT2']

euv_map=0
if(euv_map==1):
    for i in range(9):
        mwa=Tb[i]
        fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=False,figsize=(9,8))
        im1=plt.imshow(d,cmap='hot',origin='lower',extent=[X11[0][0],X11[0][-1],Y12[:,0][0],Y12[:,0][-1]],norm=LogNorm(vmin=50,vmax=8000))
        beam_centre=[-1400,-1400]
        e1 = patches.Ellipse((beam_centre[0], beam_centre[1]), bmin[i], bmax[i],angle=0, linewidth=2, fill=False, zorder=2,color='magenta')
        ax1.add_patch(e1)
        x2=np.arange(mwa.shape[0])
        y2=np.arange(mwa.shape[1])
        X2,Y2=np.meshgrid(x2,y2)
        X21=(X2-mwa.shape[0]*0.5)*50
        Y21=(Y2-mwa.shape[1]*0.5)*50
        #levels=np.array([0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])*np.max(mwa)
        levels=np.array([0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.9])*np.max(mwa[np.isfinite(mwa)])
        plt.contour(X21,Y21,mwa, levels, hold='on', colors='cyan',linewidths=2)
        plt.xlim([-1800,1800])
        plt.ylim([-1800,1800])
        plt.grid(True)
        plt.xlabel('X (arcsecs)')
        plt.ylabel('Y (arcsecs)')
        #plt.title('239 MHz')
        plt.annotate(str(freq[i])+' MHz', xy=(-300, 1400),  xycoords='data',xytext=(-300, 1400))
        plt.show()



plot_Tb=0
if(plot_Tb):
    for i in range(9):
        fig,ax=plt.subplots()
        plt.plot(np.arange(577)*0.5,np.mean(Tb_all[i,:,41:46,58:63],axis=(1,2))/1.e6,'o-')
        plt.ylabel('T$_B$ (MK)')
        plt.xlabel('Time (sec)')
        plt.title(str(freq[i])+' MHz')
        plt.show()

plot_spectrum=0
if(plot_spectrum):
    freq1=np.hstack((freq[0],freq[3:-1]))
    Tb_all1=np.concatenate((Tb_all[0].reshape(1,577,100,100),Tb_all[3:-1]))
    plt.plot(freq1,np.mean(Tb_all1[:,:,41:46,58:63],axis=(1,2,3))/1.e6,'o',color='red')
    plt.errorbar(freq1,np.mean(Tb_all1[:,:,41:46,58:63],axis=(1,2,3))/1.e6,yerr=np.std(Tb_all1[:,:,41:46,58:63],axis=(1,2,3))/1.e6,color='red',label='Active Region')
    plt.plot(freq1,np.mean(Tb_all1[:,:,40:45,42:47],axis=(1,2,3))/1.e6,'o',color='green')
    plt.errorbar(freq1,np.mean(Tb_all1[:,:,40:45,42:47],axis=(1,2,3))/1.e6,yerr=np.std(Tb_all1[:,:,40:45,42:47],axis=(1,2,3))/1.e6,color='green',label='Quiet Sun')
    plt.plot(freq1,np.mean(Tb_all1[:,:,55:61,48:54],axis=(1,2,3))/1.e6,'o',color='blue')
    plt.errorbar(freq1,np.mean(Tb_all1[:,:,55:61,48:54],axis=(1,2,3))/1.e6,yerr=np.std(Tb_all1[:,:,55:61,48:54],axis=(1,2,3))/1.e6,color='blue',label='Coronal Hole')
    plt.legend(loc=2)
    plt.ylabel('T$_B$ (MK)')
    plt.xlabel('Frequency (MHz)')
    plt.show()

plot_make_image=0
if(plot_make_image):
    for i in range(9):
        fig,ax=plt.subplots()
        im=ax.imshow(Tb[i]/1.e6,origin=True,extent=[-41.67,41.67,-41.67,41.67],vmax=0.5,vmin=0.1)
        fig.colorbar(im,label='T$_{B}$ (MK)')
        ax.contour(Tb[i]/Tb[i][np.isfinite(Tb[i])].max(),extent=[-41.67,41.67,-41.67,41.67],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],colors='k')
        ax.set_title(str(freq[i])+' MHz')
        ax.set_xlabel('X (arcmin)')
        ax.set_ylabel('Y (arcmin)')
        pt.add_beam(ax,-30,30,bmin[i]/60.,bmax[i]/60.,0)
        pt.make_circle(ax,0.5,0.5,50,100,16)
        pt.make_circle(ax,0.5,0.5,50,100,32)
        ax.grid(True)
        plt.show()

plot_spatial_distn=1
if(plot_spatial_distn):
    for i in range(9):
        sig_im=diff_udata_all[i][0,30:70,30:60].flatten()
        sig_im=sig_im[np.where(sig_im!=0)]*fact[i]
        im=udata_all[:,30:70,30:60].mean()
        res=udata_dirty_all[i][0,30:70,30:60].flatten()*fact[i]
        res=res[np.where(res!=0)]
        plt.hist(sig_im,bins=70,histtype='step',linewidth=3,normed=1,label='IMAGE')
        plt.hist(res,bins=70,histtype='step',linewidth=3,normed=1,label='RESIDUAL')
        plt.title(str(freq[i])+' MHz')
        plt.legend()
        plt.xlabel('T$_B$ (K)')
        print freq[i], np.std(res), np.std(sig_im), np.std(sig_im)-np.std(res) 
        plt.close()

plot_time_distn=1
if(plot_time_distn):
    for i in range(9):
        res=udata_dirty_all[i][:,30:70,30:60].flatten()*fact[i]
        res=res[np.where(res!=0)]
        sig_ch=diff_udata_all[i][:,55:61,48:54].flatten()
        sig_ch=sig_ch[np.where(sig_ch!=0)]
        sig_ar=diff_udata_all[i][:,41:46,58:63].flatten()
        sig_ar=sig_ar[np.where(sig_ar!=0)]
        sig_qs=diff_udata_all[i][:,40:45,42:47].flatten()
        sig_qs=sig_qs[np.where(sig_qs!=0)]
        sig_all=diff_udata_all[i][:,30:70,30:60].flatten()
        sig_all=sig_all[np.where(sig_all!=0)]*fact[i]
        #sig=diff_udata_all[i][:,30:70,30:60].flatten()
        ch=udata_all[:,55:61,48:54].mean()
        ar=udata_all[:,41:46,58:63].mean()
        qs=udata_all[:,40:45,42:47].mean()
        #plt.hist(sig_ch,bins=200,histtype='step',normed=1)
        #plt.hist(sig_ar,bins=200,histtype='step',normed=1)
        #plt.hist(sig_qs,bins=200,histtype='step',normed=1)
        plt.hist(sig_all,bins=300,histtype='step',normed=1,label='IMAGE')
        plt.hist(res,bins=300,histtype='step',normed=1,label='RESIDUAL')
        plt.title(str(freq[i])+' MHz')
        plt.legend()
        plt.xlabel('T$_B$ (K)')
        print freq[i], np.std(res), np.std(sig_all),np.std(sig_all)-np.std(res) 
        plt.close()


plot_variation=0
if(plot_variation):
    plt.plot(freq,flux_ar_rms,'o-',color='r',label='AR')
    plt.plot(freq,flux_qs_rms,'o-',color='b',label='QS')
    plt.plot(freq,flux_ch_rms,'o-',color='g',label='CH')
    plt.plot(freq,np.mean(np.array(flux_ar_res),axis=1),'o--',color='r',label='AR (RES)')
    plt.plot(freq,np.mean(np.array(flux_qs_res),axis=1),'o--',color='b',label='QS (RES)')
    plt.plot(freq,np.mean(np.array(flux_ch_res),axis=1),'o--',color='g',label='CH (RES)')
    plt.plot(freq,snoise_F,'o-',color='k',label='SELF NOISE')
    plt.plot(freq,SEFD,'o--',color='k',label='THERMAL NOISE')
    plt.legend(loc='upper center',fontsize=15)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.show()




######## FORWARD

fwd_list=sorted(glob.glob('/home/i4ds1807205/20151203/*psimas.sav'))[1:]
Tb_fwd=[0]*len(fwd_list)
Tb_fwd_con=[0]*len(fwd_list)
freq_fwd=[0]*len(fwd_list)
beam=[0]*len(fwd_list)
size=[0]*len(fwd_list)
flux_fwd=[0]*len(fwd_list)
flux_fwd_int=[0]*len(fwd_list)
for i in range(len(fwd_list)):
    freq_fwd[i]=int(fwd_list[i].split('_')[1].split('M')[0])
    fwd=readsav(fwd_list[i])
    Tb_fwd[i]=fwd['stokesstruct'][0][0]
    fwd_dx=fwd['quantmap'][0][3]*16*60 # in arcsec
    beam[i]=ut.makeGaussian(60, bmin[i]/fwd_dx, bmax[i]/fwd_dx , center=None)
    Tb_fwd_con[i] = signal.convolve(Tb_fwd[i],beam[i], mode='same')/np.sum(beam[i])
    size[i]=fwd_dx*Tb_fwd_con[i].shape[0]
    flux_fwd[i]=ut.Tb2flux(Tb_fwd[i], fwd_dx, fwd_dx, freq_fwd[i]/1000.)
    flux_fwd_int[i]=np.sum(flux_fwd[i])

Tb_fwd=np.array(Tb_fwd)
Tb_fwd_con=np.array(Tb_fwd_con)

######## MWA and FORWARD
mwa_fwd=0
if(mwa_fwd):
    mwa_coord=np.linspace(-2500,2500,100)
    forward_coord=np.linspace(-2880,2880,256)
    mwax,mway=np.meshgrid(mwa_coord, mwa_coord)
    fwdx,fwdy=np.meshgrid(forward_coord,forward_coord)
    Tb_fwd_interp=[0]*9
    for i in range(9):
        print i
        finterp = interpolate.interp2d(fwdx, fwdy, Tb_fwd_con[i], kind='linear')
        Tb_fwd_interp[i]=finterp(mwax,mway)
    Tb_fwd_interp=np.array(Tb_fwd_interp)
    pickle.dump(Tb_fwd_interp,open('Tb_fwd_interp.p','wb'))


######## PLOT
print 'Plotting ...'

plot_rms=0
if(plot_rms):
    plt.plot(freq,udata_ar_rms,'o-',color='red',label='AR')
    plt.plot(freq,udata_qs_rms,'o-',color='green',label='QS')
    plt.plot(freq,udata_ch_rms,'o-',color='blue',label='CH')
    plt.plot(freq,udata_ar_res_rms_t,'o--',color='red',label='AR (RES)')
    plt.plot(freq,udata_qs_res_rms_t,'o--',color='green',label='QS (RES)')
    plt.plot(freq,udata_ch_res_rms_t,'o--',color='blue',label='CH (RES)')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('RMS (Jy/beam)')
    plt.legend(loc=2)
    plt.show()

plot_rms_ratio=0
if(plot_rms_ratio):
    plt.plot(freq,udata_ar_rms/udata_ar_res_rms_st,'o-',color='red',label='AR')
    plt.plot(freq,udata_qs_rms/udata_qs_res_rms_st,'o-',color='green',label='QS')
    plt.plot(freq,udata_ch_rms/udata_ch_res_rms_st,'o-',color='blue',label='CH')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Ratio (CLEAN RMS/RESIDUAL RMS)')
    plt.title('MEAN of RMSs')
    plt.legend(loc=2)
    plt.show()

plot_image_rms=0
if(plot_image_rms):
    plt.plot(freq,udata_ar_res_mean,'o-',color='red',label='AR (RES)')
    plt.plot(freq,udata_qs_res_mean,'o-',color='green',label='QS (RES)')
    plt.plot(freq,udata_ch_res_mean,'o-',color='blue',label='CH (RES)')
    plt.plot(freq,sun_res_rms,'o--',color='k',label='solar (RES)')
    plt.plot(freq,sun_rms,'o-',color='k',label='solar (IMG)')
    plt.plot(freq,nonsun_rms,'o:',color='k',label='Non-solar (RES)')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('RESIDUAL RMS')
    plt.legend(loc=2)
    plt.show()

compare_clean_dirty=0
if(compare_clean_dirty):
    ratio_ar=[0]*9
    ratio_qs=[0]*9
    ratio_ch=[0]*9
    del_ar=[0]*9
    del_qs=[0]*9
    del_ch=[0]*9
    for i in range(9):
        ratio_ar[i]=np.std(Tb_ar[i])/np.std(Tb_ar_dirty[i])
        ratio_ch[i]=np.std(Tb_ch[i])/np.std(Tb_ch_dirty[i])
        ratio_qs[i]=np.std(Tb_qs[i])/np.std(Tb_qs_dirty[i])
        del_ar[i]=100*abs((np.std(Tb_ar_dirty[i])-np.std(Tb_ar[i]))/((np.std(Tb_ar_dirty[i])+np.std(Tb_ar[i]))/2))
        del_qs[i]=100*abs((np.std(Tb_qs_dirty[i])-np.std(Tb_qs[i]))/((np.std(Tb_qs_dirty[i])+np.std(Tb_qs[i]))/2))
        del_ch[i]=100*abs((np.std(Tb_ch_dirty[i])-np.std(Tb_ch[i]))/((np.std(Tb_ch_dirty[i])+np.std(Tb_ch[i]))/2))
        print freq[i],np.std(Tb_ar[i]),np.std(Tb_ar_dirty[i]),np.std(Tb_ar[i])/np.std(Tb_ar_dirty[i])
    #plt.plot(freq,ratio_ar,'o-',label='Active Region')
    #plt.plot(freq,ratio_qs,'o-',label='Quiet Sun')
    #plt.plot(freq,ratio_ch,'o-',label='Coronal Hole')
    #plt.axhline(y=1,c='k')
    plt.plot(freq,del_ar,'o-',label='Active Region')
    plt.plot(freq,del_qs,'o-',label='Quiet Sun')
    plt.plot(freq,del_ch,'o-',label='Coronal Hole')
    #plt.plot(Tb_ar[0],label='CLEAN')
    #plt.plot(Tb_ar_dirty[0],label='DIRTY')
    #plt.xlabel('Time (half-sec)')
    plt.xlabel('Frequency (MHz)')
    #plt.ylabel('1 $\sigma$ Noise Ratio (CLEAN/DIRTY)')
    plt.ylabel('1 $\sigma$ Noise percentage ($\%$)')
    plt.legend()
    plt.show()
        



movie_Tb=0
if(movie_Tb):
    for i in range(Tb_all.shape[1]):
        fig = plt.figure(figsize=(15,30))
        ax = fig.gca(projection='3d')
        X,Y=np.meshgrid(np.linspace(-2500,2500,100),np.linspace(-2500,2500,100))
        for j in range(Tb.shape[0]):
            cax=ax.contourf(X,Y,Tb_all[j,i]/1.e6,100,zdir='z', offset=freq[j],alpha=0.8,cmap='YlOrRd')
            ax.contour(X,Y,Tb_all[j,i]/np.max(Tb_all[j,i]),zdir='z',offset=freq[j],colors='k',lev=[0.4,0.6,0.8,0.9])
        ax.set_zlim(100, 250)
        ax.view_init(elev=-168, azim=53)
        ax.dist=9
        fig.colorbar(cax,label='T$_{B}$ (MK)',orientation='horizontal',norm=mpl.colors.Normalize(vmin=0.01, vmax=0.4),fraction=0.046, pad=0.01)
        ax.set_xlabel('$X (arcsec)$', fontsize=40, rotation=150)
        ax.set_ylabel('$Y (arcsec)$', fontsize=40, rotation=150)
        ax.set_zlabel('$Frequency (MHz)$', fontsize=40, rotation=0)
        ax.set_title('Time: '+str(i*0.5)+' sec')
        ax.tick_params(axis="y",direction="in", pad=-22)
        ax.tick_params(axis="x",direction="in", pad=-22)
        ax.tick_params(axis="z",direction="in", pad=-22)
        plt.savefig('plots_3d/Tb_'+str('%03d'%i)+'.png')
        plt.close()

compare_flux=0
if(compare_flux):
    plt.plot(freq,flux,'o-',color='red')
    plt.errorbar(freq,flux,yerr=flux_std,color='red',label='MWA MAPS')
    plt.plot(freq,flux_fwd_int,'o-',label='FORWARD MAPS')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend()
    plt.show()


fwd=0
if(fwd):
    fig = plt.figure(figsize=(15,30))
    ax = fig.gca(projection='3d')
    X,Y=np.meshgrid(np.linspace(-1440,1440,256),np.linspace(-1440,1440,256))
    for i in range(Tb_fwd.shape[0]):
        cax=ax.contourf(X,Y,Tb_fwd_con[i]/1.e6,100,zdir='z', offset=freq_fwd[i],alpha=0.8,vmin=1.e5/1.e6,vmax=2.e6/1.e6,cmap='YlOrRd')
        ax.contour(X,Y,Tb_fwd_con[i]/Tb_fwd[i].max(),zdir='z',offset=freq_fwd[i],colors='k',lev=[0.4,0.6,0.8,0.9])
    ax.set_zlim(100, 250)
    ax.view_init(elev=-168, azim=53)
    ax.dist=9
    fig.colorbar(cax,label='T$_{B}$ (MK)',orientation='horizontal',fraction=0.046, pad=0.04)
    ax.set_xlabel('$X (arcsec)$', fontsize=40, rotation=150)
    ax.set_ylabel('$Y (arcsec)$', fontsize=40, rotation=150)
    ax.set_zlabel('$Frequency (MHz)$', fontsize=40, rotation=0)
    ax.set_xlim([-2500,2500])
    ax.set_ylim([-2500,2500])
    ax.tick_params(axis="y",direction="in", pad=-22)
    ax.tick_params(axis="x",direction="in", pad=-22)
    ax.tick_params(axis="z",direction="in", pad=-22)
    plt.show()


cont_rms=0
if(cont_rms):
    fig = plt.figure(figsize=(15,30))
    ax = fig.gca(projection='3d')
    X,Y=np.meshgrid(np.linspace(-2500,2500,100),np.linspace(-2500,2500,100))
    for i in range(Tb.shape[0]):
        cax=ax.contourf(X,Y,Tb_rms[i],100,zdir='z', offset=freq[i],alpha=0.8,vmin=100,vmax=4000,cmap='jet')
        #ax.contour(X,Y,Tb_rms[i],zdir='z',offset=freq[i],colors='k',lev=np.linspace(0,1,100)*1.e3)
    ax.set_zlim(100, 250)
    ax.view_init(elev=-168, azim=53)
    ax.dist=9
    fig.colorbar(cax,label='T$_{B}$ (kK)',orientation='vertical',fraction=0.046, pad=0.04)
    ax.set_xlabel('$X (arcsec)$', fontsize=40, rotation=150)
    ax.set_ylabel('$Y (arcsec)$', fontsize=40, rotation=150)
    ax.set_zlabel('$Frequency (MHz)$', fontsize=40, rotation=0)
    ax.tick_params(axis="y",direction="in", pad=-22)
    ax.tick_params(axis="x",direction="in", pad=-22)
    ax.tick_params(axis="z",direction="in", pad=-22)
    plt.show()

cont_Tb=0
if(cont_Tb):
    fig = plt.figure(figsize=(15,30))
    ax = fig.gca(projection='3d')
    X,Y=np.meshgrid(np.linspace(-2500,2500,100),np.linspace(-2500,2500,100))
    for i in range(Tb.shape[0]):
        if(i!=1):
            cax=ax.contourf(X,Y,Tb[i]/1.e6,100,zdir='z', offset=freq[i],alpha=0.8,vmin=1.e4/1.e6,vmax=5.e5/1.e6,cmap='YlOrRd')
            ax.contour(X,Y,Tb[i]/Tbmax[i],zdir='z',offset=freq[i],colors='k',lev=[0.4,0.6,0.8,0.9])
    ax.set_zlim(100, 250)
    ax.view_init(elev=-168, azim=53)
    ax.dist=9
    fig.colorbar(cax,label='T$_{B}$ (MK)',orientation='horizontal',fraction=0.046, pad=0.04)
    ax.set_xlabel('$X (arcsec)$', fontsize=40, rotation=150)
    ax.set_ylabel('$Y (arcsec)$', fontsize=40, rotation=150)
    ax.set_zlabel('$Frequency (MHz)$', fontsize=40, rotation=0)
    ax.tick_params(axis="y",direction="in", pad=-22)
    ax.tick_params(axis="x",direction="in", pad=-22)
    ax.tick_params(axis="z",direction="in", pad=-22)
    plt.show()
ml=0
if(ml):
    z,x,y=np.mgrid[0:8:9j,-2500:2500:100j,-2500:2500:100j]
    iso=mlab.points3d(z,x,y,Tb,vmin=Tb.min(),vmax=Tb.max(),opacity=0.2)
    iso.contour.number_of_contours = 15
    mlab.show()


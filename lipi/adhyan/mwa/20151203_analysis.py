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
from surya.radio import main as rd
import matplotlib.cm as cm
from surya.gm import main as gm
import pywt


def fit_cubic(y):
    x=np.arange(y.shape[0])
    ynew=np.zeros((y.shape[0],y.shape[1],y.shape[2]))
    snew=np.zeros((y.shape[0],y.shape[1],y.shape[2]))
    for i in range(y.shape[1]):
        for j in range(y.shape[2]):
            f = np.polyfit(x, y[:,i,j], 5)
            p = np.poly1d(f)
            snew[:,i,j]=p(x)
            ynew[:,i,j]=y[:,i,j]-p(x)
    return ynew,snew

def robust_median(data_,n):
	mdata=-10*np.ones((data_.shape[0],data_.shape[1]-n))
	ndata=-10*np.ones((data_.shape[0],data_.shape[1]-n))
	for ch in range(data_.shape[0]):
		print 'Channel: ',ch
		for i in range(data_.shape[1]-n):
			data=data_[ch][i:n+i]
			data[np.isnan(data)] = -10
			rms=np.std(data)
			indx=np.arange(data.shape[0])
			x=(data<np.median(data)+1.2*np.std(data)) & (data>np.median(data)-1.2*np.std(data))
			datarms=data[x]
			mdata[ch,i]=np.median(data[x])
		mdata[ch][np.where(mdata[ch]==0)]=100
		ndata[ch]=data_[ch][int(n/2):-1*int(n/2)]-mdata[ch]
	return ndata,mdata

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')
#freq=[108.0,120.0,132.0,145.0,161.0,179.0,196.0,217.0,240.0]
freq=[108.0,132.0,145.0,161.0,179.0,196.0,217.0,240.0]
#freq=[108.0,145.0,161.0,179.0,196.0,217.0,240.0]
#fband=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
fband=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
#fband=['084-085','113-114','125-126','139-140','153-154','169-170','187-188']
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
Ssun_mean=[0]*len(freq)
Ssun_std=[0]*len(freq)
calS_mean=[0]*len(freq)
calS_std=[0]*len(freq)
std_ni=[0]*len(freq)
SEFD_base=[0]*len(freq)
NEFD_base=[0]*len(freq)
snoise_base=[0]*len(freq)
th_noise_base=[0]*len(freq)
th_noise_array=[0]*len(freq)
Tb_beam=[0]*len(freq)
Tsys=[0]*len(freq)
sa=[0]*len(freq)
Aeff=[0]*len(freq)
polTb=[0]*len(freq)
Ssun_mean_diff=[0]*len(freq)
Ssun_mean_cont=[0]*len(freq)
print 'Reading files...'

c=3.e8 # S.I. Unit
kb=1.38e-23 # S.I. Units
Tsky_=np.array([665,483,371,290,225,184,179,181,122,94])
Trec_=np.array([30,28,26,24,21,20,21,23,27,32])
Tgrd_=np.array([20,17,15,13,12,12,13,18,10,9])
f=np.array([103,117,131,148,167,189,213,240,272,299])
finters=interpolate.interp1d(f, Tsky_,kind='cubic')
finterr=interpolate.interp1d(f, Trec_,kind='cubic')
finterg=interpolate.interp1d(f, Tgrd_,kind='cubic')
Tsky=finters(freq)
Trec=finterr(freq)
Tgrd=finterg(freq)
baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
flux_path='/media/rohit/MWA/20151203_MWA_NEW/new_pickle/20151203_'
cal_path='/media/rohit/VLA/20151203_cal/20151203_cal/pickle/flux_V1_1133139776-%b'
gm_path='/media/rohit/VLA/20151203_MWA/pickle/gm/'
img_path='/media/rohit/VLA/20151203_MWA/images_all/1133149192-%b'
outdir='/media/rohit/VLA/20151203_MWA/Tb_new/'
centre_file='/home/i4ds1807205/20151203/20151203_nasa_horizon.dat'
ids='1133149192'

def get_solidangle(S,Tbeam,f):
    '''
    Input:  Tbeam (Beam averaged solar temp in K), S (flux in SFU), f (in MHz)
    Output: Solid Angle (in radians)
    '''
    c=3.e8
    kb=1.38e-23
    sa=(S*1.e-22)*(c/(f*1.e6))**2/(2*kb*Tbeam)
    return np.nanmean(sa)


def get_flux(flux_path,baseline_filelist,ff,i,sol):
    '''
    Output:
    mean flux in frequency and baselines, std flux in frequency and time
    '''
    Ssun_all,Tb_beam,time_str,timesec=tb.mean_flux(flux_path,fband[ff],baseline_filelist,0.5,sol)
    if(i):
        Ssun_all[:,2,:]=Ssun_all[:,1,:]
        Ssun_all[:,8,:]=Ssun_all[:,7,:]
        Tb_beam[:,2,:]=Tb_beam[:,1,:]
        Tb_beam[:,8,:]=Tb_beam[:,7,:]
    Ssun_all_=np.concatenate((Ssun_all[:,6:13],Ssun_all[:,19:26],Ssun_all[:,38:45],Ssun_all[:,51:58]),axis=1)
    Ssun_mean=np.nanmean(Ssun_all_,axis=(1))[2]
    Tb_beam=np.nanmean(Tb_beam,axis=(1,2))[2]
    Ssun_std=np.nanstd(Ssun_all,axis=(1))[2]
    return Ssun_mean,Ssun_std,Tb_beam
    
############################ Pre- Requisite ############################

get_Tb=1
if(get_Tb==1):
    print 'Starting Tb computation..'
    del_=50
    angle=-15.4 # Negative angle imples clockwise rotation and vice versa
    res=50 # arcsec
    for ff in range(len(fband)):
    #for ff in range(1):
        print str(freq[ff])+' MHz'
        ### GET FLUX ###
        Ssun_mean[ff],Ssun_std[ff],Tb_beam[ff]=get_flux(flux_path,baseline_filelist,ff,1,1)
        ### GET CAL FLUX ###
        calS_mean[ff],calS_std[ff],calT=get_flux(cal_path,baseline_filelist,ff,0,0)
        ### GM Parameters ##
        gm_filepath=gm_path+str(int(freq[ff]))+'MHz/T'
        std_ni[ff]=[0]*len(baseline_filelist)
        bb=0
        for b in baseline_filelist:
            gm_data=pickle.load(open(gm_filepath+str(b)+'_f'+fband[ff]+'_P0_C0-2_paramnorm1','r'))
            std_ni[ff][bb]=float(gm_data[0][1])
            bb=bb+1
        ### GET THERMAL NOISE ###
        sa[ff]=get_solidangle(Ssun_mean[ff],Tb_beam[ff],freq[ff])
        Aeff[ff]=(c/(freq[ff]*1.e6))**2 /sa[ff]
        Tsys[ff]=Tsky[ff]+Trec[ff]+Tgrd[ff]+Tb_beam[ff]
        SEFD_base[ff]=rd.get_SEFD(Tsys[ff],Aeff[ff])
        th_noise_base[ff]=rd.thermal_noise(SEFD_base[ff],2,4.e4,0.5,2,1)
        th_noise_array[ff]=rd.thermal_noise(SEFD_base[ff],128,2.e6,0.5,2,1)
        ### IMAGE ANALYSIS ##
        imglist=sorted(glob.glob(img_path+fband[ff]+'*.image.FITS'))[0:570]
        reslist=sorted(glob.glob(img_path+fband[ff]+'*.residual.FITS'))[0:570]
        img_time=[0]*len(imglist)
        polTb[ff]=[0]*len(imglist)
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
        flux=[0]*len(imglist)
        bmaj=[0]*len(imglist)
        bmin=[0]*len(imglist)
        bpa=[0]*len(imglist)
        ndata=[0]*len(imglist)
        offsun_mean=[0]*len(imglist)
        offsun_std=[0]*len(imglist)
        centre=[0]*len(imglist)
        j=0
        for fitsfile in imglist[:1]:
            xc,yc,img_time[j]=tb.solar_center_pixel(fitsfile,centre_file)
            centre[j]=[xc,yc]
            Tb[j],flux[j],offsun_mean[j],offsun_std[j],bmaj[j],bmin[j],bpa[j],ndata[j]=tb.compute_Tb(fitsfile,xc,yc,del_,angle,res,freq[ff]*1.e6,0.001,Ssun_mean[ff][j])
            resi_sun[j],resi_=tb.get_residuals(fitsfile,xc,yc,del_,angle)
            resi[j]=resi_[0:200,0:200]
            polTb[ff]=ut.cart2polar(Tb[j])[0]
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
        flux=np.array(flux)
        ndata=np.array(ndata)
        # do full imglist 
        print 'Mean Tb: ',np.max(Tb[0]),np.std(Tb[0]),np.max(ndata[0]),bmaj[0],bmin[0],len(np.where(ndata[0]!=0)[0])
        #pickle.dump([Tb,offsun_mean,offsun_std,bmaj,bmin,bpa,centre,flux,ndata,polTb],open(outdir+'Tb_'+str(ids)+'-'+str(fband[ff])+'.p','w'))
        #pickle.dump([resi,resi_sun,Tb_sun_resi,flux_sun_resi,Tb_fac,flux_fac],open(outdir+'res_'+str(ids)+'-'+str(fband[ff])+'.p','w'))
    polTb=np.array(polTb)
    ##### Non-Imaging Analysis #################
    N=10
    sigrand=np.random.normal(0,1,570)
    sigrand_cal=np.random.normal(0,1,210)
    nwin=[2,10,50,100]
    sigrand=[0]*4
    sigrand_smooth=[0]*4
    win=[0]*4
    nn=0
    for N in nwin:
        win[nn]=np.arange(N)/N
        sigrand_smooth[nn]=np.correlate(sigrand,win[nn]/np.sum(win[nn]),'valid')
        nn=nn+1
    sigrand_smooth=np.array(sigrand_smooth)
    win=np.array(win)
    ######
    N=2
    M=3944
    Mcal=210
    sigrand=np.random.normal(0,1,M)
    sigrand_cal=np.random.normal(0,1,Mcal)
    S_ts=[0]*len(fband)
    diff_S=[0]*len(fband)
    diff_S_smooth=[0]*len(fband)
    auto_S=[0]*len(fband)
    diff_S_std=[0]*len(fband)
    S_ts_cal=[0]*len(fband)
    diff_S_cal=[0]*len(fband)
    diff_S_smooth_cal=[0]*len(fband)
    auto_S_cal=[0]*len(fband)
    diff_S_std_cal=[0]*len(fband)
    auto_rand=[0]*len(fband)
    diff_rand_smooth=[0]*len(fband)
    kernel=(np.arange(N)/N)/np.sum(np.arange(N)/N)
    for i in range(len(fband)):
        S_ts[i]=Ssun_mean[i][0:M]#570]
        S_ts[i][np.where(S_ts[i]<0)]=0
        #S_ts[i][np.isnan(S_ts[i])]=0
        MM=S_ts[i].shape[0]
        #diff_S_=fit_cubic(S_ts[i].reshape(M,1,1))#(570,1,1))
        diff_S_=robust_median(S_ts[i].reshape(1,M,),200)#(570,1,1))
        diff_S_[0][np.where(diff_S_[0]<-1)]=0
        diff_S[i]=diff_S_[0][0]
        #diff_S_smooth[i]=np.correlate(diff_S[i],kernel,'valid')
        diff_S_std[i]=np.nanstd(diff_S[i])
        #auto_S[i]=np.correlate(diff_S_smooth[i]/diff_S_std[i],diff_S_smooth[i]/diff_S_std[i],'same')
        #auto_S[i]=np.correlate(diff_S[i]/diff_S_std[i],diff_S[i]/diff_S_std[i],'same')
        auto_S[i]=np.correlate(diff_S[i],diff_S[i],'same')
        # CAL
        #S_ts_cal[i]=calS_mean[i][0:Mcal]
        #diff_S_cal_=fit_cubic(S_ts_cal[i].reshape(Mcal,1,1))
        #diff_S_cal[i]=diff_S_cal_[0].flatten()
        #diff_S_smooth_cal[i]=np.correlate(diff_S_cal[i],kernel,'valid')
        #diff_S_std_cal[i]=np.std(diff_S_cal[i])
        #auto_S_cal[i]=np.correlate(diff_S_cal[i]/diff_S_std_cal[i],diff_S_cal[i]/diff_S_std_cal[i],'same')
        # RANDOM
        #diff_rand_smooth=np.correlate(sigrand,kernel,'valid')
        #auto_rand=np.correlate(sigrand/np.std(sigrand),sigrand/np.std(sigrand),'same')
        auto_rand=np.correlate(sigrand*th_noise_base[i]/np.std(sigrand),sigrand*th_noise_base[i]/np.std(sigrand),'same')
        #diff_rand_smooth_cal=np.correlate(sigrand_cal,kernel,'valid')
        #auto_rand_cal=np.correlate(diff_rand_smooth_cal/np.std(diff_rand_smooth_cal),diff_rand_smooth_cal/np.std(diff_rand_smooth_cal),'same')
        #auto_rand_cal=np.correlate(sigrand_cal/np.std(sigrand_cal),sigrand_cal/np.std(sigrand_cal),'same')

    lag=0.5*(np.arange(len(diff_S[0]))-int(diff_S[0].shape[0]/2))
    lagrand=0.5*(np.arange(len(sigrand))-int(sigrand.shape[0]/2))
    #lagcal=0.5*(np.arange(len(diff_S_cal[0]))-int(diff_S_cal[0].shape[0]/2))
    plot_auto=1
    if(plot_auto):
        #plt.plot(lagcal,auto_S_cal[-2],'o-',label='CAL')
        plt.plot(lagrand,auto_rand,'o',label='RANDOM')
        plt.plot(lag,auto_S[-2],'o-',label='218 MHz')
        plt.xlabel('Lag Time (sec)')
        plt.ylabel('Auto-correlation')
        plt.legend()
        plt.show()
    do_cross=0
    if(do_cross):
        cross_S=[0]*28
        fpair=[0]*28
        cross_max=[0]*28
        cross_min=[0]*28
        cross_diff=[0]*28
        m=0
        for k in range(len(fband)):
            for l in range(k):
                fpair[m]=np.array((freq[k],freq[l]))
                cross_S[m]=np.correlate(diff_S_smooth[k]/diff_S_std[k],diff_S_smooth[l]/diff_S_std[l],'same')
                cross_max[m]=np.where(cross_S[m]==np.max(cross_S[m]))[0][0]
                cross_min[m]=np.where(cross_S[m]==np.min(cross_S[m]))[0][0]
                cross_diff[m]=cross_min[m]-cross_max[m]
                m=m+1
    ## 
    plot_thermal=0
    if(plot_thermal):
        plt.plot(freq,diff_S_std,'o-',label='SOL')
        plt.plot(freq,diff_S_std_cal,'o-',label='CAL')
        plt.plot(freq,th_noise_base,'o-',label='Thermal Noise')
        plt.ylabel('Flux (SFU)')
        plt.xlabel('Frequency (MHz)')
        plt.legend(loc=2)
        plt.show()


    ## wavelets
    #widths=np.arange(1,600)
    widths=np.logspace(1.0,3.0,num=400)*0.5
    #mother='gaus1'
    mother='morl'
    cwtmatr=[0]*len(fband)
    freqs=[0]*len(fband)
    power=[0]*len(fband)
    period=[0]*len(fband)
    cwtmatr_cal=[0]*len(fband)
    freqs_cal=[0]*len(fband)
    power_cal=[0]*len(fband)
    period_cal=[0]*len(fband)
    cwtmatr_rand=[0]*len(fband)
    freqs_rand=[0]*len(fband)
    power_rand=[0]*len(fband)
    period_rand=[0]*len(fband)
    slope=[0]*len(fband)
    slope_cal=[0]*len(fband)
    for i in range(len(fband)):
        wdata=diff_S[i]/np.nanstd(diff_S[i])
        cwtmatr[i], freqs[i] = pywt.cwt(wdata[0], widths, mother)
        #cwtmatr[i]=cwtmatr[i]*diff_S_std[i]
        power[i] = (np.abs(cwtmatr[i])) ** 2
        period[i] = 1 / freqs[i]
        # random
        cwtmatr_rand[i], freqs_rand[i] = pywt.cwt(sigrand, widths, mother)
        #cwtmatr_rand[i]=cwtmatr_rand[i]*np.std(sigrand)
        power_rand[i] = (np.abs(cwtmatr_rand[i])) ** 2
        period_rand[i] = 1 / freqs_rand[i]
        # 5:50 for non variance case, 10:100 with variance
        slope[i]= ut.fit_1d(freqs[i][10:100],power[i].mean(axis=1)[10:100])
        #wdata_cal=diff_S_smooth_cal[i]/diff_S_std_cal[i]
        wdata_cal=diff_S_cal[i]/np.std(diff_S_cal[i])
        cwtmatr_cal[i], freqs_cal[i] = pywt.cwt(wdata_cal, widths, mother)
        #cwtmatr_cal[i]=cwtmatr_cal[i]*diff_S_std_cal[i]
        power_cal[i] = (np.abs(cwtmatr_cal[i])) ** 2
        period_cal[i] = 1 / freqs_cal[i]
        slope_cal[i]= ut.fit_1d(freqs_cal[i][10:100],power_cal[i].mean(axis=1)[10:100])
    slope=np.array(slope)
    slope_cal=np.array(slope_cal)
    power=np.array(power)
    power_cal=np.array(power_cal)
    power_rand=np.array(power_rand)
    ## cross-correlation of reconstructed waves
    t=20
    cross_wave=[0]*28
    fpair_wave=[0]*28
    m=0
    for k in range(len(fband)):
        for l in range(k):
            fpair_wave[m]=np.array((freq[k],freq[l]))
            cross_wave[m]=np.correlate(power[k][t]/np.std(power[k][t]),power[l][t]/np.std(power[l][t]),'same')
            m=m+1


    wavelet_power_mean=[0]*len(fband)
    wavelet_power_sum=[0]*len(fband)
    wavelet_power_mean_cal=[0]*len(fband)
    wavelet_power_sum_cal=[0]*len(fband)
    wavelet_power_mean_rand=[0]*len(fband)
    wavelet_max=[0]*len(fband)
    for i in range(len(fband)):
        wavelet_power_mean[i]=power[i].mean(axis=1)
        wavelet_power_sum[i]=power[i].sum(axis=1)
        wavelet_power_mean_cal[i]=power_cal[i].mean(axis=1)
        wavelet_power_sum_cal[i]=power_cal[i].sum(axis=1)
        wavelet_power_mean_rand[i]=power_rand[i].mean(axis=1)
        ## plot
        plt.plot(widths,wavelet_power_mean[i],'o-',label=str(freq[i])+' MHz')
        #plt.plot(widths,wavelet_power_mean_cal[i],'o-')
        #plt.plot(widths,wavelet_power_mean[i]/wavelet_power_mean_cal[i],'o-')
        wavelet_max[i]=widths[np.where(wavelet_power_mean[i]==np.max(wavelet_power_mean[i]))]
    #plt.plot(1/widths,np.array(wavelet_power_mean_cal).mean(axis=0)*170,'o-',label='CAL (scaled by 170)')
    slope_mean=ut.fit_1d(1/widths[20:50],np.array(wavelet_power_mean).mean(axis=0)[20:50])
    slope_mean_cal=ut.fit_1d(1/widths[20:50],np.array(wavelet_power_mean_cal).mean(axis=0)[20:50])
    plt.plot(1/widths,np.array(wavelet_power_mean_cal).mean(axis=0),'o-',color='blue',label='CAL')
    plt.plot(1/widths,np.array(wavelet_power_mean_rand).mean(axis=0),'o-',color='black',label='RANDOM')
    plt.plot(1/widths,np.array(wavelet_power_mean).mean(axis=0),'o-',color='red',label='SUN')
    plt.plot(1/widths,10**(np.log10(1/widths)*slope_mean[0]+slope_mean[1]),color='red',label='Sun (fit): '+str(np.round(slope_mean[0],2))+'$\pm$'+str(np.round(slope_mean[2],2)))
    plt.plot(1/widths,10**(np.log10(1/widths)*slope_mean_cal[0]+slope_mean_cal[1]),color='blue',label='CAL (fit): '+str(np.round(slope_mean_cal[0],2))+'$\pm$'+str(np.round(slope_mean[2],2)))
    plt.plot()
    plt.xlabel('Scale (sec$^-1$)')
    #plt.ylabel('Power (SFU$^2$)')
    plt.ylabel('Power')
    plt.legend()
    #plt.show()
    plt.close()

        ## Plot wavelets
    for i in range(len(fband)):
        title=str(freq[i])+' MHz'
        label='Flux'
        units='SFU'
        t=np.arange(len(diff_S[0]))*0.5
        glbl_power = power[i].mean(axis=1)
        scale_avg=power[i].mean(axis=0)
        plt.ioff()
        figprops = dict(figsize=(20, 20), dpi=72)
        fig = plt.figure(**figprops)
        # First sub-plot, the original time series anomaly and inverse wavelet
        # transform.
        ax = plt.axes([0.1, 0.75, 0.65, 0.2])
        ax.plot(t, diff_S[i], '-', linewidth=1, color=[0.5, 0.5, 0.5])
        ax.set_title('a) {}'.format(title))
        ax.grid(True)
        ax.set_ylabel(r'{} [{}]'.format(label, units))
        # Second sub-plot, the normalized wavelet power spectrum and significance
        # level contour lines and cone of influece hatched area. Note that period
        # scale is logarithmic.
        bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax)
        levels = [0.015625,0.03125,0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
        bx.contourf(t, np.log2(period[i]), np.log10(power[i]), np.log2(levels),
                            extend='both', cmap=plt.cm.viridis)
        extent = [t.min(), t.max(), 0, max(period[i])]
        #bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt,t[:1] - dt, t[:1] - dt]),np.concatenate([np.log2(coi), [1e-9], np.log2(period[i][-1:]),np.log2(period[i][-1:]), [1e-9]]),'k', alpha=0.3, hatch='x')
        bx.set_title('b) Wavelet Power Spectrum ({})'.format(label, 'MORLET'))
        bx.set_ylabel('log2(Period) (sec)')
        bx.grid(True)
        #
        #Yticks = 2 ** np.arange(np.ceil(np.log2(period[i].min())),
        #                                   np.ceil(np.log2(period[i].max())))
        #bx.set_yticks(np.log2(Yticks))
        #bx.set_yticklabels(Yticks)

        # Third sub-plot, the global wavelet and Fourier power spectra and theoretical
        # noise spectra. Note that period scale is logarithmic.
        cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
        cx.plot(glbl_power, np.log2(period[i]), 'k-', linewidth=1.5)
        cx.set_title('c) Global Wavelet Spectrum')
        cx.set_xlabel(r'Power [({})^2]'.format(units))
        #cx.set_xlim([4.e-4, glbl_power.max() + diff_S_std[i]])
        cx.set_ylim(np.log2([period[i].min(), period[i].max()]))
        cx.grid(True)
        cx.set_xscale('log')
        #cx.set_yticks(np.log2(Yticks))
        #cx.set_yticklabels(Yticks)
        #plt.setp(cx.get_yticklabels(), visible=False)

        # Fourth sub-plot, the scale averaged wavelet spectrum.
        dx = plt.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
        dx.plot(t, scale_avg, 'k-', linewidth=1.5)
        dx.set_title('d) Scale-averaged power')
        dx.set_xlabel('Time (sec)')
        dx.set_ylabel(r'Reconstructed Flux [{}]'.format(units))
        dx.grid(True)
        ax.set_xlim([t.min(), t.max()])
        #plt.savefig('wavelet_f'+str(int(freq[i]))+'.png')
        plt.close()

        ## Plotting..
        plot_ts=0
        if(plot_ts):
            for k in range(4,8):
                plt.plot(0.5*np.arange(len(diff_S_smooth[k])),diff_S_smooth[k],'-',label=str(freq[k])+' MHz')
            plt.legend()
            plt.ylabel('Flux (SFU)')
            plt.xlabel('Time (in sec)')
            plt.show()
        plot_auto=0
        if(plot_auto):
            for k in range(4,8):
                plt.plot(lag,auto_S[k],'-',label=str(freq[k])+' MHz')
            plt.legend()
            plt.ylabel('Auto-correlated Power')
            plt.xlabel('Lag (in sec)')
            plt.show()

        plot_cross=0
        if(plot_cross):
            for k in range(28):
                plt.plot(lag,cross_S[k],'o-',label=str(fpair[k])+' MHz')
                plt.legend()
                plt.ylabel('Auto-correlated Power')
                plt.xlabel('Lag (in sec)')
                plt.savefig('cross_corr_f'+str(int(fpair[k][0]))+'-'+str(int(fpair[k][1]))+'.png')
                plt.close()

azimuth_profile=0
if(azimuth_profile):
    colors = cm.rainbow(np.linspace(0, 1, 8))
    for i in range(len(fband)):
        y=np.mean(polTb[i],axis=1)/1.e6
        x=(50/60.)*np.arange(y.shape[0])
        plt.plot(x,y,'o-',color=colors[i],label=str(freq[i])+' MHz')
    plt.legend()
    plt.xlabel('Radial Distance (arcmin)')
    plt.ylabel('Brightness Temperature (MK)')
    plt.show()

sys.exit()
################################## MAIN ################################################

Tb=[0]*len(fband)
flux=[0]*len(fband)
bmaj=[0]*len(fband)
bmin=[0]*len(fband)
bpa=[0]*len(fband)
ndata=[0]*len(fband)
polTb=[0]*len(fband)
offsun_mean=[0]*len(fband)
offsun_std=[0]*len(fband)
centre=[0]*len(fband)
resi=[0]*len(fband)
resi_sun=[0]*len(fband)
Tb_sun_resi=[0]*len(fband)
Tb_fac=[0]*len(fband)
flux_sun_resi=[0]*len(fband)
flux_fac=[0]*len(fband)

for ff in range(len(fband)):
    print 'Reading..'+fband[ff]+' MHz'
    Tb[ff],offsun_mean[ff],offsun_std[ff],bmaj[ff],bmin[ff],bpa[ff],centre[ff],flux[ff],ndata[ff],polTb[ff]=pickle.load(open(outdir+'Tb_'+str(ids)+'-'+str(fband[ff])+'.p','r'))
    resi[ff],resi_sun[ff],Tb_sun_resi[ff],flux_sun_resi[ff],Tb_fac[ff],flux_fac[ff]=pickle.load(open(outdir+'res_'+str(ids)+'-'+str(fband[ff])+'.p','r'))

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

#Ssun_mean_time=np.array(Ssun_mean[0:580]).mean(axis=1)
#Ssun_rms_time=np.array(Ssun_mean[0:580]).std(axis=1)
#Ssun_rms_baseline=np.array(Ssun_std[0:580]).mean(axis=1)

######## FORWARD

fwd_list=sorted(glob.glob('/home/i4ds1807205/20151203/*psimas.sav'))[1:]
del fwd_list[1]
Tb_fwd=[0]*len(fwd_list)
Tb_fwd_con=[0]*len(fwd_list)
freq_fwd=[0]*len(fwd_list)
beam_fwd=[0]*len(fwd_list)
size=[0]*len(fwd_list)
flux_fwd=[0]*len(fwd_list)
flux_fwd_int=[0]*len(fwd_list)
for i in range(len(fwd_list)):
    freq_fwd[i]=int(fwd_list[i].split('_')[1].split('M')[0])
    fwd=readsav(fwd_list[i])
    Tb_fwd[i]=fwd['stokesstruct'][0][0]
    fwd_dx=fwd['quantmap'][0][3]*16*60 # in arcsec
    beam_fwd[i]=ut.make2DGaussian(60, bmin[i]/fwd_dx, bmax[i]/fwd_dx , center=None)
    Tb_fwd_con[i] = signal.convolve(Tb_fwd[i],beam_fwd[i], mode='same')/np.sum(beam_fwd[i])
    size[i]=fwd_dx*Tb_fwd_con[i].shape[0]
    flux_fwd[i]=ut.Tb2flux(Tb_fwd[i], fwd_dx, fwd_dx, freq_fwd[i]/1000.)
    flux_fwd_int[i]=np.sum(flux_fwd[i])

Tb_fwd=np.array(Tb_fwd)
Tb_fwd_con=np.array(Tb_fwd_con)

sys.exit()
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
########## Maps analysis ################

N=2
diff_Tb_smooth=[0]*100
auto_Tb=[0]*100
diff_Tb_std=[0]*100
for l in range(100):
    diff_Tb_smooth[l]=[0]*100
    auto_Tb[l]=[0]*100
    diff_Tb_std[l]=[0]*100
    for m in range(100):
        diff_Tb_smooth[l][m]=[0]*len(fband)
        auto_Tb[l][m]=[0]*len(fband)
        diff_Tb_std[l][m]=[0]*len(fband)
        for i in range(len(fband)):
            #diff_Tb_smooth[l][m][i]=np.correlate(diff_Tb_all[i,:,l,m],np.arange(N)/N,'same')
            diff_Tb_smooth[l][m][i]=diff_Tb_all[i,:,l,m]
            diff_Tb_std[l][m][i]=np.std(diff_Tb_smooth[l][m][i])
            auto_Tb[l][m][i]=np.correlate(diff_Tb_smooth[l][m][i]/diff_Tb_std[l][m][i],diff_Tb_smooth[l][m][i]/diff_Tb_std[l][m][i],'same')
diff_Tb_smooth=np.array(diff_Tb_smooth)
diff_Tb_std=np.array(diff_Tb_std)
auto_Tb=np.array(auto_Tb)
lagTb=0.5*(np.arange(len(diff_Tb_smooth[0,0,0]))-int(diff_Tb_smooth[0,0,0].shape[0]/2))

cross_Tb=[0]*100
for ll in range(100):
    cross_Tb[ll]=[0]*100
    for mm in range(100):
        cross_Tb[ll][mm]=[0]*28
        fpair=[0]*28
        cross_max=[0]*28
        cross_min=[0]*28
        cross_diff=[0]*28
        m=0
        for k in range(len(fband)):
            for l in range(k):
                fpair[m]=np.array((freq[k],freq[l]))
                cross_Tb[ll][mm][m]=np.correlate(diff_Tb_smooth[ll,mm,k]/diff_Tb_std[ll,mm,k],diff_Tb_smooth[ll,mm,l]/diff_Tb_std[ll,mm,l],'same')
                #cross_max[m]=np.where(cross_S[m]==np.max(cross_S[m]))[0][0]
                #cross_min[m]=np.where(cross_S[m]==np.min(cross_S[m]))[0][0]
                #cross_diff[m]=cross_min[m]-cross_max[m]
                m=m+1
cross_Tb=np.array(cross_Tb)


## wavelets
#mother='morl'
mother='gaus1'
sarrx=[[43,43,43,42,42,42,44,44,44],[58,58,58,57,57,57,59,59,59],[43,43,43,42,42,42,44,44,44]] # AR, CH, QS
sarry=[[59,60,61,59,60,61,59,60,61],[50,51,52,50,51,52,50,51,52],[44,45,46,44,45,46,44,45,46]]
cwtmatr_cal=[0]*3
cwtmatr=[0]*3
freqs=[0]*3
power=[0]*3
period=[0]*3
cwtmatr_res=[0]*3
freqs_res=[0]*3
power_res=[0]*3
period_res=[0]*3
power_flux_res=[0]*3
power_flux=[0]*3
slope_Tb=[0]*3

for ll in range(3):
    print ll
    cwtmatr[ll]=[0]*len(sarrx[ll])
    freqs[ll]=[0]*len(sarrx[ll])
    power[ll]=[0]*len(sarrx[ll])
    power_flux[ll]=[0]*len(sarrx[ll])
    period[ll]=[0]*len(sarrx[ll])
    cwtmatr_res[ll]=[0]*len(sarrx[ll])
    freqs_res[ll]=[0]*len(sarrx[ll])
    power_res[ll]=[0]*len(sarrx[ll])
    power_flux_res[ll]=[0]*len(sarrx[ll])
    period_res[ll]=[0]*len(sarrx[ll])
    slope_Tb[ll]=[0]*len(sarrx[ll])
    for mm in range(len(sarrx[ll])):
        cwtmatr[ll][mm]=[0]*len(fband)
        freqs[ll][mm]=[0]*len(fband)
        power[ll][mm]=[0]*len(fband)
        period[ll][mm]=[0]*len(fband)
        power_flux[ll][mm]=[0]*len(fband)
        cwtmatr_res[ll][mm]=[0]*len(fband)
        freqs_res[ll][mm]=[0]*len(fband)
        power_res[ll][mm]=[0]*len(fband)
        period_res[ll][mm]=[0]*len(fband)
        power_flux_res[ll][mm]=[0]*len(fband)
        slope_Tb[ll][mm]=[0]*len(fband)
        for i in range(len(fband)):
            ll_=sarrx[ll][mm]
            mm_=sarry[ll][mm]
            widths=np.arange(1,401)*0.5
            wdata=diff_Tb_smooth[ll_][mm_][i]/diff_Tb_std[ll_][mm_][i]
            #wdata=diff_Tb_all[ll_][mm_][i]/diff_Tb_std[ll_][mm_][i]
            cwtmatr[ll][mm][i], freqs[ll][mm][i] = pywt.cwt(wdata, widths, mother)
            cwtmatr[ll][mm][i]=cwtmatr[ll][mm][i]*diff_Tb_std[ll_][mm_][i]
            power[ll][mm][i] = (np.abs(cwtmatr[ll][mm][i])) ** 2
            power_flux[ll][mm][i]=ut.Tb2flux(np.sqrt(power[ll][mm][i]), 50, 50, freq[i]/1000.)**2
            period[ll][mm][i] = 1 / freqs[ll][mm][i]
            ###
            tbres=udata_dirty_all[i,:,ll_,mm_]*fact[i]
            wdata_res=tbres/np.std(tbres)
            cwtmatr_res[ll][mm][i], freqs_res[ll][mm][i] = pywt.cwt(wdata_res, widths, mother)
            cwtmatr_res[ll][mm][i]=cwtmatr_res[ll][mm][i]*np.std(tbres)
            power_res[ll][mm][i] = (np.abs(cwtmatr_res[ll][mm][i])) ** 2
            power_flux_res[ll][mm][i]=ut.Tb2flux(np.sqrt(power_res[ll][mm][i]), 50, 50, freq[i]/1000.)**2
            period_res[ll][mm][i] = 1 / freqs_res[ll][mm][i]
            ###
            slope_Tb[ll][mm][i]=ut.fit_1d(1/widths[20:50],power[ll][mm][i].mean(axis=1)[20:50])

cwtmatr=np.array(cwtmatr)
slope_Tb=np.array(slope_Tb)
power=np.array(power)
period=np.array(period)
power_flux=np.array(power_flux)
cwtmatr_res=np.array(cwtmatr_res)
power_res=np.array(power_res)
period_res=np.array(period_res)
power_flux_res=np.array(power_flux_res)
freqs=np.array(freqs)
freqs_res=np.array(freqs_res)
sar=ut.fit_1d(freqs.mean(axis=(1,2))[0][20:50],power_flux.mean(axis=(1,2,4))[0][20:50])
sch=ut.fit_1d(freqs.mean(axis=(1,2))[1][20:50],power_flux.mean(axis=(1,2,4))[1][20:50])
sqs=ut.fit_1d(freqs.mean(axis=(1,2))[2][20:50],power_flux.mean(axis=(1,2,4))[2][20:50])
lab=['AR','CH','QS']


plot_spec_region=0
if(plot_spec_region):
    plt.plot(freqs.mean(axis=(1,2))[0],power_flux.mean(axis=(1,2,4))[0],'o-',color='red',label='AR')
    plt.plot(freqs.mean(axis=(1,2))[1],power_flux.mean(axis=(1,2,4))[1],'o-',color='blue',label='CH')
    plt.plot(freqs.mean(axis=(1,2))[2],power_flux.mean(axis=(1,2,4))[2],'o-',color='green',label='QS')
    plt.plot(freqs.mean(axis=(1,2))[0],power_flux_res.mean(axis=(1,2,4))[0],'o--',color='red',label='AR (res)')
    plt.plot(freqs.mean(axis=(1,2))[1],power_flux_res.mean(axis=(1,2,4))[1],'o--',color='blue',label='CH (res)')
    plt.plot(freqs.mean(axis=(1,2))[2],power_flux_res.mean(axis=(1,2,4))[2],'o--',color='green',label='QS (res)')
    #plt.plot(freqs.mean(axis=(1,2))[0],10**(np.log10(freqs.mean(axis=(1,2))[0])*sar[0]+sar[1]),color='red',label='AR fit: '+str(np.round(sar[0],2))+'$\pm$'+str(np.round(sar[2],2)))
    #plt.plot(freqs.mean(axis=(1,2))[1],10**(np.log10(freqs.mean(axis=(1,2))[1])*sch[0]+sch[1]),color='blue',label='CH fit: '+str(np.round(sch[0],2))+'$\pm$'+str(np.round(sch[2],2)))
    #plt.plot(freqs.mean(axis=(1,2))[2],10**(np.log10(freqs.mean(axis=(1,2))[2])*sqs[0]+sqs[1]),color='green',label='QS fit: '+str(np.round(sqs[0],2))+'$\pm$'+str(np.round(sqs[2],2)))
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power (SFU$^2$)')
    plt.show()

#### EUV #####
aiafile='/home/i4ds1807205/20151203/saia_00193_fd_20151203_192129.fts'
h=fits.open(aiafile)
d=h[0].data
x1=np.arange(d.shape[0])
y1=np.arange(d.shape[1])
X1,Y1=np.meshgrid(x1,y1)
X11=(X1-d.shape[0]*0.5)*h[0].header['CDELT1']
Y12=(Y1-d.shape[1]*0.5)*h[0].header['CDELT2']

sig_quiet=[0.199,0.171,0.058,0.156,0.172,0.19,0.215,0.268,0.372]
sig_active=[0.118,0.146,0.163,0.175,0.313,0.424,0.545,0.566,0.79]

######## PLOT
print 'Plotting ...'
compare_flux=1
if(compare_flux):
    plt.plot(freq,flux,'o-',color='red')
    plt.errorbar(freq,flux,yerr=flux_std,color='red',label='MWA MAPS')
    plt.plot(freq,flux_fwd_int,'o-',label='FORWARD MAPS')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend()
    plt.show()
sys.exit()


plot_sigma=0
if(plot_sigma):
    plt.plot(flist,Ssun_mean_time,'o-',color='red',label='03-12-2015')
    plt.errorbar(flist,Ssun_mean_time,yerr=Ssun_rms_baseline,color='red')
    plt.legend(loc=2)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')


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
        plt.show()


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


import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pywt
import sunpy
from surya.utils import main as ut
import pickle
import pycwt as wavelet
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import matplotlib
import os
from sunpy.map import Map
from astropy import units as u
from matplotlib.patches import Ellipse
from scipy.io import readsav
from scipy.optimize import curve_fit 

#allmaps=pickle.load(open('/media/rohit/VLA/20160409/20160409_submap_50ms.p','rb'))
xcimax,ycimax,xci90,yci90,maxTbi,Tbi_r1,Tbi_r2,areai50,eTbi=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_i.p','rb'),encoding='latin-1')
Tbi_r1=np.array(Tbi_r1).reshape(32,2000);Tbi_r2=np.array(Tbi_r2).reshape(32,2000);max_Tbi=np.array(maxTbi).reshape(32,2000)
xcimax=np.array(xcimax).reshape(32,2000);ycimax=np.array(ycimax).reshape(32,2000);areai50=np.array(areai50).reshape(32,2000)
xcvmax,ycvmax,xcv90,ycv90,maxTbv,Tbv_r1,Tbv_r2,eTbv=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_v.p','rb'),encoding='latin-1')
xclmax,yclmax,xcl90,ycl90,maxTbl,Tbl_r1,Tbl_r2,eTbl=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_l.p','rb'),encoding='latin-1')
xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r.p','rb'),encoding='latin-1')
qsx_0,qsy_0,qsxcr90_0,qsycr90_0,qsmaxTbr_0,qsTbr_r1_0,qsTbr_r2_0,qsarear50_0,qstimevla_0=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs.p','rb'),encoding='latin-1')
qsx_0,qsy_0,qsmaxTbr_0,qsTbr_r1_0,qsTbr_r2_0,qsarear50_0=np.array(qsx_0),np.array(qsy_0),np.array(qsmaxTbr_0),np.array(qsTbr_r1_0),np.array(qsTbr_r2_0),np.array(qsarear50_0)
qsx_3,qsy_3,qsxcr90_3,qsycr90_3,qsmaxTbr_3,qsTbr_r1_3,qsTbr_r2_3,qsarear50_3,qstimevla_3=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs_3.p','rb'),encoding='latin-1')
qsx_3,qsy_3,qsmaxTbr_3,qsTbr_r1_3,qsTbr_r2_3,qsarear50_3=np.array(qsx_3),np.array(qsy_3),np.array(qsmaxTbr_3),np.array(qsTbr_r1_3),np.array(qsTbr_r2_3),np.array(qsarear50_3)
qsx_5,qsy_5,qsxcr90_5,qsycr90_5,qsmaxTbr_5,qsTbr_r1_5,qsTbr_r2_5,qsarear50_5,qstimevla_5=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs_5.p','rb'),encoding='latin-1')
qsx_5,qsy_5,qsmaxTbr_5,qsTbr_r1_5,qsTbr_r2_5,qsarear50_5=np.array(qsx_5),np.array(qsy_5),np.array(qsmaxTbr_5),np.array(qsTbr_r1_5),np.array(qsTbr_r2_5),np.array(qsarear50_5)
qsmaxTbr_0[2132:2170]=2.5e7;qsmaxTbr_0[4531:4569]=2.5e7;qsmaxTbr_3[2132:2170]=0.9e7;qsmaxTbr_3[4531:4569]=0.9e7;qsmaxTbr_5[2132:2170]=0.65e7;qsmaxTbr_5[4531:4569]=0.65e7
ds_=aa=np.load('/media/rohit/VLA/20160409/sun_L_20160409.1s.ms.dspec.npz')
#timevla=allmaps['vla']['timevla']
#time94=allmaps['aia94']['time94'];time131=allmaps['aia131']['time131'];time335=allmaps['aia335']['time335']
#time1600=allmaps['aia1600']['time1600'];time1700=allmaps['aia1700']['time1700']
freq=np.round(np.linspace(0.997,2.005,32),3)
############################################
listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_4/*spw.*0-15*.FITS'))[0:2000]
mapp=Map(sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.*0-15*.FITS'))[800]);mapp.peek();plt.show()
l=0;bmin=[0]*32;bmaj=[0]*32;spwlist=['0-15','16-31','32-47','48-63']
for i in range(8):
    for j in range(4):
        ll=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_'+str(i)+'/*spw.*'+str(spwlist[j])+'*.FITS'))[0]
        mm=Map(ll);bmaj[l]=mm.meta['bmin']*3600;bmin[l]=mm.meta['bmaj']*3600
        l=l+1
bmaj=np.array(bmaj);bmin=np.array(bmin)
############################################
h=np.arange(200)*1400/1.e3 # in Mm
h1=np.arange(200)*500/1.e3 # in Mm
sh=np.linspace(10,500,100)
ux=np.linspace(0,50,50)
dem_arr=[0]*100;ne=[0]*100;ne_dem=[0]*100;wpe_dem=[0]*100;omegap_dem=[0]*100
tem=2.6e27
for i in range(100):
    dem_arr[i]=[0]*200;ne[i]=[0]*200
    #ne[i]=1.16e17*(1+(h*1000/(10*sh[i])))**(-10+1)
    ne[i]=1.16e17*(np.exp(-1*h1*1000/sh[i]))+9.5e8*np.exp(-1*h1*1000/31.e3)#+4.2e4 *10**(4.32/((695700.+h*1000)/695700.))
    #ne[i]=4.2e4 *10**(4.32/((695700.+h*1000)/695700.))+1.16e17*(np.exp(-1*h*1000/sh[i]))
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+1.e6
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+(3.09e8*(1-0.5*np.sin(th))/((695700.+h*1000)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+h*1000)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+h*1000)/695700.)**2.5)
    for j in range(200):
        dem_arr[i][j]=np.sum(ne[i][j:]*ne[i][j:])*(h[1]-h[0])*1.e8
    idx=ut.find_nearest(dem_arr[i],tem)[0]
    ne_dem[i]=ne[i][idx];wpe_dem[i]=9000*np.sqrt(ne_dem[i])/1.e9

dem_arr=np.array(dem_arr);ne=np.array(ne);ne_dem=np.array(ne_dem);wpe_dem=np.array(wpe_dem)
expolB_=readsav('/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav')
expolB=expolB_['box']
bx,by,bz=expolB['bx'][0],expolB['by'][0],expolB['bz'][0]
babs=np.sqrt(bx*bx+by*by+bz*bz)
dr=expolB['dr'] # 3 values of extrapolated box in units Rsun
b1=280;b2=170
idxb1=np.where((babs<b1+1) & (babs>b1-1))
idxb2=np.where((babs<b2+1) & (babs>b2-1))

hh2_0=np.array([43.4,40.5,40.1,39.2,38.5,37.7,37.1,36.3,36,35.6,35,34.1,33.5,33.1,32.6,32.1,31.6,31.1,30.8,30.3,29.8,29.1,28.8,28.3,27.9,27.5,27.2,27,26.5,26.2,25.8])
babs_h2=freq/0.0028/2;hh2=[0]*len(babs_h2);ne2=[0]*len(babs_h2)
for i in range(len(babs_h2)):
    hh2[i]=ut.find_nearest(babs[:,160,172],babs_h2[i])[0]*1400/1000.
    ne2[i]=ne[28][ut.find_nearest(babs[:,160,172],babs_h2[i])[0]]
ne2=np.array(ne2);babs_h2=np.array(babs_h2)
vA=2.18e6/np.sqrt(ne2)*babs_h2/1.e3 #in Mm/s


############################################
get_qp=0
if(get_qp):
    listqp=sorted(glob.glob('/media/rohit/VLA/20160409/blob/qb/*.p'))
    #listqp=sorted(glob.glob('/media/rohit/VLA1/20160409/blob/b_new/*spw_0_16-31*.p'),key=len)
    x0=np.zeros(len(listqp));y0=np.zeros(len(listqp));sigx0=np.zeros(len(listqp));sigy0=np.zeros(len(listqp));rot0=np.zeros(len(listqp));tb0=np.zeros(len(listqp))
    x1=np.zeros(len(listqp));y1=np.zeros(len(listqp));sigx1=np.zeros(len(listqp));sigy1=np.zeros(len(listqp));rot1=np.zeros(len(listqp));tb1=np.zeros(len(listqp));s0=[0]*len(listqp);s1=[0]*len(listqp);s2=[0]*len(listqp)
    for i in range(len(listqp)):
        aa=pickle.load(open(listqp[i],'rb'))
        s0[i]=aa[3];s1[i]=aa[4];s2[i]=aa[5];x0[i],y0[i],sigx0[i],sigy0[i],rot0[i],tb0[i]=aa[0];x1[i],y1[i],sigx1[i],sigy1[i],rot1[i],tb1[i]=aa[1]

get_gauss=0
if(get_gauss):
    spwlist=['0-15','16-31','32-47','48-63']
    #spwl=[0,1,2,3,4,5,6,7];m=len(spwlist);n=len(spwl)
    #spwlist=['0-15']#,'16-31','32-47','48-63']
    spwl=[0,1,2,3,4,5,6,7];m=len(spwlist);n=len(spwl);tt=2000
    x0=[0]*n;y0=[0]*n;sigx0=[0]*n;sigy0=[0]*n;rot0=[0]*n;tb0=[0]*n;x1=[0]*n;y1=[0]*n;sigx1=[0]*n;sigy1=[0]*n;rot1=[0]*n;tb1=[0]*n;s0=[0]*n;s1=[0]*n;s2=[0]*n;w=[0]*n
    for s in range(len(spwl)):
        print(spwl[s])
        j=0;x0[s]=[0]*m;y0[s]=[0]*m;sigx0[s]=[0]*m;sigy0[s]=[0]*m;rot0[s]=[0]*m;tb0[s]=[0]*m;x1[s]=[0]*m;y1[s]=[0]*m;sigx1[s]=[0]*m;sigy1[s]=[0]*m;rot1[s]=[0]*m;tb1[s]=[0]*m;s0[s]=[0]*m;s1[s]=[0]*m;s2[s]=[0]*m;w[s]=[0]*m
        for ss in spwlist:
            x0[s][j]=np.zeros(tt);y0[s][j]=np.zeros(tt);sigx0[s][j]=np.zeros(tt);sigy0[s][j]=np.zeros(tt);rot0[s][j]=np.zeros(tt);tb0[s][j]=np.zeros(tt);w[s][j]=np.zeros(tt)
            x1[s][j]=np.zeros(tt);y1[s][j]=np.zeros(tt);sigx1[s][j]=np.zeros(tt);sigy1[s][j]=np.zeros(tt);rot1[s][j]=np.zeros(tt);tb1[s][j]=np.zeros(tt);s0[s][j]=[0]*tt;s1[s][j]=[0]*tt;s2[s][j]=[0]*tt
            for k in range(tt):
                kk="%04d" %k;f='/media/rohit/VLA/20160409/blob/b_new/blob_param_spw_'+str(spwl[s])+'_'+str(ss)+'_'+str(kk)+'.p'
                if(os.path.isfile(f)==False):
                    k=k-1;kk="%04d" %k
                    s0_=aa[3];s1_=aa[4];s2_=aa[5];x0[s][j][k],y0[s][j][k],sigx0[s][j][k],sigy0[s][j][k],rot0[s][j][k],tb0[s][j][k]=aa[0];x1[s][j][k],y1[s][j][k],sigx1[s][j][k],sigy1[s][j][k],rot1[s][j][k],tb1[s][j][k]=aa[1];w[s][j][k]=aa[2]
                else:
                    ff=open(f,'rb');aa=pickle.load(ff)
                    #s0[s][j][k]=aa[3];s1[s][j][k]=aa[4];s2[s][j][k]=aa[5];x0[s][j][k],y0[s][j][k],sigx0[s][j][k],sigy0[s][j][k],rot0[s][j][k],tb0[s][j][k]=aa[0];x1[s][j][k],y1[s][j][k],sigx1[s][j][k],sigy1[s][j][k],rot1[s][j][k],tb1[s][j][k]=aa[1]
                    s0_=aa[3];s1_=aa[4];s2_=aa[5];x0[s][j][k],y0[s][j][k],sigx0[s][j][k],sigy0[s][j][k],rot0[s][j][k],tb0[s][j][k]=aa[0];x1[s][j][k],y1[s][j][k],sigx1[s][j][k],sigy1[s][j][k],rot1[s][j][k],tb1[s][j][k]=aa[1];w[s][j][k]=aa[2]
            j=j+1
    x0=np.array(x0).reshape(n*m,tt);y0=np.array(y0).reshape(n*m,tt);sigx0=np.array(sigx0).reshape(n*m,tt);sigy0=np.array(sigy0).reshape(n*m,tt);rot0=np.array(rot0).reshape(n*m,tt);tb0=np.array(tb0).reshape(n*m,tt)
    x1=np.array(x1).reshape(n*m,tt);y1=np.array(y1).reshape(n*m,tt);sigx1=np.array(sigx1).reshape(n*m,tt);sigy1=np.array(sigy1).reshape(n*m,tt);rot1=np.array(rot1).reshape(n*m,tt);tb1=np.array(tb1).reshape(n*m,tt)
    sys.exit()
    pickle.dump([x0,y0,sigx0,sigy0,rot0,tb0,x1,y1,sigx1,sigy1,rot1,tb1],open('/media/rohit/VLA/20160409/blob/all_params.p','wb'))
x0,y0,sigx0,sigy0,rot0,tb0,x1,y1,sigx1,sigy1,rot1,tb1=pickle.load(open('/media/rohit/VLA/20160409/blob/all_params.p','rb'),encoding='latin-1')
#[[qx0,qy0,qsigx0,qsigy0,qrot0,qtb0],[qx1,qy1,qsigx1,qsigy1,qrot1,qtb1]]=pickle.load(open('/media/rohit/VLA/20160409/blob/blob_qall.p','rb'))
#[[x0,y0,sigx0,sigy0,rot0,tb0],[x1,y1,sigx1,sigy1,rot1,tb1]]=pickle.load(open('/media/rohit/VLA/20160409/blob/blob_all.p','rb'))


plot_tb_line=1
if(plot_tb_line):
    f,ax=plt.subplots(2,1);nn=7
    for i in range(nn):
        ax[0].plot(np.arange(len(tb0[i]))*0.05,tb0[i]/1.e6,'-',label=str(freq[4*i])+' MHz')
        ax[1].plot(np.arange(len(tb1[i]))*0.05,tb1[i]/1.e6,'-',label=str(freq[4*i])+' MHz')
    ax[0].legend(loc=2);ax[1].legend(loc=2);ax[0].set_ylabel('T$_B$ (MK)');ax[1].set_ylabel('T$_B$ (MK)');ax[1].set_xlabel('Time (sec)')
    ax[0].set_title('Source Source');ax[1].set_title('North Source')
    plt.show()

plot_area=1
if(plot_area):
    plt.plot(freq[0:18],np.sqrt(areai50[:,858][0:18]/np.pi),'o-',label='Source Area (50%)')
    plt.plot(freq[0:18],np.sqrt(bmaj[0:18]*bmin[0:18]/4.),'-',label='PSF')
    plt.xlabel('Frequency (GHz)');plt.ylabel('Effective Radius (arcsec)');plt.legend();plt.ylim(5,20)
    plt.show()

plot_tb_ns=1
if(plot_tb_ns):
    f,ax=plt.subplots(2,1);n1=0;n2=4
    ax[0].plot(np.arange(len(tb0[n1]))*0.05,tb1[n1]/1.e6,'-',label='North')
    ax[0].plot(np.arange(len(tb1[n1]))*0.05,tb0[n1]/1.e6,'-',label='South')
    ax[1].plot(np.arange(len(tb0[n2]))*0.05,tb1[n2]/1.e6,'-',label='North')
    ax[1].plot(np.arange(len(tb1[n2]))*0.05,tb0[n2]/1.e6,'-',label='South')
    ax[0].legend(loc=2);ax[1].legend(loc=2);ax[0].set_ylabel('T$_B$ (MK)');ax[1].set_ylabel('T$_B$ (MK)');ax[1].set_xlabel('Time (sec)')
    ax[0].set_title(str(freq[4*n1])+' MHz');ax[1].set_title(str(freq[4*n2])+' MHz')
    plt.show()

plot_ns_ts=0
if(plot_ns_ts):
    f,ax=plt.subplots(2,1)
    ax[0].plot(np.arange(len(qsmaxTbr)+2000)*0.05,np.hstack((qsmaxTbr,maxTbr[0]))/1.e6,label='Max $T_B$')
    ax[1].plot(np.arange(len(tb0)+len(qtb0))*0.05,np.hstack((qtb0,tb0))/1.e6,label='North Source')
    ax[1].plot(np.arange(len(tb1)+len(qtb1))*0.05,np.hstack((qtb1,tb1))/1.e6,label='South Source')
    #ax[1].plot(np.arange(len(tb1)+len(qtb1))*0.05,(np.hstack((qtb1,tb1))+np.hstack((qtb0,tb0)))/1.e6,label='South Source')
    ax[0].legend();ax[1].legend()
    ax[1].set_xlabel('Time (sec)');ax[0].set_ylabel('T$_B$ (MK)')
    plt.show()


def do_wavelet(Tbr, N):
    #datawave_=(Tbr-Tbr[0])/1.e6
    datawave_=Tbr/1.e6
    #datawave_=np.convolve(datawave_, np.ones(N)/N, mode='valid')
    datawave=datawave_[int(N/2):int(-1*N/2+1)]-np.convolve(datawave_, np.ones(N)/N, mode='valid')
    std=np.std(datawave);datawave_std=datawave/std
    timewave=np.arange(len(datawave))*0.05
    mother=wavelet.Morlet(6)#mother=wavelet.MexicanHat()
    dt=0.05;s0 = 2 * dt;dj = 1. / 12 # Lowest scale s0
    J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
    alpha, _, _ = wavelet.ar1(datawave_std)
    wave, scales, freqs, coi1, fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
    power = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power /= scales[:, None]
    return datawave_std,power,iwave,coi1,period

def plot_wavelet(Tb,t,dt,period,power,coi1,fil,tfreq):
    plt.ioff();label='T$_B$';units='(MK)'
    figprops = dict(figsize=(20, 20), dpi=72)
    fig = plt.figure(**figprops)
    ax = plt.axes([0.1, 0.68, 0.67, 0.2])
    #ax.plot(t, iwave[i], '-', linewidth=1, color=[0.5, 0.5, 0.5])
    ax.plot(t,Tb, 'k', linewidth=1.5)
    ax.grid(True)
    ax.set_ylabel(r'{} [{}]'.format(label, units))
    bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax);levels = np.linspace(1,100,90)
    #bx.contourf(t, periodds, np.log2(powerds), np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
    bx.contourf(t, period, np.log2(power),np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
    extent = [t.min(), t.max(), 0, max(period)]
    bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([coi1, [1e-9], period[-1:],period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
    #bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[i])+' GHz({})'.format(label, 'MORLET'))
    bx.set_ylabel('Period (sec)');bx.set_xlabel('Time (sec)');bx.set_yscale('log')
    bx.grid(True)
    cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
    cx.plot(power.mean(axis=1), period, 'k-', linewidth=1.5)
    cx.set_title('Global Wavelet Spectrum')
    cx.set_xlabel(r'Power [({})^2]'.format(units))
    #cx.set_ylim(([periodds.min(), periodds.max()]))
    cx.grid(True)
    cx.set_xscale('log');plt.title('Frequency: '+str(tfreq)+' GHz')
    plt.savefig(fil,dpi=100)
    plt.close()

if_qs=1
if(if_qs):
    N=400
    #Tbrs=np.hstack((qtb1,tb1));Tbrn=np.hstack((qtb0,tb0))
    #Tbrs=np.hstack((tb1));Tbrn=np.hstack((tb0))
    #Tbrs_wave,s_power,s_iwave,s_coi1,period=do_wavelet(Tbrs, N)
    #Tbrn_wave,n_power,n_iwave,n_coi1,period=do_wavelet(Tbrn, N)
    qTbmax_wave0,s_power0,s_iwave0,s_coi10,period=do_wavelet(np.array(qsmaxTbr_0), N)
    qTbmax_wave3,s_power3,s_iwave3,s_coi13,period=do_wavelet(np.array(qsmaxTbr_3), N)
    qTbmax_wave5,s_power5,s_iwave5,s_coi15,period=do_wavelet(np.array(qsmaxTbr_5), N)
    dt=0.05;tq=np.arange(len(qTbmax_wave0))*dt
    #plot_wavelet(qTbmax_wave0,t,dt,period,s_power0,s_coi10,'power0.png','')
    #plot_wavelet(Tbrn_wave,t,dt,period,s_power3,n_coi1)

N=400
Tbrs_wave=[0]*x0.shape[0];Tbrn_wave=[0]*x0.shape[0];s_power=[0]*x0.shape[0];n_power=[0]*x0.shape[0]
s_coi1=[0]*x0.shape[0];n_coi1=[0]*x0.shape[0];ns_power=[0]*x0.shape[0];nn_power=[0]*x0.shape[0];Tbmax_wave=[0]*x0.shape[0];max_power=[0]*x0.shape[0];max_coi=[0]*x0.shape[0]
for i in range(x0.shape[0]):
    Tbrs_wave[i],s_power[i],s_iwave,s_coi1[i],period=do_wavelet(tb0[i], N);ns_power[i]=s_power[i]/np.nanmax(s_power[i])
    Tbrn_wave[i],n_power[i],n_iwave,n_coi1[i],period=do_wavelet(tb1[i], N);nn_power[i]=n_power[i]/np.nanmax(n_power[i])
    Tbmax_wave[i],max_power[i],max_iwave,max_coi[i],period=do_wavelet(maxTbr[i], N)
    dt=0.05;t=np.arange(len(Tbrs_wave[i]))*dt
    #plot_wavelet(Tbrs_wave[i],t,dt,period,s_power[i],s_coi1[i],'/home/i4ds1807205/Dropbox/20160409/wavelets/blob/south_wavelet_'+str(i)+'.png',freq[i])
    #plot_wavelet(Tbrn_wave[i],t,dt,period,n_power[i],n_coi1[i],'/home/i4ds1807205/Dropbox/20160409/wavelets/blob/north_wavelet_'+str(i)+'.png',freq[i])

s_power=np.array(s_power);n_power=np.array(n_power);glbl_spower=np.mean(s_power,axis=2);glbl_npower=np.mean(n_power,axis=2);glbl_maxpower=np.mean(max_power,axis=2)
ns_power=np.array(ns_power);nn_power=np.array(nn_power);nglbl_spower=np.mean(ns_power,axis=2);nglbl_npower=np.mean(nn_power,axis=2)
#------------ Plot Wavelet Paper

f,ax=plt.subplots(1,1)
ax.imshow(max_power[0],aspect='auto',origin='lower')
ax.set_ylabel('Wavelet Scales (sec)')
ax.set_xlabel('Time (sec)')
plt.show()

i=0
plt.ioff();label='T$_B$';units='MK'
figprops = dict(figsize=(20, 20), dpi=72)
fig = plt.figure(**figprops)
ax = plt.axes([0.1, 0.68, 0.61, 0.2])
ax.plot(t,Tbmax_wave[i], 'k', linewidth=1.5)
ax.grid(True)
ax.set_ylabel(r'{} [{}]'.format(label, units))
bx = plt.axes([0.1, 0.17, 0.61, 0.48], sharex=ax);levels = np.linspace(1,100,90)
bx.contourf(t, period, np.log2(max_power[0]),np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
extent = [t.min(), t.max(), 0, max(period)]
bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([max_coi[i], [1e-9], period[-1:],period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
#bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[i])+' GHz({})'.format(label, 'MORLET'))
bx.set_ylabel('Period (sec)');bx.set_xlabel('Time (sec)');bx.set_yscale('log')
bx.grid(True);bx.set_ylim(0.1,27);bx.set_xlim(0,80)
cx = plt.axes([0.75, 0.17, 0.15, 0.48], sharey=bx)
cx.plot(max_power[i].mean(axis=1), period, 'k-', linewidth=1.5)
cx.set_title('Global Wavelet Spectrum');cx.set_xlabel(r'Power (${}^2$)'.format(units))
cx.grid(True);cx.set_xscale('log');plt.title('Frequency: '+str(freq[i])+' GHz')
plt.show()

plt.ioff();label='T$_B$';units='MK'
figprops = dict(figsize=(20, 20), dpi=72)
fig = plt.figure(**figprops)
ax = plt.axes([0.1, 0.68, 0.61, 0.2])
ax.plot(tq,qTbmax_wave0, 'k', linewidth=1.5)
ax.grid(True)
ax.set_ylabel(r'{} [{}]'.format(label, units))
bx = plt.axes([0.1, 0.17, 0.61, 0.48], sharex=ax);levels = np.linspace(1,100,90)
bx.contourf(tq, period, np.log2(s_power0),np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
extent = [tq.min(), tq.max(), 0, max(period)]
bx.fill(np.concatenate([tq, tq[-1:] + dt, tq[-1:] + dt, tq[:1] - dt, tq[:1] - dt]),np.concatenate([s_coi10, [1e-9], period[-1:],period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
bx.set_ylabel('Period (sec)');bx.set_xlabel('Time (sec)');bx.set_yscale('log')
bx.grid(True);bx.set_ylim(0.1,27);bx.set_xlim(0,220)
cx = plt.axes([0.75, 0.17, 0.15, 0.48], sharey=bx)
cx.plot(s_power0.mean(axis=1), period, 'k-', linewidth=1.5)
cx.set_title('Global Wavelet Spectrum');cx.set_xlabel(r'Power (${}^2$)'.format(units))
cx.grid(True);cx.set_xscale('linear');plt.title('Frequency: '+str(freq[i])+' GHz')
plt.show()

#------------ Wait-time distribution
def fit_poisson(x,lam):
    return lam*np.exp(-lam*x)

def fit_linear(x,m,c):
    return m*x+c

Pwtime=[0]*3;ts1=[0]*3;numf=[0]*3;wtime=[0]*3;Pwtime_fit=[0]*3;popt=[0]*3;wtime_int=[0]*3;epopt=[0]*3
qTbmax_wave=[qTbmax_wave0,qTbmax_wave3,qTbmax_wave5]
for j in range(3):
    lin=np.concatenate((qTbmax_wave[j],Tbmax_wave[j*15]*10))
    twt=np.arange(len(lin))*0.05
    idx=np.where(lin>0*np.std(lin[0:1000]));lin_=lin[idx];Pwtime[j]=[0]*len(lin_);tw=twt[idx]
    diffs = np.diff(idx)[0] != 1;indexes = np.nonzero(diffs)[0] + 1;ts1[j] = np.split(idx[0], indexes);numf[j]=len(ts1[j])
    wtime[j]=[0]*numf[j];wtime[j][0]=0
    for i in range(1,numf[j]):
        wtime[j][i]=(ts1[j][i][0]-ts1[j][i-1][-1])*0.05
    Pwtime[j]=np.histogram(wtime[j],bins=10)[0];wtime_int[j]=np.histogram(wtime[j],bins=10)[1]
    #Pwtime[j]=[0]*numf[j]
    #for i in range(numf[j]):
    #    Pwtime[j][i]=np.exp(-wtime[j][i]/np.mean(wtime[j][0:i]))/np.mean(wtime[j][0:i])
    #popt[j], _ = curve_fit(fit_function,np.array(wtime[j])[2:], np.array(Pwtime[j])[2:]);popt[j]=popt[j][0]
    #Pwtime_fit[j]=fit_function(np.array(wtime[j])[2:], popt[j])
    popt[j], cov = curve_fit(fit_linear,np.log10(wtime_int[j][1:]),np.log10(Pwtime[j]+1));epopt[j]=np.sqrt(cov[0][0])
    Pwtime_fit[j]=10**fit_linear(np.log10(np.linspace(0.1,10,100)), popt[j][0],popt[j][1])

plt.plot(np.array(wtime_int[0])[1:],Pwtime[0],'o',label='1 GHz',color='r')
plt.plot(np.linspace(0.1,10,100),Pwtime_fit[0],'-',color='r',label='$\lambda$='+str(np.round(popt[0][0],2))+'$\pm$'+str(np.round(epopt[0],2)))
plt.plot(np.array(wtime_int[1])[1:],Pwtime[1],'o',label='1.5 GHz',color='g')
plt.plot(np.linspace(0.1,10,100),Pwtime_fit[1],'-',color='g',label='$\lambda$='+str(np.round(popt[1][0],2))+'$\pm$'+str(np.round(epopt[1],2)))
plt.plot(np.array(wtime_int[2])[1:],Pwtime[2],'o',label='2 GHz',color='b')
plt.plot(np.linspace(0.1,10,100),Pwtime_fit[2],'-',color='b',label='$\lambda$='+str(np.round(popt[2][0],2))+'$\pm$'+str(np.round(epopt[2],2)))
plt.legend();plt.xlabel('Wait Time (sec)');plt.ylabel('Occupancy (arbitrary)')
plt.show()



#------------ Spectral density
from scipy.signal import periodogram as ped
sd_ped=[0]*32;sdfit=[0]*32;yfit=[0]*32;errfit=[0]*32
for i in range(32):
    f, sd_ped[i] = ped(Tbmax_wave[i],20,scaling='density')
    xfit=np.log10(f[200:800]);p=np.polyfit(xfit,np.log10(sd_ped[i][200:800]),1,cov=True);sdfit[i]=p[0][0];yfit[i]=10**(p[0][0]*np.log10(f[50:1000])+p[0][1]);errfit[i]=3*np.sqrt(np.diag(p[1]))[0]
sd_ped=np.array(sd_ped);sdfit=np.array(sdfit);yfit=np.array(yfit);errfit=np.array(errfit)
qf, qsd_ped0 = ped(qTbmax_wave0,20,scaling='density');qxfit0=np.log10(qf[200:800]);p=np.polyfit(qxfit0,np.log10(qsd_ped0[200:800]),1,cov=True);qsdfit0=p[0][0];qyfit0=10**(p[0][0]*qxfit0+p[0][1]);qerrfit0=3*np.sqrt(np.diag(p[1]))[0]
qf, qsd_ped3 = ped(qTbmax_wave3,20,scaling='density');qxfit3=np.log10(qf[200:800]);p=np.polyfit(qxfit3,np.log10(qsd_ped3[200:800]),1,cov=True);qsdfit3=p[0][0];qyfit3=10**(p[0][0]*qxfit3+p[0][1]);qerrfit3=3*np.sqrt(np.diag(p[1]))[0]
qf, qsd_ped5 = ped(qTbmax_wave5,20,scaling='density');qxfit5=np.log10(qf[200:800]);p=np.polyfit(qxfit5,np.log10(qsd_ped5[200:800]),1,cov=True);qsdfit5=p[0][0];qyfit5=10**(p[0][0]*qxfit5+p[0][1]);qerrfit5=3*np.sqrt(np.diag(p[1]))[0]

plt.plot(f/vA[6],sd_ped[6],label=str(freq[6])+' GHz',color='red')
plt.plot(f/vA[30],sd_ped[30],label=str(freq[30])+' GHz',color='green')
plt.plot(f[50:1000]/vA[6],yfit[6],'--',color='b',label='Fit: '+str(np.round(sdfit[6],3))+'$\pm$'+str(np.round(errfit[6],3)))
plt.plot(f[50:1000]/vA[30],yfit[30],'--',color='k',label='Fit: '+str(np.round(sdfit[30],3))+'$\pm$'+str(np.round(errfit[30],3)))
plt.xlabel('Wavenumber (Mm$^{-1}$)');plt.xscale('log');plt.yscale('log');plt.ylabel('PSD');plt.ylim(10**(-5),10**2);plt.legend(loc=3);plt.show()

plt.plot(freq,sdfit*-1,'o');plt.errorbar(freq,sdfit*-1,yerr=errfit,color='b',label='bursts');plt.xlabel('Frequency (GHz)');plt.axhline(y=1.67,linestyle='--',label='Kolmogorov index',color='k')
plt.plot([freq[0],freq[12],freq[22]],[qsdfit0*-1,qsdfit3*-1,qsdfit5*-1],'o',color='r');plt.errorbar([freq[0],freq[12],freq[22]],[qsdfit0*-1,qsdfit3*-1,qsdfit5*-1],yerr=[qerrfit0,qerrfit3,qerrfit5],color='r',label='pre-bursts')
plt.legend(loc=3);plt.ylim(-1,3);plt.ylabel('PSD index');plt.show()

print('Frequency (GHz) & B (G) & $n_e$ ($\\times 10$^{8}$ cm^{-3}$) & v$_A$ (Mm/s) \\\\ \\hline')
for i in range(8):
    print(np.round(freq[i*4+1],3),'&',np.round(babs_h2[i*4+1],1),'&',np.round(ne2[i*4+1]/1.e8,1),'&',np.round(vA[i*4+1],1),'\\\\')
#------------ test correlate
xx=np.linspace(0,10,100);x=np.linspace(-5,5,100);y1=np.exp(-x**2)*np.sin(x);y2=np.exp((-(x-0)**2))*np.sin(x-0)+np.exp((-(x-2)**2)/2.0)*np.sin(x-2)*0.5
co=np.correlate(y2,y1, "same")
f,ax=plt.subplots(2,1);ax0=ax[0];ax1=ax[1]
ax0.plot(xx,y1,label='series 1');ax0.plot(xx,y2,label='series 2');ax0.legend();ax1.plot(x,co)
plt.show()

#------------ Correlate

ccqs03=np.correlate(qTbmax_wave0,qTbmax_wave3, "same");ccqs53=np.correlate(qTbmax_wave5,qTbmax_wave3, "same");ccqs33=np.correlate(qTbmax_wave3,qTbmax_wave3, "same")
t1=ut.find_nearest(ccqs03[2000:2100],np.max(ccqs03[2000:2200])*0.5)[0];t2=ut.find_nearest(ccqs03[2100:2200][::-1],np.max(ccqs03[2000:2200])*0.5)[0];del03=-1*(t2-t1+1)*0.05
t1=ut.find_nearest(ccqs53[2000:2100],np.max(ccqs53[2000:2200])*0.5)[0];t2=ut.find_nearest(ccqs53[2100:2200][::-1],np.max(ccqs53[2000:2200])*0.5)[0];del53=-1*(t2-t1+1)*0.05
plt.plot((np.arange(4200)-2100)*0.05,ccqs03/np.max(ccqs03),'o-',label=str(freq[0])+' GHz');plt.plot((np.arange(4200)-2100)*0.05,ccqs33/np.max(ccqs33),'o-',label=str(freq[10])+' GHz');plt.plot((np.arange(4200)-2100)*0.05,ccqs53/np.max(ccqs53),'o-',label=str(freq[20])+' GHz')
plt.xlabel('Lag Time (s)');plt.ylabel('Cross-Correlation');plt.legend();plt.xlim(-5,5);plt.show()

cc0_s=[0]*32;cc9_s=[0]*32;auto_s=[0]*32
cc0_n=[0]*32;cc9_n=[0]*32;auto_n=[0]*32
pcc0_s=[0]*32;pcc9_s=[0]*32;pauto_s=[0]*32
pcc0_n=[0]*32;pcc9_n=[0]*32;pauto_n=[0]*32;delt=[0]*32;pdelt=[0]*32
for i in range(32):
    #ts0=Tbmax_wave[11];ts1=Tbmax_wave[i]
    ts0=Tbmax_wave[11][300:500];ts1=Tbmax_wave[i][300:500]
    ts9=Tbi_r1[9];tn0=Tbrn_wave[0];tn1=Tbrn_wave[i];tn9=Tbrn_wave[9]
    pts0=glbl_spower[0];pts1=glbl_spower[i];pts9=glbl_spower[9];ptn0=glbl_npower[0];ptn1=glbl_npower[i];ptn9=glbl_npower[9]
    cc0_s[i]=np.correlate(ts1,ts0, "same");cc9_s[i]=np.correlate(ts9,ts1, "same");auto_s[i]=np.correlate(ts0,ts0, "same")
    cc0_n[i]=np.correlate(tn1,tn0, "same");cc9_n[i]=np.correlate(tn9,tn1, "same");auto_n[i]=np.correlate(tn0,tn0, "same")
    pcc0_s[i]=np.correlate(pts1,pts0, "same");pcc9_s[i]=np.correlate(pts9,pts1, "same");pauto_s[i]=np.correlate(pts0,pts0, "same")
    pcc0_n[i]=np.correlate(ptn1,ptn0, "same");pcc9_n[i]=np.correlate(ptn9,ptn1, "same");pauto_n[i]=np.correlate(ptn0,ptn0, "same")
    #t1=ut.find_nearest(cc0_s[i][650:700],np.max(cc0_s[i][650:750])*0.4)[0];t2=ut.find_nearest(cc0_s[i][700:750][::-1],np.max(cc0_s[i][650:750])*0.4)[0]
    t1=ut.find_nearest(cc0_s[i][750:800],np.max(cc0_s[i][750:850])*0.5)[0];t2=ut.find_nearest(cc0_s[i][800:850][::-1],np.max(cc0_s[i][750:850])*0.5)[0]
    #t1=ut.find_nearest(cc0_s[i][135:150],np.max(cc0_s[i][135:165])*0.5)[0];t2=ut.find_nearest(cc0_s[i][150:165][::-1],np.max(cc0_s[i][135:165])*0.5)[0]
    #t1=ut.find_nearest(cc0_s[i][80:100],np.max(cc0_s[i][80:120])*0.5)[0];t2=ut.find_nearest(cc0_s[i][100:120][::-1],np.max(cc0_s[i][80:120])*0.5)[0]
    pt1=ut.find_nearest(pcc0_s[i][0:48],np.max(pcc0_s[i])*0.5)[0];pt2=ut.find_nearest(pcc0_s[i][48:96][::-1],np.max(pcc0_s[i])*0.5)[0]
    delt[i]=-1*(t2-t1+1)*0.05;pdelt[i]=period[pt1]-period[pt2];cc0_s[i]=cc0_s[i]/np.max(cc0_s[i])
cc0_s=np.array(cc0_s);cc9_s=np.array(cc9_s);auto_s=np.array(auto_s);delt=np.array(delt);pdelt=np.array(pdelt)
cc0_n=np.array(cc0_n);cc9_n=np.array(cc9_n);auto_n=np.array(auto_n)
pcc0_s=np.array(pcc0_s);pcc9_s=np.array(pcc9_s);pauto_s=np.array(pauto_s)
pcc0_n=np.array(pcc0_n);pcc9_n=np.array(pcc9_n);pauto_n=np.array(pauto_n)
plt.plot(np.arange(len(cc01[0]))-80,cc00);plt.plot(np.arange(len(cc01[0]))-80,cc01);plt.show()
plt.plot((np.arange(len(cc0_s[0]))-800)*0.05,cc0_s[5],'-',label=str(freq[5])+' GHz');plt.plot((np.arange(len(cc0_s[0]))-800)*0.05,cc0_s[15],'-',label=str(freq[15])+' GHz');plt.plot((np.arange(len(cc0_s[0]))-800)*0.05,cc0_s[11],'-',label=str(freq[11])+' GHz');plt.legend();plt.xlim(-2,5);plt.xlabel('Lag (sec)');plt.ylabel('Cross-Correlation');plt.show()
plt.plot(np.arange(2000)*0.05,maxTbr[0]/1.e6,label=str(freq[0])+' GHz');plt.plot(np.arange(2000)*0.05,maxTbr[11]/1.e6,label=str(freq[11])+' GHz');plt.plot(np.arange(2000)*0.05,maxTbr[15]/1.e6,label=str(freq[15])+' GHz');plt.legend();plt.xlabel('Time (sec)');plt.ylabel('$T_B$ (MK)');plt.show()
plt.plot(freq[0:16],delt[0:16],'o-');plt.plot([freq[0],freq[11],freq[20]],[del03,0,del53],'s-',color='k');plt.ylabel('Lag (sec)');plt.xlabel('Frequency (GHz)');plt.show()
#--------------------------
cc0_s=np.load('/media/rohit/VLA/20160409/correlation_s_spw0.npz.npy')
do_cc=1
if(do_cc):
    #cc0_s=np.load('/media/rohit/VLA/20160409/correlation_s_spw2.npz.npy');cc0_t=[0]*32
    cc0_s=np.load('/media/rohit/VLA/20160409/correlation_s_0.npz.npy');cc0_t=[0]*32
    for i in range(32):
        cc0_t[i]=[0]*100
        for k1 in range(100):
            cc0_t[i][k1]=[0]*100
            for k2 in range(100):
                t1=ut.find_nearest(cc0_s[i][k1][k2][750:800],np.max(cc0_s[i][k1][k2][750:850])*0.5)[0];t2=ut.find_nearest(cc0_s[i][k1][k2][800:850][::-1],np.max(cc0_s[i][k1][k2][750:850])*0.5)[0]
                cc0_t[i][k1][k2]=t2-t1
    cc0_t=np.array(cc0_t)
mapp=Map(sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.*0-15*.FITS'))[860]);mapp.peek();plt.show()
dd=mapp.data[78:178,78:178];dd[np.isnan(dd)]=0;idxdd=np.where(dd<np.nanmax(dd)*0.5)
cc0_t=np.load('/media/rohit/VLA/20160409/correlation_t_spw2.npz.npy')
ccmax=[0]*32;ccmin=[0]*32;ccmean0=[0]*32
for i in range(32):
    cc=-1*(cc0_t[i]+1)*0.05;cc[idxdd]=0;dd0=mapp;ccmax[i]=np.nanmax(cc);ccmin[i]=np.nanmin(cc);ccmean0[i]=np.nanmean(cc[40:47,43:50])
    xlvla=dd0.center.Tx.value-2.0*int(cc.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(cc.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(cc.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(cc.shape[0]/2)
    plt.imshow(cc,origin=0,vmin=-1.5,vmax=1.5,extent=[xlvla,xrvla,ylvla,yrvla],interpolation='None',cmap='seismic')
    xpi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Tx.value;ypi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Ty.value
    plt.plot([xpi],[ypi],'o',color='k');plt.xlim([-825,-725]);plt.ylim([200,300])
    plt.text(-815,215,'Maximum = '+str(np.round(ccmax[i],3))+'s Minimum ='+str(np.round(ccmin[i],3))+'s',fontsize=15)
    plt.title(str(freq[i])+' GHz');plt.colorbar(label='Lag Time (s)');plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
    ax = plt.gca();ellipse = Ellipse(xy=(-810, 270), width=bmaj[i], height=bmin[i], edgecolor='k', fc='None', lw=2);ax.add_patch(ellipse)
    plt.savefig('/media/rohit/VLA/20160409/corr/spw2_'+str(i)+'.png')
    plt.close()

cc0_t=np.load('/media/rohit/VLA/20160409/correlation_t_0.npz.npy')
ccmax=[0]*32;ccmin=[0]*32
for i in range(32):
    cc=-1*(cc0_t[i]+1)*0.05;cc[idxdd]=0;dd0=mapp;ccmax[i]=np.nanmax(cc);ccmin[i]=np.nanmin(cc)
    xlvla=dd0.center.Tx.value-2.0*int(cc.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(cc.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(cc.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(cc.shape[0]/2)
    plt.imshow(cc,origin=0,vmin=-1.5,vmax=1.5,extent=[xlvla,xrvla,ylvla,yrvla],interpolation='None',cmap='seismic')
    xpi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Tx.value;ypi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Ty.value
    plt.plot([xpi],[ypi],'o',color='k');plt.xlim([-825,-725]);plt.ylim([200,300])
    plt.text(-815,215,'Maximum = '+str(np.round(ccmax[i],3))+'s Minimum ='+str(np.round(ccmin[i],3))+'s',fontsize=15)
    plt.title(str(freq[i])+' GHz');plt.colorbar(label='Lag Time (s)');plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
    ax = plt.gca();ellipse = Ellipse(xy=(-810, 270), width=bmaj[i], height=bmin[i], edgecolor='k', fc='None', lw=2);ax.add_patch(ellipse)
    plt.savefig('/media/rohit/VLA/20160409/corr/'+str(i)+'.png')
    plt.close()


cc0_qt=np.load('/media/rohit/VLA/20160409/correlation_t_q.npz.npy');ccmax=[0]*3;ccmin=[0]*3;ccmean0=[0]*3
for i in range(3):
    cc=-1*(cc0_qt[i]+1)*0.05;cc[idxdd]=0;dd0=mapp;ccmax[i]=np.nanmax(cc);ccmin[i]=np.nanmin(cc);ccmean0[i]=np.nanmean(cc[40:47,43:50])
    xlvla=dd0.center.Tx.value-2.0*int(cc.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(cc.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(cc.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(cc.shape[0]/2)
    plt.imshow(cc,origin=0,vmin=-1.5,vmax=1.5,extent=[xlvla,xrvla,ylvla,yrvla],interpolation='None',cmap='seismic')
    xpi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Tx.value;ypi=dd0.pixel_to_world(124*u.pixel,135*u.pixel).Ty.value
    plt.plot([xpi],[ypi],'o',color='k');plt.xlim([-825,-725]);plt.ylim([200,300])
    plt.text(-815,215,'Maximum = '+str(np.round(ccmax[i],3))+'s Minimum ='+str(np.round(ccmin[i],3))+'s',fontsize=15)
    plt.title(str(freq[i])+' GHz');plt.colorbar(label='Lag Time (s)');plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
    ax = plt.gca();ellipse = Ellipse(xy=(-810, 270), width=bmaj[i], height=bmin[i], edgecolor='k', fc='None', lw=2);ax.add_patch(ellipse)
    plt.savefig('/media/rohit/VLA/20160409/corr/q_'+str(i)+'.png')
    plt.close()

hh2=np.array([43.4, 40.5, 40.1, 39.2, 38.5, 37.7, 37.1, 36.3, 36. , 35.6, 35. ,34.1, 33.5, 33.1, 32.6, 32.1, 31.6, 31.1, 30.8, 30.3, 29.8, 29.1, 28.8, 28.3, 27.9, 27.5, 27.2, 27. , 26.5, 26.2, 25.8, 25.4])

plt.plot(freq,delt*0.05,'o-');plt.ylabel('Lag (sec)');plt.xlabel('Frequency (GHz)');plt.show()

plt.plot((np.arange(1601)-800)*0.05,cc0_s[11,57,46],'-',label=str(freq[11])+' GHz')
plt.plot((np.arange(1601)-800)*0.05,cc0_s[14,57,46],'-',label=str(freq[14])+' GHz')
plt.plot((np.arange(1601)-800)*0.05,cc0_s[5,57,46],'-',label=str(freq[5])+' GHz')
plt.xlabel('Lag (sec)');plt.ylabel('Cross-Correlation')
plt.legend();plt.show()

f,ax=plt.subplots(1,1)
ax.plot(freq,ccmax,'o-',label='Maximum Lag')
ax.plot(freq,ccmin,'o-',label='Minimum Lag')
ax.plot(freq,ccmean0,'o-',label='Region');plt.legend(loc=2)
ax1=ax.twiny();new_tick_locations = hh2[::4];ax1.set_xlim(ax.get_xlim())
ax.set_xticks(freq[::4]);ax1.set_xticks(freq[::4])
ax1.set_xticklabels(new_tick_locations)
ax.set_xlabel('Frequency (GHz)');ax1.set_ylabel('Coronal Height (Mm)');ax.set_ylabel('Lag time (s)');plt.show()

plt.imshow(np.log10(cc0_s[:,650:950]),aspect='auto',origin=0)
plt.contour(cc0_s[:,650:950]/np.nanmax(cc0_s[:,650:950]),levels=[0.05,0.08,0.2],colors='k')
plt.yticks(np.arange(32)[::4],freq[::4]);plt.xlabel('Time (sec)');plt.ylabel('Frequency (GHz)')
plt.xticks(np.arange(300)[::50],0.05*(np.arange(300)[::50]-150))
plt.show()

plt.imshow(pcc0_s/np.max(pcc0_s),aspect='auto',origin=0)
plt.contour(pcc0_s/np.max(pcc0_s),levels=[0.6,0.7,0.8],colors='k')
plt.yticks(np.arange(32)[::4],freq[::4]);plt.xlabel('Time (sec)');plt.ylabel('Frequency (GHz)')
plt.xticks(np.arange(97)[::16],np.round(0.05*(np.arange(97)[::16]-48),4))
plt.show()

f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
ax0.imshow(np.array(Tbmax_wave),extent=[200,1800,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=-2,vmax=2)
ax1.imshow(np.array(maxTbr)/1.e6,extent=[0,2000,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=20,vmax=80)
ax0.axhline(y=16,color='k');ax1.axhline(y=16,color='k')
ax0.set_yticks(np.arange(32)[::4]);ax0.set_yticklabels(freq[::4]);ax1.set_yticks(np.arange(32)[::4]);ax1.set_yticklabels(freq[::4])
ax1.set_xlabel('Time (sec)');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Continuum Subtracted');ax1.set_title('Maximum Tb')
ax1.set_xticks(np.arange(2000)[::100]);ax1.set_xticklabels(np.round(0.05*(np.arange(2000)[::100]),3))
plt.show()

f,ax=plt.subplots(2,2,sharey=True);ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
ax00.imshow(np.array(Tbmax_wave)[:,650:680],extent=[650,680,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=0,vmax=7)
ax00.contour(np.array(Tbmax_wave)[:,650:680],levels=[1.5],colors='k',extent=[650,680,0,32],aspect='auto',origin='lower')
ax01.imshow(np.array(Tbmax_wave)[:,700:720],extent=[700,720,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=0,vmax=6)
ax01.contour(np.array(Tbmax_wave)[:,700:720],levels=[1.5],colors='k',extent=[700,720,0,32],aspect='auto',origin='lower')
ax10.imshow(np.array(Tbmax_wave)[:,330:370],extent=[330,370,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=0,vmax=3)
ax10.contour(np.array(Tbmax_wave)[:,330:370],levels=[0.7],colors='k',extent=[330,370,0,32],aspect='auto',origin='lower')
ax11.imshow(np.array(Tbmax_wave)[:,570:600],extent=[570,600,0,32],aspect='auto',origin=0,interpolation='None',cmap='coolwarm',vmin=0,vmax=5)
ax11.contour(np.array(Tbmax_wave)[:,570:600],levels=[1.5],colors='k',extent=[570,600,0,32],aspect='auto',origin='lower')
ax00.set_yticks(np.arange(32)[::4]);ax00.set_yticklabels(freq[::4]);ax01.set_yticks(np.arange(32)[::4]);ax01.set_yticklabels(freq[::4])
ax10.set_yticks(np.arange(32)[::4]);ax10.set_yticklabels(freq[::4]);ax11.set_yticks(np.arange(32)[::4]);ax11.set_yticklabels(freq[::4])
ax10.set_xticks(np.arange(330,370,4));ax10.set_xticklabels(np.round(0.05*(np.arange(330,370,4)),3))
ax11.set_xticks(np.arange(570,600,4));ax11.set_xticklabels(np.round(0.05*(np.arange(570,600,4)),3))
ax01.set_xticks(np.arange(700,720,4));ax01.set_xticklabels(np.round(0.05*(np.arange(700,720,4)),3))
ax00.set_xticks(np.arange(650,680,4));ax00.set_xticklabels(np.round(0.05*(np.arange(650,680,4)),3))
ax10.set_xlabel('Time (sec)');ax00.set_ylabel('Frequency (GHz)');ax11.set_xlabel('Time (sec)')
plt.show()

f,ax=plt.subplots(1,1,figsize=(18,6))
ax.imshow(np.array(Tbmax_wave)[:,550:850],extent=[750,1050,0,32],origin=0,aspect='auto',cmap='coolwarm',interpolation='None',vmin=-5,vmax=5)
#ax.contour(np.array(Tbmax_wave)[:,550:850]/np.array(Tbmax_wave)[:,550:850].max(),levels=[0.2,0.4,0.6],extent=[750,1050,0,32],origin='lower',colors='k')
#ax.imshow(np.array(ycrmax)[:,550:850],extent=[750,1050,0,32],origin=0,aspect='auto',cmap='coolwarm',interpolation='None',vmin=230,vmax=270)
ax.arrow(780, 31, 0, -5, head_width=5, head_length=1, fc='g', ec='g')
ax.arrow(830, 1.5, 20, 0, head_width=1, head_length=5, fc='b', ec='b')
ax.arrow(930, 20, -12, 0, head_width=1, head_length=5, fc='b', ec='b')
ax.arrow(920, 8, 20, 0, head_width=1, head_length=5, fc='g', ec='g')
ax.axvline(x=858,color='k',linestyle='--')
ax.set_xticks(np.linspace(750,1050,300)[::30]);ax.set_xticklabels(np.round(0.05*(np.arange(300)[::30]),3))
ax.set_yticks(np.arange(32)[::4]);ax.set_yticklabels(freq[::4]);ax.set_yticks(np.arange(32)[::4]);ax.set_yticklabels(freq[::4])
ax.set_ylabel('Frequency (GHz)');ax.set_xlabel('Time (s)')
ax1=ax.twinx();ax1.set_yticks(np.arange(32)[::4]);ax1.set_yticklabels(hh2[::4]);ax1.set_ylabel('Coronal Height (Mm)')
plt.show()

f,ax=plt.subplots(1,1,figsize=(18,6))
ax.imshow(np.array(Tbmax_wave),extent=[0,2000,0,32],origin=0,aspect='auto',cmap='coolwarm',interpolation='None',vmin=-1,vmax=1)
#ax.set_xticks(np.linspace(750,1050,300)[::30]);ax.set_xticklabels(np.round(0.05*(np.arange(300)[::30]),3))
ax.set_yticks(np.arange(32)[::4]);ax.set_yticklabels(freq[::4]);ax.set_yticks(np.arange(32)[::4]);ax.set_yticklabels(freq[::4])
ax.set_ylabel('Frequency (GHz)');ax.set_xlabel('Time (s)')
ax1=ax.twinx();ax1.set_yticks(np.arange(32)[::4]);ax1.set_yticklabels(hh2[::4]);ax1.set_ylabel('Coronal Height (Mm)')
plt.show()



f,ax=plt.subplots(1,1,sharex=True)
#ax.plot(np.array(qTbmax_wave0),label=str(freq[0])+' GHz')
#ax.plot(np.array(qTbmax_wave3),label=str(freq[10])+' GHz')
#ax.plot(np.array(qTbmax_wave5),label=str(freq[20])+' GHz')
ax.plot(np.array(qsmaxTbr_0),label=str(freq[0])+' GHz')
ax.plot(np.array(qsmaxTbr_3),label=str(freq[10])+' GHz')
ax.plot(np.array(qsmaxTbr_5),label=str(freq[20])+' GHz')
ax.set_xlabel('Time (sec)');ax0.set_ylabel('Frequency (GHz)')
ax1.set_xticks(np.arange(4200)[::100]);ax1.set_xticklabels(np.round(0.05*(np.arange(4200)[::100]),3))
plt.show()

plt.plot((np.arange(1601)-800)*0.05,cc0_s[0],'o-',label=str(freq[0])+' GHz')
plt.plot((np.arange(1601)-800)*0.05,cc0_s[5],'o-',label=str(freq[5])+' GHz')
plt.plot((np.arange(1601)-800)*0.05,cc0_s[10],'o-',label=str(freq[10])+' GHz')
plt.plot((np.arange(1601)-800)*0.05,cc0_s[14],'o-',label=str(freq[14])+' GHz')
plt.legend();plt.ylabel('Cross Correlation');plt.xlabel('Time (sec)');plt.show()

f,ax=plt.subplots(3,1);ax0=ax[0];ax1=ax[1];ax3=ax[2]
im0=ax0.imshow(glbl_npower/np.nanmax(glbl_npower),aspect='auto',origin='lower',cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
ax0.set_xscale('log');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('North Source')
im1=ax1.imshow(glbl_spower/np.nanmax(glbl_npower),aspect='auto',origin='lower',cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
#ax0.plot(periodqs,0.9+glbl_powerqs/np.max(glbl_powerqs)*0.5,'-',color='k',linewidth=4);ax0.set_xlim(period[0],period[-1]);ax0.set_ylim(freq[0],freq[-1])
ax1.set_xscale('log');ax1.set_xlabel('Period (s)');ax1.set_ylabel('Frequency (GHz)');ax1.set_title('South Source')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im1, cax=cax, orientation='vertical')
im3=ax3.imshow(tb0,aspect='auto',origin='lower',cmap='jet',interpolation='None',extent=[0,100,freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax3.set_xlabel('Time (s)');ax3.set_ylabel('Frequency (GHz)');ax3.set_title('Maximum $T_B$')
divider = make_axes_locatable(ax3);cax = divider.append_axes('right', size='5%', pad=0.05);cax.set_title('(K)');f.colorbar(im3, cax=cax, orientation='vertical')
plt.show()

f,ax0=plt.subplots(1,1)
im0=ax0.imshow(glbl_maxpower/np.nanmax(glbl_maxpower),aspect='auto',origin='lower',cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
ax0.set_xscale('log');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Maximum $T_B$')
ax0.axhline(y=1.48,color='k')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
plt.show()

f,ax0=plt.subplots(1,1)
im0=ax0.imshow(maxTbr/1.e6,aspect='auto',origin='lower',cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
ax0.set_xscale('log');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Maximum $T_B$')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
plt.show()

fs = 10.,noisegen = pyplnoise.RedNoise(fs, 1e-3, fs/2.)

from scipy import fft,ifft
whitenoise = np.random.uniform(0,1,2000)
timewavewn=np.arange(2000)*0.05;mother=wavelet.Morlet(6)
fouriertransformed = np.fft.fftshift(fft.fft(whitenoise))
pinktransformed = np.reciprocal(fouriertransformed)
pinknoise = ifft(np.fft.ifftshift(pinktransformed)).real
std=np.std(whitenoise);data_stdwn=whitenoise/std
dt=0.05;s0 = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdwn)
wavewn, scaleswn, freqswn, coiwn, fftwn, fftfreqswn = wavelet.cwt(data_stdwn, dt, dj, s0, J,mother)
iwavewn = wavelet.icwt(wavewn, scaleswn, dt, dj, mother) * std
powerwn = (np.abs(wavewn)) ** 2;periodwn = 1 / freqswn;powerwn /= scaleswn[:, None]
glbl_powerwn=np.mean(powerwn,axis=1)
std=np.std(pinknoise);data_stdpn=pinknoise/std
dt=0.05;s_ = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdpn)
wavepn, scalespn, freqspn, coipn, fftpn, fftfreqspn = wavelet.cwt(data_stdpn, dt, dj, s_, J,mother)
iwavepn = wavelet.icwt(wavepn, scalespn, dt, dj, mother) * std
powerpn = (np.abs(wavepn)) ** 2;periodpn = 1 / freqspn;powerpn /= scalespn[:, None]
glbl_powerpn=np.mean(powerpn,axis=1)
plot_wavelet(pinknoise,timewavewn,dt,periodpn,powerpn,coipn,'/home/i4ds1807205/Dropbox/20160409/wavelets/blob/pink_'+str(i)+'.png','')
plot_wavelet(whitenoise,timewavewn,dt,periodwn,powerwn,coiwn,'/home/i4ds1807205/Dropbox/20160409/wavelets/blob/white_'+str(i)+'.png','')
plot_wavelet(qTbmax_wave0,tq,dt,period,s_power0,s_coi10,'/home/i4ds1807205/Dropbox/20160409/wavelets/qt_0.png',freq[0])
plot_wavelet(qTbmax_wave3,tq,dt,period,s_power3,s_coi13,'/home/i4ds1807205/Dropbox/20160409/wavelets/qt_3.png',freq[10])
plot_wavelet(qTbmax_wave5,tq,dt,period,s_power5,s_coi15,'/home/i4ds1807205/Dropbox/20160409/wavelets/qt_5.png',freq[20])

f,ax0=plt.subplots(1,1)
ax0.plot(period,glbl_maxpower[0],'-',label=str(freq[0])+' GHz',color='orange')
ax0.plot(period,glbl_maxpower[4],'-',label=str(freq[4])+' GHz',color='brown')
ax0.plot(period,glbl_maxpower[10],'-',label=str(freq[10])+' GHz',color='cyan')
ax0.plot(period,glbl_maxpower[14],'-',label=str(freq[14])+' GHz',color='y')
ax0.plot(period,s_power0.mean(axis=1),'o-',label=str(freq[0])+' GHz (Pre-burst)',color='r')
ax0.plot(period,s_power3.mean(axis=1),'o-',label=str(freq[10])+' GHz (Pre-burst)',color='g')
ax0.plot(period,s_power5.mean(axis=1),'o-',label=str(freq[20])+' GHz (Pre-burst)',color='b')
ax0.plot(periodwn,glbl_powerwn,'-',label='White Noise',color='gray')
#ax0.plot(periodpn,glbl_powerpn,'-',label='Pink Noise',color='m')
ax0.legend(loc=1, prop={'size': 11});ax0.set_xscale('log');ax0.set_xlabel('Period (sec)');ax0.set_ylabel('Wavelet Power')
plt.show()

sys.exit()
from scipy import fft,ifft
whitenoise = np.random.uniform(0,1,len(Tbrs_wave))
timewavewn=np.arange(len(Tbrs_wave))*0.05;mother=wavelet.Morlet(6)
fouriertransformed = np.fft.fftshift(fft(whitenoise))
pinktransformed = np.reciprocal(fouriertransformed)
pinknoise = ifft(np.fft.ifftshift(pinktransformed)).real
std=np.std(whitenoise);data_stdwn=whitenoise/std
dt=0.05;s0 = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdwn)
wavewn, scaleswn, freqswn, coiwn, fftwn, fftfreqswn = wavelet.cwt(data_stdwn, dt, dj, s0, J,mother)
iwavewn = wavelet.icwt(wavewn, scaleswn, dt, dj, mother) * std
powerwn = (np.abs(wavewn)) ** 2;periodwn = 1 / freqswn;powerwn /= scaleswn[:, None]
glbl_powerwn=np.mean(powerwn,axis=1)
std=np.std(pinknoise);data_stdpn=pinknoise/std
dt=0.05;s_ = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdpn)
wavepn, scalespn, freqspn, coipn, fftpn, fftfreqspn = wavelet.cwt(data_stdpn, dt, dj, s_, J,mother)
iwavepn = wavelet.icwt(wavepn, scalespn, dt, dj, mother) * std
powerpn = (np.abs(wavepn)) ** 2;periodpn = 1 / freqspn;powerpn /= scalespn[:, None]
glbl_powerpn=np.mean(powerpn,axis=1)

plt.plot(period,s_power.mean(axis=1),'o-',label='Loop (A, South)')
plt.plot(period,n_power.mean(axis=1),'o-',label='Loop (B, North)')
plt.plot(periodwn,glbl_powerwn,'o-',label='White Noise')
plt.plot(periodpn,glbl_powerpn,'o-',label='Pink Noise (Real)')
plt.xlabel('Time (s)');plt.ylabel('Power');plt.legend();plt.show()

################# Analysis #################
#Tb=np.array(allmaps['vla']['datavla']).reshape(32,119,150,200)
#Tb=np.array(allmaps['vla']['datavla']).reshape(1,2399,150,201)
#Tbmax_ds=Tb.max(axis=(2,3))
#Tb_mean_r1=Tb[:,:,30:80,115:145].mean(axis=(2,3))
#m=allmaps['vla']['timevla'].shape[0]
#tidx94=[0]*m;tidx131=[0]*m
#tidx335=[0]*m;tidx1600=[0]*m;tidx1700=[0]*m
#for i in range(allmaps['vla']['timevla'].shape[0]):
#    tidx94[i]=ut.find_predecessor(allmaps['aia94']['time94'],allmaps['vla']['timevla'][i])[0]
#    tidx131[i]=ut.find_predecessor(allmaps['aia131']['time131'],allmaps['vla']['timevla'][i])[0]
#    tidx335[i]=ut.find_predecessor(allmaps['aia335']['time335'],allmaps['vla']['timevla'][i])[0]
#    tidx1600[i]=ut.find_predecessor(allmaps['aia1600']['time1600'],allmaps['vla']['timevla'][i])[0]
#    tidx1700[i]=ut.find_predecessor(allmaps['aia1700']['time1700'],allmaps['vla']['timevla'][i])[0]

#vlafreqidx=np.array(list(np.arange(32))*119).reshape(119,32).swapaxes(0,1).flatten()
vlafreqidx=np.array(list(np.arange(1))*2399).reshape(2399,1).swapaxes(0,1).flatten();timediff=[0]*32;power2=[0]*32;power1=[0]*32;powern=[0]*32;powers=[0]*32
xcrmaxn=xcrmax*1.0;xcrmaxs=xcrmax*1.0;ycrmaxn=ycrmax*1.0;ycrmaxs=ycrmax*1.0;maxTbrn=maxTbr*1.0;maxTbrs=maxTbr*1.0;coi1=[0]*32;coi2=[0]*32;coin=[0]*32;cois=[0]*32
ycrmaxn[ycrmaxn<255]=0;ycrmaxs[ycrmaxs>255]=0;xcrmaxn[ycrmaxn<255]=0;xcrmaxs[ycrmaxs>255]=0;maxTbrn[ycrmax<255]=0;maxTbrs[ycrmax>255]=0;N=500
glbl_power1=[0]*32;glbl_power2=[0]*32;glbl_powermaxn=[0]*32;glbl_powermaxs=[0]*32;scale_avg=[0]*32;max_period=[0]*32;power=[0]*32;iwave=[0]*32;Tbr1_diff=[0]*32;Tbr1_sm=[0]*32
for i in range(32):
    datawave_=(Tbr_r2[i]-Tbr_r2[i][0])/1.e6
    Tbr1_sm=np.convolve(datawave_, np.ones(N)/N, mode='valid')
    datawave=datawave_[int(N/2):int(-1*N/2+1)]-np.convolve(datawave_, np.ones(N)/N, mode='valid');datawave=datawave/np.std(datawave);Tbr1_diff[i]=datawave
    timewave=np.arange(len(datawave))*0.05;timediff[i]=timewave#timevla[:2000]-timevla[0]
    mother=wavelet.Morlet(6)#mother=wavelet.MexicanHat()
    std=1;datawave_std=datawave/std
    dt=0.05;s0 = 2 * dt;dj = 1. / 12 # Lowest scale s0
    J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
    alpha, _, _ = wavelet.ar1(datawave)
    wave, scales, freqs, coi1[i], fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
    iwave[i] = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
    power2[i] = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power[i] /= scales[:, None]
    #widths=np.logspace(1,3,num=400)*0.05#cwtmatr, freqs = pywt.cwt(wdata, widths, mother)#power = (np.abs(cwtmatr)) ** 2
    t=np.arange(len(datawave))*0.05
    glbl_power2[i] = power2[i].mean(axis=1)
    scale_avg[i]=power[i].mean(axis=0)
    max_period[i]=period[np.where(glbl_power1[i]==np.max(glbl_power1[i]))][0]
    ############################
    datawave_=(Tbr_r1[i]-Tbr_r1[i][0])/1.e6
    Tbr1_sm=np.convolve(datawave_, np.ones(N)/N, mode='valid')
    datawave=datawave_[int(N/2):int(-1*N/2+1)]-np.convolve(datawave_, np.ones(N)/N, mode='valid');datawave=datawave/np.std(datawave);Tbr1_diff[i]=datawave
    timewave=np.arange(len(datawave))*0.05#timevla[:2000]-timevla[0]
    mother=wavelet.Morlet(6)#mother=wavelet.MexicanHat()
    std=1;datawave_std=datawave/std
    dt=0.05;s0 = 2 * dt;dj = 1. / 12 # Lowest scale s0
    J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
    alpha, _, _ = wavelet.ar1(datawave)
    wave, scales, freqs, coi2[i], fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
    iwave[i] = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
    power1[i] = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power[i] /= scales[:, None]
    #widths=np.logspace(1,3,num=400)*0.05#cwtmatr, freqs = pywt.cwt(wdata, widths, mother)#power = (np.abs(cwtmatr)) ** 2
    glbl_power1[i] = power1[i].mean(axis=1)
    ############################
    if(maxTbrn[i].max()!=0):
        datawave_=(maxTbrn[i])/1.e6;datawave=datawave_/np.std(datawave_)
        timewave=np.arange(len(datawave))*0.05#timevla[:2000]-timevla[0]
        mother=wavelet.Morlet(6)#mother=wavelet.MexicanHat()
        std=1;datawave_std=datawave/std
        dt=0.05;s0 = 2 * dt;dj = 1. / 12; J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
        alpha, _, _ = wavelet.ar1(datawave)
        wave, scales, freqs, coin[i], fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
        iwave[i] = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
        powern[i] = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power[i] /= scales[:, None]
        #widths=np.logspace(1,3,num=400)*0.05#cwtmatr, freqs = pywt.cwt(wdata, widths, mother)#power = (np.abs(cwtmatr)) ** 2
        glbl_powermaxn[i] = powern[i].mean(axis=1)
    else:
        glbl_powermaxn[i]=np.zeros(97)
    ############################
    if(maxTbrs[i].max()!=0):
        datawave_=(maxTbrs[i])/1.e6;datawave=datawave_/np.std(datawave_)
        timewave=np.arange(len(datawave))*0.05#timevla[:2000]-timevla[0]
        mother=wavelet.Morlet(6)#mother=wavelet.MexicanHat()
        std=1;datawave_std=datawave/std
        dt=0.05;s0 = 2 * dt;dj = 1. / 12; J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
        alpha, _, _ = wavelet.ar1(datawave)
        wave, scales, freqs, cois[i], fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
        iwave[i] = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
        powers[i] = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power[i] /= scales[:, None]
        #widths=np.logspace(1,3,num=400)*0.05#cwtmatr, freqs = pywt.cwt(wdata, widths, mother)#power = (np.abs(cwtmatr)) ** 2
        glbl_powermaxs[i] = powers[i].mean(axis=1)
    else:
        glbl_powermaxs[i]=np.zeros(97)
glbl_powermaxn=np.array(glbl_powermaxn);glbl_powermaxs=np.array(glbl_powermaxs)
glbl_power1=np.array(glbl_power1);glbl_power2=np.array(glbl_power2);scale_avg=np.array(scale_avg);max_period=np.array(max_period)




N=500
qsTbr_r1_1=np.hstack((qsTbr_r1[:2131],qsTbr_r1[2170:4531],qsTbr_r1[4420:4570]))
qsTbr_r1_sm=np.convolve(qsTbr_r1_1, np.ones(N)/N, mode='valid')
qsTbr_r1_diff=qsTbr_r1_1[int(N/2):int(-1*N/2+1)]-qsTbr_r1_sm
f,ax=plt.subplots(2,1)
ax[0].plot(qsTbr_r1_1[int(N/2):int(-N/2)],'o-')
ax[0].plot(qsTbr_r1_sm,'-')
ax[1].plot(qsTbr_r1_diff,'o-')
plt.show()

dataqs=qsTbr_r1_diff/1.e6
timewaveqs=np.arange(len(dataqs))*0.05#;mother=wavelet.MexicanHat()
mother=wavelet.Morlet(6)
std=np.std(dataqs);data_stdqs=dataqs/std
dt=0.05;s0 = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(dataqs)
waveqs, scalesqs, freqsqs, coiqs, fftqs, fftfreqsqs = wavelet.cwt(data_stdqs, dt, dj, s0, J,mother)
iwaveqs = wavelet.icwt(waveqs, scalesqs, dt, dj, mother) * std
powerqs = (np.abs(waveqs)) ** 2;periodqs = 1 / freqsqs;powerqs /= scalesqs[:, None]
glbl_powerqs=np.mean(powerqs,axis=1)

i=24;plt.ioff();label='TB';units='K'
figprops = dict(figsize=(20, 20), dpi=72)
fig = plt.figure(**figprops)
ax = plt.axes([0.1, 0.68, 0.67, 0.2])
#ax.plot(t, iwave[i], '-', linewidth=1, color=[0.5, 0.5, 0.5])
t=timediff[i]
ax.plot(timediff[i],Tbr1_diff[i], 'k', linewidth=1.5)
ax.grid(True)
ax.set_ylabel(r'{} [{}]'.format(label, units))
bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax);levels = np.linspace(1.e-5,450,90)
#bx.contourf(t, periodds, np.log2(powerds), np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
bx.contourf(t, period, np.log2(power1[i]),extend='both', cmap=plt.cm.jet)
extent = [t.min(), t.max(), 0, max(period)]
bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([coi1[i], [1e-9], period[-1:],period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[i])+' GHz({})'.format(label, 'MORLET'))
bx.set_ylabel('Period (sec)');bx.set_xlabel('Time (sec)');bx.set_yscale('log')
bx.grid(True)
cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
cx.plot(glbl_power1[i], period, 'k-', linewidth=1.5)
cx.set_title('c) Global Wavelet Spectrum')
cx.set_xlabel(r'Power [({})^2]'.format(units))
#cx.set_ylim(([periodds.min(), periodds.max()]))
cx.grid(True)
cx.set_xscale('log')
plt.show()
############

from scipy import fft,ifft

whitenoise = np.random.uniform(0,1,len(dataqs))
timewavewn=np.arange(len(dataqs))*0.05#;mother=wavelet.MexicanHat()
fouriertransformed = np.fft.fftshift(fft(whitenoise))
pinktransformed = np.reciprocal(fouriertransformed)
pinknoise = ifft(np.fft.ifftshift(pinktransformed)).real
std=np.std(whitenoise);data_stdwn=whitenoise/std
dt=0.05;s0 = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdwn)
wavewn, scaleswn, freqswn, coiwn, fftwn, fftfreqswn = wavelet.cwt(data_stdwn, dt, dj, s0, J,mother)
iwavewn = wavelet.icwt(wavewn, scaleswn, dt, dj, mother) * std
powerwn = (np.abs(wavewn)) ** 2;periodwn = 1 / freqswn;powerwn /= scaleswn[:, None]
glbl_powerwn=np.mean(powerwn,axis=1)
std=np.std(pinknoise);data_stdpn=pinknoise/std
dt=0.05;s0 = 2 * dt;dj = 1. / 12;J = 8. / dj;alpha, _, _ = wavelet.ar1(data_stdpn)
wavepn, scalespn, freqspn, coipn, fftpn, fftfreqspn = wavelet.cwt(data_stdpn, dt, dj, s0, J,mother)
iwavepn = wavelet.icwt(wavepn, scalespn, dt, dj, mother) * std
powerpn = (np.abs(wavepn)) ** 2;periodpn = 1 / freqspn;powerpn /= scalespn[:, None]
glbl_powerpn=np.mean(powerpn,axis=1)


f,ax=plt.subplots(1,1)
ax.plot(freqswn,glbl_powerwn,'o-',label='White Noise')
ax.plot(freqsqs,glbl_powerqs,'o-',label='Pre-burst')
ax.plot(freqspn,glbl_powerpn,'o-',label='Pink Noise')
ax.plot(freqs,glbl_power1[0],'o-',label='Burst')
ax.legend();ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Power')
plt.show()


sys.exit()

#######3####

ds=ds_['spec'].mean(axis=(0,1))
qds=ds[:,900:2400];qds1=qds[0:20].mean(axis=0)
for i in range(12):
    qds1[111+i*120]=qds1[110+i*120]

t=np.arange(len(dataqs))*0.05
mother=wavelet.Morlet(6)
std=1;datawave_std=dataqs/std;dt=0.05;s0 = 8 * dt;dj = 1. / 12;J = 8. / dj
alpha, _, _ = wavelet.ar1(dataqs);waveds, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dataqs, dt, dj, s0, J,mother)
iwaveds= wavelet.icwt(waveds, scales, dt, dj, mother) * std
powerds = (np.abs(waveds)) ** 2;fft_powerds = np.abs(fft) ** 2;periodds = 1 / freqs;powerds /= scales[:, None]
glbl_powerds = powerds.mean(axis=1)
scale_avgds = powerds.mean(axis=0)
plt.ioff()
figprops = dict(figsize=(20, 20), dpi=72)
fig = plt.figure(**figprops)
ax = plt.axes([0.1, 0.75, 0.65, 0.2])
ax.plot(t, iwaveds, '-', linewidth=1, color=[0.5, 0.5, 0.5])
ax.plot(t, dataqs, 'k', linewidth=1.5)
ax.grid(True)
ax.set_ylabel(r'{} [{}]'.format(label, units))
bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax);levels = np.linspace(1.e-5,450,90)
#bx.contourf(t, periodds, np.log2(powerds), np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
bx.contourf(t, periodds, np.log2(powerds),extend='both', cmap=plt.cm.YlOrRd)
extent = [t.min(), t.max(), 0, max(periodds)]
bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([coi, [1e-9], periodds[-1:],periodds[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[0])+' GHz({})'.format(label, 'MORLET'))
bx.set_ylabel('Period (sec)');bx.set_yscale('log')
bx.grid(True)
cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
cx.plot(glbl_powerds, periodds, 'k-', linewidth=1.5)
cx.set_title('c) Global Wavelet Spectrum')
cx.set_xlabel(r'Power [({})^2]'.format(units))
#cx.set_ylim(([periodds.min(), periodds.max()]))
cx.grid(True)
cx.set_xscale('log')
# Fourth sub-plot, the scale averaged wavelet spectrum.
dx = plt.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
dx.plot(t, scale_avgds, 'k-', linewidth=1.5)
dx.set_title('d) Scale-averaged power')
dx.set_xlabel('Time (sec)')
dx.set_ylabel(r'Reconstructed $T_B$ [{}]'.format(units))
dx.grid(True)
plt.show()


############

plot_wavelet=1
if(plot_wavelet):
    k=0
    for k in range(32):
        #datawave=(Tbi_r1[k]-Tbi_r1[k][0])/1.e6
        datawave=Tbr1_diff[k]/1.e6
        plt.ioff();figprops = dict(figsize=(20, 20), dpi=72)
        fig = plt.figure(**figprops)
        # First sub-plot, the original time series anomaly and inverse wavelet transform.
        ax = plt.axes([0.1, 0.75, 0.65, 0.2])
        ax.plot(t, iwave[k], '-', linewidth=1, color=[0.5, 0.5, 0.5])
        ax.plot(t, datawave, 'k', linewidth=1.5)
        ax.set_title('a) {}'.format(title))
        ax.grid(True)
        ax.set_ylabel(r'{} [{}]'.format(label, units))
        # Second sub-plot, the normalized wavelet power spectrum and significance level contour lines and cone of influece hatched area. Note that period scale is logarithmic.
        bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax)
        levels = np.linspace(1.e-5,450,90)#[0.015625,0.03125,0.0625, 0.125, 0.25, 0.5, 1,2 ]
        bx.contourf(t, period, np.log2(power[k]), np.log2(levels),
                            extend='both', cmap=plt.cm.YlOrRd)
        extent = [t.min(), t.max(), 0, max(period)]
        #bx.contour(t, numpy.log2(period), sig95, [-99, 1], colors='k', linewidths=2,extent=extent)
        bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([coi, [1e-9], period[-1:],period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
        bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[k])+' GHz({})'.format(label, 'MORLET'))
        bx.set_ylabel('Period (sec)');bx.set_yscale('log')
        bx.grid(True)
        #Yticks = 2 ** np.arange(np.ceil(np.log2(period[i].min())),
        #bx.set_yticks(np.log2(Yticks))
        #bx.set_yticklabels(Yticks)
        # Third sub-plot, the global wavelet and Fourier power spectra and theoretical noise spectra. Note that period scale is logarithmic.
        cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
        cx.plot(glbl_power[k], period, 'k-', linewidth=1.5)
        cx.set_title('c) Global Wavelet Spectrum')
        cx.set_xlabel(r'Power [({})^2]'.format(units))
        #cx.set_xlim([4.e-4, glbl_power.max() + diff_S_std[i]])
        cx.set_ylim(([period.min(), period.max()]))
        cx.grid(True)
        cx.set_xscale('log')
        #cx.set_yticks(np.log2(Yticks))
        #cx.set_yticklabels(Yticks)
        #plt.setp(cx.get_yticklabels(), visible=False)
        # Fourth sub-plot, the scale averaged wavelet spectrum.
        dx = plt.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
        dx.plot(t, scale_avg[k], 'k-', linewidth=1.5)
        dx.set_title('d) Scale-averaged power')
        dx.set_xlabel('Time (sec)')
        dx.set_ylabel(r'Reconstructed $T_B$ [{}]'.format(units))
        dx.grid(True)
        ax.set_xlim([t.min(), t.max()])
        plt.savefig('/home/i4ds1807205/Dropbox/20160409/wavelet_diff_'+str(k)+'.png')
        plt.close()

########### 1 -> North; 2-> South #####################

glbl_power=np.array(glbl_power)
plot_powerspec_all=1
if(plot_powerspec_all):
    for i in [0,5,10,15,20,25]:
        plt.plot(period,glbl_power1[i],'o-',label=str(freq[i])+' GHz')
    plt.plot(periodqs,powerqs.mean(axis=1),'*-',color='k',label='Before Burst')
    plt.xlabel('Time (s)')
    plt.ylabel('Global Power Spectrum')
    plt.legend()
    plt.show()

plot_time_all=1
if(plot_time_all):
    for i in [0,5,10,15,20,25]:
        plt.plot(np.arange(len(Tbi_r1[i]))*0.05,Tbi_r1[i]/1.e6,'-',label=str(freq[i])+' GHz')
    plt.xlabel('Time (s)')
    plt.ylabel('T$_B$ (MK)')
    plt.legend()
    plt.show()

plot_powerspec_all=1
if(plot_powerspec_all):
    for i in range(32):
        plt.plot(period,glbl_power[i],'o-',label=str(freq[i])+' GHz')
        plt.xlabel('Time (s)')
        plt.ylabel('Global Power Spectrum')
        plt.xscale('log');plt.yscale('linear');plt.legend(loc=4);plt.xlim(0.5,12);plt.ylim(1.e-2,200)
        plt.savefig('/home/i4ds1807205/Dropbox/20160409/wavelet_power_'+str(i)+'.png')
        plt.close()

plt.plot(np.arange(len(Tbi_r1[0]))*0.05,Tbi_r1[0]/Tbi_r1[0].max(),'-',label=str(freq[0])+' GHz')
plt.plot(np.arange(len(Tbi_r1[15]))*0.05,Tbi_r1[15]/Tbi_r1[15].max(),'-',linewidth=2,label=str(freq[15])+' GHz')
plt.plot(np.arange(len(Tbi_r1[20]))*0.05,Tbi_r1[20]/Tbi_r1[20].max(),'-',linewidth=2,label=str(freq[20])+' GHz')
plt.legend()
plt.show()

plt.plot(np.arange(len(Tbr1_diff[0]))*0.05,Tbr1_diff[0]/Tbr1_diff[0].max(),'-',label=str(freq[0])+' GHz')
plt.plot(np.arange(len(Tbr1_diff[15]))*0.05,Tbr1_diff[15]/Tbr1_diff[15],'-',linewidth=2,label=str(freq[15])+' GHz')
plt.plot(np.arange(len(Tbr1_diff[20]))*0.05,Tbr1_diff[20]/Tbr1_diff[20],'-',linewidth=2,label=str(freq[20])+' GHz')
plt.legend()
plt.show()

f,ax=plt.subplots(3,1);ax0=ax[0];ax1=ax[1];ax3=ax[2]
im0=ax0.imshow(glbl_power1/np.max(glbl_power1),aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
ax0.set_xscale('log');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Northern Source')
im1=ax1.imshow(glbl_power2/np.max(glbl_power2),aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]])
ax1.plot(periodqs,1+glbl_powerqs/np.max(glbl_powerqs)*0.5,'-',color='white',linewidth=4);ax1.set_xlim(period[0],period[-1]);ax1.set_ylim(freq[0],freq[-1])
ax1.set_xscale('log');ax1.set_xlabel('Period (s)');ax1.set_ylabel('Freqeuncy (GHz)');ax1.set_title('Southern Source')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im1, cax=cax, orientation='vertical')
im3=ax3.imshow(Tbr_r2,aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[0,100,freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax3.set_xlabel('Time (s)');ax3.set_ylabel('Frequency (GHz)');ax3.set_title('Southern Source')
divider = make_axes_locatable(ax3);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im3, cax=cax, orientation='vertical')
plt.show()

f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
im0=ax0.imshow(tb0,aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]],vmin=8.e6,vmax=6.e7)
ax0.set_xscale('linear');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Northern Source')
im1=ax1.imshow(tb1,aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]],vmin=4.e6,vmax=4.e7)
ax1.set_xscale('linear');ax1.set_xlabel('Time (s)');ax1.set_ylabel('Freqeuncy (GHz)');ax1.set_title('Southern Source')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im1, cax=cax, orientation='vertical')
plt.show()

f,ax=plt.subplots(3,1);ax0=ax[0];ax1=ax[1];ax3=ax[2]
im0=ax0.imshow(glbl_powermaxn/np.max(glbl_powermaxn),aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax0.set_xscale('log');ax0.set_ylabel('Frequency (GHz)');ax0.set_title('Northern Source')
im1=ax1.imshow(glbl_powermaxs/np.max(glbl_powermaxs),aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[period[0],period[-1],freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax0.plot(periodqs,0.9+glbl_powerqs/np.max(glbl_powerqs)*0.5,'-',color='k',linewidth=4);ax0.set_xlim(period[0],period[-1]);ax0.set_ylim(freq[0],freq[-1])
ax1.set_xscale('log');ax1.set_xlabel('Period (s)');ax1.set_ylabel('Freqeuncy (GHz)');ax1.set_title('Southern Source')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im1, cax=cax, orientation='vertical')
im3=ax3.imshow(maxTbr,aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[0,100,freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax3.set_xlabel('Time (s)');ax3.set_ylabel('Freqeuncy (GHz)');ax3.set_title('Maximum $T_B$')
divider = make_axes_locatable(ax3);cax = divider.append_axes('right', size='5%', pad=0.05);cax.set_title('(K)');f.colorbar(im3, cax=cax, orientation='vertical')
plt.show()

f,ax=plt.subplots(3,1,sharex=True);ax0=ax[0];ax1=ax[1];ax3=ax[2]
im0=ax0.imshow(xcrmax,aspect='auto',origin=0,cmap='jet',interpolation='None',vmin=-780,vmax=-740,extent=[0,100,freq[0],freq[-1]])
ax0.set_ylabel('Frequency (GHz)')
im1=ax1.imshow(ycrmax,aspect='auto',origin=0,cmap='jet',interpolation='None',vmin=230,vmax=250,extent=[0,100,freq[0],freq[-1]])
ax1.set_xlabel('Period (s)');ax1.set_ylabel('Frequency (GHz)');ax1.set_title('Y-coordinate');ax0.set_title('X-coordinate')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);cax.set_title('(arcsec)');f.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);cax.set_title('(arcsec)');f.colorbar(im1, cax=cax, orientation='vertical')
im3=ax3.imshow(maxTbr,aspect='auto',origin=0,cmap='jet',interpolation='None',extent=[0,100,freq[0],freq[-1]],norm=matplotlib.colors.LogNorm())
ax3.set_xlabel('Time (s)');ax3.set_ylabel('Frequency (GHz)');ax3.set_title('Southern Source')
divider = make_axes_locatable(ax3);cax = divider.append_axes('right', size='5%', pad=0.05);cax.set_title('(K)');f.colorbar(im3, cax=cax, orientation='vertical')
plt.show()


plt.plot(freq,max_period,'o-')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Period of max Power (s)')
plt.show()


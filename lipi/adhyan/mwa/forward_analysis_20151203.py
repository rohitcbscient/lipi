import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
import pickle
from scipy import signal
from surya.utils import main as ut
from surya.utils import model as md
import glob
import matplotlib.cm as cm


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

dir_='/media/rohit/VLA/20151203_MWA/Tb_new/'
data_100MHz=pickle.load(open(dir_+'Tb_1133149192-084-085.p','r'))

data_160MHz=pickle.load(open(dir_+'Tb_1133149192-125-126.p','r'))
data_240MHz=pickle.load(open(dir_+'Tb_1133149192-187-188.p','r'))


mwa_100MHz=np.mean(np.array(data_100MHz[0][0:577]),axis=0)
mwa_160MHz=np.mean(np.array(data_160MHz[0][0:577]),axis=0)
mwa_240MHz=np.mean(np.array(data_240MHz[0][0:577]),axis=0)

bmin_100MHz=data_100MHz[3][0]*60 # In arcmin
bmax_100MHz=data_100MHz[4][0]*60 # In arcmin
bmin_160MHz=data_160MHz[3][0]*60 # In arcmin
bmax_160MHz=data_160MHz[4][0]*60 # In arcmin
bmin_240MHz=data_240MHz[3][0]*60 # In arcmin
bmax_240MHz=data_240MHz[4][0]*60 # In arcmin

#flist=[108,120,132,145,160,179,196,217,240]
#flabel=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
flist=[108,132,145,160,179,196,217,240]
flabel=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
base=['000-008','000-009','000-010','008-009','008-010','009-010']
mwa_flux=[0]*len(base)
i=0
for b in base:
    print b
    mwa_flux[i]=[0]*len(flabel)
    j=0
    for f in flabel:
        d=pickle.load(open('/home/i4ds1807205/20151203/pickle/flux_V1_1133149192-%b'+str(f)+'_T'+str(b)+'.p','rb'))
        d_=d[17][3][0][np.isfinite(d[17][3][0])]
        mwa_flux[i][j]=np.mean(d_)
        #print mwa_flux[i][j]
        j=j+1
    i=i+1

emwa_flux=np.std(np.array(mwa_flux),axis=0)
mwa_flux=np.mean(np.array(mwa_flux),axis=0)

Tb_all=[0]*len(flabel)
bmin=[0]*len(flabel)
bmax=[0]*len(flabel)
Tb_fwd=[0]*len(flabel)
tau=[0]*len(flabel)
h=[0]*len(flabel)
hkm=[0]*len(flabel)
Tb_convolved=[0]*len(flabel)
flux_fwd=[0]*len(flabel)
flux_mwa=[0]*len(flabel)
tempobs=[0]*len(flabel)
taur=[0]*len(flabel)
taul=[0]*len(flabel)
dtaur=[0]*len(flabel)
dtaul=[0]*len(flabel)
tau_fwd=[0]*len(flabel)
tau_convolved=[0]*len(flabel)

k=0
for f in flabel:
    print flist[k],' MHz'
    data=pickle.load(open(dir_+'Tb_1133149192-'+f+'.p','r'))
    Tb_all[k]=data[0][0:577]
    bmin[k]=data[3][0]*60
    bmax[k]=data[4][0]*60
    forward=readsav('/home/i4ds1807205/20151203/20151203_'+str(flist[k])+'MHz_psimas.sav')
    rt=readsav('/home/i4ds1807205/20151203/RT_params_'+str(flist[k])+'MHz.sav')
    tempobs=rt['tempall']
    taur[k]=rt['taur'][:,:,-1]
    taul[k]=rt['taul'][:,:,-1]
    tau_fwd[k]=taul[k]+taur[k]
    #dtaur[k]=rt['dtaur']
    #dtaul[k]=rt['dtaul']
    Tb_fwd[k]=forward['stokesstruct'][0][0]
    ### CONVOLUTION
    forward_dx=forward['quantmap'][0][3]*16*60 # In arcsec
    forward_dy=forward['quantmap'][0][4]*16*60 # In arcsec; originally in Rsun
    fwdx=forward['quantmap'][0][3]*16*60 # Rsun to arcsec
    fwdy=forward['quantmap'][0][4]*16*60 
    fwhmx=int(bmin[k]*60/forward_dx)
    fwhmy=int(bmax[k]*60/forward_dy)
    beam=ut.makeGaussian(60, fwhmx, fwhmy, center=None)
    Tb_convolved[k] = signal.convolve(Tb_fwd[k],beam, mode='same')/np.sum(beam)
    tau_convolved[k] = signal.convolve(tau_fwd[k],beam, mode='same')/np.sum(beam)
    size=fwdx#*Tb_fwd[k].shape[0]
    flux_fwd[k]=ut.Tb2flux(Tb_convolved[k], size, size, flist[k]/1000.)
    flux_mwa[k]=np.mean(Tb_all[k],axis=0)*mwa_flux[k]/np.sum(Tb_all[k][0])#ut.Tb2flux(Tb_all[k][0], 50, 50, flist[k]/1000.)
    ###
    tau[k]=np.log((Tb_fwd[k].max()/np.mean(Tb_all[k],axis=(1,2)))-1)
    #print flist[k],np.max(Tb_fwd[k]),np.max(Tb_all[k]),np.std(Tb_all[k])#,np.mean(tau[k]),np.std(tau[k])
    h[k]=md.nk_freq2r(flist[k],1)[0]
    hkm[k]=(h[k]-1)*3.e5
    k=k+1

Tb_fwd=np.array(Tb_fwd)
Tb_all=np.array(Tb_all)

lev=np.linspace(0.1,0.9,10)
Tb_bimage=[0]*8
Tb_fwd_bimage=[0]*8
flux_bimage=[0]*8
fw_bimage=[0]*8
fw_bsize=[0]*8
fw_tau_bimage=[0]*8
mwa_bsize=[0]*8
mwab=[0]*8
fwdb=[0]*8
for i in range(8):
    Tb_bimage[i]=[0]*10
    flux_bimage[i]=[0]*10
    fw_bimage[i]=[0]*10
    fw_tau_bimage[i]=[0]*10
    Tb_fwd_bimage[i]=[0]*10
    fw_bsize[i]=[0]*10
    mwa_bsize[i]=[0]*10
    mwab[i]=[0]*10
    fwdb[i]=[0]*10
    for j in range(10):
        mwab[i][j]=ut.get_bimage(flux_mwa[i],lev[j])
        fwdb[i][j]=ut.get_bimage(flux_fwd[i],lev[j]) # flux_fwd is convolved
        flux_bimage[i][j]=flux_mwa[i]*ut.get_bimage(flux_mwa[i],lev[j])
        Tb_bimage[i][j]=np.mean(Tb_all[i],axis=0)*ut.get_bimage(flux_mwa[i],lev[j])
        fw_bimage[i][j]=flux_fwd[i]*ut.get_bimage(flux_fwd[i],lev[j])
        Tb_fwd_bimage[i][j]=Tb_convolved[i]*ut.get_bimage(flux_fwd[i],lev[j])
        fw_tau_bimage[i][j]=tau_convolved[i]*ut.get_bimage(flux_fwd[i],lev[j])
        mwa_bsize[i][j]=50*50*np.sum(ut.get_bimage(flux_mwa[i],lev[j]))
        fw_bsize[i][j]=22.5*22.5*np.sum(ut.get_bimage(flux_fwd[i],lev[j]))


fw_bsize=np.array(fw_bsize)
mwa_bsize=np.array(mwa_bsize)
flux_bimage=np.array(flux_bimage)
Tb_bimage=np.array(Tb_bimage)
Tb_fwd_bimage=np.array(Tb_fwd_bimage)
fw_bimage=np.array(fw_bimage)
fw_tau_bimage=np.array(fw_tau_bimage)
flux_bimage_sum=np.sum(flux_bimage,axis=(2,3))
fw_bimage_sum=np.sum(fw_bimage,axis=(2,3))

tau=[0]*len(lev)
ep=[0]*len(lev)
fwd_tau_mean=[0]*len(lev)
mwa_Tb_mean=[0]*len(lev)
fwd_Tb_mean=[0]*len(lev)
fwd_size1=[0]*len(lev)
mwa_size1=[0]*len(lev)
print 'Freq','Tb','Tb fwd','tau','tau fwd'
for c in range(len(lev)):
    tau[c]=[0]*len(flist)
    ep[c]=[0]*len(flist)
    fwd_tau_mean[c]=[0]*len(flist)
    mwa_Tb_mean[c]=[0]*len(flist)
    fwd_Tb_mean[c]=[0]*len(flist)
    fwd_size1[c]=[0]*len(flist)
    mwa_size1[c]=[0]*len(flist)
    k=0
    for f in range(len(flist)):
        mwa_sum=np.sum(flux_bimage[f,c])
        idx,arridx=ut.find_nearest(fw_bimage_sum[f],mwa_sum)
        mwa_size1[c][k]=mwa_bsize[k,c]
        fwd_size1[c][k]=fw_bsize[f,idx]
        fwd_tau_mean[c][k]=fw_tau_bimage[f,idx][np.nonzero(Tb_fwd_bimage[f,idx])].mean()
        mwa_Tb_mean[c][k]=Tb_bimage[f,c][np.nonzero(Tb_bimage[f,c])].mean()
        fwd_Tb_mean[c][k]=Tb_fwd_bimage[f,idx][np.nonzero(Tb_fwd_bimage[f,idx])].mean()
        tau[c][k]=np.log(fwd_Tb_mean[c][k]/(fwd_Tb_mean[c][k]-mwa_Tb_mean[c][k]))
        ep[c][k]=np.sqrt(tau[c][k]/fwd_tau_mean[c][k])
        #print flist[k],mwa_Tb_mean[k],fwd_Tb_mean[k],tau[k],fwd_tau_mean[k],ep[k]
        k=k+1
mwa_size1=np.array(mwa_size1)
fwd_size1=np.array(fwd_size1)
fwd_tau_mean=np.array(fwd_tau_mean)
tau=np.array(tau)
mwa_Tb_mean=np.array(mwa_Tb_mean)
fwd_Tb_mean=np.array(fwd_Tb_mean)
ep=np.array(ep)

cc=9
for k in range(len(flist)):
    print flist[k],mwa_Tb_mean[cc][k],fwd_Tb_mean[cc][k],tau[cc][k],fwd_tau_mean[cc][k],ep[cc][k]

diffpos=[0]*len(flabel)    
xpos_fwd=[0]*len(flabel)    
ypos_fwd=[0]*len(flabel)    
xpos_mwa=[0]*len(flabel)    
ypos_mwa=[0]*len(flabel)    
for i in range(8):
    p0=np.array(np.where(Tb_fwd[i]==np.max(Tb_fwd[i])))*22.5
    p1=np.array(np.where(Tb_all[i][0]==np.max(Tb_all[i][0])))*50
    xpos_fwd[i]=p0[0][0]-128*22.5    
    ypos_fwd[i]=p0[1][0]-128*22.5
    xpos_mwa[i]=p1[0][0]-50*50.
    ypos_mwa[i]=p1[1][0]-50*50.
    diffpos[i]=np.sqrt((xpos_mwa[i]-xpos_fwd[i])**2 + (ypos_mwa[i]-ypos_fwd[i])**2)/60.
    
#mwa_flux=[2.52,3.45,4.50,5.49,6.59,8.32,11.52,16.0,18.64]
#emwa_flux=[0.34,0.38,1.30,0.57,0.68,0.86,1.19,1.64,1.93]
############### PLOT
plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')

plot_tau=0
if(plot_tau==1):        
    for j in range(len(epsbyh)):
        plt.plot(r,tau[j])
    plt.axhline(y=1,color='k')
    plt.xlabel('Heliocentric distance (R$_{sun}$)')
    plt.ylabel('$\\tau$')
    plt.title('For 40 MHz')
    #plt.xlim([1.7,2.7])
    plt.show()
plot_ratio=0
if(plot_ratio):
    plt.plot(flist,np.mean(Tb_fwd,axis=(1,2))/np.mean(Tb_all,axis=(1,2,3)),'o-')
    plt.errorbar(flist,np.mean(Tb_fwd,axis=(1,2))/np.mean(Tb_all,axis=(1,2,3)),yerr=np.std(Tb_all,axis=(1,2,3))/np.max(Tb_all)+np.std(Tb_fwd,axis=(1,2))/np.max(Tb_fwd))
    plt.title('Ratio of Mean Forward T$_{B}$ and Observed T$_{B}$')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Ratio')
    plt.show()

plot_deviation=0
if(plot_deviation):
    plt.plot(flist,diffpos,'o-')
    plt.errorbar(flist,diffpos,yerr=bmax)
    plt.title('')
    plt.ylabel('Radial distance (arcmin)')
    plt.xlabel('Frequency (MHz)')
    plt.show()
plot_fluxes=0
if(plot_fluxes==1):
    plt.plot(fall,np.array(fluxall),'o-',color='red',label='FORWARD')
    mwa_Tb_sum=Tb_bimage[f,idx][np.nonzero(Tb_bimage[f,idx])].mean()
    mwa_Tb_sum=Tb_bimage[f,idx][np.nonzero(Tb_bimage[f,idx])].mean()
    plt.errorbar(flist,mwa_flux,emwa_flux,color='blue')
    plt.plot(flist,mwa_flux,'o',color='blue',label='MWA')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    plt.show()


plot_contour_maps=0
if(plot_contour_maps):
    levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    for c in range(len(lev)):
        cont=np.round(lev[c]*100)
        for f in range(len(flist)):
            fig, ax = plt.subplots(1, 1, figsize=(10,10),sharex=True)
            ax.contour(fw_bimage[f,c]/np.max(fw_bimage[f,c]),levels=levels_,extent=[-2500,2500,-2500,2500],colors='blue')
            ax.contour(flux_bimage[f,c]/np.max(flux_bimage[f,c]),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
            ax.set_title('Contour level '+str(cont)+'% '+'Frequency: '+str(int(flist[f]))+' MHz')
            plt.savefig('/home/i4ds1807205/Dropbox/20151203/contours/cont_'+str(int(cont))+'_freq_'+str(int(flist[f]))+'.png')
            plt.close()
            
plot_source_size=0
if(plot_source_size):
    colors = cm.rainbow(np.linspace(0, 1, len(flist)))
    for i in range(8):
        plt.plot(mwa_bsize[i,:]/3600.,fw_bsize[i,:]/3600.,'o-',color=colors[i],label=str(flist[i])+' MHz')
    plt.plot(np.arange(2000),np.arange(2000),'-',color='k',linewidth=4)
    plt.legend(loc=2)
    plt.xlabel('MWA source size (arcmin$^2$)')
    plt.ylabel('FORWARD source size (arcmin$^2$)')
    plt.show()

plot_source_size_analysed=0
if(plot_source_size_analysed):
    colors = cm.rainbow(np.linspace(0, 1, len(flist)))
    for i in range(8):
        y=[0]*len(lev)
        for c in range(len(lev)):
            mwa_sum=np.sum(flux_bimage[i,c])
            idx,arridx=ut.find_nearest(fwd_sum[i],mwa_sum)
            y[c]=fw_bsize[i,idx]
        y=np.array(y)
        plt.plot(mwa_bsize[i,:]/3600.,y/3600.,'o-',color=colors[i],label=str(flist[i])+' MHz')
    plt.plot(np.arange(2000),np.arange(2000),'-',color='k',linewidth=4)
    plt.legend(loc=2)
    plt.xlabel('MWA source size (arcmin$^2$)')
    plt.ylabel('FORWARD source size (arcmin$^2$)')
    plt.show()

plot_contour_maps_analysed=0
if(plot_contour_maps_analysed):
    levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    fwd_sum=np.sum(fw_bimage,axis=(2,3))
    for c in range(len(lev)):
    #for c in range(1):
        cont=np.round(lev[c]*100)
        for f in range(len(flist)):
        #for f in range(1):
            mwa_sum=np.sum(flux_bimage[f,c])
            idx,arridx=ut.find_nearest(fwd_sum[f],mwa_sum)
            mwa_size=np.round(len(np.where(flux_bimage[f,c]!=0)[0])*(50./60)**2)
            fw_size=np.round(len(np.where(fw_bimage[f,idx]!=0)[0])*(22.5/60)**2)
            fig, ax = plt.subplots(1, 1, figsize=(10,10),sharex=True)
            ax.contour(fw_bimage[f,idx]/np.max(fw_bimage[f,idx]),levels=levels_,extent=[-2880,2880,-2880,2880],colors='blue')
            ax.contour(flux_bimage[f,c]/np.max(flux_bimage[f,c]),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
            ax.text(-2100,2000,'FORWARD contour:'+str(np.round(lev[idx]*100))+'% | flux: '+str(np.round(np.sum(fwd_sum[f,idx]),1))+' SFU |'+' Area: '+str(fw_size)+' arcmin$^2$',{'color':'blue','fontsize': 15})
            ax.text(-2100,2200,'MWA flux: '+str(np.round(mwa_sum,1))+' SFU |'+' Area: '+str(mwa_size)+' arcmin$^2$',{'color':'red','fontsize': 15})
            ax.set_title('Contour level '+str(cont)+'% '+'Frequency: '+str(int(flist[f]))+' MHz')
            plt.savefig('/home/i4ds1807205/Dropbox/20151203/contours_analysed/cont_'+str(int(cont))+'_freq_'+str(int(flist[f]))+'.png')
            plt.close()

plot_flux_resolve=0
if(plot_flux_resolve):
    for i in range(8):
        plt.plot(np.linspace(0.1,0.9,10),flux_bimage_sum[i],'o-',color='red',label='MWA')
        plt.plot(np.linspace(0.1,0.9,10),fw_bimage_sum[i],'o-',label='FORWARD')
        plt.errorbar(np.linspace(0.1,0.9,10),flux_bimage_sum[i],yerr=np.ones(10)*emwa_flux[i],color='red')
        plt.legend()
        plt.xlabel('Contours levels (%)')
        plt.ylabel('Flux (SFU)')
        plt.title(str(flist[i])+' MHz')
        plt.show()

plot_composite=0
if(plot_composite==1):
    levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    f, ax = plt.subplots(3, 3, figsize=(20,10),sharex=True)
    im00=ax[0,0].imshow(mwa_100MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im00,ax=ax[0,0],label='T$_{B}$(MK)')
    ax[0,0].contour(mwa_100MHz/np.max(mwa_100MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,0].set_title('100 MHz')
    ax[0,0].grid(True)
    ax[0,0].set_ylabel('arcsec')
    im10=ax[1,0].imshow(I_100MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im10,ax=ax[1,0],label='T$_{B}$(MK)')
    #ax[1,0].contour(I_100MHz/np.max(I_100MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,0].set_xlim([-2500,2500])
    ax[1,0].set_ylim([-2500,2500])
    ax[1,0].grid(True)
    ax[1,0].set_xlabel('arcsec')
    ax[1,0].set_ylabel('arcsec')
    im20=ax[2,0].imshow(I_100MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im20,ax=ax[2,0],label='T$_{B}$(MK)')
    ax[2,0].contour(I_100MHz_convolved/np.max(I_100MHz_convolved),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[2,0].set_xlim([-2500,2500])
    ax[2,0].set_ylim([-2500,2500])
    ax[2,0].grid(True)
    ax[2,0].set_xlabel('arcsec')
    ax[2,0].set_ylabel('arcsec')
    ##
    im01=ax[0,1].imshow(mwa_160MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im01,ax=ax[0,1],label='T$_{B}$(MK)')
    ax[0,1].contour(mwa_160MHz/np.max(mwa_160MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,1].set_title('160 MHz')
    ax[0,1].grid(True)
    im11=ax[1,1].imshow(I_160MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im11,ax=ax[1,1],label='T$_{B}$(MK)')
    #ax[1,1].contour(I_160MHz/np.max(I_160MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,1].set_xlim([-2500,2500])
    ax[1,1].set_ylim([-2500,2500])
    ax[1,1].grid(True)
    ax[1,1].set_xlabel('arcsec')
    im21=ax[2,1].imshow(I_160MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im21,ax=ax[2,1],label='T$_{B}$(MK)')
    ax[2,1].contour(I_160MHz_convolved/np.max(I_160MHz_convolved),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[2,1].set_xlim([-2500,2500])
    ax[2,1].set_ylim([-2500,2500])
    ax[2,1].grid(True)
    ax[2,1].set_xlabel('arcsec')
    ##
    im02=ax[0,2].imshow(mwa_240MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im02,ax=ax[0,2],label='T$_{B}$(MK)')
    ax[0,2].contour(mwa_240MHz/np.max(mwa_240MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,2].set_title('240 MHz')
    ax[0,2].grid(True)
    im12=ax[1,2].imshow(I_240MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im12,ax=ax[1,2],label='T$_{B}$(MK)')
    #ax[1,2].contour(I_240MHz/np.max(I_240MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,2].set_xlim([-2500,2500])
    ax[1,2].set_ylim([-2500,2500])
    ax[1,2].grid(True)
    ax[1,2].set_xlabel('arcsec')
    im22=ax[2,2].imshow(I_240MHz_convolved/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im22,ax=ax[2,2],label='T$_{B}$(MK)')
    ax[2,2].contour(I_240MHz_convolved/np.max(I_240MHz_convolved),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[2,2].set_xlim([-2500,2500])
    ax[2,2].set_ylim([-2500,2500])
    ax[2,2].grid(True)
    ax[2,2].set_xlabel('arcsec')
    f.tight_layout()
    plt.show()

plot_composite_raw=0
if(plot_composite_raw==1):
    levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    f, ax = plt.subplots(2, 3, figsize=(20,10),sharex=True)
    im00=ax[0,0].imshow(mwa_100MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True)
    f.colorbar(im00,ax=ax[0,0],label='T$_{B}$(MK)')
    ax[0,0].contour(mwa_100MHz/np.max(mwa_100MHz),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,0].set_title('100 MHz')
    im10=ax[1,0].imshow(I_100MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im10,ax=ax[1,0],label='T$_{B}$(MK)')
    ax[1,0].contour(I_100MHz/np.max(I_100MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
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
    im10=ax[1,1].imshow(I_160MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im10,ax=ax[1,1],label='T$_{B}$(MK)')
    ax[1,1].contour(I_160MHz/np.max(I_160MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
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
    im12=ax[1,2].imshow(I_240MHz/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True)
    f.colorbar(im12,ax=ax[1,2],label='T$_{B}$(MK)')
    ax[1,2].contour(I_240MHz/np.max(I_240MHz),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,2].set_xlim([-2500,2500])
    ax[1,2].set_ylim([-2500,2500])
    ax[0,2].grid(True)
    ax[1,2].grid(True)
    ax[1,2].set_xlabel('arcsec')
    f.tight_layout()
    plt.show()



#############################################################################################

list_=sorted(glob.glob('/home/i4ds1807205/20151203/*psimas.sav'))
fall,fluxall,Tball=get_forward(list_)

forward_100MHz=readsav('/home/i4ds1807205/20151203/20151203_100MHz_psimas.sav')
forward_160MHz=readsav('/home/i4ds1807205/20151203/20151203_160MHz_psimas.sav')
forward_240MHz=readsav('/home/i4ds1807205/20151203/20151203_240MHz_psimas.sav')

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

#beam_101MHz=beam_101MHz/np.sum(beam_101MHz)
#beam_160MHz=beam_160MHz/np.sum(beam_160MHz)
#beam_240MHz=beam_240MHz/np.sum(beam_240MHz)

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


r=np.linspace(1,215,10000) # Units: Rsun
ne_nk=(4.2*10**(4+(4.32/r))) # Unit: cm^{-3}
ne_nk4=4*ne_nk
ne_st=(1.36e12*(1/r**2.14)+1.68e14*(1/r**6.13))*1.0e-6 # Unit: cm^{-3}
ne_st10=10*ne_st
Rsun=6.957e5 # km

f=240 # in MHz
f=f*1.e6
e=1.9e-19
eps0=8.85e-12
me=9.1e-31
#ompe=np.sqrt((e*e*ne_nk*1.e6)/(eps0*me))
#fpe=ompe/(2*np.pi)
fpe=9000*np.sqrt(ne_nk)
frac=fpe**4 /(f*f-fpe*fpe)**2

dr=r[1]-r[0]
epsbyh=np.linspace(4.5,7.0,50)*1.e-5*(Rsun)
tau=[0]*len(epsbyh)
for j in range(len(epsbyh)):
    tau[j]=[0]*len(r)
    for i in range(len(r)):
        tau[j][i]=np.sum(frac[i:]*(dr)*epsbyh[j])*(0.5*np.sqrt(np.pi))
        #tau[j][i]=np.sum((dr*(len(r)-i)))*(0.5*np.sqrt(np.pi))



plot_mlab=0
if(plot_mlab==1):
    z,x,y=np.mgrid[0:2:3j,-2500:2500:100j,-2500:2500:100j]
    iso=mlab.contour3d(z,x,y,mwa,vmin=mwa.min(),vmax=mwa.max(),opacity=0.3)
    iso.contour.number_of_contours = 15
    mlab.show()

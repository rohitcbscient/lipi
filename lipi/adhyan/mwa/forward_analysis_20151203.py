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
#flist=[108,120,132,145,160,179,196,217,240]
#flabel=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
flist=[108,132,145,160,179,196,217,240]
flabel=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
base=['000-008','000-009','000-010','008-009','008-010','009-010']
mwa_flux_=[0]*len(base)
corr_fact=[0]*len(base)
fringe_fact=[0]*len(base)
i=0
for b in base:
    print b
    mwa_flux_[i]=[0]*len(flabel)
    corr_fact[i]=[0]*len(flabel)
    fringe_fact[i]=[0]*len(flabel)
    j=0
    for f in flabel:
        d=pickle.load(open('/home/i4ds1807205/20151203/pickle/flux_pickle/flux_V1_1133149192-%b'+str(f)+'_T'+str(b)+'.p','rb'))
        d[17][3][0][np.isnan(d[17][3][0])]=0
        mwa_flux_[i][j]=d[17][3][0]
        corr_fact[i][j]=np.mean(d[17][2][0][d[17][2][0]!=0])
        fringe_fact[i][j]=np.mean(d[17][7][0][d[17][7][0]!=0])
        #print mwa_flux[i][j]
        j=j+1
    i=i+1
corr_fact=np.array(corr_fact)
fringe_fact=np.array(fringe_fact)
plot_factor=1
if(plot_factor):
    for i in range(len(base)):
        plt.plot(flist,100.0/fringe_fact[i],label='Tile:'+base[i])
    plt.plot(flist,100.0/corr_fact[0],'o-',label='XX')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('factor')
    plt.legend(loc=3)
    plt.show()

mwa_flux=[0]*len(flist)
emwa_flux=[0]*len(flist)
for i in range(len(flist)):
    emwa_flux[i]=np.std(np.array(mwa_flux_)[:,i,:,:][np.array(mwa_flux_)[:,i,:,:]!=0])
    mwa_flux[i]=np.mean(np.array(mwa_flux_)[:,i,:,:][np.array(mwa_flux_)[:,i,:,:]!=0])
mwa_flux=np.array(mwa_flux) # Average over time, baselines, channels..
emwa_flux=np.array(emwa_flux)
emwa_flux[1]=emwa_flux[1]/2
list_=sorted(glob.glob('/home/i4ds1807205/20151203/*psimas.sav'))
fall,fluxall,Tball=get_forward(list_)

Tb_all=[0]*len(flabel)
bmin=[0]*len(flabel)
bmax=[0]*len(flabel)
Tb_fwd=[0]*len(flabel)
tau=[0]*len(flabel)
h=[0]*len(flabel)
hkm=[0]*len(flabel)
Tb_convolved=[0]*len(flabel)
flux_fwd=[0]*len(flabel)
flux_fwd_decon=[0]*len(flabel)
flux_mwa=[0]*len(flabel)
tempobs=[0]*len(flabel)
taur=[0]*len(flabel)
taul=[0]*len(flabel)
dtaur=[0]*len(flabel)
dtaul=[0]*len(flabel)
tau_fwd=[0]*len(flabel)
tau_convolved=[0]*len(flabel)
eTb_frac=[0]*len(flabel)
rtt=1
k=0
for f in flabel:
    print flist[k],' MHz'
    data=pickle.load(open(dir_+'Tb_1133149192-'+f+'.p','r'))
    data_res=pickle.load(open(dir_+'res_1133149192-'+f+'.p','r'))
    Tb_all[k]=data[0][0:577]
    bmin[k]=data[3][0]*60
    bmax[k]=data[4][0]*60
    eTb_frac[k]=np.std(data_res[0][0])/data[8][0].max()
    forward=readsav('/home/i4ds1807205/20151203/20151203_'+str(flist[k])+'MHz_psimas.sav')
    Tb_fwd[k]=forward['stokesstruct'][0][0]
    ### CONVOLUTION
    forward_dx=forward['quantmap'][0][3]*16*60 # In arcsec
    forward_dy=forward['quantmap'][0][4]*16*60 # In arcsec; originally in Rsun
    fwdx=forward['quantmap'][0][3]*16*60 # Rsun to arcsec
    fwdy=forward['quantmap'][0][4]*16*60 
    fwhmx=int(bmin[k]*60/forward_dx)
    fwhmy=int(bmax[k]*60/forward_dy)
    beam=ut.make2DGaussian(60, fwhmx, fwhmy, center=None)
    Tb_convolved[k] = signal.convolve(Tb_fwd[k],beam, mode='same')/np.sum(beam)
    size=fwdx#*Tb_fwd[k].shape[0]
    flux_fwd[k]=ut.Tb2flux(Tb_convolved[k], size, size, flist[k]/1000.)
    flux_fwd_decon[k]=ut.Tb2flux(Tb_fwd[k], size, size, flist[k]/1000.)
    flux_mwa[k]=np.mean(Tb_all[k],axis=0)*mwa_flux[k]/np.sum(Tb_all[k][0])#ut.Tb2flux(Tb_all[k][0], 50, 50, flist[k]/1000.)
    ###
    if(rtt):
        rt=readsav('/home/i4ds1807205/20151203/RT_params_'+str(flist[k])+'MHz.sav')
        tempobs=rt['tempall']
        densobs=rt['densall']
        brobs=rt['brall']
        bthobs=rt['bthall']
        bphobs=rt['bphall']
        taur[k]=rt['taur'][:,:,-1]
        taul[k]=rt['taul'][:,:,-1]
        tau_fwd[k]=taul[k]+taur[k]
        dtaur[k]=rt['dtaur']
        dtaul[k]=rt['dtaul']
        tau_convolved[k] = signal.convolve(tau_fwd[k],beam, mode='same')/np.sum(beam)
    tau[k]=np.log((0.2*Tb_fwd[k].max()/np.mean(Tb_all[k],axis=(1,2)))-1)
    #print flist[k],np.max(Tb_fwd[k]),np.max(Tb_all[k]),np.std(Tb_all[k]),np.mean(tau[k]),np.std(tau[k])
    h[k]=md.nk_freq2r(flist[k],1)[0]
    hkm[k]=(h[k]-1)*3.e5
    k=k+1
Tb_fwd=np.array(Tb_fwd)
Tb_all=np.array(Tb_all)

nlev=30
lev=np.linspace(0.1,0.9,nlev)
Tb_bimage=[0]*8
Tb_fwd_bimage=[0]*8
flux_bimage=[0]*8
fw_bimage=[0]*8
fw_bimage_decon=[0]*8
fw_bsize=[0]*8
fw_bsize_decon=[0]*8
fw_tau_bimage=[0]*8
mwa_bsize=[0]*8
mwa_full_size=[0]*8
fwd_full_size=[0]*8
mwa_limit=[0]*8
mwab=[0]*8
fwdb=[0]*8
for i in range(8):
    Tb_bimage[i]=[0]*nlev
    flux_bimage[i]=[0]*nlev
    fw_bimage[i]=[0]*nlev
    fw_bimage_decon[i]=[0]*nlev
    fw_tau_bimage[i]=[0]*nlev
    Tb_fwd_bimage[i]=[0]*nlev
    fw_bsize[i]=[0]*nlev
    fw_bsize_decon[i]=[0]*nlev
    mwa_bsize[i]=[0]*nlev
    mwab[i]=[0]*nlev
    fwdb[i]=[0]*nlev
    for j in range(nlev):
        mwab[i][j]=ut.get_bimage(flux_mwa[i],lev[j])
        fwdb[i][j]=ut.get_bimage(flux_fwd[i],lev[j]) # flux_fwd is convolved
        flux_bimage[i][j]=flux_mwa[i]*ut.get_bimage(flux_mwa[i],lev[j])
        Tb_bimage[i][j]=np.mean(Tb_all[i],axis=0)*ut.get_bimage(flux_mwa[i],lev[j])
        fw_bimage[i][j]=flux_fwd[i]*ut.get_bimage(flux_fwd[i],lev[j])
        fw_bimage_decon[i][j]=flux_fwd_decon[i]*ut.get_bimage(flux_fwd_decon[i],lev[j])
        Tb_fwd_bimage[i][j]=Tb_convolved[i]*ut.get_bimage(flux_fwd[i],lev[j])
        fw_tau_bimage[i][j]=tau_convolved[i]*ut.get_bimage(flux_fwd[i],lev[j])
        mwa_bsize[i][j]=50*50*np.sum(ut.get_bimage(flux_mwa[i],lev[j]))
        fw_bsize[i][j]=22.5*22.5*np.sum(ut.get_bimage(flux_fwd[i],lev[j]))
        fw_bsize_decon[i][j]=22.5*22.5*np.sum(ut.get_bimage(flux_fwd_decon[i],lev[j]))
    mwa_flux_array=flux_mwa[i][flux_mwa[i]!=0]
    mwa_Tb_array=Tb_all[i][0][Tb_all[i][0]!=0]
    mwa_full_size[i]=50*50*len(mwa_flux_array)/3600.
    mwa_limit_flux=np.min(mwa_flux_array)
    mwa_limit[i]=np.min(mwa_Tb_array)
    fwd_flux_array=flux_fwd[i][flux_fwd[i]>mwa_limit_flux]
    fwd_full_size[i]=22.5*22.5*len(fwd_flux_array)/3600.


mwa_limit=np.array(mwa_limit)
fwd_full_size=np.array(fwd_full_size)
mwa_full_size=np.array(mwa_full_size)
fw_bsize=np.array(fw_bsize)
fw_bsize_decon=np.array(fw_bsize_decon)
mwa_bsize=np.array(mwa_bsize)
flux_bimage=np.array(flux_bimage)
Tb_bimage=np.array(Tb_bimage)
Tb_fwd_bimage=np.array(Tb_fwd_bimage)
fw_bimage=np.array(fw_bimage)
fw_bimage_decon=np.array(fw_bimage_decon)
fw_tau_bimage=np.array(fw_tau_bimage)
flux_bimage_sum=np.sum(flux_bimage,axis=(2,3))
fw_bimage_sum=np.sum(fw_bimage,axis=(2,3))
fw_bimage_sum_decon=np.sum(fw_bimage_decon,axis=(2,3))

tau=[0]*len(lev)
ep=[0]*len(lev)
fwd_tau_mean=[0]*len(lev)
mwa_Tb_mean=[0]*len(lev)
fwd_Tb_mean=[0]*len(lev)
fwd_size1=[0]*len(lev)
mwa_size1=[0]*len(lev)
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
        tau[c][k]=np.log(0.2*fwd_Tb_mean[c][k]/(0.2*fwd_Tb_mean[c][k]-mwa_Tb_mean[c][k]))
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
thermal_noise_array=np.array([1.44e-5,1.375e-5,1.374e-5,1.376e-5,1.49e-5,1.799e-5,2.319e-5,2.716e-5])

print_size=1
if(print_size):
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(mwa_full_size[k]/1.e3,3),'$\pm$',np.round(np.sqrt(mwa_full_size[k])/1.e3,3),' & ',np.round(mwa_limit[k]/1.e6,3),'\\\\'

print_Tb=1
if(print_Tb):
    cc=0
    print 'Frequency (MHz)','&','T$_{B,MWA}$ (MK)','&','T$_{B,FWD}$ (MK)','&','$\\tau_{MWA}$','$\\tau_{FWD}$'
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(mwa_Tb_mean[cc][k]/1.e6,3),'$\pm$',np.round(mwa_Tb_mean[cc][k]*(emwa_flux[k]/mwa_flux[k])/1.e6,3),' & ',np.round(fwd_Tb_mean[cc][k]*0.2/1.e6,2),' & ',np.round(tau[cc][k],4),'$\pm$',np.round(tau[cc][k]*eTb_frac[k]*tau[cc][k],4),' & ',np.round(fwd_tau_mean[cc][k],3),'\\\\'#,ep[cc][k]

levels_=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]
fwd_sum=np.sum(fw_bimage,axis=(2,3))
fwd_sum_decon=np.sum(fw_bimage_decon,axis=(2,3))
mwa_size=[0]*len(lev)
fw_size=[0]*len(lev)
efw_size=[0]*len(lev)
fw_size_decon=[0]*len(lev)
mwa_sum=[0]*len(lev)
fwd_sum_mwa=[0]*len(lev)
diff_flux=[0]*len(lev)
idx=[0]*len(lev)
for c in range(len(lev)):
#for c in range(1):
    cont=np.round(lev[c]*100)
    mwa_size[c]=[0]*len(flist)
    mwa_sum[c]=[0]*len(flist)
    fwd_sum_mwa[c]=[0]*len(flist)
    fw_size[c]=[0]*len(flist)
    efw_size[c]=[0]*len(flist)
    fw_size_decon[c]=[0]*len(flist)
    diff_flux[c]=[0]*len(flist)
    idx[c]=[0]*len(flist)
    for f in range(len(flist)):
    #for f in range(1):
        mwa_sum[c][f]=np.sum(flux_bimage[f,c])
        idx[c][f],arridx=ut.find_nearest(fwd_sum[f],mwa_sum[c][f])
        fwd_sum_mwa[c][f]=np.sum(fwd_sum[f][idx[c][f]])
        mwa_size[c][f]=np.round(len(np.where(flux_bimage[f,c]!=0)[0])*(50./60)**2)
        fw_size[c][f]=np.round(len(np.where(fw_bimage[f,idx[c][f]]!=0)[0])*(22.5/60)**2)
        fw_size_decon[c][f]=np.round(len(np.where(fw_bimage_decon[f,idx[c][f]]!=0)[0])*(22.5/60)**2)
        ### Upper limit ###
        idx_up,arridx_up=ut.find_nearest(fwd_sum[f],mwa_sum[c][f]+emwa_flux[f]*0.5)
        fw_size_up=np.round(len(np.where(fw_bimage[f,idx_up]!=0)[0])*(22.5/60)**2)
        efw_size[c][f]=fw_size_up-fw_size[c][f]
        diff_flux[c][f]=2*abs(mwa_sum[c][f]-fwd_sum_mwa[c][f])/(mwa_sum[c][f]+fwd_sum_mwa[c][f])

mwa_size=np.array(mwa_size)
fwd_sum_mwa=np.array(fwd_sum_mwa)
mwa_sum=np.array(mwa_sum)
fw_size=np.array(fw_size)
fw_size_decon=np.array(fw_size_decon)
idx=np.array(idx)
diff_flux=np.array(diff_flux)

## Radius of the disk
fmiss_size=mwa_size[0]-fw_size[0]
rmwa=np.sqrt(mwa_size[0]/np.pi)
reff=(rmwa/np.pi)*(np.sqrt(1+(np.pi*fmiss_size/rmwa**2))-1)
plot_reff=1
if(plot_reff):
    plt.plot(flist,reff,'o-',color='blue')
    plt.errorbar(flist,reff,yerr=reff*((emwa_flux/mwa_flux)+(efw_size[0]/fw_size[0])),color='blue')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Additional Radius (in arcmin)')
    plt.xlim(100,250)
    plt.show()

plot_size_comparison=0
if(plot_size_comparison):
    cont_list=[0]
    for cc in cont_list:
        plt.plot(flist,mwa_size[cc],'o',color='blue')
        plt.errorbar(flist,mwa_size[cc],yerr=mwa_size[cc]*emwa_flux/mwa_flux,label='MWA',color='blue')
        plt.plot(flist,fw_size[cc],'o-',color='red')
        plt.errorbar(flist,fw_size[cc],yerr=efw_size[cc],color='red',label='FORWARD')
    plt.legend(loc=1)
    plt.xlim(100,250)
    plt.ylim(300,2400)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Size (arcmin$^{2}$)')
    plt.show()
plot_size_comparison_diff=0
if(plot_size_comparison_diff):
    cont_list=[0]
    for cc in cont_list:
        plt.plot(flist,mwa_size[cc]-fw_size[cc],'o',color='k')
        plt.errorbar(flist,mwa_size[cc]-fw_size[cc],yerr=(mwa_size[cc]-fw_size[cc])*((emwa_flux/mwa_flux)+(efw_size[cc]/fw_size[cc])),label='MWA',color='k')
    plt.legend(loc=1)
    plt.xlim(100,250)
    plt.ylim(-10,2400)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Size (arcmin$^{2}$)')
    plt.show()

Tb_mean=np.nanmean(Tb_all,axis=1)
Tb_mwa_ch=np.mean(Tb_mean[:,55:61,48:54],axis=(1,2))/1.e6
Tb_mwa_qs=np.mean(Tb_mean[:,40:45,42:47],axis=(1,2))/1.e6
Tb_mwa_ar=np.mean(Tb_mean[:,41:46,58:63],axis=(1,2))/1.e6
Tb_fwd_ch=np.mean(Tb_fwd[:,150:165,110:130],axis=(1,2))/1.e6
Tb_fwd_ar=np.mean(Tb_fwd[:,125:140,160:172],axis=(1,2))/1.e6
Tb_fwd_qs=np.mean(Tb_fwd[:,106:116,100:115],axis=(1,2))/1.e6
plot_Tb_regions=0
if(plot_Tb_regions):
    plt.plot(freq,Tb_mwa_ch,'o',color='blue')
    plt.errorbar(freq,Tb_mwa_ch,yerr=Tb_mwa_ch*(emwa_flux/mwa_flux),color='blue',label='MWA CH')
    plt.plot(freq,Tb_fwd_ch,'o--',color='blue',label='FWD CH')
    plt.plot(freq,Tb_mwa_ar,'o',color='red')
    plt.errorbar(freq,Tb_mwa_ar,yerr=Tb_mwa_ar*(emwa_flux/mwa_flux),color='red',label='MWA AR')
    plt.plot(freq,Tb_fwd_ar,'o--',color='red',label='FWD AR')
    plt.plot(freq,Tb_mwa_qs,'o',color='green')
    plt.errorbar(freq,Tb_mwa_qs,yerr=Tb_mwa_qs*(emwa_flux/mwa_flux),color='green',label='MWA QS')
    plt.plot(freq,Tb_fwd_qs,'o--',color='green',label='FWD QS')
    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('T$_{B}$ (MK)')
    plt.show()

sys.exit()
powlaw=[0]*len(flist)
epowlaw=[0]*len(flist)
mwa_total_area=[0]*len(flist)
for i in range(len(flist)):
    powlaw[i],intcpt,epowlaw[i],eintcpt=ut.fit_1d(mwa_size[:,i],fw_size_decon[:,i])
    mwa_total_area[i]=(1600)*powlaw[i]+intcpt


##########################################

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
plot_ratio=1
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
plot_fluxes=1
if(plot_fluxes==1):
    plt.plot(fall,np.array(fluxall),'o-',color='red',label='FORWARD')
    #mwa_Tb_sum=Tb_bimage[f,idx][np.nonzero(Tb_bimage[f,idx])].mean()
    plt.errorbar(flist,mwa_flux,emwa_flux,color='blue',label='MWA')
    plt.plot(flist,mwa_flux,'o',color='blue')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    plt.show()

plot_fluxes_diff=1
if(plot_fluxes_diff==1):
    del fluxall[0]
    del fluxall[1]
    plt.errorbar(flist,np.array(fluxall)-mwa_flux,emwa_flux,color='k')
    plt.plot(flist,np.array(fluxall)-mwa_flux,'o',color='k')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('$\Delta S$ (SFU)')
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


plot_contour_maps_analysed=1
if(plot_contour_maps_analysed):
            fig, ax = plt.subplots(1, 1, figsize=(10,10),sharex=True)
            ax.contour(fw_bimage[f,idx]/np.max(fw_bimage[f,idx]),levels=levels_,extent=[-2880,2880,-2880,2880],colors='blue')
            ax.contour(flux_bimage[f,c]/np.max(flux_bimage[f,c]),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
            ax.text(-2100,2000,'FORWARD contour:'+str(np.round(lev[idx]*100))+'% | flux: '+str(np.round(np.sum(fwd_sum[f,idx]),1))+' SFU |'+' Area: '+str(fw_size)+' arcmin$^2$',{'color':'blue','fontsize': 15})
            ax.text(-2100,2200,'MWA flux: '+str(np.round(mwa_sum,1))+' SFU |'+' Area: '+str(mwa_size)+' arcmin$^2$',{'color':'red','fontsize': 15})
            ax.set_title('Contour level '+str(cont)+'% '+'Frequency: '+str(int(flist[f]))+' MHz')
            #plt.savefig('/home/i4ds1807205/Dropbox/20151203/contours_analysed/cont_'+str(int(cont))+'_freq_'+str(int(flist[f]))+'.png')
            plt.close()

image_fwd_mwa=1
if(image_fwd_mwa):
    c=0
    fr1=7
    fr2=6
    f, ax = plt.subplots(2, 2, figsize=(20,10),sharex=True)
    im00=ax[0,0].imshow(Tb_bimage[fr1,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.1,vmax=0.4)
    f.colorbar(im00,ax=ax[0,0],label='T$_{B}$(MK)')
    ax[0,0].contour(Tb_bimage[fr1,c]/Tb_bimage[fr1,c].max(),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,0].set_title('240 MHz')
    ax[0,0].set_ylabel('arcsec')
    ax[0,0].grid(True) # 0.25 is due to difference in resolution 
    im10=ax[1,0].imshow(0.25*Tb_fwd_bimage[fr1,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True,vmin=0.1,vmax=0.4)
    f.colorbar(im10,ax=ax[1,0],label='T$_{B}$(MK)')
    ax[1,0].contour(Tb_fwd_bimage[fr1,c]/np.max(Tb_fwd_bimage[fr1,c]),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,0].set_xlim([-2500,2500])
    ax[1,0].grid(True)
    ax[1,0].set_ylabel('arcsec')
    ax[1,0].set_xlabel('arcsec')
    im01=ax[0,1].imshow(Tb_bimage[fr2,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.1,vmax=0.4)
    f.colorbar(im01,ax=ax[0,1],label='T$_{B}$(MK)')
    ax[0,1].contour(Tb_bimage[fr2,c]/Tb_bimage[fr2,c].max(),levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,1].set_title('217 MHz')
    ax[0,1].grid(True)
    im11=ax[1,1].imshow(0.25*Tb_fwd_bimage[fr2,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-1440*2,1440*2,-1440*2,1440*2],origin=True,vmin=0.1,vmax=0.4)
    f.colorbar(im11,ax=ax[1,1],label='T$_{B}$(MK)')
    ax[1,1].contour(Tb_fwd_bimage[fr2,c]/np.max(Tb_fwd_bimage[fr2,c]),levels=levels_,extent=[-1440*2,1440*2,-1440*2,1440*2],colors='red')
    ax[1,1].set_xlim([-2500,2500])
    ax[1,1].grid(True)
    ax[1,1].set_xlabel('arcsec')
    plt.show()
    

plot_source_size_analysed=1
if(plot_source_size_analysed):
    i_=[0,1,6,7]
    colors = cm.rainbow(np.linspace(0, 1, len(i_)))
    for i in i_:
        y=[0]*len(lev)
        y_decon=[0]*len(lev)
        for c in range(len(lev)):
            mwa_sum=np.sum(flux_bimage[i,c])
            idx_decon,arridx_decon=ut.find_nearest(fwd_sum_decon[i],mwa_sum)
            idx,arridx=ut.find_nearest(fwd_sum[i],mwa_sum)
            y[c]=fw_bsize[i,idx]
            y_decon[c]=fw_bsize_decon[i,idx_decon]
        y=np.array(y)
        y_decon=np.array(y_decon)
        plt.plot(mwa_bsize[i,:]/3600.,y/3600.,'o-',color=colors[i_.index(i)],label=str(flist[i])+' MHz')
        #plt.plot(mwa_bsize[i,:]/3600.,y_decon/3600.,'o-',color=colors[i],label=str(flist[i])+' MHz')
    plt.plot(np.arange(2000),np.arange(2000),'-',color='k',linewidth=4)
    plt.legend(loc=2)
    plt.xlabel('MWA source size (arcmin$^2$)')
    plt.ylabel('FORWARD source size (arcmin$^2$)')
    plt.show()
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

beam_101MHz=ut.make2DGaussian(60, fwhmx_100MHz, fwhmy_100MHz , center=None)
beam_160MHz=ut.make2DGaussian(60, fwhmx_160MHz, fwhmy_160MHz , center=None)
beam_240MHz=ut.make2DGaussian(60, fwhmx_240MHz, fwhmy_240MHz , center=None)

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

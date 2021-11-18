import numpy as np
import matplotlib.pyplot as plt
import glob
import pickle
from surya.utils import main as ut
from surya.plot import main as spl
from surya.utils import main as ut
from scipy.io import readsav
from scipy import signal
from surya.utils import model as md



dir_='/media/rohit/VLA/20151203_MWA/Tb_new/'
flabel=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
flist=[108,132,145,160,179,196,217,240]
bmin=[3.63,3.26,3.28,3.09,3.09,3.11,3.09,3.37]
bmax=[3.75,3.42,3.35,3.30,3.21,3.14,3.13,3.48]
Tb_fwd=[0]*len(flabel)
h=[0]*len(flabel)
hkm=[0]*len(flabel)
Tb_convolved=[0]*len(flabel)
Tb_convolved_los=[0]*len(flabel)
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
flux_los_mwa=[0]*len(flabel)
flux_los_fpe=[0]*len(flabel)
flux_los_r=[0]*len(flabel)
flux_los=[0]*len(flabel)
rflux_los=[0]*len(flabel)
Tb_fwd_los=[0]*len(flabel)
rtt=1
k=0
for f in flabel:
    print flist[k],' MHz'
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
    ###
    if(rtt):
        rt=readsav('/home/i4ds1807205/20151203/RT_params_'+str(flist[k])+'MHz.sav')
        r=rt['r3dall']
        tempobs=rt['tempall']
        densobs=rt['densall']
        brobs=rt['brall']
        bthobs=rt['bthall']
        bphobs=rt['bphall']
        bobs=np.sqrt(brobs*brobs+bthobs*bthobs+bphobs*bphobs)
        taur[k]=rt['taur']
        taul[k]=rt['taul']
        tau_fwd[k]=(taul[k]+taur[k])*0.5
        dtaur[k]=rt['dtaur']
        dtaul[k]=rt['dtaul']
        fpe=9000.*np.sqrt(densobs)
        flux_los[k]=[0]*(tempobs.shape[2]-1)
        Tb_fwd_los[k]=[0]*(tempobs.shape[2])
        Tb_fwd_los[k][0]=np.zeros((100,100))+3.0
        for l in range(tempobs.shape[2]-1):
            tempobs_los=signal.convolve(tempobs[:,:,l],beam,mode='same')/np.sum(beam)
            dtaul_los=signal.convolve(dtaul[k][:,:,l],beam,mode='same')/np.sum(beam)
            dtaur_los=signal.convolve(dtaur[k][:,:,l],beam,mode='same')/np.sum(beam)
            Tb_fwd_los[k][l+1]=(Tb_fwd_los[k][l]*np.exp(-1*dtaul[k][:,:,l])+Tb_fwd_los[k][l]*np.exp(-1*dtaur[k][:,:,l])+tempobs[:,:,l]*(1.0-np.exp(-1*dtaul[k][:,:,l]))+tempobs[:,:,l]*(1.0-np.exp(-1*dtaur[k][:,:,l])))/2.0
            #Tb_fwd_los=np.sum((tempobs[:,:,:l]*(1.0-np.exp(-1*dtaul[k][:,:,:l]))+tempobs[:,:,:l]*(1.0-np.exp(-1*dtaur[k][:,:,:l])))/2.0,axis=2)
            flux_los_=ut.Tb2flux(Tb_fwd_los[k][l+1], size, size, flist[k]/1000.)
            flux_los[k][l]=signal.convolve(flux_los_,beam,mode='same')/np.sum(beam)
        tau_convolved[k] = signal.convolve(tau_fwd[k][:,:,-1],beam, mode='same')/np.sum(beam)
    print np.array(flux_los[k])[-1].sum()
    flux_los_fpe[k]=flux_los[k][-1].sum(axis=(0,1))-np.array(flux_los[k]).sum(axis=(1,2))[ut.find_nearest(fpe[50,50],flist[k]*1.e6)[0]]
    flux_los_r[k]=r[50,50][ut.find_nearest(fpe[50,50],flist[k]*1.e6)[0]]
    h[k]=md.nk_freq2r(flist[k],1)[0]
    hkm[k]=(h[k]-1)*3.e5
    k=k+1
Tb_fwd=np.array(Tb_fwd)
flux_los=np.array(flux_los)
Tb_fwd_los=np.array(Tb_fwd_los)

## Polar
polTb_fwd=[0]*len(flist)
for j in range(len(flist)):
    polTb_fwd[j],r_fwd,th_fwd=ut.cart2polar(Tb_fwd[j])
polTb_fwd=np.array(polTb_fwd)


sys.exit()
####
baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
freq_filelist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
freq=[108.0,120.0,133.0,145.0,161.0,179.0,197.0,217.0,240.0]
freq_filelist=[freq_filelist[0],freq_filelist[1],freq_filelist[7],freq_filelist[8]]
idx=[0,1,7,8]
Tb_sub_max=[0]*4
Tb_sub_hot=[0]*4
Tb_sub_mean=[0]*4
Tb_sub_cool=[0]*4
for k in range(len(freq_filelist)):
    f=freq_filelist[k]
    j=idx[k]
    Tb_path='/nas08-data02/rohit/20151203_MWA/Tb_new/'
    Tb_sub_path='/nas08-data02/rohit/20151203_sub/Tb/'
    Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(Tb_path+'/Tb_20151203_'+f+'.p','rb'))
    Tb_sub,Tb_sub_std,xc,yc,time_string,time_sec,bmaj,bmin,bpa=pickle.load(open(Tb_sub_path+'/Tb_20151203_'+f+'_sub.p','rb'))
    Tb_sub_max[k]=[0]*2000#*len(Tb_sub)
    Tb_sub_hot[k]=[0]*2000#*len(Tb_sub)
    Tb_sub_cool[k]=[0]*2000#*len(Tb_sub)
    Tb_sub_mean[k]=[0]*2000#*len(Tb_sub)
    #n=(Tb_sub_std[99]*15.0)/Tb_sub[99].max()
    n=(Tb_sub_std[99]*5.0)/Tb_sub[99].max()
    #for i in range(len(Tb_sub)):
    for i in range(2000):
        Tb_sub[i][np.where(Tb_sub[i]==0)]=1.e-16
        Tb_sub_max[k][i]=np.nanmax(Tb_sub[i])
        bimage_hot=ut.get_bimage(Tb_sub[i],n)
        bimage_mean=ut.get_bimage(abs(Tb_sub[i]),n)
        Tb_sub_mean[k][i]=np.nanmean(Tb_sub[i][bimage_mean.astype(bool)])
        Tb_sub_hot[k][i]=np.sum(Tb_sub[i][bimage_hot.astype(bool)])/np.sum(bimage_hot)
        Tb_neg=Tb_sub[i]*1.0
        Tb_neg[Tb_neg>0]=1.e-16
        bimage_cool=ut.get_bimage(abs(Tb_neg),n)
        Tb_sub_cool[k][i]=np.sum(Tb_neg[bimage_cool.astype(bool)])/np.sum(bimage_cool)

Tb_sub_max=np.array(Tb_sub_max)
Tb_sub_mean=np.array(Tb_sub_mean)
Tb_sub_cool=np.array(Tb_sub_cool)
Tb_sub_hot=np.array(Tb_sub_hot)
Tb_sub_mean[:,822:1300]=np.nan
Tb_sub_max[:,822:1300]=np.nan
Tb_sub_hot[:,822:1300]=np.nan
Tb_sub_cool[:,822:1300]=np.nan

fig=plt.figure()
ax1=fig.add_subplot(211,aspect='auto')
for k in range(4):
    ax1.plot(-1*Tb_sub_cool[k],'o',label=str(freq[idx[k]])+' MHz')
#ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('$\Delta T_B$ (K)')
ax1.legend(loc=2)
ax2=fig.add_subplot(212,aspect='auto',sharex=ax1,sharey=ax1)
for k in range(4):
    ax2.plot(Tb_sub_hot[k],'o',label=str(freq[idx[k]])+' MHz')
ax2.set_xlabel('Time (sec)')
ax2.set_ylabel('$\Delta T_B$ (K)')
ax2.legend(loc=2)
#ax2.set_xlim([0,800])
ax2.set_ylim([0,2200])
plt.show()

fig=plt.figure()
ax1=fig.add_subplot(211,aspect='auto')
ax1.plot(-1*Tb_sub_cool[0],'o',label=str(freq[idx[0]])+' MHz')
ax1.plot(Tb_sub_hot[0],'o',label=str(freq[idx[0]])+' MHz')
#ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('$\Delta T_B$ (K)')
ax2=fig.add_subplot(212,aspect='auto',sharex=ax1,sharey=ax1)
ax2.plot(-1*Tb_sub_cool[3],'o',label=str(freq[idx[3]])+' MHz')
ax2.plot(Tb_sub_hot[3],'o',label=str(freq[idx[3]])+' MHz')
ax2.set_xlabel('Time (sec)')
ax2.set_ylabel('$\Delta T_B$ (K)')
ax1.legend(loc=2)
ax2.legend(loc=2)
#ax2.set_xlim([0,800])
ax2.set_ylim([0,2200])
plt.show()


sys.exit()



plt.plot(abs(np.array(Tb_sub_mean)),'o-',label='MEAN')
plt.plot(abs(np.array(Tb_sub_hot)),'o-',label='HOT') 
plt.plot(abs(np.array(Tb_sub_cool)),'o-',label='COOL')
plt.plot(abs(np.array(Tb_sub_max)),'o-',label='MAX')
plt.ylabel('$T_B$ (K)')
plt.xlabel('Time (sec)')
plt.legend(loc=2)
plt.savefig('pngs_'+str(int(freq[j]))+'MHz/timeseries_'+str(f)+'.png')
plt.close()

plot_Tbsub=1
if(plot_Tbsub):
    for i in range(len(Tb_sub)):
        aa=Tb_sub[i]/1.e3
        levels=np.array([-30,-20,-10,-5,5,10,20,30])*Tb_sub_std[99]/1.e3
        plt.rcParams["contour.negative_linestyle"] = 'dashed'
        stri="%04d"%i
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=-4,vmax=4,cmap='coolwarm')
        #CS=ax.contour(aa, levels,extent=[-2500,2500,-2500,2500],  colors='k',linewidths=2)
        im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-5000,5000,-5000,5000],origin='lower',vmin=-4,vmax=4,cmap='coolwarm')
        CS=ax.contour(aa, levels,extent=[-5000,5000,-5000,5000],  colors='k',linewidths=2)
        #ax.contour(aa, levels,linestyles='--',extent=[-2500,2500,-2500,2500], colors='k',linewidths=2)
        spl.add_beam(ax,-2000, -2000,bmaj*3600,bmin*3600,bpa)
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
        #ax.set_xlim(-2500,2500)
        #ax.set_ylim(-2500,2500)
        ax.set_xlim(-5000,5000)
        ax.set_ylim(-5000,5000)
        ax.set_title(str(freq[j])+' MHz Time:'+str(time_string[i])+' UT')
        r1=16.*60
        r2=32.*60
        circ1=plt.Circle((0.5,0.5), radius=r1, color='brown', linewidth=4,fill=False)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='brown', linewidth=4,fill=False)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.grid(True)
        fig.colorbar(im,label='(kK)')
        fig.savefig('pngs_'+str(int(freq[j]))+'MHz/contour_'+str(stri)+'.png',dpi=100)
        plt.close()


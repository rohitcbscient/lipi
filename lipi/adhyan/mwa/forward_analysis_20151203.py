import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from scipy.io import readsav
import pickle
from scipy import signal
from surya.utils import main as ut
from surya.utils import model as md
import glob
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import proj3d
from matplotlib import patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')


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

# central_freq,chan_list,auto_t1,auto_t2,cross,ncross,phase_ncross,u,v,w,azi_pointing,ele_pointing,ph_az,ph_el,start_time,mid_time,end_time,[0,0,corr_factor,S_sun,T_sun,Un_Tbeam_Sun,Temp_beam_sun,fringe_factor,T_baseline,Tsky_integrated]
corr_fact=np.array(corr_fact)
fringe_fact=np.array(fringe_fact)
plot_factor=0
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
    data=pickle.load(open(dir_+'Tb_1133149192-'+f+'.p','r'))
    data_res=pickle.load(open(dir_+'res_1133149192-'+f+'.p','r'))
    Tb_all[k]=np.array(data[0][0:577])*np.pi # Due to extra pi in the compute_Tb in omega 
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
    rflux_los[k]=r[:,:,ut.find_nearest(flux_los[k][-1].sum(axis=(0,1))-np.array(flux_los[k]).sum(axis=(1,2)),mwa_flux[k])[0]] # location in r which matches MWA flux
    flux_los_mwa[k]=flux_los[k][-1].sum(axis=(0,1))-np.array(flux_los[k]).sum(axis=(1,2))[ut.find_nearest(flux_los[k][-1].sum(axis=(0,1))-np.array(flux_los[k]).sum(axis=(1,2)),mwa_flux[k])[0]]
    flux_los_fpe[k]=flux_los[k][-1].sum(axis=(0,1))-np.array(flux_los[k]).sum(axis=(1,2))[ut.find_nearest(fpe[50,50],flist[k]*1.e6)[0]]  
    flux_los_r[k]=r[50,50][ut.find_nearest(fpe[50,50],flist[k]*1.e6)[0]]  
    tau[k]=np.log((0.2*Tb_fwd[k].max()/np.mean(Tb_all[k],axis=(1,2)))-1)
    #print flist[k],np.max(Tb_fwd[k]),np.max(Tb_all[k]),np.std(Tb_all[k]),np.mean(tau[k]),np.std(tau[k])
    h[k]=md.nk_freq2r(flist[k],1)[0]
    hkm[k]=(h[k]-1)*3.e5
    k=k+1
Tb_fwd=np.array(Tb_fwd)
Tb_all=np.array(Tb_all)
flux_los=np.array(flux_los)
Tb_fwd_los=np.array(Tb_fwd_los)

############################## RLOS #######################
rlos=np.linspace(1.01,1.2,501)
rlos=np.linspace(1.01,2.2,501)
#rlos=np.linspace(1.01,2.0,101)
flux_los_r1=[0]*8
flux_los_r2=[0]*8
rflux_los_mwa=[0]*8
rflux_los_fpe=[0]*8
flux_los_fpe1=[0]*8
fpe_mean=[0]*(len(rlos)-1)
fpe_std=[0]*(len(rlos)-1)
fpe_hist=[0]*(len(rlos)-1)
densfpe_hist=[0]*(len(rlos)-1)
rfpe_mean=[0]*8
rfpe_std=[0]*8
for i in range(8):
    print str(flist[i])+' MHz'
    flux_los_r1[i]=[0]*(len(rlos)-1)
    flux_los_r2[i]=[0]*(len(rlos)-1)
    for j in range(len(rlos)-1):
        los_idx_for_dens=np.where((rlos[j]<r[:,:,:-1]) & (r[:,:,:-1]<rlos[j+1]))
        flux_los_r1[i][j]=[0]*100
        for l in range(100):
            flux_los_r1[i][j][l]=[0]*100
            for m in range(100):
                los_idx=np.array(np.where(r[l,m,:-2]>rlos[j]))
                if(len(los_idx[0])>1):
                    flux_los_r1[i][j][l][m]=flux_los[i,-1,l,m]-flux_los[i,los_idx[0][0],l,m]
                else:
                    flux_los_r1[i][j][l][m]=0
        #flux_los_r1[i][j]=flux_los[i][-1].sum()-flux_los[i][los_idx[2],los_idx[1],los_idx[0]].sum()
        flux_los_r2[i][j]=np.sum(np.array(flux_los_r1[i][j]))
        fpe_hist[j]=fpe[los_idx_for_dens][fpe[los_idx_for_dens]>1.e6]
        densfpe_hist[j]=densobs[los_idx_for_dens][fpe[los_idx_for_dens]>1.e6]
        fpe_mean[j]=np.median(fpe_hist[j])
        fpe_std[j]=np.std(fpe_hist[j])
        #flux_los_r1[i][j]=np.sum(flux_los[i].swapaxes(0,2)[los_idx])
        #denom[i][j]=(np.sqrt(np.pi)/2)*((fpe_mean[j]**4)/(freq**2 -fpe_mean[j,finite_r[5:-1]]**2)**2)
    rfpe_mean[i]=rlos[ut.find_nearest(fpe_mean[2:],int(flist[i])*1.e6)[0]]
    rfpe_std[i]=1.0*(-rlos[ut.find_nearest(np.array(fpe_mean)[2:],int(flist[i])*1.e6)[0]]+rlos[ut.find_nearest(np.array(fpe_mean)[2:]+1.0*np.array(fpe_std)[2:],int(flist[i])*1.e6)[0]])
    rflux_los_mwa[i]=rlos[ut.find_nearest(flux_los_r2[i],mwa_flux[i])[0]]
    rflux_los_fpe[i]=rlos[ut.find_nearest(np.array(fpe_mean)[2:],int(flist[i])*1.e6)[0]-1]
    flux_los_fpe1[i]=flux_los_r2[i][ut.find_nearest(np.array(fpe_mean)[2:],int(flist[i])*1.e6)[0]-1]
flux_los_r2=np.array(flux_los_r2)
fpe_mean=np.array(fpe_mean)
fpe_std=np.array(fpe_std)

plot_hist_fpe=1
if(plot_hist_fpe):
    plt.hist(fpe_hist[10]/1.e6, linewidth=3, histtype='step',label=str(np.round(rlos[10],2))+' R$_{\odot}$')
    plt.hist(fpe_hist[70]/1.e6, linewidth=3, histtype='step',label=str(np.round(rlos[70],2))+' R$_{\odot}$')
    plt.hist(fpe_hist[250]/1.e6, linewidth=3, histtype='step',label=str(np.round(rlos[250],2))+' R$_{\odot}$')
    plt.hist(fpe_hist[399]/1.e6, linewidth=3, histtype='step',label=str(np.round(rlos[399],2))+' R$_{\odot}$')
    plt.legend()
    plt.xlabel('Plasma Frequency (MHz)')
    plt.show()

plot_hist_fpe=0
if(plot_hist_fpe):
    plt.hist(densfpe_hist[0], histtype='step',label=str(rlos[0])+' R$_{\odot}$')
    plt.hist(densfpe_hist[10], histtype='step',label=str(rlos[10])+' R$_{\odot}$')
    plt.hist(densfpe_hist[100], histtype='step',label=str(rlos[100])+' R$_{\odot}$')
    plt.hist(densfpe_hist[-1], histtype='step',label=str(rlos[-1])+' R$_{\odot}$')
    plt.legend()
    plt.xlabel('Plasma Density (cm$^{-3}$)')
    plt.show()

plot_flux_los=1
if(plot_flux_los):
    colors = pl.cm.jet(np.linspace(0,1,8))[::-1]
    for k in range(8):
        #plt.plot(r[50,50],flux_los[k][-1].sum(axis=(0,1))-flux_los[k,:].sum(axis=(1,2)),color=colors[k],label=str(flist[k])+' MHz')
        plt.plot(rflux_los_mwa[k],mwa_flux[k],'s',color=colors[k],markersize=20)
        plt.errorbar(rflux_los_mwa[k],flux_los_mwa[k],yerr=emwa_flux[k],color=colors[k])
        plt.plot(rflux_los_fpe[k],flux_los_fpe1[k],'*',color=colors[k],markersize=20)
        plt.errorbar(rflux_los_fpe[k],flux_los_fpe1[k],xerr=rfpe_std[k],color=colors[k])
        plt.plot(rlos[:-1],flux_los_r2[k],'-',markersize=5,color=colors[k],label=str(flist[k])+' MHz')
    plt.xlabel('Radial Coordinate (R$_{\odot}$)')
    plt.ylabel('Flux density (SFU)')
    plt.legend()
    plt.ylim([0,21])
    #plt.xlim([1,1.1])
    plt.show()

plot_rflux_los=0
if(plot_rflux_los):
    plt.plot(flist,np.array(rflux_los)[:,50,50],'o-',color='k')
    #plt.errorbar(flist,np.array(rflux_los)[:,50,50],yerr=np.array(rflux_los)[:,50,50]*emwa_flux/mwa_flux,color='k')
    plt.ylim([1.0,1.05])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Radial Coordinate (R$_{\odot}$)')
    plt.show()

plot_image_radial=0
if(plot_image_radial):
    rfpe_240=signal.convolve(np.array(rfpe)[-1],beam,mode='same')/np.sum(beam)
    rfpe_108=signal.convolve(np.array(rfpe)[0],beam,mode='same')/np.sum(beam)
    #plt.imshow(np.array(rflux_los)[0],aspect='auto',extent=[-2880,2880,-2880,2880])
    #plt.imshow(np.array(rflux_los)[0]-r[:,:,0],aspect='auto',extent=[-2880,2880,-2880,2880],vmin=0,vmax=0.02)
    fig,ax=plt.subplots(1)
    im=ax.imshow(rfpe_108-rfpe_240,origin=0,aspect='equal',extent=[-2500,2500,-2500,2500],vmin=0,vmax=0.15,cmap='YlGnBu')
    r1=(50*16)
    r2=(50*32)
    circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
    fig.colorbar(im,label='R$_{\odot}$')
    ax.add_patch(circ1)
    ax.add_patch(circ2)
    plt.xlabel('X (arcsec)')
    plt.ylabel('Y (arcsec)')
    plt.show()

## Scattering parameters

fpe=9000.*np.sqrt(densobs)
rsun2km=6.95e5 # in km
en=20
ep2byh=np.linspace(1,500,en)*1.0e-5 # in km^-1
r=rt['r3dall']
th=rt['theta3dall']
phi=rt['phi3dall']
do_los=0
if(do_los):
    tau_sc=[0]*len(flist)
    tau_sc_convolved=[0]*len(flist)
    tau_sc_sum=[0]*len(flist)
    ep2byh_freq=[0]*len(flist)
    sc_size=[0]*len(flist)
    rfpe=[0]*len(flist)
    bobs_fpe=[0]*len(flist)
    tempobs_fpe=[0]*len(flist)
    densobs_fpe=[0]*len(flist)
    beta=[0]*len(flist)
    ep2byh_freq_mean=[0]*len(flist)
    ave_fpe=[0]*len(flist)
    r_taueq1=[0]*len(flist)
    flux_fwd_fpe=[0]*len(flist)
    Tb_fwd_fpe=[0]*len(flist)
    flux_fwd_tau=[0]*len(flist)
    Tb_fwd_tau=[0]*len(flist)
    flux_fwd_fpe_conv=[0]*len(flist)
    Tb_fwd_fpe_conv=[0]*len(flist)
    Tb_fwd_fpe_mean=[0]*len(flist)
    for k in range(len(flist)):
        freq=flist[k]*1.e6
        print flist[k],' MHz'
        tau_sc[k]=[0]*r.shape[0]
        shne=[0]*r.shape[0]
        rfpe[k]=[0]*r.shape[0]
        bobs_fpe[k]=[0]*r.shape[0]
        tempobs_fpe[k]=[0]*r.shape[0]
        densobs_fpe[k]=[0]*r.shape[0]
        beta[k]=[0]*r.shape[0]
        ave_fpe[k]=[0]*r.shape[0]
        r_taueq1[k]=[0]*r.shape[0]
        flux_fwd_fpe[k]=[0]*r.shape[0]
        Tb_fwd_fpe[k]=[0]*r.shape[0]
        flux_fwd_tau[k]=[0]*r.shape[0]
        Tb_fwd_tau[k]=[0]*r.shape[0]
        Tb_fwd_fpe_mean[k]=[0]*r.shape[0]
        for i in range(r.shape[0]):
            tau_sc[k][i]=[0]*r.shape[1]
            shne[i]=[0]*r.shape[0]
            rfpe[k][i]=[0]*r.shape[0]
            bobs_fpe[k][i]=[0]*r.shape[0]
            tempobs_fpe[k][i]=[0]*r.shape[0]
            densobs_fpe[k][i]=[0]*r.shape[0]
            beta[k][i]=[0]*r.shape[0]
            ave_fpe[k][i]=[0]*r.shape[0]
            r_taueq1[k][i]=[0]*r.shape[0]
            flux_fwd_fpe[k][i]=[0]*r.shape[0]
            Tb_fwd_fpe[k][i]=[0]*r.shape[0]
            flux_fwd_tau[k][i]=[0]*r.shape[0]
            Tb_fwd_tau[k][i]=[0]*r.shape[0]
            Tb_fwd_fpe_mean[k][i]=[0]*r.shape[0]
            for j in range(r.shape[1]):
                r_taueq1[k][i][j]=r[i][j][ut.find_nearest(tau_fwd[k][i,j],1)[0]]
                ii=np.where(fpe[i,j]<freq)[0];dii=-ii[:-1]+ii[1:]
                jj=np.where(np.array(tau_fwd)[k,i,j]>1)[0];djj=-jj[:-1]+jj[1:]
                if(len(np.where(dii>1)[0])!=0):
                    ini=np.where(dii>1)[0][0]+1
                    finite_r=np.array(np.where(fpe[i,j]<freq)[0])[ini:]
                else:
                    ini=0
                    finite_r=np.array(np.where(fpe[i,j]<freq)[0])[ini:]
                templos=[0]*len(finite_r)
                #templos[0]=tempobs[i,j,finite_r[0]]
                templos[0]=0#Tb_fwd_los[k][finite_r[0]][i][j]
                for l in range(1,len(finite_r)):
                    ll=finite_r[l]
                    templos[l]=(templos[l-1]*np.exp(-1*dtaul[k][i,j,ll])+templos[l-1]*np.exp(-1*dtaur[k][i,j,ll])+tempobs[i,j,l]*(1.0-np.exp(-1*dtaul[k][i,j,ll]))+tempobs[i,j,l]*(1.0-np.exp(-1*dtaur[k][i,j,ll])))/2.0
                    #templos[l]=templos[l-1]*np.exp(-0.5*(dtaul[k][i,j,l]+dtaur[k][i,j,l]))+tempobs[i,j,l]*(1.0-np.exp(-0.5*(dtaul[k][i,j,l]+dtaur[k][i,j,l])))
                Tb_fwd_fpe[k][i][j]=templos[-1]
                #flux_fwd_fpe[k][i][j]=np.sum((flux_los[k].swapaxes(0,2)[j,i,1:]-flux_los[k].swapaxes(0,2)[j,i,:-1])[np.array(finite_r)[0]:])
                #Tb_fwd_fpe[k][i][j]=np.sum((Tb_fwd_los[k].swapaxes(0,2)[j,i,1:]-Tb_fwd_los[k].swapaxes(0,2)[j,i,:-1])[np.array(finite_r)[0]:])
                #Tb_fwd_fpe[k][i][j]=np.array(Tb_fwd_los)[k,np.array(finite_r)[0],i,j]
                #flux_fwd_tau[k][i][j]=np.sum((flux_los[k].swapaxes(0,2)[i,j,1:]-flux_los[k].swapaxes(0,2)[i,j,:-1])[np.array(finite_r_tau)[0]:])
                #Tb_fwd_tau[k][i][j]=np.array(Tb_fwd_los)[k,np.array(finite_r_tau)[0],i,j]
                dr=abs(r[i,j][:-1]-r[i,j][1:])[finite_r[5:-1]]
                tau_sc_=(np.sqrt(np.pi)/2)*((fpe[i,j,finite_r[5:-1]]**4)/(freq**2 -fpe[i,j,finite_r[5:-1]]**2)**2)
                tau_sc[k][i][j]=np.sum(tau_sc_*dr*rsun2km)
                rfpe[k][i][j]=r[i][j][finite_r]
                densobs_fpe[k][i][j]=densobs[i][j][finite_r]
                bobs_fpe[k][i][j]=bobs[i][j][finite_r]
                tempobs_fpe[k][i][j]=tempobs[i][j][finite_r]
                beta[k][i][j]=(densobs_fpe[k][i][j]*1.38e-16*tempobs_fpe[k][i][j])/(bobs_fpe[k][i][j]**2/(8*np.pi))
                shne[i][j]=(r[i][j][ut.find_nearest(densobs_fpe[k][i][j],np.max(densobs_fpe[k][i][j])/np.e)[0]]-r[i][j][ut.find_nearest(densobs_fpe[k][i][j],np.max(densobs_fpe[k][i][j]))[0]])*6.95e5
                #flux_fwd_rfpe[k][i][j]=flux_los[k,-1,i,j]-flux_los[k,ut.find_nearest(fpe[i,j],freq)[0],i,j]
                #flux_fwd_fpe[k][i][j]=flux_los[k,-1,i,j]-flux_los[k,np.array(finite_r)[0],i,j]
                #ave_fpe[i][j]=np.mean(fpe[i][j][ut.find_nearest(densarray,np.max(densarray)/np.e)[0]:ut.find_nearest(densarray,np.max(densarray))[0]])
        flux_fwd_fpe[k]=ut.Tb2flux(np.array(Tb_fwd_fpe[k]), size, size, flist[k]/1000.)
        print np.array(flux_fwd_fpe[k]).sum(),tau_sc[k][50][50]/1.e5#, np.array(flux_fwd_tau[k]).sum()
        fwhmx=int(bmin[k]*60/forward_dx)
        fwhmy=int(bmax[k]*60/forward_dy)
        beam=ut.make2DGaussian(60, fwhmx, fwhmy, center=None)
        Tb_fwd_fpe_conv[k]=signal.convolve(np.array(Tb_fwd_fpe[k]),beam,mode='same')/np.sum(beam)
    pickle.dump([rfpe,flux_fwd_fpe,Tb_fwd_fpe,Tb_fwd_fpe_conv,beta,shne,tau_sc],open('/media/rohit/MWA/20151203_MWA_NEW/20151203_flux_los.p','wb'))
        #Tb_fwd_tau_conv[k]=signal.convolve(np.array(Tb_fwd_tau[k]),beam,mode='same')/np.sum(beam)
        ##ep2byh_freq[k]=tau[0][k]/np.array(tau_sc[k])
        #flux_fwd_fpe_conv[k]=signal.convolve(np.array(flux_fwd_fpe[k]),beam,mode='same')/np.sum(beam)
        #Tb_fwd_fpe_conv[k]=signal.convolve(np.array(Tb_fwd_fpe_mean[k]),beam,mode='same')/np.sum(beam)
        #tau_sc_convolved[k]=10**(signal.convolve(np.log10(tau_sc[k]),beam, mode='same')/np.sum(beam))
        #ep2byh_freq[k]=tau_1MK[0][k]/np.array(tau_sc[k])
        #ep2byh_freq_mean[k]=np.mean(ep2byh_freq[k][35:63,38:60])
        #sc_size[k]=5.8e-9*(dr.max()*ep2byh_freq[k]*densobs[:,:,-1].mean()**2)/((np.sqrt(1-(fpe.mean()/freq)**2))*(freq/1.e6)**4)
    pickle.dump([rfpe,flux_fwd_fpe,Tb_fwd_fpe,Tb_fwd_fpe_conv,beta,shne,tau_sc],open('/media/rohit/MWA/20151203_MWA_NEW/20151203_flux_los.p','wb'))

Tb_fwd_fpe_conv=np.array(Tb_fwd_fpe_conv)

rfpe_full,flux_fwd_fpe_full,Tb_fwd_fpe_full,Tb_fwd_fpe_conv_full,beta_full,shne_full,tau_sc_full=pickle.load(open('/media/rohit/MWA/20151203_MWA_NEW/20151203_flux_los_full_los.p','rb'))
rfpe_zero,flux_fwd_fpe_zero,Tb_fwd_fpe_zero,Tb_fwd_fpe_conv_zero,beta_zero,shne_zero,tau_sc_zero=pickle.load(open('/media/rohit/MWA/20151203_MWA_NEW/20151203_flux_los_zero_los.p','rb'))
#rfpe,flux_fwd_fpe,Tb_fwd_fpe,flux_fwd_tau,Tb_fwd_tau,Tb_fwd_fpe_conv,beta,shne=pickle.load(open('/media/rohit/MWA/20151203_MWA_NEW/20151203_flux_los.p','rb'))

mwa_Tb_mean=[0]*8
fwd_Tb_mean=[0]*8
tau=[0]*8
fwd_tau_mean=[0]*8
mwa_size_tb=[0]*8
for i in range(8):
    mwa_Tb_mean[i]=Tb_all[i][0][Tb_all[i][0]!=0].mean()
    fwd_Tb_mean[i]=np.array(Tb_convolved)[i][np.array(Tb_convolved)[i]>2.e5].mean()
    tau[i]=np.log(fwd_Tb_mean[i]/(fwd_Tb_mean[i]-mwa_Tb_mean[i]))
    fwd_tau_mean[i]=np.array(tau_convolved)[i][np.array(Tb_convolved)[i]>2.e5].mean()


mwa_size_tb=[0]*len(flist)
fwd_size_tb=[0]*len(flist)
for i in range(len(flist)):
    Tb_all[i,0][np.isnan(Tb_all[i,0])]=0
    mwa_min=np.nanmin(Tb_all[i,0][Tb_all[i,0]!=0])
    mwa_max=np.nanmax(Tb_all[i,0][Tb_all[i,0]!=0])
    mwa_size_tb[i]=np.nansum(ut.get_bimage(Tb_all[i,0],2.e5/mwa_max))
    fwd_size_tb[i]=np.nansum(ut.get_bimage(Tb_fwd[i],2.e5/mwa_max))
    print mwa_max,mwa_min,mwa_size_tb[i],2.e5/mwa_max

print_Tb_new=1
if(print_Tb_new):
    print 'Frequency', '&', 'T$_{B,MWA}$',  '&', 'T$_{B,FWD}$', '&', '$\\tau_{MWA}$', '&', '$\\tau_{FWD}$', '&', '$\Phi_{MWA}$',  '&','r$_{eff}$',' & ','$\Phi_{FWD}$','&','h','\\\\'
    print  '(MHz)', ' &', ' (MK)', ' &', ' (MK)', ' &', ' &', ' &',  '($10^3$ arcmin$^2$)', ' &',  '(arcmin)',  '&',  '($10^3$ arcmin$^2$)','&','(arcmin)','\\\\'
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(mwa_Tb_mean[k]/1.e6,2),'$\pm$',np.round(mwa_Tb_mean[k]*(emwa_flux[k]/mwa_flux[k])/1.e6,2),' & ',np.round(fwd_Tb_mean[k]/1.e6,2),' & ',np.round(tau[k],3),'$\pm$',np.round(tau[k]*eTb_frac[k]*tau[k],3),' & ',np.round(fwd_tau_mean[k],2),' & ',np.round(mwa_size_tb[k]/1.e3,2),'$\pm$',np.round(np.sqrt(mwa_size_tb[k])/1.e3,2),' & ',np.round(np.sqrt((mwa_size_tb[k])/np.pi),2),'$\pm$',np.round(2*np.sqrt(mwa_size_tb[k])/1.e3,2),' & ',np.round(fwd_size_tb[k]/1.e3,2),'&',np.round(np.sqrt(np.array(mwa_size_tb[k])/np.pi)-np.sqrt(np.array(fwd_size_tb[k])/np.pi),2),'\\\\'

print_Tb_regions=1
if(print_Tb_regions):
    fqs=(np.array(Tb_qs_mwa)[6:]-np.array(Tb_qs_fwd)[6:])/np.array(Tb_qs_mwa)[6:]*100;fch=80*(np.array(Tb_qs_mwa)[6:]-np.array(Tb_ch_fwd)[6:])/np.array(Tb_ch_mwa)[6:];far=100*(np.array(Tb_ar1_mwa)[6:]-np.array(Tb_ar1_fwd)[6:])/np.array(Tb_ar1_mwa)[6:]
    print 'Parameter', '&','217 MHz',' & ','240 MHz','\\\\'
    print 'T$_{B,MWA,AR}$',' & ',np.round(np.array(Tb_ar1_mwa)[6]/1.e6,2),'$\pm$',np.round(np.array(Tb_ar1_mwa)[7]*(emwa_flux[7]/mwa_flux[7])/1.e6,2),' & ',np.round(np.array(Tb_ar1_mwa)[7]/1.e6,2),'$\pm$',np.round(np.array(Tb_ar1_mwa)[7]*(emwa_flux[7]/mwa_flux[7]/1.e6),2),'\\\\'
    print 'T$_{B,FWD,AR}$',' & ',np.round(np.array(Tb_ar1_fwd)[6]/1.e6,2),' & ',np.round(np.array(Tb_ar1_fwd)[7]/1.e6,2),'\\\\'
    print 'T$_{B,MWA,QS}$',' & ',np.round(np.array(Tb_qs_mwa)[6]/1.e6,2),'$\pm$',np.round(np.array(Tb_qs_mwa)[7]*(emwa_flux[7]/mwa_flux[7])/1.e6,2),' & ',np.round(np.array(Tb_qs_mwa)[7]/1.e6,2),'$\pm$',np.round(np.array(Tb_qs_mwa)[7]*(emwa_flux[7]/mwa_flux[7]/1.e6),2),'\\\\'
    print 'T$_{B,FWD,QS}$',' & ',np.round(np.array(Tb_qs_fwd)[6]/1.e6,2),' & ',np.round(np.array(Tb_qs_fwd)[7]/1.e6,2),'\\\\'
    print 'T$_{B,MWA,CH}$',' & ',np.round(np.array(Tb_ch_mwa)[6]/1.e6,2),'$\pm$',np.round(np.array(Tb_ch_mwa)[7]*(emwa_flux[7]/mwa_flux[7])/1.e6,2),' & ',np.round(np.array(Tb_ch_mwa)[7]/1.e6,2),'$\pm$',np.round(np.array(Tb_ch_mwa)[7]*(emwa_flux[7]/mwa_flux[7]/1.e6),2),'\\\\'
    print 'T$_{B,FWD,CH}$',' & ',np.round(np.array(Tb_ch_fwd)[6]/1.e6,2),' & ',np.round(np.array(Tb_ch_fwd)[7]/1.e6,2),'\\\\'
    print 'f_{SC,QS}',' & ',np.round(fqs[0],2),'$\pm$',np.round(fqs[0]*(emwa_flux[6]/mwa_flux[6]),2),' & ',np.round(fqs[1],2),'$\pm$',np.round(fqs[1]*(emwa_flux[7]/mwa_flux[7]),2),'\\\\'
    print 'f_{SC,CH}',' & ',np.round(fch[0],2),'$\pm$',np.round(fch[0]*(emwa_flux[6]/mwa_flux[6]),2),' & ',np.round(fch[1],2),'$\pm$',np.round(fch[1]*(emwa_flux[7]/mwa_flux[7]),2),'\\\\'
    print 'f_{SC,AR}',' & ',np.round(far[0],2),'$\pm$',np.round(far[0]*(emwa_flux[6]/mwa_flux[6]),2),' & ',np.round(far[1],2),'$\pm$',np.round(far[1]*(emwa_flux[7]/mwa_flux[7]),2),'\\\\'

sys.exit()

nlev=300
flux_mwa_bimage=[0]*8
flux_mwa_bimage_=[0]*8
flux_fwd_bimage=[0]*8
Tb_mwa_bimage=[0]*8
Tb_mwa_bimage_=[0]*8
Tb_fwd_bimage=[0]*8
fw_bsize=[0]*8
mwa_bsize=[0]*8
mwab=[0]*8
fwdb=[0]*8
for i in range(8):
    lev=np.linspace(0.80,0.99,nlev)
    flux_mwa_bimage[i]=[0]*nlev
    flux_mwa_bimage_[i]=[0]*nlev
    flux_fwd_bimage[i]=[0]*nlev
    Tb_mwa_bimage[i]=[0]*nlev
    Tb_mwa_bimage_[i]=[0]*nlev
    Tb_fwd_bimage[i]=[0]*nlev
    fw_bsize[i]=[0]*nlev
    mwa_bsize[i]=[0]*nlev
    mwab[i]=[0]*nlev
    fwdb[i]=[0]*nlev
    lev_mwa=[0]*nlev
    for j in range(nlev):
        fwdb[i][j]=ut.get_bimage(np.array(flux_fwd[i]),lev[j]) # flux_fwd is convolved
        flux_fwd_bimage[i][j]=np.array(flux_fwd[i])*fwdb[i][j]
        Tb_fwd_bimage[i][j]=np.array(Tb_convolved[i])*fwdb[i][j]
        mwab[i][j]=ut.get_bimage(flux_mwa[i],lev[j])
        fw_bsize[i][j]=50*50*np.sum(fwdb[i][j])
        flux_mwa_bimage_[i][j]=flux_mwa[i]*mwab[i][j]
        Tb_mwa_bimage_[i][j]=np.mean(Tb_all[i],axis=0)*mwab[i][j]
    flux_fwd_total=np.array(flux_fwd_bimage[i]).sum(axis=(1,2))
    flux_mwa_total=np.array(flux_mwa_bimage_[i]).sum(axis=(1,2))
    for k in range(nlev):
        nn=ut.find_nearest(flux_mwa_total,flux_fwd_total[k])[0]
        flux_mwa_bimage[i][k]=flux_mwa_bimage_[i][nn]
        Tb_mwa_bimage[i][k]=Tb_mwa_bimage_[i][nn]
        mwab[i][k]=ut.get_bimage(flux_mwa[i],lev[nn])
        mwa_bsize[i][k]=50*50*np.sum(mwab[i][k])

lev=np.array(lev)
lev_mwa=np.array(lev_mwa)
mwab=np.array(mwab)
fwdb=np.array(fwdb)
fw_bsize=np.array(fw_bsize)
mwa_bsize=np.array(mwa_bsize)
Tb_fwd_bimage=np.array(Tb_fwd_bimage)
Tb_mwa_bimage=np.array(Tb_mwa_bimage)
flux_mwa_bimage=np.array(flux_mwa_bimage)
flux_fwd_bimage=np.array(flux_fwd_bimage)


plot_10percent=0
if(plot_10percent):
    plt.plot(flist,mwa_bsize[:,0]/3600.,'o-',color='blue',label='MWA')
    plt.errorbar(flist,mwa_bsize[:,0]/3600.,yerr=np.sqrt(mwa_bsize[:,0]/3600.),color='blue')
    plt.plot(flist,fw_bsize_decon[:,0]/3600.,'o-',color='red',label='FORWARD')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Area (arcmin$^2$)')
    plt.legend()
    plt.xlim(100,250)
    plt.show()





tau=[0]*len(lev)
tau_1MK=[0]*len(lev)
ep=[0]*len(lev)
fwd_tau_mean=[0]*len(lev)
mwa_Tb_mean=[0]*len(lev)
fwd_Tb_mean=[0]*len(lev)
fwd_size1=[0]*len(lev)
mwa_size1=[0]*len(lev)
for c in range(len(lev)):
    tau[c]=[0]*len(flist)
    tau_1MK[c]=[0]*len(flist)
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
        tau_1MK[c][k]=np.log(1.12e6/(1.12e6-mwa_Tb_mean[c][k]))
        ep[c][k]=np.sqrt(tau[c][k]/fwd_tau_mean[c][k])
        #print flist[k],mwa_Tb_mean[k],fwd_Tb_mean[k],tau[k],fwd_tau_mean[k],ep[k]
        k=k+1
mwa_size1=np.array(mwa_size1)
fwd_tau_mean=np.array(fwd_tau_mean)
tau=np.array(tau)
tau_1MK=np.array(tau_1MK)
mwa_Tb_mean=np.array(mwa_Tb_mean)
fwd_Tb_mean=np.array(fwd_Tb_mean)
ep=np.array(ep)
thermal_noise_array=np.array([1.44e-5,1.375e-5,1.374e-5,1.376e-5,1.49e-5,1.799e-5,2.319e-5,2.716e-5])

sys.exit()

Tb_fwd_fpe_mean_conv=[0]*len(flist)
for k in range(len(flist)):
    fwhmx=int(bmin[k]*60/forward_dx)
    fwhmy=int(bmax[k]*60/forward_dy)
    beam=ut.make2DGaussian(60, fwhmx, fwhmy, center=None)
    Tb_fwd_fpe_mean_conv[k]=signal.convolve(np.array(Tb_fwd_fpe_mean[k]),beam,mode='same')/np.sum(beam)
Tb_fwd_fpe=np.array(Tb_fwd_fpe)
Tb_fwd_fpe_conv=np.array(Tb_fwd_fpe_conv)
ep2byh_freq=np.array(ep2byh_freq)
r_taueq1=np.array(r_taueq1)
rfpe=np.array(rfpe)
densfpe=np.array(densfpe)
bobs_rfpe=np.array(bobs_rfpe)
tempobs_rfpe=np.array(tempobs_rfpe)
densobs_rfpe=np.array(densobs_rfpe)
beta=np.array(beta)
ave_fpe=np.array(ave_fpe)
flux_fwd_rfpe=np.array(flux_fwd_rfpe)
flux_fwd_rfpe_conv=np.array(flux_fwd_rfpe_conv)

################################################################### AR1, AR2, AR3, QS, CH
rgn=['AR1','AR2','AR3','QS','CH']
# MWA
mwa_ar1xl,mwa_ar1xr,mwa_ar1yl,mwa_ar1yr=57,62,40,48
mwa_ar2xl,mwa_ar2xr,mwa_ar2yl,mwa_ar2yr=43,53,31,41
mwa_ar3xl,mwa_ar3xr,mwa_ar3yl,mwa_ar3yr=33,43,41,50
mwa_qsxl,mwa_qsxr,mwa_qsyl,mwa_qsyr=49,53,42,47
mwa_chxl,mwa_chxr,mwa_chyl,mwa_chyr=43,53,52,58
# FORWARD 
fwd_ar1xl,fwd_ar1xr,fwd_ar1yl,fwd_ar1yr=58,68,48,58
fwd_ar2xl,fwd_ar2xr,fwd_ar2yl,fwd_ar2yr=34,42,34,44
fwd_ar3xl,fwd_ar3xr,fwd_ar3yl,fwd_ar3yr=30,38,50,56
fwd_qsxl,fwd_qsxr,fwd_qsyl,fwd_qsyr=41,48,45,51
fwd_chxl,fwd_chxr,fwd_chyl,fwd_chyr=44,52,60,68

Tb_ar1_mwa=[0]*len(flist)
Tb_ar2_mwa=[0]*len(flist)
Tb_ar3_mwa=[0]*len(flist)
Tb_qs_mwa=[0]*len(flist)
Tb_ch_mwa=[0]*len(flist)
Tb_ar1_mwam=[0]*len(flist)
Tb_ar2_mwam=[0]*len(flist)
Tb_ar3_mwam=[0]*len(flist)
Tb_qs_mwam=[0]*len(flist)
Tb_ch_mwam=[0]*len(flist)
Tb_ar1_fwd=[0]*len(flist)
Tb_ar2_fwd=[0]*len(flist)
Tb_ar3_fwd=[0]*len(flist)
Tb_qs_fwd=[0]*len(flist)
Tb_ch_fwd=[0]*len(flist)
Tb_ar1_fwd_nbr=[0]*len(flist)
Tb_qs_fwd_nbr=[0]*len(flist)
Tb_ch_fwd_nbr=[0]*len(flist)
rto_ar1=[0]*len(flist)
rto_ar2=[0]*len(flist)
rto_ar3=[0]*len(flist)
rto_qs=[0]*len(flist)
rto_ch=[0]*len(flist)
tau_ar1=[0]*len(flist)
tau_ar2=[0]*len(flist)
tau_ar3=[0]*len(flist)
tau_qs=[0]*len(flist)
tau_ch=[0]*len(flist)
tau_sc_param_ar1=[0]*len(flist)
tau_sc_param_ar2=[0]*len(flist)
tau_sc_param_ar3=[0]*len(flist)
tau_sc_param_qs=[0]*len(flist)
tau_sc_param_ch=[0]*len(flist)
tau_sc_convolved=[0]*len(flist)
epbyh_ar=[0]*len(flist)
epbyh_qs=[0]*len(flist)
epbyh_ch=[0]*len(flist)
diff_max=[0]*len(flist);alpha_ar1=[0]*len(flist);alpha_qs=[0]*len(flist);alpha_ch=[0]*len(flist);Tb_mwam=[0]*len(flist)
shift=[0]*len(flist);alpha=[0]*len(flist)
ep2byh_shift=[0]*len(flist)
fwdbymwa=np.array([2.03,1.57,1.52,1.50,1.44,1.21,1.03,1.02])
Tb_fwd_fpe_conv=np.array(Tb_fwd_fpe_conv)
mm=10
for i in range(len(flist)):
    fwhmx=int(bmin[i]*60/forward_dx)
    fwhmy=int(bmax[i]*60/forward_dy)
    beam=ut.make2DGaussian(60, fwhmx, fwhmy, center=None)
    tau_sc_convolved[i]=signal.convolve(np.array(tau_sc_full[i]),beam,mode='same')/np.sum(beam)
    #####
    Tb_ar1_mwa[i]=Tb_all[i,0][mwa_ar1yl:mwa_ar1yr,mwa_ar1xl:mwa_ar1xr].mean()
    Tb_ar2_mwa[i]=Tb_all[i,0][mwa_ar2yl:mwa_ar2yr,mwa_ar2xl:mwa_ar2xr].mean()
    Tb_ar3_mwa[i]=Tb_all[i,0][mwa_ar3yl:mwa_ar3yr,mwa_ar3xl:mwa_ar3xr].mean()
    Tb_qs_mwa[i]=Tb_all[i,0][mwa_qsyl:mwa_qsyr,mwa_qsxl:mwa_qsxr].mean()
    Tb_ch_mwa[i]=Tb_all[i,0][mwa_chyl:mwa_chyr,mwa_chxl:mwa_chxr].mean()
    Tb_ar1_mwam[i]=np.nanmean(Tb_all[i,0][mwa_ar1yl-13:mwa_ar1yr+13,mwa_ar1xl-13:mwa_ar1xr+13])
    Tb_ar2_mwam[i]=np.nanmean(Tb_all[i,0][mwa_ar2yl-13:mwa_ar2yr+13,mwa_ar2xl-13:mwa_ar2xr+13])
    Tb_ar3_mwam[i]=np.nanmean(Tb_all[i,0][mwa_ar3yl-13:mwa_ar3yr+13,mwa_ar3xl-13:mwa_ar3xr+13])
    Tb_qs_mwam[i]=np.nanmean(Tb_all[i,0][mwa_qsyl-13:mwa_qsyr+13,mwa_qsxl-13:mwa_qsxr+13])
    Tb_ch_mwam[i]=np.nanmean(Tb_all[i,0][mwa_chyl-13:mwa_chyr+13,mwa_chxl-13:mwa_chxr+13])
    #####
    Tb_ar1_fwd[i]=Tb_convolved[i][fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr].max()#/fwdbymwa[i]
    Tb_ar2_fwd[i]=Tb_convolved[i][fwd_ar2yl:fwd_ar2yr,fwd_ar2xl:fwd_ar2xr].mean()#/fwdbymwa[i]
    Tb_ar3_fwd[i]=Tb_convolved[i][fwd_ar3yl:fwd_ar3yr,fwd_ar3xl:fwd_ar1xr].mean()#/fwdbymwa[i]
    Tb_qs_fwd[i]=Tb_convolved[i][fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr].mean()#/fwdbymwa[i]
    Tb_ch_fwd[i]=Tb_convolved[i][fwd_chyl:fwd_chyr,fwd_chxl:fwd_chxr].mean()#/fwdbymwa[i]
    ######
    #alpha_ar1[i]=(-Tb_ar1_mwa[i]+Tb_ar1_fwd[i])/(Tb_ar1_fwd[i])
    #alpha_qs[i]=(Tb_qs_mwa[i]-Tb_qs_fwd[i])/(Tb_qs_mwam[i]-Tb_qs_fwd[i])
    #alpha_ch[i]=(Tb_ch_mwa[i]-Tb_ch_fwd[i])/(Tb_ch_mwam[i]-Tb_qs_fwd[i])
    Tb_mwam[i]=Tb_all[i,0][np.where((Tb_all[i,0]/Tb_all[i,0].max())>0.6)].mean()
    alpha_ar1[i]=(Tb_ar1_mwa[i]-Tb_ar1_fwd[i])/(Tb_mwam[i]-Tb_ar1_fwd[i])
    alpha_qs[i]=(Tb_qs_mwa[i]-Tb_qs_fwd[i])/(Tb_mwam[i]-Tb_qs_fwd[i])
    alpha_ch[i]=(Tb_ch_mwa[i]-Tb_ch_fwd[i])/(Tb_mwam[i]-Tb_qs_fwd[i])
    ######
    rto_ar1[i]=Tb_ar1_mwa[i]/(Tb_ar1_fwd[i])#*fwdbymwa[i])
    rto_ar2[i]=Tb_ar2_mwa[i]/(Tb_ar2_fwd[i])#*fwdbymwa[i])
    rto_ar3[i]=Tb_ar3_mwa[i]/(Tb_ar3_fwd[i])#*fwdbymwa[i])
    rto_qs[i]=Tb_qs_mwa[i]/(Tb_qs_fwd[i])#*fwdbymwa[i])
    rto_ch[i]=Tb_ch_mwa[i]/(Tb_ch_fwd[i])#*fwdbymwa[i])
    ######
    Tb_ar1_fwd_nbr[i]=(Tb_convolved[i][fwd_ar1yl-mm:fwd_ar1yr+mm,fwd_ar1xl-mm:fwd_ar1xl].mean()+Tb_convolved[i][fwd_ar1yl-mm:fwd_ar1yl,fwd_ar1xl-mm:fwd_ar1xr].mean()+Tb_convolved[i][fwd_ar1yl-mm:fwd_ar1yr+mm,fwd_ar1xr:fwd_ar1xr+mm].mean())/3#/fwdbymwa[i]
    Tb_qs_fwd_nbr[i]=(Tb_convolved[i][fwd_qsyl-mm:fwd_qsyr+mm,fwd_qsxl-mm:fwd_qsxl].mean()+Tb_convolved[i][fwd_qsyl-mm:fwd_qsyl,fwd_qsxl-mm:fwd_qsxr].mean()+Tb_convolved[i][fwd_qsyl-mm:fwd_qsyr+mm,fwd_qsxr:fwd_qsxr+mm].mean())/3#/fwdbymwa[i]
    Tb_ch_fwd_nbr[i]=(Tb_convolved[i][fwd_chyl-mm:fwd_chyr+mm,fwd_chxl-mm:fwd_chxl].mean()+Tb_convolved[i][fwd_chyl-mm:fwd_chyl,fwd_chxl-mm:fwd_chxr].mean()+Tb_convolved[i][fwd_chyl-mm:fwd_chyr+mm,fwd_chxr:fwd_chxr+mm].mean())/3#/fwdbymwa[i]
    ######
    tau_ar1[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar1_fwd[i]-Tb_ar1_mwa[i]))
    tau_ar2[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar2_fwd[i]-Tb_ar2_mwa[i]))
    tau_ar3[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar3_fwd[i]-Tb_ar3_mwa[i]))
    tau_qs[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar1_fwd[i]-Tb_qs_mwa[i]))
    tau_ch[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar1_fwd[i]-Tb_ch_mwa[i]))
    #####
    tau_sc_param_ar1[i]=tau_sc_convolved[i][fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr].mean()
    tau_sc_param_ar2[i]=tau_sc_convolved[i][fwd_ar2yl:fwd_ar2yr,fwd_ar2xl:fwd_ar2xr].mean()
    tau_sc_param_ar3[i]=tau_sc_convolved[i][fwd_ar3yl:fwd_ar3yr,fwd_ar3xl:fwd_ar3xr].mean()
    tau_sc_param_qs[i]=tau_sc_convolved[i][fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr].mean()
    #tau_sc_param_qs[i]=tau_sc_convolved[i][50][50]
    tau_sc_param_ch[i]=tau_sc_convolved[i][fwd_chyl:fwd_chyr,fwd_chxl:fwd_chxr].mean()
    #####
    epbyh_ar[i]=np.log(Tb_ar1_fwd[i]/Tb_ar1_mwa[i])/tau_sc_param_ar1[i]
    epbyh_ch[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar1_fwd[i]-Tb_ch_mwa[i]))/tau_sc_param_ch[i]
    if(Tb_qs_fwd[i]<Tb_qs_mwa[i]):
        epbyh_qs[i]=np.log(Tb_ar1_fwd[i]/(Tb_ar1_fwd[i]-Tb_qs_mwa[i]))/tau_sc_param_qs[i]
    if(Tb_qs_fwd[i]>Tb_qs_mwa[i]):
        epbyh_qs[i]=np.log(Tb_ar1_fwd[i]/Tb_qs_mwa[i])/tau_sc_param_qs[i]
    #####
    (mwax,mway)=np.where(Tb_all[i][0]==np.nanmax(Tb_all[i][0]))
    (fwdx,fwdy)=np.where(Tb_fwd[i]==np.nanmax(Tb_fwd[i]))
    diff_max[i]=np.sqrt((mwax-fwdx)**2 + (mway-fwdy)**2)[0]*50./60
    mu=np.sqrt(1-(fpe[fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr].mean()/1.e6/flist[i])**2)
    fpe_ar=np.nanmean(fpe[fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr])
    delS=np.nanmean(np.array(shne_full)[fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr])# in km
    #Ne=np.nanmean(densobs_fpe[i][fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr])
    #shift[i]=5.8e-9*Ne*Ne*delS/(mu*flist[i])**4
    #ep2byh_shift[i]=(diff_max[i]*(np.pi/(60.*180)))**2*mu**4*flist[i]**4/(5.8e-9*Ne**2*delS)



plot_epbyh=1
if(plot_epbyh):
    plt.plot(flist,np.array(epbyh_ar)/5.e-7,'o-',label='AR')
    plt.plot(flist,np.array(epbyh_qs)/5.e-7,'o-',label='QS')
    plt.plot(flist,np.array(epbyh_ch)/5.e-7,'o-',label='CH')
    plt.plot(flist,np.array((np.array(epbyh_ar),np.array(epbyh_qs),np.array(epbyh_ch))).mean(axis=0)/5.e-7,'o-',label='Mean')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('$\Delta N/N$ (%)')
    plt.legend(loc=2)
    plt.show()


Tb_ar1_mwa=np.array(Tb_ar1_mwa)
Tb_ar2_mwa=np.array(Tb_ar2_mwa)
Tb_ar3_mwa=np.array(Tb_ar3_mwa)
Tb_qs_mwa=np.array(Tb_qs_mwa)
Tb_ch_mwa=np.array(Tb_ch_mwa)
Tb_ar1_fwd=np.array(Tb_ar1_fwd)
Tb_ar2_fwd=np.array(Tb_ar2_fwd)
Tb_ar3_fwd=np.array(Tb_ar3_fwd)
Tb_qs_fwd=np.array(Tb_qs_fwd)
Tb_ch_fwd=np.array(Tb_ch_fwd)

tau_r1=np.linspace(0,5,200)
tau_r2=np.linspace(0,5,200)
simTb_ar1=np.zeros((8,200,200))
simTb_qs=np.zeros((8,200,200))
simTb_ch=np.zeros((8,200,200))
tau_ar_lmt=np.zeros(8)
etau_ar_lmt=np.zeros(8)
tau_qs_lmt=np.zeros(8)
etau_qs_lmt=np.zeros(8)
tau_ch_lmt=np.zeros(8)
etau_ch_lmt=np.zeros(8)
ratio_Tb_sc_ar=np.zeros(8)
ratio_Tb_sc_qs=np.zeros(8)
ratio_Tb_sc_ch=np.zeros(8)
ep2byh_ar=np.zeros(8)
ep2byh_ch=np.zeros(8)
ep2byh_qs=np.zeros(8)
eep2byh_ar=np.zeros(8)
eep2byh_ch=np.zeros(8)
eep2byh_qs=np.zeros(8)
for k in range(8):
    for i in range(200):
        for j in range(200):
            simTb_ar1[k,i,j]=Tb_ar1_fwd[k]*(1-np.e**(-1*tau_r1[i]))+np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_r2[j]))
            simTb_qs[k,i,j]=Tb_qs_fwd[k]*(1-np.e**(-1*tau_r1[i]))+np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_r2[j]))
            simTb_ch[k,i,j]=Tb_ch_fwd[k]*(1-np.e**(-1*tau_r1[i]))+np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_r2[j]))
    tau_ar_lmt[k]=tau_r1[ut.find_nearest(simTb_ar1[k,40],Tb_ar1_mwa[k])[0]]
    etau_ar_lmt[k]=tau_ar_lmt[k]-tau_r1[ut.find_nearest(simTb_ar1[k,40],Tb_ar1_mwa[k]-0.5*Tb_ar1_mwa[k]*emwa_flux[k]/mwa_flux[k])[0]]
    tau_qs_lmt[k]=tau_r1[ut.find_nearest(simTb_qs[k,40],Tb_qs_mwa[k])[0]]
    etau_qs_lmt[k]=tau_qs_lmt[k]-tau_r1[ut.find_nearest(simTb_qs[k,40],Tb_qs_mwa[k]-0.5*Tb_qs_mwa[k]*emwa_flux[k]/mwa_flux[k])[0]] 
    tau_ch_lmt[k]=tau_r1[ut.find_nearest(simTb_ch[k,40],Tb_ch_mwa[k])[0]]
    etau_ch_lmt[k]=tau_ch_lmt[k]-tau_r1[ut.find_nearest(simTb_ch[k,40],Tb_ch_mwa[k]-0.5*Tb_ch_mwa[k]*emwa_flux[k]/mwa_flux[k])[0]]
    ratio_Tb_sc_ar[k]=np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_ar_lmt[k]))/Tb_ar1_mwa[k]
    ratio_Tb_sc_qs[k]=np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_qs_lmt[k]))/Tb_qs_mwa[k]
    ratio_Tb_sc_ch[k]=np.array(Tb_qs_fwd_nbr)[k]*(1-np.e**(-1*tau_ch_lmt[k]))/Tb_ch_mwa[k]
    ep2byh_ar[k]=tau_qs_lmt[k]/np.array(tau_sc_param_ar1)[k]
    eep2byh_ar[k]=etau_qs_lmt[k]/np.array(tau_sc_param_ar1)[k]
    ep2byh_qs[k]=tau_qs_lmt[k]/np.array(tau_sc_param_qs)[k]
    eep2byh_qs[k]=etau_qs_lmt[k]/np.array(tau_sc_param_qs)[k]
    ep2byh_ch[k]=tau_qs_lmt[k]/np.array(tau_sc_param_ch)[k]
    eep2byh_ch[k]=etau_qs_lmt[k]/np.array(tau_sc_param_ch)[k]




print 'Frequency (MHz)',' & ','$\\tau_{SC,AR}$',' & ','$\\tau_{SC,QS}$',' & ','\\tau_{SC,CH}$','\\\\'
for k in [0,1,2,3,4,5,6,7]:
    print flist[k],' & ',np.round(tau_ar_lmt[k],2),' $\pm$ ',np.round(etau_ar_lmt[k],2),' & ',np.round(tau_qs_lmt[k],2),' $\pm$ ',np.round(etau_qs_lmt[k],2),' & ',np.round(tau_ch_lmt[k],2),' $\pm$ ',np.round(etau_ch_lmt[k],2),'\\\\'


print 'Frequency (MHz)',' & ','T$_{B,MWA,AR}$ (MK)',' & ','T$_{B,MWA,QS}$ (MK)',' & ','T$_{B,MWA,CH}$ (MK)',' & ','R$_{AR1}$ ',' & ','R$_{QS}$',' & ','R$_{CH}$',' & ','Source Shift(\')','\\\\'
for k in [0,1,2,3,4,5,6,7]:
    print flist[k],' & ',np.round(Tb_ar1_mwa[k]/1.e6,2),'$\pm$',np.round(Tb_ar1_mwa[k]*(emwa_flux[k]/mwa_flux[k])/1.e6,2),' & ',np.round(Tb_qs_mwa[k]/1.e6,2),'$\pm$',np.round(Tb_qs_mwa[k]*(emwa_flux[k]/mwa_flux[k])/1.e6,2),' & ',np.round(Tb_ch_mwa[k]/1.e6,2),'$\pm$',np.round(Tb_ch_mwa[k]*(emwa_flux[k]/mwa_flux[k])/1.e6,2),' & ',np.round(rto_ar1[k],2),'$\pm$',np.round(rto_ar1[k]*(emwa_flux[k]/mwa_flux[k]),2),' & ',np.round(rto_qs[k],2),'$\pm$',np.round(rto_qs[k]*(emwa_flux[k]/mwa_flux[k]),2),' & ',np.round(rto_ch[k],2),'$\pm$',np.round(rto_ch[k]*(emwa_flux[k]/mwa_flux[k]),2),' & ',np.round(diff_max[k],2),'$\pm$',np.round(bmax[k],2),'\\\\'

print 'Frequency (MHz)',' & ','Mean $\epsilon^{2}/h$ ($\\times 10^{-5}$ km$^{-1}$)',' & ','1 $\sigma$ $\epsilon^{2} /h$ ($\\times 10^{-5}$ km$^{-1}$)',' \\\\'
for k in range(8):
    ep2h=tau_ar1[k]/np.array(tau_sc[k])[fwd_ar1yl:fwd_ar1yr,fwd_ar1xl:fwd_ar1xr].mean()
    print flist[k],' & ',np.round(ep2byh_freq[k][35:63,38:60].mean()*1.e5,2),' & ',np.round(ep2byh_freq[k][35:63,38:60].std()*1.e5,2),' \\\\'

print 'Frequency (MHz)',' & ','Mean $\epsilon^{2}/h$ ($\\times 10^{-5}$ km$^{-1}$)',' & ','1 $\sigma$ $\epsilon^{2} /h$ ($\\times 10^{-5}$ km$^{-1}$)',' \\\\'
for k in range(8):
    print flist[k],' & ',np.round(ep2byh_freq[k][35:63,38:60].mean()*1.e5,2),' & ',np.round(ep2byh_freq[k][35:63,38:60].std()*1.e5,2),' \\\\'

print 'Frequency',' & ','FORWARD', ' & ','MWA',' & ','Ratio',' & ','Shift' ,'\\\\'
print ' (MHz) ',' & ',' (in PSFs) ',' & ',' (in PSFs) ',' & ',' & ','(arcmin)','\\\\'
for k in range(8):
    print flist[k],' & ',np.round(fw_bsize[k,250]/(25*np.pi*2500/4),2),' & ',np.round(mwa_bsize[k,250]/(25*np.pi*2500/4),2), ' & ',np.round(mwa_bsize[k,250]/fw_bsize[k,250],2),' & ',np.round(diff_max[k],2),'$\pm$',np.round(bmax[k],2)

plot_ratio_paper=1
if(plot_ratio_paper):
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax1.plot(flist,np.round(np.array(Tb_ar1_mwa)/1.e6,2),'o',color='blue')
    ax1.errorbar(flist,np.array(Tb_ar1_mwa)/1.e6,yerr=np.round(np.array(Tb_ar1_mwa)*(emwa_flux/mwa_flux)/1.e6,2),color='blue',label='AR')
    ax1.plot(flist,np.round(np.array(Tb_qs_mwa)/1.e6,2),'o',color='green')
    ax1.plot(flist,np.round(np.array(Tb_ch_mwa)/1.e6,2),'o',color='red')
    ax1.errorbar(flist,np.array(Tb_qs_mwa)/1.e6,yerr=np.round(np.array(Tb_qs_mwa)*(emwa_flux/mwa_flux)/1.e6,2),color='green',label='QS')
    ax1.errorbar(flist,np.array(Tb_ch_mwa)/1.e6,yerr=np.round(np.array(Tb_ch_mwa)*(emwa_flux/mwa_flux)/1.e6,2),color='red',label='CH')
    ax2=fig.add_subplot(212,sharex=ax1)
    ax2.plot(flist,np.array(rto_ar1),'o',color='blue')
    ax2.errorbar(flist,np.array(rto_ar1),yerr=np.array(rto_ar1)*(emwa_flux/mwa_flux),color='blue',label='AR')
    ax2.plot(flist,np.array(rto_qs),'o',color='green')
    ax2.plot(flist,np.array(rto_ch),'o',color='red')
    ax2.errorbar(flist,np.array(rto_qs),yerr=np.array(rto_qs)*(emwa_flux/mwa_flux),color='green',label='QS')
    ax2.errorbar(flist,np.array(rto_ch),yerr=np.array(rto_ch)*(emwa_flux/mwa_flux),color='red',label='CH')
    ax1.legend(loc=2)
    ax2.legend(loc=2)
    ax2.set_xlabel('Frequency (MHz)')
    ax1.set_ylabel('$T_{B,MWA}$ (MK)')
    ax2.set_ylabel('$T_{B,MWA}/T_{B,FWD}$')
    ax2.set_xlim([100,250])
    plt.show()

plot_ratio=1
if(plot_ratio):
    plt.plot(flist,ratio_Tb_sc_ar,'o-',color='red',label='Active Region')
    plt.plot(flist,ratio_Tb_sc_qs,'o-',color='green',label='Quiet Sun')
    plt.plot(flist,ratio_Tb_sc_ch,'o-',color='blue',label='Coronal Hole')
    plt.errorbar(flist,ratio_Tb_sc_ar,yerr=ratio_Tb_sc_ar*(emwa_flux/mwa_flux),color='red')
    plt.errorbar(flist,ratio_Tb_sc_qs,yerr=ratio_Tb_sc_qs*(emwa_flux/mwa_flux),color='green')
    plt.errorbar(flist,ratio_Tb_sc_ch,yerr=ratio_Tb_sc_ch*(emwa_flux/mwa_flux),color='blue')
    plt.ylabel('$T_{B,SC}/T_{B,MWA}$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=3)
    plt.xlim([100,250])
    plt.show()

plot_tau_max=1
if(plot_tau_max):
    tau_max=np.log10(1/(1-(np.max(Tb_all[:,0],axis=(1,2))/Tb_ar1_fwd)))
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.plot(flist,tau_max,'o',color='blue')
    ax1.errorbar(flist,tau_max,yerr=np.round(tau_max*(emwa_flux/mwa_flux),2),color='blue',label='')
    ax1.set_xlabel('Frequency (MHz)')
    ax1.set_ylabel('$\\tau$')
    plt.xlim([100,250])
    plt.show()
    

plot_tau_old=0
if(plot_tau_old):
    k=5
    #aa=(simTb_qs[k]/1.e6)
    #aa[(aa>1.02) & (aa<1.12)]=np.nan
    aa=(simTb_ar1[k]/1.e6)
    #aa[(aa>1.16) & (aa<1.30)]=np.nan
    #aa=(simTb_ch[k]/1.e6)
    up=Tb_ar1_mwa[k]/1.e6 + 0.5*emwa_flux[k]/mwa_flux[k]*Tb_ar1_mwa[k]/1.e6
    down=Tb_ar1_mwa[k]/1.e6 - 0.5*emwa_flux[k]/mwa_flux[k]*Tb_ar1_mwa[k]/1.e6
    #up=Tb_qs_mwa[k]/1.e6 + 0.5*emwa_flux[k]/mwa_flux[k]*Tb_qs_mwa[k]/1.e6
    #down=Tb_qs_mwa[k]/1.e6 - 0.5*emwa_flux[k]/mwa_flux[k]*Tb_qs_mwa[k]/1.e6
    #up=Tb_ch_mwa[k]/1.e6 + 0.5*emwa_flux[k]/mwa_flux[k]*Tb_ch_mwa[k]/1.e6
    #down=Tb_ch_mwa[k]/1.e6 - 0.5*emwa_flux[k]/mwa_flux[k]*Tb_ch_mwa[k]/1.e6
    aa[(aa>down) & (aa<up)]=np.nan
    plt.imshow(aa,origin=True,extent=[0,5,0,5],interpolation='None',aspect='auto')
    plt.colorbar(label='$T_B (MK)$')
    plt.ylabel('$\\tau_{R1}$')
    plt.xlabel('$\\tau_{SC}$')
    plt.show()

############## NEW TAU

te_max=tempobs[52,67]
te_cent=tempobs[50,50]
fpe_cent=fpe[50,50]
dtau_l_max=np.array(dtaul)[:,52,67,:]
dtau_r_max=np.array(dtaur)[:,52,67,:]
dtau_l_cent=np.array(dtaul)[:,50,50,:]
dtau_r_cent=np.array(dtaur)[:,50,50,:]
tau_cent=(dtau_l_cent+dtau_r_cent)*0.5
r_cent=r[50,50];dr_cent=r_cent[1:]-r_cent[:-1]
Temean=np.mean(Tb_fwd_los,axis=(2,3))
ep_list=np.linspace(1.e-6,1.e-0,400)
ep_list=np.logspace(-6, 0, num=400)
mean_lev=[0.2,0.2,0.2,0.2,0.2,0.5,0.9,0.9]
mean_lev=[2,3,4,4,5,10,15,17]
Tb_cent=[0]*8
Tb_max=[0]*8
dtausc=[0]*8
ep_cent=[0]*8
ep_cent_up=[0]*8
for k in range(8):
    print k
    Tb_cent[k]=[0]*len(ep_list)
    dtausc[k]=[0]*len(ep_list)
    for e in range(len(ep_list)):
        Tb_cent[k][e]=np.zeros(200)
        Tb_cent[k][e][0]=te_cent[0]
        Tb_max[k]=np.zeros(200)
        Tb_max[k][0]=te_max[0]
        dtausc[k][e]=np.zeros(200)
        for l in range(1,200):
            dtausc[k][e][l]=(np.sqrt(np.pi)/2.0)*(fpe_cent[l]**4)/(2*((flist[k]*1.e6)**2 - fpe_cent[l]**2)**2)*dr_cent[l-1]*7.e5*ep_list[e]
            #dtau_r_cent[k][l]=dtau_r_cent[k][l]+dtausc[k][l];dtau_l_cent[k][l]=dtau_l_cent[k][l]+dtausc[k][l]
            Tb_max[k][l]=(Tb_max[k][l-1]*np.exp(-1*dtau_l_max[k][l])+Tb_max[k][l-1]*np.exp(-1*dtau_r_max[k][l])+ \
                    te_max[l]*(1.0-np.exp(-1*dtau_l_max[k][l]))+te_max[l]*(1.0-np.exp(-1*dtau_r_max[k][l])))/2.0 
            Tb_cent[k][e][l]=(Tb_cent[k][e][l-1]*np.exp(-1*dtau_l_cent[k][l])+Tb_cent[k][e][l-1]*np.exp(-1*dtau_r_cent[k][l])+ \
                    te_cent[l]*(1.0-np.exp(-1*dtau_l_cent[k][l]))+te_cent[l]*(1.0-np.exp(-1*dtau_r_cent[k][l])))/2.0 
            #Tb_cent[k][e][l]=Tb_cent[k][e][l]*(np.exp(-1*dtausc[k][e][l]))+mean_lev[k]*Tb_max[k][e][l]*(1.0-np.exp(-1*dtausc[k][e][l]))
            Tb_cent[k][e][l]=Tb_cent[k][e][l]*(np.exp(-1*dtausc[k][e][l]))+mean_lev[k]*Temean[k][l]*(1.0-np.exp(-1*dtausc[k][e][l]))
            #Tb_cent[k][e][l]=Tb_cent[k][e][l]*(np.exp(-1*dtausc[k][e][l]))+tempobs[40:60,40:60].mean(axis=(0,1))[l]*(1.0-np.exp(-1*dtausc[k][e][l]))
    ep_cent[k]=ep_list[ut.find_nearest(np.array(Tb_cent[k])[:,-1][:250],Tb_all[k,0,50,50])[0]]
    ep_cent_up[k]=ep_list[ut.find_nearest(np.array(Tb_cent[k])[:,-1][:250],Tb_all[k,0,50,50]+0.5*Tb_all[k,0,50,50]*(emwa_flux[k]/mwa_flux[k]))[0]]
ep_cent=np.array(ep_cent)
Tb_cent=np.array(Tb_cent)
Tb_max=np.array(Tb_max)
dtausc=np.array(dtausc)
ep_error=(ep_cent_up-ep_cent)
#ep_error[np.where(ep_error==0)]=1e-5
ep_cent_up=np.array(list(ep_cent[:5])+list(ep_cent[6:]))
ep_error_up=np.array(list(ep_error[:5])+list(ep_error[6:]))
flist_up=np.array(list(flist[:5])+list(flist[6:]))
z=np.polyfit(np.log10(flist), np.log10(ep_cent), 1, full=True)
z1=np.polyfit(np.log10(flist), np.log10(ep_cent), 1)
p = np.poly1d(z1);yfit=10**p(np.log10(np.arange(300)+100))

plot_epsilon=1
if(plot_epsilon):
    colors = cm.rainbow(np.linspace(0, 1, len(flist)))[::-1]
    for i in range(8):
        plt.plot(ep_list/1.e-5,Tb_cent[i][:,-1]/1.e6,color=colors[i],linewidth=3,label=str(flist[i])+' MHz')
        plt.errorbar(ep_cent[i]/1.e-5,Tb_all[i,0,50,50]/1.e6,yerr=0.5*Tb_all[i,0,50,50]*(emwa_flux[i]/mwa_flux[i])/1.e6,fmt='o',color=colors[i])
    plt.legend(loc=2,prop={'size':12})
    plt.xlabel('$\epsilon^2/h (\\times10^{-5})$ $km^{-1}$')
    plt.ylabel('$T_B$ (MK)')
    plt.xscale('log')
    plt.show()

plot_Tb=1
if(plot_Tb):
    colors = cm.rainbow(np.linspace(0, 1, len(flist)))[::-1]
    plt.plot(r[50,50],te_cent/1.e6,color='k',linewidth=3,label='T$_{e,LOS}$')
    for i in range(8):
        plt.plot(r[50,50],0.60*mean_lev[i]*Temean[i]/1.e6,color=colors[i],linewidth=3,label=str(flist[i])+' MHz')
    plt.legend(loc='upper center',ncol=5,prop={'size':14})
    plt.xlabel('$r_{LOS} (R_{\odot})$')
    plt.ylabel('Temperature (MK)')
    plt.xlim([1,1.5]);plt.ylim([0,2])
    #plt.xscale('log')
    plt.show()

plot_epsilon1=1
if(plot_epsilon1):
    plt.plot(flist,ep_cent/1.e-5,'-',color='k')
    plt.errorbar(flist,ep_cent/1.e-5,yerr=ep_error/1.e-5,fmt='o',color='k')
    #plt.plot(flist_up,ep_cent_up/1.e-5,'-',color='k')
    #plt.errorbar(flist_up,ep_cent_up/1.e-5,yerr=ep_error_up/1.e-5,fmt='o',color='k')
    plt.plot(np.arange(300)+100,yfit/1.e-5,'--',label='$\\epsilon^{2}/h = (5\\times10^{-9}) f^{(1.74\\pm1.41)}$ $km^{-1}$')
    plt.ylabel('$\epsilon^2/h (\\times10^{-5})$ $km^{-1}$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    plt.ylim([0.1,11])
    plt.yscale('log')
    plt.show()

plot_epsilon1=1
if(plot_epsilon1):
    ep=np.sqrt(ep_cent*40);eep=np.sqrt(ep_cent*40)*2*ep_error/ep_cent
    plt.plot(flist,ep,'-',color='k')
    plt.errorbar(flist,ep,yerr=eep,fmt='o',color='k')
    plt.plot(np.arange(300)+100,np.sqrt(yfit*40),'--',label='$\\epsilon = (4.4\\times10^{-4}) f^{(0.87\\pm0.69)}$')
    plt.ylabel('$\epsilon$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    #plt.ylim([0.1,11])
    plt.yscale('log')
    plt.show()

tausc_qs=np.log((Tb_qs_fwd-np.array(Tb_ar1_fwd))/(Tb_qs_mwa-np.array(Tb_ar1_fwd)))
tausc_ch=np.log((Tb_ch_fwd-np.array(Tb_ar1_fwd))/(Tb_ch_mwa-np.array(Tb_ar1_fwd)))
epbyh_qs=tausc_qs/np.array(tau_sc_param_qs)
epbyh_ch=tausc_qs/np.array(tau_sc_param_ch)
plot_tau=1
if(plot_tau):
    plt.plot(flist,tausc_qs,'o-',color='green')
    plt.plot(flist,tausc_ch,'o-',color='blue')
    plt.errorbar(flist,tausc_qs,yerr=tausc_qs*(emwa_flux/mwa_flux),color='green',label='QS')
    plt.errorbar(flist,tausc_ch,yerr=tausc_ch*(emwa_flux/mwa_flux),color='blue',label='CH')
    plt.ylabel('$\\tau_{SC}$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=1)
    plt.xlim([100,250])
    plt.show()

plot_frac_sc=1
if(plot_frac_sc):
    frac_qs=np.array(Tb_ar1_fwd)*(1-np.e**(-1*np.array(tausc_qs)))/Tb_qs_mwa
    frac_ch=np.array(Tb_ar1_fwd)*(1-np.e**(-1*np.array(tausc_ch)))/Tb_ch_mwa
    plt.plot(flist,frac_qs,'o-',color='green')
    plt.plot(flist,frac_ch,'o-',color='blue')
    plt.errorbar(flist,frac_qs,yerr=frac_qs*(emwa_flux/mwa_flux),color='green',label='QS')
    plt.errorbar(flist,frac_ch,yerr=frac_ch*(emwa_flux/mwa_flux),color='blue',label='CH')
    plt.ylabel('$T_{B,SC}/T_{B}$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=1)
    plt.xlim([100,250])
    plt.show()
    

frac_qs=np.array(Tb_ar1_fwd)*(1-np.e**(-1*tausc_qs))/Tb_qs_mwa
frac_ch=np.array(Tb_ar1_fwd)*(1-np.e**(-1*tausc_ch))/Tb_qs_mwa
print 'Frequency (MHz)',' & ','$<\epsilon^{2}/h>_{AR}$ ($\\times 10^{-5}$ km$^{-1}$)',' & ','$<\epsilon^{2}/h>_{QS}$ ($\\times 10^{-5}$ km$^{-1}$)',' & ','$<\epsilon^{2}/h>_{CH}$ ($\\times 10^{-5}$ km$^{-1}$)','\\\\'
for k in [6,7]:
    print flist[k],' & ',np.round(ep2byh_ar[k]/1.e-5,3),' $\pm$ ',np.round(eep2byh_ar[k]/1.e-5,3),' & ',np.round(ep2byh_qs[k]/1.e-5,3),' $\pm$ ',np.round(eep2byh_qs[k]/1.e-5,3),' & ',np.round(ep2byh_ch[k]/1.e-5,1),' $\pm$ ',np.round(eep2byh_ch[k]/1.e-5,1),'\\\\'

print 'Frequency (MHz)',' & ','217 MHz',' & ','240 MHz',' \\\\'
print '\\tau_{SC}',' & ',np.round(tausc_qs[-2],2),' & ',np.round(tausc_qs[-1],2),' \\\\'
print '$T_{B,SC}/T_{B}$',' & ',np.round(frac_qs[-2],2),' $\pm$ ',np.round(frac_qs[-2]*emwa_flux[-2]/mwa_flux[-2],2),' & ',np.round(frac_qs[-1],2),' $\pm$ ',np.round(frac_qs[-1]*emwa_flux[-1]/mwa_flux[-1],2),' \\\\'
print '$<\epsilon^{2}/h>_{QS}$ ($\\times 10^{-5}$ km$^{-1}$)',' & ',np.round(ep2byh_qs[-2]/1.e-5,3),' $\pm$ ',np.round(eep2byh_qs[-2]/1.e-5,3),' & ',np.round(ep2byh_qs[-1]/1.e-5,3),' $\pm$ ',np.round(eep2byh_qs[-1]/1.e-5,3),' \\\\'

plot_epbyh=0
if(plot_epbyh):
    plt.plot(flist,epbyh_qs,'o-',color='green')
    plt.plot(flist,epbyh_ch,'o-',color='blue')
    plt.errorbar(flist,epbyh_qs,yerr=epbyh_qs*(emwa_flux/mwa_flux),color='green',label='QS')
    plt.errorbar(flist,epbyh_ch,yerr=epbyh_ch*(emwa_flux/mwa_flux),color='blue',label='CH')
    plt.ylabel('$T_{B,SC}/T_{B}$')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=1)
    plt.xlim([100,250])
    plt.show()

plot_Tb_mwa_contour=1
from surya.plot import main as spl
if(plot_Tb_mwa_contour):
    res=50
    m=4
    for i in [0,3,6]:
        levels=(np.array([0.7,0.75,0.8,0.85,0.9,0.95]))*np.nanmax(Tb_all[i][0:30].mean(axis=0))
        aa=Tb_all[i][0:30].mean(axis=0)
        aa=np.vstack((np.zeros((m,100)),aa))[int(m/2):-int(m/2)]
        aa[np.where(aa==0)]=np.nan
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-res,res,-res,res],origin='lower',vmin=mn,vmax=mx)
        im=ax.imshow(aa/1.e6,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=1.5)
        ax.contour(aa, levels,extent=[-2500,2500,-2500,2500], hold='on', colors='k',linewidths=2)
        spl.add_beam(ax,-2000, -2000,bmax[i]*60,bmin[i]*60,0)
        ax.set_xlabel('Solar X (arcsec)')
        ax.set_ylabel('Solar Y (arcsec)')
        ax.set_xlim(-2500,2500)
        ax.set_ylim(-2500,2500)
        flist[3]='162'
        ax.set_title(str(flist[i])+' MHz')
        r1=16.*60
        r2=32.*60
        circ1=plt.Circle((0.5,0.5), radius=r1, color='brown', linewidth=4,fill=False)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='brown', linewidth=4,fill=False)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.grid(True)
        fig.colorbar(im,label='(MK)')
        plt.show()

plot_Tb_fwd_contour=1
from surya.plot import main as spl
if(plot_Tb_fwd_contour):
    res=50
    m=4
    for i in [0,3,6]:
        levels=(np.array([0.7,0.75,0.8,0.85,0.9,0.95]))*np.nanmax(Tb_convolved[i])#*np.nanmax(Tb_all[i][0:30].mean(axis=0))
        aa=Tb_convolved[i]
        aa[np.where(aa<0.2e6)]=np.nan
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-res,res,-res,res],origin='lower',vmin=mn,vmax=mx)
        im=ax.imshow(aa/1.e6,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=1.5)
        ax.contour(aa, levels,extent=[-2500,2500,-2500,2500], hold='on', colors='k',linewidths=2)
        spl.add_beam(ax,-2000, -2000,bmax[i]*60,bmin[i]*60,0)
        ax.set_xlabel('Solar X (arcsec)')
        ax.set_ylabel('Solar Y (arcsec)')
        ax.set_xlim(-2500,2500)
        ax.set_ylim(-2500,2500)
        flist[3]='162'
        ax.set_title(str(flist[i])+' MHz')
        r1=16.*60
        r2=32.*60
        circ1=plt.Circle((0.5,0.5), radius=r1, color='brown', linewidth=4,fill=False)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='brown', linewidth=4,fill=False)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.grid(True)
        fig.colorbar(im,label='(MK)')
        plt.show()

plot_Tb_fwd_rfpe_contour=1
from surya.plot import main as spl
if(plot_Tb_fwd_rfpe_contour):
    levels=(np.array([0.7,0.75,0.8,0.85,0.9,0.95]))
    res=50
    m=4
    for i in [0,3,6]:
        levels=(np.array([0.7,0.75,0.8,0.85,0.9,0.95]))*np.nanmax(Tb_convolved[i])
        aa=np.array(Tb_fwd_fpe_conv_full)[i]
        aa[np.where(aa<0.2e6)]=np.nan
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-res,res,-res,res],origin='lower',vmin=mn,vmax=mx)
        im=ax.imshow(aa/1.e6,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=1.5)
        ax.contour(aa/np.nanmax(aa), levels,extent=[-2500,2500,-2500,2500], hold='on', colors='k',linewidths=2)
        spl.add_beam(ax,-2000, -2000,bmax[i]*60,bmin[i]*60,0)
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
        ax.set_xlim(-2500,2500)
        ax.set_ylim(-2500,2500)
        flist[3]='162'
        ax.set_title(str(flist[i])+' MHz')
        r1=16.*60
        r2=32.*60
        circ1=plt.Circle((0.5,0.5), radius=r1, color='brown', linewidth=4,fill=False)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='brown', linewidth=4,fill=False)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.grid(True)
        fig.colorbar(im,label='(MK)')
        plt.show()

old_loop=0
if(old_loop):
    for k in range(len(flist)):
        freq=flist[k]*1.e6
        print flist[k],' MHz'
        tau_sc[k]=[0]*en
        tau_sc_sum[k]=[0]*en
        for l in range(en):
            tau_sc[k][l]=[0]*r.shape[0]
            shne=[0]*r.shape[0]
            for i in range(r.shape[0]):
                tau_sc[k][l][i]=[0]*r.shape[1]
                shne[i]=[0]*r.shape[0]
                for j in range(r.shape[1]):
                    finite_r=np.where(fpe[i,j]<freq)[0][2:] # Ignoring first two points
                    dr=abs(r[i,j][:-1]-r[i,j][1:])[finite_r[:-1]]
                    tau_sc_=(np.sqrt(np.pi)/2)*((fpe[i,j,finite_r[:-1]]**4)/(freq**2 -fpe[i,j,finite_r[:-1]]**2)**2)
                    tau_sc[k][l][i][j]=np.sum(tau_sc_*ep2byh[l]*dr*rsun2km)
                    densarray=densobs[i,j][fpe[i,j]<freq]
                    shne[i][j]=(r[i][j][ut.find_nearest(densarray,np.max(densarray)/np.e)[0]]-r[i][j][ut.find_nearest(densarray,np.max(densarray))[0]])*6.95e5
            tau_sc_sum[k][l]=np.mean(np.array(tau_sc[k][l]))
        ep2byh_freq[k]=ep2byh[ut.find_nearest(tau_sc_sum[k],tau[0][k])[0]]
        print ep2byh_freq[k]
        sc_size[k]=5.8e-9*(dr.max()*ep2byh_freq[k]*densobs[:,:,-1].mean()**2)/((np.sqrt(1-(fpe.mean()/freq)**2))*(freq/1.e6)**4)


print_Tb=1
if(print_Tb):
    cc=0
    print 'Frequency (MHz)','&','T$_{B,MWA}$ (MK)','&','T$_{B,FWD}$ (MK)','&','$\\tau_{MWA}$','&','$\\tau_{FWD}$','&','$\\epsilon^{2}/h$','\\\\'
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(mwa_Tb_mean[cc][k]/1.e6,3),'$\pm$',np.round(mwa_Tb_mean[cc][k]*(emwa_flux[k]/mwa_flux[k])/1.e6,3),' & ',np.round(fwd_Tb_mean[cc][k]/1.e6,2),' & ',np.round(tau[cc][k],4),'$\pm$',np.round(tau[cc][k]*eTb_frac[k]*tau[cc][k],4),' & ',np.round(fwd_tau_mean[cc][k],3),' & ',np.round(ep2byh_freq[k]/1.e-5,2),'$\\times 10^{-5}$','\\\\'
print_size=1
if(print_size):
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(mwa_full_size[k]/1.e3,3),'$\pm$',np.round(np.sqrt(mwa_full_size[k])/1.e3,3),' & ',np.round(np.sqrt((mwa_full_size[k])/np.pi),3),'$\pm$',np.round(2*np.sqrt(mwa_full_size[k])/1.e3,3),'\\\\'


plot_ep2byh=1
if(plot_ep2byh):
    plt.imshow(np.arrayS(tau_sc_sum).swapaxes(0,1),origin=True,extent=[0,7,1,500],interpolation=None,aspect='auto')
    plt.plot(np.arange(6)+0.5,np.array(ep2byh_freq)[:-2]/1.0e-7,'o-',color='white')
    plt.xticks([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0],flist)
    plt.yticks([1,100,200,300,400,500],['1$\\times 10^{-7}$','1$\\times 10^{-5}$','2$\\times 10^{-5}$','3$\\times 10^{-5}$','4$\\times 10^{-5}$','5$\\times 10^{-5}$'])
    plt.colorbar(label='$\\tau$')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('$\epsilon^2$/h (km$^{-1}$)')
    plt.show()

plot_height_diff=1
if(plot_height_diff):
    qs_fwd_tau=rt['taur'][fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr,:].mean(axis=(0,1))+rt['taul'][fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr,:].mean(axis=(0,1))
    qs_fwd_Te=tempobs[fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr,:].mean(axis=(0,1))
    qs_fwd_Tb_sum=(1-np.e**(+qs_fwd_tau[:-1]-qs_fwd_tau[1:]))*qs_fwd_Te[1:]
    qs_fwd_Tb_cumsum=np.cumsum(qs_fwd_Tb_sum[::-1])
    qs_fwd_Tb=np.array(Tb_convolved)[:,fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr].mean(axis=(1,2))
    qs_fwd_r=r[fwd_qsyl:fwd_qsyr,fwd_qsxl:fwd_qsxr,:].mean(axis=(0,1))
    


tau_sc=np.array(tau_sc)
tau_sc_sum=np.array(tau_sc_sum)
ep2byh_freq=np.array(ep2byh_freq)



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
        fw_size[c][f]=np.round(len(np.where(fw_bimage[f,idx[c][f]]!=0)[0])*(50./60)**2)
        fw_size_decon[c][f]=np.round(len(np.where(fw_bimage_decon[f,idx[c][f]]!=0)[0])*(50./60)**2)
        ### Upper limit ###
        idx_up,arridx_up=ut.find_nearest(fwd_sum[f],mwa_sum[c][f]+emwa_flux[f]*0.5)
        fw_size_up=np.round(len(np.where(fw_bimage[f,idx_up]!=0)[0])*(50./60)**2)
        efw_size[c][f]=fw_size_up-fw_size[c][f]
        diff_flux[c][f]=2*abs(mwa_sum[c][f]-fwd_sum_mwa[c][f])/(mwa_sum[c][f]+fwd_sum_mwa[c][f])

mwa_size=np.array(mwa_size)
fwd_sum_mwa=np.array(fwd_sum_mwa)
mwa_sum=np.array(mwa_sum)
fw_size=np.array(fw_size)
fw_size_decon=np.array(fw_size_decon)
idx=np.array(idx)
diff_flux=np.array(diff_flux)

### Polar maps
polTb_mwa=[0]*len(flist)
polTb_fwd=[0]*len(flist)
for j in range(len(flist)):
    polTb_mwa[j],r_mwa,th_mwa=ut.cart2polar(Tb_all[j,0])
    polTb_fwd[j],r_fwd,th_fwd=ut.cart2polar(Tb_fwd[j])
    
polTb_mwa=np.array(polTb_mwa)
polTb_fwd=np.array(polTb_fwd)

plot_limb=1
if(plot_limb):
    rmwa=r_mwa[:,0]
    thmwa=th_mwa[0]*180/np.pi
    colors = cm.rainbow(np.linspace(0, 1, len(flist)))[::-1]
    f, (ax0,ax1) = plt.subplots(2, 1, figsize=(6,12))
    im0=ax0.imshow(polTb_mwa[0]/1.e6,extent=[thmwa[0],thmwa[-1],0,rmwa[-1]*50./60],aspect='auto',origin='lower',cmap='jet',vmin=0,vmax=1.0)
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im0,label='$T_B$ (MK)',cax=cax)
    ax0.axvline(x=thmwa[46],color='orange')
    ax0.set_xlabel('$\\theta$ (degrees)')
    ax0.set_ylabel('r (arcmin)')
    for i in range(len(flist)):
        ax1.plot(polTb_mwa[i][:,46]/1.e6,color=colors[i],label=str(flist[i])+' MHz')
        ax1.plot(polTb_fwd[i][:,46]/1.e6,'--',color=colors[i])
        ax1.set_xlabel('r (arcmin)')
        ax1.set_ylabel('$T_B$ (MK)')
        ax1.legend()
    plt.show()

plot_radius_240=1
if(plot_radius_240):
    plt.plot(r_mwa[:,0]*(50/60.),polTb_mwa[-1].mean(axis=1)/1.e6,'o-',color='blue',label='MWA')
    plt.plot(r_mwa[:,0]*(50/60.),polTb_fwd[-1].mean(axis=1)/1.e6,'o-',color='red',label='FORWARD')
    plt.title('240 MHz')
    plt.xlabel('Radial Coordinate (arcmin)')
    plt.ylabel('$T_B$ (MK)')
    plt.legend()
    plt.show()

plot_density=1
if(plot_density):
    plt.plot(r[43,60,2:],densobs[43,60,2:]/1.e9,'-',color='b',label='AR')
    plt.plot(rfpe[:,43,60],densfpe[:,43,60]/1.e9,'*',markersize=15,color='b')
    plt.plot(r[42,45,4:],densobs[42,45,4:]/1.e9,'-',color='g',label='QS')
    plt.plot(rfpe[:,42,45],densfpe[:,42,45]/1.e9,'*',markersize=15,color='g')
    plt.plot(r[58,52,8:],densobs[58,52,8:]/1.e9,'-',color='r',label='CH')
    plt.plot(rfpe[:,58,52],densfpe[:,58,52]/1.e9,'*',markersize=15,color='r')
    plt.plot(md.nk_freq2r(np.array(flist),1)[0],md.nk_freq2r(np.array(flist),1)[1]/1.e9,'*-',color='k',markersize=15,label='Newkirk')
    plt.xlabel('Radial distance ($R_{\odot}$)')
    plt.ylabel('Coronal Density ($\\times 10^{9}cm^{-3}$)')
    plt.legend()
    plt.show()

plot_taueq1=1
if(plot_taueq1):
    plt.plot(flist,r_taueq1[:,43,60],'o-',color='b',label='AR ($\\tau_{ff}=1$)')
    plt.plot(flist,rfpe[:,43,60],'s--',color='b',label='AR ($\\omega=\\omega_{pe}$)')
    plt.plot(flist,r_taueq1[:,42,45],'o-',color='g',label='QS ($\\tau_{ff}=1$)')
    plt.plot(flist,rfpe[:,42,45],'s--',color='g',label='QS ($\\omega=\\omega_{pe}$)')
    plt.plot(flist,r_taueq1[:,58,52],'o-',color='r',label='CH ($\\tau_{ff}=1$)')
    plt.plot(flist,rfpe[:,58,52],'s--',color='r',label='CH ($\\omega=\\omega_{pe}$)')
    plt.ylabel('Radial distance ($R_{\odot}$)')
    plt.xlabel('Frequency (MHz)')
    plt.legend(loc=2)
    plt.show()


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

flux_from_size=(2*1.23e-23*mwa_limit*(fmiss_size*((3.14159)**2/(180*60)**2))/(3.e2/np.array(flist))**2)*1.e22
diff_S=(np.array(fluxall)-mwa_flux*fringe_fact.mean(axis=0))
inc_factor=diff_S/newflux
reff_4m_flux=(rmwa/np.pi)*(np.sqrt(1+(np.pi*inc_factor*fmiss_size/rmwa**2))-1)
print_radius=1
if(print_radius):
    print 'Frequency (MHz)','&','T$_{B,limit}$','&','r$_{eff}$','&','r$_{eff,limit}$',' \\\\'
    new_radius=[0]*len(flist)
    for k in range(len(flist)):
        new_radius[k]=int(reff_4m_flux[k]+np.sqrt((mwa_full_size[k])/np.pi))
        print np.round(flist[k],1),' & ', np.round(mwa_limit[k]/1.e6,3), ' & ',np.round(reff[k],3),' $\pm$ ', np.round(np.sqrt(mwa_full_size[k])/1.e3,3)*2,' & ',np.round(reff_4m_flux[k],3),' $\pm$ ',np.round(np.sqrt((mwa_full_size[k])/1.e3)*2+np.round(reff_4m_flux[k]*eTb_frac[k],3),3),'\\\\'

##### Power law in Sizes
pw=[0]*len(flist)
epw=[0]*len(flist)
for i in range(len(flist)):
    pw_=np.polyfit(mwa_size[:,i], fw_size[:,i], 1,cov=True)
    pw[i]=pw_[0][0]
    epw[i]=np.sqrt(np.diag(pw_[1]))[0]

print_pw=1
if(print_pw):
    print 'Frequency (MHz)','&','$\\beta$','\\\\'
    for k in range(len(flist)):
        print np.round(flist[k],1),' & ',np.round(pw[k],3),' $\pm$ ', np.round(epw[k],3),'\\\\'

plot_radius_108=0
if(plot_radius_108):
    plt.plot(r_mwa[:,0],polTb_mwa[0].mean(axis=1)/1.e6,'o-',color='blue',label='MWA')
    plt.plot(r_mwa[:,0],polTb_fwd[0].mean(axis=1)/1.e6,'o-',color='red',label='FORWARD')
    plt.axhline(0.2,label='Image Limit',color='green')
    #plt.axvline(new_radius[0],label='r$_{aux,limit}$',color='orange')
    plt.title('108 MHz')
    plt.xlabel('Radial Distance (arcmin)')
    plt.ylabel('$T_B$ (MK)')
    plt.legend()
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

plot_source_size_analysed=0
if(plot_source_size_analysed):
    i_=[0,1,6,7]
    colors = cm.rainbow(np.linspace(0, 1, len(i_)))[::-1]
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
    p0=np.array(np.where(Tb_fwd[i]==np.max(Tb_fwd[i])))*50
    p1=np.array(np.where(Tb_all[i][0]==np.max(Tb_all[i][0])))*50
    xpos_fwd[i]=p0[0][0]-128*50.  
    ypos_fwd[i]=p0[1][0]-128*50.
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

plot_fluxes=1
if(plot_fluxes==1):
    plt.plot(fall,np.array(fluxall),'o-',color='red',label='FORWARD (full integration)')
    plt.plot(flist,np.array(flux_fwd_fpe_full).sum(axis=(1,2)),'s--',color='red',label='FORWARD ($f_{pe}$)')
    plt.plot(flist,np.array(flux_fwd_fpe_zero).sum(axis=(1,2)),'s-',color='red',label='FORWARD ($f_{pe}$)')
    #mwa_Tb_sum=Tb_bimage[f,idx][np.nonzero(Tb_bimage[f,idx])].mean()
    plt.errorbar(flist,mwa_flux*fringe_fact.mean(axis=0),emwa_flux,color='blue',label='MWA')
    plt.plot(flist,mwa_flux*fringe_fact.mean(axis=0),'o',color='blue')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (SFU)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    plt.show()

plot_fluxes_only=1
if(plot_fluxes_only==1):
    plt.plot(fall,np.array(fluxall),'o-',color='red',label='FORWARD')
    plt.errorbar(flist,mwa_flux*fringe_fact.mean(axis=0),emwa_flux,color='blue',label='MWA')
    plt.plot(flist,mwa_flux*fringe_fact.mean(axis=0),'o',color='blue')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux density (SFU)')
    plt.legend(loc=2)
    plt.xlim([100,250])
    plt.show()

plot_fluxes_diff=0
if(plot_fluxes_diff==1):
    del fluxall[0]
    del fluxall[1]
    plt.errorbar(flist,np.array(fluxall)-mwa_flux*fringe_fact.mean(axis=0),emwa_flux,color='k')
    plt.plot(flist,np.array(fluxall)-mwa_flux*fringe_fact.mean(axis=0),'o',color='k')
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
        import matplotlib.gridspec as gridspec
        Tb_fwd_bimage[np.isnan(Tb_fwd_bimage)]=0 
        Tb_mwa_bimage[np.isnan(Tb_mwa_bimage)]=0 
        f=7;idx=250;levels_=[0.0]
        #colors = cm.rainbow(np.linspace(0, 1, len(flist)))
        colors_=['b','cornflowerblue','g','r','c','purple','m','y']
        fig=plt.figure(figsize=(20,8))#,dpi=200)
        spec2 = gridspec.GridSpec(ncols=4, nrows=2)
        f2_ax00 = fig.add_subplot(spec2[0, 0])
        f2_ax01 = fig.add_subplot(spec2[0, 1],sharey=f2_ax00,sharex=f2_ax00);f2_ax02 = fig.add_subplot(spec2[0, 2],sharey=f2_ax00,sharex=f2_ax00);f2_ax03 = fig.add_subplot(spec2[0, 3],sharey=f2_ax00,sharex=f2_ax00)
        f2_ax10 = fig.add_subplot(spec2[1, 0],sharey=f2_ax00,sharex=f2_ax00);f2_ax11 = fig.add_subplot(spec2[1, 1],sharey=f2_ax00,sharex=f2_ax00);f2_ax12 = fig.add_subplot(spec2[1, 2],sharey=f2_ax00,sharex=f2_ax00);f2_ax13 = fig.add_subplot(spec2[1, 3],sharey=f2_ax00,sharex=f2_ax00)
        f2_ax00.contour(Tb_fwd_bimage[0,idx]/np.max(Tb_fwd_bimage[0,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax00.contour(Tb_mwa_bimage[0,idx]/np.max(Tb_mwa_bimage[0,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        #f2_ax00.text(-2100,2000,'FORWARD contour:'+str(np.round(lev[idx]*100))+'% | flux: '+str(np.round(np.sum(flux_fwd_bimage[0,idx]),3))+' SFU |'+' Area: '+str(fw_bsize[0,idx]/2500.)+' arcmin$^2$',{'color':'blue','fontsize': 15})
        #f2_ax00.text(-2100,2200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[0,idx]),3))+' SFU |'+' Area: '+str(mwa_bsize[0,idx]/2500.)+' arcmin$^2$',{'color':'red','fontsize': 15})
        f2_ax00.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[0,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax00.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[0,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax00.text(-700,0,'(A) '+str(flist[0])+' MHz',{'fontsize': 20})
        #f2_ax00.set_xlim([-1000,1400]);f2_ax00.set_ylim([-1000,600])
        f2_ax00.set_xlim([-1000,1000]);f2_ax00.set_ylim([-1000,1000])
        f2_ax00.set_ylabel('Solar Y (arcsec)')
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[0,idx]==np.nanmax(Tb_mwa_bimage[0,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[0,idx]==np.nanmax(Tb_mwa_bimage[0,idx]))[1]];f2_ax00.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[0,idx]==np.nanmax(Tb_fwd_bimage[0,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[0,idx]==np.nanmax(Tb_fwd_bimage[0,idx]))[1]];f2_ax00.scatter(yf,xf,marker='o',s=300)
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax00.add_patch(circ1)
        ############
        f2_ax01.contour(Tb_fwd_bimage[1,idx]/np.max(Tb_fwd_bimage[1,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax01.contour(Tb_mwa_bimage[1,idx]/np.max(Tb_mwa_bimage[1,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[1,idx]==np.nanmax(Tb_mwa_bimage[1,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[1,idx]==np.nanmax(Tb_mwa_bimage[1,idx]))[1]];f2_ax01.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[1,idx]==np.nanmax(Tb_fwd_bimage[1,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[1,idx]==np.nanmax(Tb_fwd_bimage[1,idx]))[1]];f2_ax01.scatter(yf,xf,marker='o',s=300)
        f2_ax01.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[1,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax01.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[1,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax01.text(-700,0,'(B) '+str(flist[1])+' MHz',{'fontsize': 20})
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax01.add_patch(circ1)
        ############
        f2_ax02.contour(Tb_fwd_bimage[2,idx]/np.max(Tb_fwd_bimage[2,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax02.contour(Tb_mwa_bimage[2,idx]/np.max(Tb_mwa_bimage[2,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[2,idx]==np.nanmax(Tb_mwa_bimage[2,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[2,idx]==np.nanmax(Tb_mwa_bimage[2,idx]))[1]];f2_ax02.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[2,idx]==np.nanmax(Tb_fwd_bimage[2,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[2,idx]==np.nanmax(Tb_fwd_bimage[2,idx]))[1]];f2_ax02.scatter(yf,xf,marker='o',s=300)
        f2_ax02.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[2,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax02.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[2,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax02.text(-700,0,'(C) '+str(flist[2])+' MHz',{'fontsize': 20})
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax02.add_patch(circ1)
        ############
        f2_ax03.contour(Tb_fwd_bimage[3,idx]/np.max(Tb_fwd_bimage[3,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax03.contour(Tb_mwa_bimage[3,idx]/np.max(Tb_mwa_bimage[3,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[3,idx]==np.nanmax(Tb_mwa_bimage[3,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[3,idx]==np.nanmax(Tb_mwa_bimage[3,idx]))[1]];f2_ax03.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[3,idx]==np.nanmax(Tb_fwd_bimage[3,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[3,idx]==np.nanmax(Tb_fwd_bimage[3,idx]))[1]];f2_ax03.scatter(yf,xf,marker='o',s=300)
        f2_ax03.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[3,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax03.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[3,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax03.text(-700,0,'(D) '+str(flist[3])+' MHz',{'fontsize': 20})
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax03.add_patch(circ1)
        ############
        f2_ax10.contour(Tb_fwd_bimage[4,idx]/np.max(Tb_fwd_bimage[4,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax10.contour(Tb_mwa_bimage[4,idx]/np.max(Tb_mwa_bimage[4,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[4,idx]==np.nanmax(Tb_mwa_bimage[4,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[4,idx]==np.nanmax(Tb_mwa_bimage[4,idx]))[1]];f2_ax10.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[4,idx]==np.nanmax(Tb_fwd_bimage[4,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[4,idx]==np.nanmax(Tb_fwd_bimage[4,idx]))[1]];f2_ax10.scatter(yf,xf,marker='o',s=300)
        f2_ax10.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[4,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax10.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[4,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax10.text(-700,0,'(E) '+str(flist[4])+' MHz',{'fontsize': 20})
        f2_ax10.set_ylabel('Solar Y (arcsec)')
        f2_ax10.set_xlabel('Solar X (arcsec)')
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax10.add_patch(circ1)
        ############
        f2_ax11.contour(Tb_fwd_bimage[5,idx]/np.max(Tb_fwd_bimage[5,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax11.contour(Tb_mwa_bimage[5,idx]/np.max(Tb_mwa_bimage[5,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[5,idx]==np.nanmax(Tb_mwa_bimage[5,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[5,idx]==np.nanmax(Tb_mwa_bimage[5,idx]))[1]];f2_ax11.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[5,idx]==np.nanmax(Tb_fwd_bimage[5,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[5,idx]==np.nanmax(Tb_fwd_bimage[5,idx]))[1]];f2_ax11.scatter(yf,xf,marker='o',s=300)
        f2_ax11.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[5,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax11.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[5,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax11.text(-700,0,'(F) '+str(flist[5])+' MHz',{'fontsize': 20})
        f2_ax11.set_xlabel('Solar X (arcsec)')
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax11.add_patch(circ1)
        ############
        f2_ax12.contour(Tb_fwd_bimage[6,idx]/np.max(Tb_fwd_bimage[6,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax12.contour(Tb_mwa_bimage[6,idx]/np.max(Tb_mwa_bimage[6,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[6,idx]==np.nanmax(Tb_mwa_bimage[6,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[6,idx]==np.nanmax(Tb_mwa_bimage[6,idx]))[1]];f2_ax12.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[6,idx]==np.nanmax(Tb_fwd_bimage[6,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[6,idx]==np.nanmax(Tb_fwd_bimage[6,idx]))[1]];f2_ax12.scatter(yf,xf,marker='o',s=300)
        f2_ax12.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[6,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax12.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[6,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax12.text(-700,0,'(G) '+str(flist[6])+' MHz',{'fontsize': 20})
        f2_ax12.set_xlabel('Solar X (arcsec)')
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax12.add_patch(circ1)
        ############
        f2_ax13.contour(Tb_fwd_bimage[7,idx]/np.max(Tb_fwd_bimage[7,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='-',linewidths=5)
        f2_ax13.contour(Tb_mwa_bimage[7,idx]/np.max(Tb_mwa_bimage[7,idx]),levels=levels_,extent=[-2500,2500,-2500,2500],linestyles='--',linewidths=5)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[7,idx]==np.nanmax(Tb_mwa_bimage[7,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_mwa_bimage[7,idx]==np.nanmax(Tb_mwa_bimage[7,idx]))[1]];f2_ax13.scatter(yf,xf,marker='*',s=300)
        xf,yf=np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[7,idx]==np.nanmax(Tb_fwd_bimage[7,idx]))[0]],np.linspace(-2500,2500,100)[np.where(Tb_fwd_bimage[7,idx]==np.nanmax(Tb_fwd_bimage[7,idx]))[1]];f2_ax13.scatter(yf,xf,marker='o',s=300)
        f2_ax13.text(-700,-200,'MWA flux: '+str(np.round(np.sum(flux_mwa_bimage[7,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax13.text(-700,-100,'FORWARD flux: '+str(np.round(np.sum(flux_fwd_bimage[7,idx]),3))+' SFU',{'fontsize': 15})
        f2_ax13.text(-700,0,'(H) '+str(flist[7])+' MHz',{'fontsize': 20})
        f2_ax13.set_xlabel('Solar X (arcsec)')
        r1=(50*16);circ1=plt.Circle((0.0,0.0), radius=r1, color='black',linestyle='-', linewidth=2,fill=False);f2_ax13.add_patch(circ1)
        ############
        plt.show()

image_dens=1
if(image_dens):
    plt.imshow(np.log10(densobs[:,:,-1]),origin=True,extent=[-2500,2500,-2500,2500],aspect='equal')
    plt.title('(a) Density')
    plt.colorbar(label='log(n$_e$) (cm$^{-3}$)')
    plt.grid(True)
    plt.xlabel('Solar X (arcsec)')
    plt.ylabel('Solar Y (arcsec)')
    plt.show()

image_temp=1
if(image_temp):
    plt.imshow(np.log10(tempobs[:,:,-1]),origin=True,extent=[-2500,2500,-2500,2500],aspect='equal')
    plt.title('(b) Temperature')
    plt.colorbar(label='log(T$_e$) (K)')
    plt.grid(True)
    plt.xlabel('Solar X (arcsec)')
    #plt.ylabel('Solar Y (arcsec)')
    plt.show()

image_bfield=1
if(image_bfield):
    plt.imshow(bobs[:,:,-1],origin=True,extent=[-2500,2500,-2500,2500],aspect='equal',vmin=0.001,vmax=3)
    plt.title('(c) Absolute Magnetic field')
    plt.colorbar(label='B$_{total}$ (G)')
    plt.grid(True)
    plt.xlabel('Solar X (arcsec)')
    #plt.ylabel('Solar Y (arcsec)')
    plt.show()
image_tau=1
if(image_tau):
    fig,ax=plt.subplots(1)
    im=ax.imshow(np.log10(tau_convolved[-1]),origin=True,extent=[-2500,2500,-2500,2500],aspect='equal')
    fig.colorbar(im,label='log$_{10}(\\tau)$')
    x,y=np.meshgrid(np.linspace(-2500,2500,100),np.linspace(-2500,2500,100))
    ax.contour(x,y,np.log10(tau_convolved[-1])/np.max(np.log10(tau_convolved[-1])), [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], hold='on', colors='k',linewidths=2,extent=[-2500,2500,-2500,2500])
    ax.set_title('(d) $\\tau$ (240 MHz)')
    ax.grid(True)
    plt.xlabel('Solar X (arcsec)')
    #plt.ylabel('Solar Y (arcsec)')
    r1=(50*16)
    r2=(50*32)
    circ1=plt.Circle((0.5,0.5), radius=r1, color='white',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='white',linestyle='--', linewidth=4,fill=False)
    ax.add_patch(circ1)
    ax.add_patch(circ2)
    plt.show()

Tb_3D=1
if(Tb_3D):
    fig = plt.figure(figsize=(8,30))
    ax = fig.gca(projection='3d')
    X,Y=np.meshgrid(np.linspace(2500,-2500,100),np.linspace(2500,-2500,100))
    for j in range(Tb_all.shape[0]):
        TT=Tb_all[j,0]
        TT[TT==0]=np.nan
        cax=ax.contourf(X,Y,TT[:,::-1]/1.e6,100,zdir='z', offset=flist[j],alpha=0.8,cmap='YlOrRd',vmin=0.6, vmax=1.2)
        ax.contour(X,Y,Tb_all[j,i][:,::-1]/np.max(Tb_all[j,i]),zdir='z',offset=flist[j],colors='k',lev=[0.4,0.6,0.8,0.9])
    ax.set_zlim(100, 250)
    ax.view_init(elev=-168, azim=53)
    ax.dist=9
    fig.colorbar(cax,label='T$_{B}$ (MK)',orientation='horizontal',norm=mpl.colors.Normalize(vmin=0.6, vmax=1.5),fraction=0.046, pad=0.05,ticks=[0.4,0.6,0.8,1.0,1.2,1.5])
    ax.set_xlabel('$Solar X (arcsec)$', fontsize=20, rotation=90)
    ax.set_ylabel('$Solar Y (arcsec)$', fontsize=20, rotation=90)
    ax.set_zlabel('$Frequency (MHz)$', fontsize=20, rotation=0)
    ax.tick_params(axis="y",direction="in", pad=-22)
    ax.tick_params(axis="x",direction="in", pad=-22)
    ax.tick_params(axis="z",direction="in", pad=-22)
    #plt.savefig('plots_3d/Tb_'+str('%03d'%i)+'.png')
    plt.show()


image_fwd_mwa=1
if(image_fwd_mwa):
    c=0
    fr1=7
    fr2=6
    levels_=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])*1.5e6
    f, ax = plt.subplots(2, 2, figsize=(20,10))
    im00=ax[0,0].imshow(Tb_bimage[fr1,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.01,vmax=1.5)
    #f.colorbar(im00,ax=ax[0,0],label='T$_{B}$(MK)')
    ax[0,0].contour(Tb_bimage[fr1,c],levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,0].set_title('240 MHz')
    ax[0,0].set_ylabel('arcsec')
    ax[0,0].grid(True)  
    r1=(50*16)
    r2=(50*32)
    circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
    ax00=ax[0,0]
    ax00.add_patch(circ1)
    ax00.add_patch(circ2)
    ax00.set_xlim([-2500,2500])
    ax00.set_ylim([-2500,2500])
    im10=ax[1,0].imshow(Tb_fwd_bimage[fr1,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.01,vmax=1.5)
    #f.colorbar(im10,ax=ax[1,0],label='T$_{B}$(MK)')
    ax[1,0].contour(Tb_fwd_bimage[fr1,c],levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[1,0].set_xlim([-2500,2500])
    ax[1,0].grid(True)
    ax[1,0].set_ylabel('arcsec')
    ax[1,0].set_xlabel('arcsec')
    circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
    ax10=ax[1,0]
    ax10.add_patch(circ2)
    ax10.add_patch(circ1)
    ax10.set_xlim([-2500,2500])
    ax10.set_ylim([-2500,2500])
    im01=ax[0,1].imshow(Tb_bimage[fr2,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.01,vmax=1.5)
    f.colorbar(im01,ax=ax[0,1],label='T$_{B}$(MK)')
    ax[0,1].contour(Tb_bimage[fr2,c],levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[0,1].set_title('217 MHz')
    ax[0,1].grid(True)
    circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
    ax01=ax[0,1]
    ax01.add_patch(circ2)
    ax01.add_patch(circ1)
    ax01.set_xlim([-2500,2500])
    ax01.set_ylim([-2500,2500])
    im11=ax[1,1].imshow(Tb_fwd_bimage[fr2,c]/1.e6,aspect='equal',cmap='YlGnBu',extent=[-2500,2500,-2500,2500],origin=True,vmin=0.01,vmax=1.5)
    f.colorbar(im11,ax=ax[1,1],label='T$_{B}$(MK)')
    ax[1,1].contour(Tb_fwd_bimage[fr2,c],levels=levels_,extent=[-2500,2500,-2500,2500],colors='red')
    ax[1,1].set_xlim([-2500,2500])
    ax[1,1].grid(True)
    ax[1,1].set_xlabel('arcsec')
    circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
    ax11=ax[1,1]
    ax11.add_patch(circ2)
    ax11.add_patch(circ1)
    ax11.set_xlim([-2500,2500])
    ax11.set_ylim([-2500,2500])
    plt.show()
    
def add_beam(ax,xcenter,ycenter,width, height,angle):
        theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
        x = 0.5 * width * np.cos(theta)
        y = 0.5 * height * np.sin(theta)
        rtheta = np.radians(angle)
        R = np.array([[np.cos(rtheta), -np.sin(rtheta)],[np.sin(rtheta),np.cos(rtheta)],])
        x, y = np.dot(R, np.array([x, y]))
        x += xcenter
        y += ycenter
        ax.fill(x, y, alpha=0.8, facecolor='yellow', edgecolor='yellow', linewidth=0.5, zorder=1)

        e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)
        ax.add_patch(e1)

def plot_mwa_regions(aa,bmaj,bmin,bpa,t,mn,mx,filename,lev_max,res):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        aa[np.where(aa==0)]=np.nan
        im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=mn,vmax=mx)
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-res,res,-res,res],origin='lower')
        n=(2*res/aa.shape[0])
        levels=(np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]))*lev_max
        x,y=np.meshgrid(np.linspace(-2500,2500,aa.shape[0]),np.linspace(-2500,2500,aa.shape[1]))
        ax.contour(x,y,aa, levels, hold='on', colors='k',linewidths=2,extent=[-2500,2500,-2500,2500])
        add_beam(ax,-40*res, -40*res,bmaj*60,bmin*60,bpa)
        ax.annotate('Active Region', xy=((64-50)*50., (45-50)*50), xytext=(20*50, 30*50),arrowprops=dict(facecolor='blue', shrink=0.05))
        ax.annotate('Quiet Sun', xy=((46-50)*50, (45-50)*50), xytext=(0*50, -40*50),arrowprops=dict(facecolor='blue', shrink=0.05))
        ax.annotate('Coronal Hole', xy=((48-50)*50., (58-50)*50), xytext=(0*res*50, 40*50),arrowprops=dict(facecolor='blue', shrink=0.05))
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
        rect1 = patches.Rectangle((56*50-2500,40*50-2500),50*8,50*8,linewidth=3,edgecolor='orange',facecolor='none')
        rect2 = patches.Rectangle((42*50-2500,42*50-2500),50*8,50*8,linewidth=3,edgecolor='orange',facecolor='none')
        rect3 = patches.Rectangle((43*50-2500,53*50-2500),50*10,50*10,linewidth=3,edgecolor='orange',facecolor='none')
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.add_patch(rect3)
        #ax.set_xlim(-40,40)
        #ax.set_ylim(-40,40)
        ax.set_title(t)
        r1=(50*16)
        r2=(50*32)
        circ1=plt.Circle((0.5,0.5), radius=r1, color='black',linestyle='--', linewidth=4,fill=False)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='black',linestyle='--', linewidth=4,fill=False)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.grid(True)
        fig.colorbar(im,label='(MK)')
        fig.savefig(filename)
        plt.show()
i=7
plot_mwa_regions(Tb_all[i][0]/1.e6,bmax[i],bmin[i],0,'',0,1.,'',np.nanmax(Tb_all[i][0])/1.e6,50)

def plot_fwd_regions(aa,bmaj,bmin,bpa,t,mn,mx,filename,lev_max,res):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='auto')
        aa[np.where(aa==0)]=np.nan
        im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=mn,vmax=mx)
        #im=ax.imshow(aa,aspect='equal',interpolation='none',extent=[-res,res,-res,res],origin='lower')
        n=(2*res/aa.shape[0])
        levels=(np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]))*lev_max
        x,y=np.meshgrid(np.linspace(-2500,2500,aa.shape[0]),np.linspace(-2500,2500,aa.shape[1]))
        ax.contour(x,y,aa, levels, hold='on', colors='k',linewidths=2,extent=[-2500,2500,-2500,2500])
        add_beam(ax,-40*res, -40*res,bmaj*60,bmin*60,bpa)
        ax.annotate('Active Region', xy=((64-50)*50., (49-50)*50), xytext=(20*50, 30*50),color='white',arrowprops=dict(facecolor='white', shrink=0.05))
        ax.annotate('Quiet Sun', xy=((46-50)*50, (45-50)*50), xytext=(0*50, -40*50),color='white',arrowprops=dict(facecolor='white', shrink=0.05))
        ax.annotate('Coronal Hole', xy=((48-50)*50., (64-50)*50), xytext=(0*res*50, 40*50),color='white',arrowprops=dict(facecolor='white', shrink=0.05))
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
        rect1 = patches.Rectangle((60*50-2500,48*50-2500),50*8,50*8,linewidth=3,edgecolor='orange',facecolor='none')
        rect2 = patches.Rectangle((42*50-2500,45*50-2500),50*8,50*8,linewidth=3,edgecolor='orange',facecolor='none')
        rect3 = patches.Rectangle((43*50-2500,58*50-2500),50*10,50*10,linewidth=3,edgecolor='orange',facecolor='none')
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.add_patch(rect3)
        #ax.set_xlim(-40,40)
        #ax.set_ylim(-40,40)
        ax.set_title(t)
        r1=(50*16)
        r2=(50*32);r3=(50*5.5)
        circ1=plt.Circle((0.5,0.5), radius=r1, color='white', linestyle='--',linewidth=4,fill=False)#,transform=ax.transAxes)
        circ2=plt.Circle((0.5,0.5), radius=r2, color='white',linestyle='--', linewidth=4,fill=False)#,transform=ax.transAxes)
        circ3=plt.Circle((150,0.5), radius=r3, color='white', linestyle='-',linewidth=4,fill=False)#,transform=ax.transAxes)
        ax.add_patch(circ2)
        ax.add_patch(circ1)
        ax.add_patch(circ3)
        ax.grid(True)
        #fig.colorbar(im,label='(MK)')
        fig.savefig(filename)
        plt.show()
i=7
plot_fwd_regions(Tb_convolved[i]/1.e6,bmax[i],bmin[i],0,'',0,1.,'',np.nanmax(Tb_fwd[i])/1.e6,50)

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

plot_composite=1
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

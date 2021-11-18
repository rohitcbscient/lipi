import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from surya.plot import main as pl
from surya.utils import main as ut
import matplotlib as mpl
import matplotlib.cm as cm

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')

def read_rhessi():
    rr=open('/home/i4ds1807205/vla_data/analysis/20120225_rhessi_corrected.txt','r')
    r=rr.readlines()
    chan=[3,6,12,25,50,100,300,800,7000,20000]
    rts=10
    rtn=len(r)-rts
    counts=[0]*(len(chan)-1)
    for i in range(len(chan)-1):
        counts[i]=[0]*rtn
        rtime=[0]*rtn
        for j in range(rtn):
            counts[i][j]=float([x for x in r[rts+j].split('    ') if x != ''][i+1])
            rtime[j]=ut.hms2sec_c(' '+[x for x in r[rts+j].split(' ') if x != ''][0])
    counts=np.array(counts)
    return rtime,counts

def read_goes():
    ff=open('/home/i4ds1807205/vla_data/analysis/20120225_goes.txt','r')
    f=ff.readlines()
    gs=36249
    ge=36850
    ng=ge-gs
    gtime=[0]*ng
    flux_l=[0]*ng
    flux_h=[0]*ng
    temp=[0]*ng
    em=[0]*ng
    for i in range(ng):
        time=f[gs+i].split('     ')[0].split(' ')[1]
        flux_l[i]= float(f[gs+i].split('     ')[1])
        flux_h[i]= float(f[gs+i].split('     ')[2])
        temp[i]= float([x for x in f[gs+i].split('    ') if x != ''][4])
        em[i]= float(f[gs+i].split('     ')[3])
        gtime[i]=ut.hms2sec_c(' '+time)
    return gtime,flux_l,flux_h


def read_ds():
    spec=[0]*3
    tim=[0]*3
    freq_=[0]*3
    bl=['7_16','7_16','16_19']
    for i in range(3):
            data=np.load('/home/i4ds1807205/vla_data/2050.1s.cal.ms.bl'+bl[i]+'.spec.npz')
            spec[i]=data['spec'][0][0]
    tim=data['tim']-data['tim'][0]
    freq_=data['freq']
    spec=np.array(spec)
    spec_ave=np.mean(spec,axis=0)
    spec_ave_freq=np.mean(spec_ave[516:-1],axis=0)
    # Flag bad data
    spec_ave_freq[np.where(spec_ave_freq<300)]=np.nan
    svtime=ut.hms2sec_c(' 20:40:00')
    vtime=svtime+np.arange(tim.shape[0])
    ####
    spec_im=[0]*spec_ave.shape[0]
    for i in range(spec_ave.shape[0]):
        vmean=np.mean(spec_ave[i,0:80])
        spec_im[i]=spec_ave[i]-vmean
        #spec_im[i][np.where(spec_im[i])<-200]=np.nan
    spec_im=np.array(spec_im)
    return vtime,svtime,spec_im,spec_ave_freq

def euv_vla(cmap,vlasubmap,xl,xr,yl,yr):
    lev=np.linspace(0.35,20.95,10)*7
    plt.imshow(cmap[40],origin=True,extent=[xl,xr,yl,yr],cmap='hot')
    plt.contour(vlasubmap[::-1,:],extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='white')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    #plt.title(t)
    plt.grid(True)
    plt.show()

def plot_centroids(xc,yc):
    xc[93]=xc[92]
    plt.figure(figsize=(18,9))
    plt.plot(xc,'o-',color='blue')
    #plt.errorbar(np.arange(len(xc)),xc,yerr=np.ones(len(xc))*2.0,color='blue')
    plt.xticks([0,40,80,120,160,200],['20:46:00','20:46:40','20:47:20','20:48:00','20:48:40','20:49:20'])
    plt.tick_params(axis='y',colors='blue')
    plt.ylabel('Centroid location X (arcsec)')
    plt.text(60+1,465,'A',color='black')
    plt.text(71,465,'B',color='black')
    plt.text(79,465,'C',color='black')
    plt.text(106,466,'D',color='black')
    plt.text(128,459,'E',color='black')
    #plt.text(152,459,'F',color='black')
    plt.xlabel('Time (HH:MM:SS UT)')
    plt.twinx()
    plt.plot(yc,'o-',color='red')
    #plt.errorbar(np.arange(len(yc)),yc,yerr=np.ones(len(yc))*2.0,color='red')
    plt.tick_params(colors='red')
    plt.ylabel('Centroid location Y (arcsec)')
    plt.grid()
    plt.show()

def euv_vla_qs_centroids(ccmap,map_qs,freq,xc,yc,xl,xr,yl,yr,bmaj,bmin,angle,t):
    lev=np.linspace(map_qs.min(),map_qs.max(),10)
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='BuGn')
    #plt.contour(map_qs,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='white')
    ss=ax.scatter(xc,yc,c=freq,s=80,cmap='spring')
    plt.colorbar(ss,label='Frequency (MHz)')
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    plt.title(t)
    ax.grid(True)

def plot_ds():
    spec_im=read_ds()[2]
    spec_im[np.where(spec_im<-50)]=np.nan
    fig,ax=plt.subplots()
    im=ax.imshow(spec_im[516:-1],cmap='jet',origin=True,vmax=95,vmin=20)
    ax.set_xlabel('Time (HH:MM UT)')
    ax.set_ylabel('Frequency (MHz)')
    ax.set_yticks(np.array([0,100,200,300,400,500]))
    ax.set_yticklabels(['1516','1616','1715','1815','1915','2015'])
    ax.set_xticks(np.array([0,300,600,900,1200]))
    ax.set_xticklabels(['20:40','20:45','20:50','20:55','21:00'])
    #plt.colorbar()
    plt.show()

def goes_vla_line_plot():
    gtime,flux_l,flux_h=read_goes()
    rtime,counts=read_rhessi()
    vtime,svtime,spec,spec_ave_freq=read_ds()
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(nrows=4,sharex=True,figsize=(14,14))
    ax0.plot(gtime,flux_l,'o-',color='blue',label='GOES X-ray (1.0-8.0 $\AA$)')
    ax0.set_ylabel('Flux (W/m$^{2}$)')
    #ax0.set_yscale('log')
    ax0.grid(True)
    ax0.legend()
    ax1.plot(gtime,flux_h,'o-',color='red',label='GOES X-ray (0.5-4.0 $\AA$)')
    ax1.set_ylabel('Flux (W/m$^{2}$)')
    #ax1.set_yscale('log')
    ax1.grid(True)
    ax1.legend()
    #ax1.plot(gtime,temp,'o-',label='Temperature')
    #ax1.set_ylabel('Temperature (MK)')
    ax2.plot(rtime,counts[0,:],'o-',label='RHESSI (3-6 keV)')
    ax2.plot(rtime,counts[1,:],'o-',label='RHESSI (6-12 keV)')
    ax2.set_ylabel('Corrected Counts')
    ax2.grid(True)
    ax2.legend()
    ax3.plot(vtime,spec_ave_freq,'o-',color='k',label='VLA (1.5-2.0 GHz)')
    ax3.set_xlabel('Time (HH:MM UT)')
    ax3.set_ylabel('Amplitude')
    ax3.grid(True)
    ax3.legend()
    ax3.set_xticks(np.array([0,300,600,900,1200])+svtime)
    ax3.set_xticklabels(['20:40','20:45','20:50','20:55','21:00'])
    fig.show()


def rhessi_vla_line_plot():
    gtime,flux_l,flux_h=read_goes()
    rtime,counts=read_rhessi()
    vtime,svtime,spec,spec_ave_freq=read_ds()
    fig, (ax0) = plt.subplots(nrows=1,sharex=True,figsize=(14,14))
    #ax1.set_yscale('log')
    ax0.grid(True)
    ax0.legend()
    #ax1.plot(gtime,temp,'o-',label='Temperature')
    #ax1.set_ylabel('Temperature (MK)')
    ax0.plot(rtime,counts[0,:]/counts[0,20],'-',label='RHESSI (3-6 KeV)')
    ax0.plot(rtime,counts[1,:]/counts[1,20],'-',label='RHESSI (6-12 KeV)')
    ax0.set_ylabel('Corrected Counts')
    ax0.grid(True)
    ax0.legend()
    ax1=ax0.twinx()
    ax1.plot(vtime,spec_ave_freq/spec_ave_freq[20],'-',color='k',label='VLA (1.5-2.0 GHz)')
    ax0.set_xlabel('Time (HH:MM UT)')
    ax0.set_ylabel('Amplitude')
    ax0.grid(True)
    ax1.legend(loc=2)
    ax1.set_xticks(np.array([0,300,600,900,1200])+svtime)
    ax1.set_xticklabels(['20:40','20:45','20:50','20:55','21:00'])
    ax1.set_yticklabels([''])
    fig.show()


def euv_vla_rhessi_qs_centroids(ccmap,map_qs,xc_err,yc_err,rhmap_l,rhmap_h,lev,freq,xc,yc,xl,xr,yl,yr,bmaj,bmin,angle,t,at,vt,rt):
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='Greys',vmin=0.01,vmax=40)
    ax.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=3,colors='blue')
    ax.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=3,colors='magenta')
    ax.text(xl+11, yr-30, '6-10 keV', color='magenta',bbox=dict(facecolor='white', alpha=0.))
    ax.text(xl+11, yr-35, '10-18 keV', color='blue',bbox=dict(facecolor='white', alpha=0.))
    ax.text(xl+34, yr-25, 'AIA: '+str(at)+' UT', color='maroon',fontweight='bold',fontsize=20,bbox=dict(facecolor='white', alpha=0.))
    ax.text(xl+34, yr-28, 'VLA: '+str(vt)+' UT', color='maroon',fontweight='bold',fontsize=20,bbox=dict(facecolor='white', alpha=0.))
    ax.text(xl+34, yr-31, 'RHESSI: '+str(rt)+' UT', color='maroon',fontweight='bold',fontsize=20,bbox=dict(facecolor='white', alpha=0.))
    cmap = mpl.cm.get_cmap('jet')
    norm = mpl.colors.Normalize(vmin=min(freq), vmax=max(freq),clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='jet')
    time_color = np.array([(mapper.to_rgba(v)) for v in freq])
    ss=plt.scatter(xc,yc,c=freq,s=90,cmap=cmap)
    for x, y, ex,ey, color in zip(xc, yc, xc_err, yc_err, time_color):
        ax.plot(x, y, 'o', color=color)
        ax.errorbar(x, y, xerr=ex, yerr=ey, lw=1, capsize=3, color=color)
    ax.set_xlabel('Solar X (arcsec)',fontsize=20)
    ax.set_ylabel('Solar Y (arcsec)',fontsize=20)
    ax.annotate('N',xy=(xl+60,yl+58),color='k')
    ax.annotate('W',xy=(xl+65,yl+52),color='k')
    ax.arrow(xl+60,yl+52,4.5,0,head_width=1.0,fc="k", ec="k", head_length=1)
    ax.arrow(xl+60,yl+52,0,4.5,head_width=1.0,fc="k", ec="k", head_length=1)
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    #ax.set_ylim(yl,yr)
    #ax.set_xlim(xl,xr)
    ax.set_ylim(330,370)
    ax.set_xlim(460,500)
    #pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    #ax.set_title(t,fontsize=20)
    ax.grid(True)
    cbar=fig.colorbar(ss,fraction=0.056, pad=0.04)
    cbar.ax.tick_params(labelsize=20);cbar.set_label(label='Frequency (MHz)',size=20)

def euv_vla_rhessi_contour(ccmap,map_qs,lev,xl,xr,yl,yr,bmaj,bmin,angle,t):
    fig=plt.figure(figsize=(8,8))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='BuGn',vmin=0.01,vmax=40)
    ax.contour(map_qs/np.max(map_qs),extent=[xl,xr,yl,yr],levels=lev,linewidths=3,colors='red')
    ax.set_xlabel('Solar X (arcsec)')
    ax.set_ylabel('Solar Y (arcsec)')
    ax.annotate('N',xy=(xl+60,yl+58),color='k')
    ax.annotate('W',xy=(xl+65,yl+52),color='k')
    ax.arrow(xl+60,yl+52,4.5,0,head_width=1.0,fc="k", ec="k", head_length=1)
    ax.arrow(xl+60,yl+52,0,4.5,head_width=1.0,fc="k", ec="k", head_length=1)
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    ax.set_title(t,fontsize=20)
    ax.grid(True)

def euv(ccmap,xl,xr,yl,yr,t):
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='sdoaia94',vmin=0.01,vmax=20)# 400 for 304; 2000 for 171;
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    plt.title(t,fontsize=20)
    ax.grid(True)

def hmi_vla_rhessi_qs_centroids(hmi,ccmap,map_qs,rhmap_l,rhmap_h,lev,lev_1,freq,xc,yc,xl,xr,yl,yr,bmaj,bmin,angle,t):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(113,aspect='auto')
    im=ax[0].imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray')
    plt.contourf(ccmap,extent=[xl,xr,yl,yr],levels=lev_1,alpha=0.4,cmap='YlGn')
    plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    plt.colorbar(ss,label='Frequency (MHz)')
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    plt.title(t)
    ax.grid(True)

def plot_spec_movie(Tb,eTb,freq):
    Tb=Tb/1.e6
    eTb=eTb/1.e6
    for i in range(Tb.shape[0]):
        ii="%03d" %i
        fig=plt.figure(1,figsize=(12,12))
        gridspec.GridSpec(2,2)
        plt.subplot2grid((2,2),(0,0),colspan=2,rowspan=1)
        plt.locator_params(axis='x',nbins=5)
        plt.locator_params(axis='y',nbins=5)
        plt.plot(Tb.mean(axis=1),'o-')
        plt.axvline(x=i,color='k')
        plt.ylabel('Temperature (MK)')
        plt.xticks([0,40,80,120,160,200],['20:46:00','20:46:40','20:47:20','20:48:00','20:48:40','20:49:20'])
        plt.xlabel('Time (HH:MM:SS)')
        plt.grid(True)
        plt.subplot2grid((2,2),(1,0),colspan=2,rowspan=2)
        plt.locator_params(axis='x',nbins=5)
        plt.locator_params(axis='y',nbins=5)
        plt.errorbar(np.arange(Tb.shape[1]),Tb[i],yerr=eTb,xerr=None)
        plt.xticks(np.arange(96)[::18],freq[::18])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Temperature (MK)')
        plt.ylim([1,8])
        plt.grid(True)
        plt.savefig('/media/rohit/VLA/20120225_sub/spec/'+'spec'+'_'+str(ii)+'.png')
        plt.close()
    



def plot_spec(Tb,eTb,freq):
    plt.rcParams['ytick.labelsize']=16;plt.rcParams['xtick.labelsize']=16;plt.rcParams['axes.labelsize']=16
    Tb=Tb/1.e6
    eTb=eTb/1.e6
    a,b,c,d,e,f=68,77,79,103+8,125+8,159+8
    fig=plt.figure(1,figsize=(9,6),constrained_layout=False)
    gs=gridspec.GridSpec(3,3)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.plot(Tb[:,90],'o-')
    ax0.set_ylabel('Temperature (MK)')
    r1,r2,r3,r4,r5,r6=60,85,95,135,145,180
    #ax0.fill_betweenx(np.linspace(1,5,10),r1,r2,facecolor='yellow',alpha=0.4)
    #ax0.fill_betweenx(np.linspace(1,5,10),r3,r4,facecolor='yellow',alpha=0.4)
    #ax0.fill_betweenx(np.linspace(1,5,10),r5,r6,facecolor='yellow',alpha=0.4)
    ax0.text(r1-1,5.6,'A',color='black')
    ax0.text(73,5.3,'B',color='black')
    ax0.text(80,4.5,'C',color='black')
    ax0.text(103,4.2,'D',color='black')
    ax0.text(124,3.4,'E',color='black')
    ax0.text(157,3.4,'F',color='black')
    ax0.set_xticks([0,40,80,120,160,200])
    ax0.set_xticklabels(['20:46:00','20:46:40','20:47:20','20:48:00','20:48:40','20:49:20'])
    ax0.set_xlabel('Time (HH:MM:SS)')
    ax10=fig.add_subplot(gs[1,0])
    ax10.errorbar(np.arange(Tb.shape[1]),Tb[a],yerr=eTb,xerr=None,label='A')
    ax10.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax10.set_ylabel('Temperature (MK)')
    ax10.set_ylim([2,8])
    ax10.text(10,7,'Burst A',fontsize=16)
    ax11=fig.add_subplot(gs[1,1])
    ax11.errorbar(np.arange(Tb.shape[1]),Tb[b],yerr=eTb,xerr=None,label='B')
    ax11.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax11.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    ax11.set_ylim([2,8])
    ax11.text(10,7,'Burst B',fontsize=16)
    ax12=fig.add_subplot(gs[1,2])
    ax12.errorbar(np.arange(Tb.shape[1]),Tb[c],yerr=eTb,xerr=None,label='C')
    ax12.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax12.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    ax12.set_ylim([2,8])
    ax12.text(10,7,'Burst C',fontsize=16)
    ax20=fig.add_subplot(gs[2,0],sharex=ax10)
    ax20.errorbar(np.arange(Tb.shape[1]),Tb[d],yerr=eTb,xerr=None,label='D')
    ax20.set_xticks(np.arange(96)[::18])
    ax20.set_xticklabels(list(np.round(freq[::18]/1000,2)))
    ax20.set_xlabel('Frequency (GHz)')
    ax20.set_ylabel('Temperature (MK)')
    ax20.set_ylim([2,8])
    ax20.text(10,7,'Burst D',fontsize=16)
    ax21=fig.add_subplot(gs[2,1],sharex=ax11) 
    ax21.errorbar(np.arange(Tb.shape[1]),Tb[e],yerr=eTb,xerr=None,label='E')
    ax21.set_xticks(np.arange(96)[::18])
    ax21.set_xticklabels(list(np.round(freq[::18]/1000,2)))
    ax21.set_xlabel('Frequency (GHz)')
    ax21.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    ax21.set_ylim([2,8])
    ax21.text(10,7,'Burst E',fontsize=16)
    ax22=fig.add_subplot(gs[2,2],sharex=ax12)
    ax22.errorbar(np.arange(Tb.shape[1]),Tb[f],yerr=eTb,xerr=None,label='F')
    ax22.set_xticks(np.arange(96)[::18])
    ax22.set_xticklabels(list(np.round(freq[::18]/1000,2)))
    ax22.set_xlabel('Frequency (GHz)')
    ax22.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    ax22.set_ylim([2,8])
    ax22.text(10,7,'Burst F',fontsize=16)
    plt.show()

def Tb(Tb,rc_mean,Tb1,freq):
    r1,r2,r3,r4,r5,r6=60,85,95,135,145,180
    fig=plt.figure(1,figsize=(35,30))
    gridspec.GridSpec(3,3)
    plt.subplot2grid((3,3),(0,0),colspan=3,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.plot(Tb/1.e6,'o-')
    plt.ylabel('Temperature (MK)')
    plt.fill_betweenx(np.linspace(1,5,10),r1,r2,facecolor='yellow',alpha=0.4)
    plt.fill_betweenx(np.linspace(1,5,10),r3,r4,facecolor='yellow',alpha=0.4)
    plt.fill_betweenx(np.linspace(1,5,10),r5,r6,facecolor='yellow',alpha=0.4)
    plt.text(r1+5,4.6,'A',color='black')
    plt.text(74,4.3,'B',color='black')
    plt.text(79,4,'C',color='black')
    plt.text(103,3.9,'D',color='black')
    plt.text(124,3.6,'E',color='black')
    plt.text(157,3.6,'F',color='black')
    plt.xticks([0,40,80,120,160,200],['20:46:00','20:46:40','20:47:20','20:48:00','20:48:40','20:49:20'])
    plt.xlabel('Time (HH:MM:SS)')
    plt.grid(True)
    ####
    plt.subplot2grid((3,3),(1,0),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(Tb1[:,r1:r2],aspect='auto',vmin=1,vmax=5)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,8,16,24],['20:47:01','20:47:09','20:47:17','20:47:25'])
    plt.colorbar(label='Temperature (MK)')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    plt.subplot2grid((3,3),(1,1),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(Tb1[:,r3:r4],aspect='auto',vmin=1,vmax=4.5)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,10,20,30],['20:47:35','20:47:45','20:47:55','20:48:05'])
    plt.colorbar(label='Temperature (MK)')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    plt.subplot2grid((3,3),(1,2),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(Tb1[:,r5:r6],aspect='auto',vmin=1,vmax=4.5)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,10,20,30],['20:48:25','20:48:35','20:48:45','20:48:55'])
    plt.colorbar(label='Temperature (MK)')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    ####
    plt.subplot2grid((3,3),(2,0),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(rc_mean[:,r1:r2],aspect='auto',vmin=8,vmax=13)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,8,16,24],['20:47:01','20:47:09','20:47:17','20:47:25'])
    plt.colorbar(label='arcsec')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    plt.subplot2grid((3,3),(2,1),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(rc_mean[:,r3:r4],aspect='auto',vmin=8,vmax=13)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,10,20,30],['20:48:25','20:48:35','20:48:45','20:48:55'])
    plt.colorbar(label='arcsec')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    plt.subplot2grid((3,3),(2,2),colspan=1,rowspan=1)
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.imshow(rc_mean[:,r5:r6],aspect='auto',vmin=8,vmax=13)
    plt.yticks(np.arange(96)[::18],freq[::18])
    plt.xticks([0,10,20,30],['20:48:25','20:48:35','20:48:45','20:48:55'])
    plt.colorbar(label='arcsec')
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Time (HH:MM:SS)')
    plt.show()

def composite_map(hmi,ccmap,lev_1,xl,xr,yl,yr):
    r1,r2,r3,r4,r5,r6=60,85,95,135,145,180
    i=0
    euv94=ccmap[0][i]
    euv131=ccmap[1][i]
    euv171=ccmap[2][i]
    euv193=ccmap[3][i]
    euv211=ccmap[4][i]
    euv304=ccmap[5][i]
    euv335=ccmap[6][i]
    euv1600=ccmap[7][i]
    euv1700=ccmap[8][i]
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    plt.contour(euv94/np.max(euv94),extent=[xl,xr,yl,yr],levels=np.linspace(0.4,0.8,2),colors='blue')
    plt.contour(euv131/np.max(euv131),extent=[xl,xr,yl,yr],levels=lev_1,colors='green')
    #plt.contour(euv171/np.max(euv171),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    plt.contour(euv304/np.max(euv304),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    plt.text(xl+3, yl+5, '94 $\AA$', color='blue',bbox=dict(facecolor='white', alpha=0.1))
    plt.text(xl+3, yl+10, '131 $\AA$', color='green',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '171 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    plt.text(xl+3, yl+15, '304 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.contour(euv193/np.max(euv193),extent=[xl,xr,yl,yr],levels=lev_1,colors='orange')
    #plt.contour(euv211/np.max(euv211),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv1600/np.max(euv1600),extent=[xl,xr,yl,yr],levels=lev_1,colors='black')
    #plt.contour(euv1700/np.max(euv1700),extent=[xl,xr,yl,yr],levels=lev_1,colors='blue')
    #plt.contourf(euv/np.max(euv),extent=[xl,xr,yl,yr],levels=lev_1,alpha=0.2,cmap='YlGn')
    #plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    #plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()


def centroid_map(hmi,ccmap,xc,yc,lev_1,xl,xr,yl,yr):
    a,b,c,d,e,f=68,77,79,103,125,159
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    #plt.contour(euv94/np.max(euv94),extent=[xl,xr,yl,yr],levels=np.linspace(0.4,0.8,2),colors='blue')
    #plt.contour(euv131/np.max(euv131),extent=[xl,xr,yl,yr],levels=lev_1,colors='green')
    #plt.contour(euv171/np.max(euv171),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv304/np.max(euv304),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.text(xl+3, yl+5, '94 $\AA$', color='blue',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+10, '131 $\AA$', color='green',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '171 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '304 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.contour(euv193/np.max(euv193),extent=[xl,xr,yl,yr],levels=lev_1,colors='orange')
    #plt.contour(euv211/np.max(euv211),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv1600/np.max(euv1600),extent=[xl,xr,yl,yr],levels=lev_1,colors='black')
    #plt.contour(euv1700/np.max(euv1700),extent=[xl,xr,yl,yr],levels=lev_1,colors='blue')
    plt.contourf(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],alpha=0.4,cmap='YlGn')
    #plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    #plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    plt.errorbar(xc[a].mean(),yc[a].mean(),xerr=xc[a].std(),yerr=yc[a].std(),elinewidth=5,c='b',label='A')
    plt.errorbar(xc[b].mean(),yc[b].mean(),xerr=xc[b].std(),yerr=yc[b].std(),elinewidth=5,c='c',label='B')
    plt.errorbar(xc[c].mean(),yc[c].mean(),xerr=xc[c].std(),yerr=yc[c].std(),elinewidth=5,c='y',label='C')
    plt.errorbar(xc[d].mean(),yc[d].mean(),xerr=xc[d].std(),yerr=yc[d].std(),elinewidth=5,c='m',label='D')
    plt.errorbar(xc[e].mean(),yc[e].mean(),xerr=xc[e].std(),yerr=yc[e].std(),elinewidth=5,c='g',label='E')
    plt.errorbar(xc[f].mean(),yc[f].mean(),xerr=xc[f].std(),yerr=yc[f].std(),elinewidth=5,c='r',label='F')
    plt.legend()
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

if __name__=='__main__':
    main();
else:
    #print 'plotting module for 20120225....'
    print('plotting module for 20120225....')



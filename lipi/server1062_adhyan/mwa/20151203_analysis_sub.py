import numpy as np
import matplotlib.pyplot as plt
import glob
import pickle
from surya.plot import main as spl
from surya.utils import main as ut


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


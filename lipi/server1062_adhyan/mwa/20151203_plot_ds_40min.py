import numpy as np
import pickle
import os
import subprocess
import glob
import itertools
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
import sys

plt.style.use('/nas08-data02/rohit/scripts/general/plt_style.py')

def ds_plot(Tsun_plot,Tn,filename,ytick_list,flist_list,x_space,tick,xlabel,e1,e2,e3,e4,vmn,vmx,cmtitle):
        f = plt.gcf()
        plt.rcParams['axes.linewidth'] = 1
        f, ax2=plt.subplots(1, 1, sharey=True)
        im2=ax2.imshow(Tsun_plot[::-1,:]/Tn,interpolation='none',aspect='auto',cmap='jet',extent=[e1,e2,e3,e4],norm=LogNorm(vmin=vmn,vmax=vmx))
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        cbar2 = plt.colorbar(im2, cax=cax2, format="%.2f",ticks=tick)
        cbar2.set_label(cmtitle,fontweight='bold')
        #ax2.set_xlim([d1,d2])
        xa = ax2.get_xaxis()
        ax2.set_xlabel(xlabel,fontweight='bold')
        ax2.set_ylabel("Frequency (MHz)",fontweight='bold')
        ax2.tick_params('both', length=5, width=2, which='major')
        ax2.tick_params('both', length=5, width=2, which='minor')
        ax2.set_yticks(ytick_list)
        ax2.set_xticks(xtick_list)
        ax2.set_yticklabels(flist_space)
        ax2.set_xticklabels(x_space)
        ax2.set_xlim(-4,116)
        plt.rc('font',weight='bold')
        #plt.savefig(filename)
        #xa.set_major_formatter(FuncFormatter(lambda x, pos:str(start + datetime.timedelta(seconds=x))))
        f.show()

def put_spaces(amp,Tbeam,flux,n):
        Tbeam_space=Tbeam
        flux_space=flux
        amp_space=amp
        for i in range(n-1):
                Tbeam_space = np.insert(Tbeam_space, 28*(i+1)+14*i, np.zeros((14,Tbeam.shape[1])), 0)
                #amp1 = np.insert(amp1, 28*(i+1)+14*i, np.zeros((14,amp1.shape[1])), 0)
                flux_space = np.insert(flux_space, 28*(i+1)+14*i, np.zeros((14,flux.shape[1])), 0)
                amp_space = np.insert(amp_space, 28*(i+1)+14*i, np.zeros((14,amp.shape[1])), 0)
        Tbeam_space[0:14,:]=Tbeam_space[14:28,:]
        Tbeam_space[84:84+14,:]=Tbeam_space[84+14:84+28,:]
        Tbeam_space[294:294+14,:]=Tbeam_space[294+14:294+28,:]
        Tbeam_space[336:336+14,:]=Tbeam_space[336+14:336+28,:]
        flux_space[0:14,:]=flux_space[14:28,:]
        flux_space[84:84+14,:]=flux_space[84+14:84+28,:]
        flux_space[294:294+14,:]=flux_space[294+14:294+28,:]
        flux_space[336:336+14,:]=flux_space[336+14:336+28,:]
        flux_space[42*5:42*5+14,:]=flux_space[42*5+14:42*5+28,:]
        return  amp_space,Tbeam_space,flux_space

def line_plot(Tbeam_space):
	x=np.arange(Tbeam_space.shape[1])*0.5
	linestyle_=['-','-','-','-','-','-','-','--','--']
	for i in range(len(freq_list)):
		plt.plot(x,Tbeam_space[(28+14)*i+1],linestyle=linestyle_[i],label=str(int(freq_list[i]))+' MHz')
	plt.xlabel('Time (sec)')
	plt.ylabel('Flux (SFU)')
	plt.legend(loc=9,ncol=5,prop={'size':11})
	plt.show()

def hist_plot(Tbeam_space,n,b):
	for i in n:
		plt.hist(Tbeam_space[(28+14)*i:(28+14)*(i+1),:].flatten(),bins=b,label=str(int(freq_list[i]))+' MHz')
	plt.xlabel('Flux (SFU)')
	plt.legend(loc=1)
	plt.show()

baseline_list=['000-008','000-009','000-010','008-009','008-010','009-010']
baseline_list_label=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
freq_list=[108.0,120.0,132.0,145.0,161.0,179.0,196.0,217.0,240.0]
flist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
file_='/nas08-data02/rohit/20151203_MWA/new_pickle/20151203_'
nt=4144
nf=28*len(flist)
central_freq,TSky,Tb,Trec,Tpickup,solidangle_str,rnmean,rnrms,Tsmean,Tsrms,smean,srms = np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list)),np.zeros(len(freq_list))
data=[0]*len(baseline_list)
sxx=[0]*len(baseline_list)
Tsxx=[0]*len(baseline_list)
rnxx=[0]*len(baseline_list)
syy=[0]*len(baseline_list)
Tsyy=[0]*len(baseline_list)
rnyy=[0]*len(baseline_list)
for j in range(len(baseline_list)):
        print baseline_list[j]
        data[j]=[0]*len(flist)
        sxx[j]=[0]*len(flist)
        Tsxx[j]=[0]*len(flist)
        rnxx[j]=[0]*len(flist)
        syy[j]=[0]*len(flist)
        Tsyy[j]=[0]*len(flist)
        rnyy[j]=[0]*len(flist)
        for i in range(len(flist)):
                data[j][i]=pickle.load(open(str(file_)+str(flist[i])+'_T'+str(baseline_list[j])+'.p','r'))
                data[j][i][5][0][np.isnan(data[j][i][5][0])] = 0
                rnxx[j][i]=data[j][i][5][0]
                rnxx[j][i]=np.vstack((rnxx[j][i][6:13],rnxx[j][i][19:26],rnxx[j][i][38:45],rnxx[j][i][51:58]))
                data[j][i][17][3][0][np.isnan(data[j][i][17][3][0])] = 0
                sxx[j][i]=data[j][i][17][3][0]
                sxx[j][i]=np.vstack((sxx[j][i][6:13],sxx[j][i][19:26],sxx[j][i][38:45],sxx[j][i][51:58]))
                data[j][i][17][6][0][np.isnan(data[j][i][17][6][0])] = 0
                Tsxx[j][i]=data[j][i][17][6][0]
                Tsxx[j][i]=np.vstack((Tsxx[j][i][6:13],Tsxx[j][i][19:26],Tsxx[j][i][38:45],Tsxx[j][i][51:58]))
                ###########################################################             
                data[j][i][5][3][np.isnan(data[j][i][5][3])] = 0
                rnyy[j][i]=data[j][i][5][3]
                rnyy[j][i]=np.vstack((rnyy[j][i][6:13],rnyy[j][i][19:26],rnyy[j][i][38:45],rnyy[j][i][51:58]))
                data[j][i][17][3][3][np.isnan(data[j][i][17][3][3])] = 0
                syy[j][i]=data[j][i][17][3][3]
                syy[j][i]=np.vstack((syy[j][i][6:13],syy[j][i][19:26],syy[j][i][38:45],syy[j][i][51:58]))
                data[j][i][17][6][3][np.isnan(data[j][i][17][6][3])] = 0
                Tsyy[j][i]=data[j][i][17][6][3]
                Tsyy[j][i]=np.vstack((Tsyy[j][i][6:13],Tsyy[j][i][19:26],Tsyy[j][i][38:45],Tsyy[j][i][51:58]))
        rnxx[j]=np.array(rnxx[j]).reshape((nf,nt))
        sxx[j]=np.array(sxx[j]).reshape((nf,nt))
        Tsxx[j]=np.array(Tsxx[j]).reshape((nf,nt))
        rnyy[j]=np.array(rnyy[j]).reshape((nf,nt))
        syy[j]=np.array(syy[j]).reshape((nf,nt))
        Tsyy[j]=np.array(Tsyy[j]).reshape((nf,nt))


amp_space,Tbeam_space,flux_space=put_spaces(rnxx[1],Tsxx[1],sxx[1],len(freq_list))
#0-11,8-21,9-22,10-23, baseline_list=['000-008','000-009','000-010','008-009','008-010','009-010']
#ytick_list=[0,28,42,70,84,112,126,154,168,196,210,238,252,280,294,322,336,364]
#flist_space=[102.25,103.50,115.00,117.50,129.00,131.50,145.00,147.50,165.00,167.50,186.00,188.50,211.00,213.50,239.25,240.50,271.25,272.50,296.00,298.50]
#flist_space=[107.52,108.8,119.04,120.32,131.84,133.12,144.64,145.92,160.00,161.28,177.92,179.2,195.84,197.12,216.32,217.6,239.36,240.64]
#x_space=['03:40','03:50','04:00','04:10','04:20','04:30','04:40']
#x_space=['04:00:40','04:01:15','04:01:50','04:02:25','04:03:00','04:03:35','04:04:10','04:04:48']
x_space=['03:24:33','03:40:18','03:41:00','03:41:42','03:42:24','03:43:06','03:43:48','03:44:32']
#x_space=['03:56:32','03:57:07','03:57:42','03:58:17','03:58:52','03:59:27','04:00:02','04:00:40']
#x_space=['03:56:32','04:00:40','04:04:48','04:08:56','04:13:04','04:17:21','04:21:20']
tick=[1,50,100,150,200,250,300,400,500]
#tick=[0.01,0.03,0.05,0.10,0.30,0.50]
xlabel='03-December-2015 (UT)'

###################### 20120903 #########################
# GO TO /Data/rohit/20130903_MWA/pickle/
dd=1
if(dd==1):
	ytick_list=[0,28,42,70,84,112,126,154,168,196,210,238,252,280,294,322,336,364]
        flist_space=[107.52,108.8,119.04,120.32,131.84,133.12,144.64,145.92,160.00,161.28,177.92,179.2,195.84,197.12,216.32,217.6,239.36,240.64]
	xtick_list=[0,16,32,48,64,80,96,112]
	x_space=['03:24:33','03:29:33','03:34:33','03:39:33','03:44:33','03:49:33','03:54:33','03:59:33']
	tick=[0.1,0.3,1,3,5,10,20]
        xlabel='03-December-2015 (HH:MM:SS UT)'
	e1,e2,e3,e4,vmn,vmx=0,xtick_list[-1],0,ytick_list[-1],1,40
        flux_space[:,592*3:592*4]=np.nan
	ds_plot(flux_space,1,'flux',ytick_list,flist_space,x_space,tick,xlabel,e1,e2,e3,e4,vmn,vmx,'(SFU)') 
        ################################# FIGURE 1 ##############################################
	#line_plot(flux_space)
	#hist_plot(flux_space,[4],600)




lines=1
if(lines):
    x1=1776;x2=2368
    fl=np.nanmean(np.array(sxx),axis=0).reshape(9,28,4144).mean(axis=1);fl[np.where(fl==0)]=np.nan
    fl[2]=np.nanmean(np.array(sxx),axis=0).reshape(9,28,4144)[2,12:25,:].mean(axis=0);fl[np.where(fl==0)]=np.nan
    fl[0][580]=np.nan;fl[0][581]=np.nan;fl[0][4132]=np.nan;fl[0][4133]=np.nan;fl[3][4132]=np.nan;fl[3][4133]=np.nan;fl[6][2948]=np.nan;fl[6][2949]=np.nan
    fl1=fl*1.0
    fl1[8][3548:]=fl1[8][3548:]*1.01;fl1[7][3548:]=fl1[7][3548:]*1.01;fl1[6][3548:]=fl1[6][3548:]*1.01
    f,ax=plt.subplots(8,2,figsize=(12,12))
    ax[0,0].plot(fl[0][0:x1],'o-',markersize=2,label='108 MHz');ax[0,1].plot(fl[0][x2:],'o-',markersize=2,label='108 MHz');ax[0,0].set_ylim(1,3);ax[0,1].set_ylim(1,3);ax[0,0].legend(loc=3,fontsize=10);ax[0,1].legend(loc=3,fontsize=10)
    ax[1,0].plot(fl[1][0:x1],'o-',markersize=2,label='120 MHz');ax[1,1].plot(fl[1][x2:],'o-',markersize=2,label='120 MHz');ax[1,0].set_ylim(2,4);ax[1,1].set_ylim(2,4);ax[1,0].legend(loc=3,fontsize=10);ax[1,1].legend(loc=3,fontsize=10)
    ax[2,0].plot(fl[2][0:x1],'o-',markersize=2,label='133 MHz');ax[2,1].plot(fl[2][x2:],'o-',markersize=2,label='133 MHz');ax[2,0].set_ylim(3,5);ax[2,1].set_ylim(3,5);ax[2,0].legend(loc=3,fontsize=10);ax[2,1].legend(loc=3,fontsize=10)
    ax[3,0].plot(fl[3][0:x1],'o-',markersize=2,label='145 MHz');ax[3,1].plot(fl[3][x2:],'o-',markersize=2,label='145 MHz');ax[3,0].set_ylim(4,6);ax[3,1].set_ylim(4,6);ax[3,0].legend(loc=3,fontsize=10);ax[3,1].legend(loc=2,fontsize=10)
    ax[4,0].plot(fl[5][0:x1],'o-',markersize=2,label='179 MHz');ax[4,1].plot(fl[5][x2:],'o-',markersize=2,label='179 MHz');ax[4,0].set_ylim(7,9);ax[4,1].set_ylim(7,9);ax[4,0].legend(loc=2,fontsize=10);ax[4,1].legend(loc=2,fontsize=10)
    ax[5,0].plot(fl1[6][0:x1],'o-',markersize=2,label='197 MHz');ax[5,1].plot(fl1[6][x2:],'o-',markersize=2,label='197 MHz');ax[5,0].set_ylim(10,12);ax[5,1].set_ylim(10,12);ax[5,0].legend(loc=2,fontsize=10);ax[5,1].legend(loc=2,fontsize=10)
    ax[6,0].plot(fl1[7][0:x1],'o-',markersize=2,label='217 MHz');ax[6,1].plot(fl1[7][x2:],'o-',markersize=2,label='217 MHz');ax[6,0].set_ylim(14,17);ax[6,1].set_ylim(14,17);ax[6,0].legend(loc=2,fontsize=10);ax[6,1].legend(loc=2,fontsize=10)
    ax[7,0].plot(fl1[8][0:x1],'o-',markersize=2,label='240 MHz');ax[7,1].plot(fl1[8][x2:],'o-',markersize=2,label='240 MHz');ax[7,0].set_ylim(16,21);ax[7,1].set_ylim(16,21);ax[7,0].legend(loc=2,fontsize=10);ax[7,1].legend(loc=2,fontsize=10)
    ax[4,0].set_ylabel('Flux Density (SFU)');ax[7,0].set_xlabel('Time (HH:MM)');ax[7,1].set_xlabel('Time (HH:MM)')
    ax[7,0].set_xticks([0,500,1000,1500]);ax[7,0].set_xticklabels(['03:24','03:28','03:32','03:36'])
    ax[7,1].set_xticks([0,500,1000,1500]);ax[7,1].set_xticklabels(['03:44','03:48','03:52','03:56'])
    ax[0,0].set_xticks([]);ax[1,0].set_xticks([]);ax[2,0].set_xticks([]);ax[3,0].set_xticks([]);ax[4,0].set_xticks([]);ax[5,0].set_xticks([]);ax[6,0].set_xticks([])
    ax[0,1].set_xticks([]);ax[1,1].set_xticks([]);ax[2,1].set_xticks([]);ax[3,1].set_xticks([]);ax[4,1].set_xticks([]);ax[5,1].set_xticks([]);ax[6,1].set_xticks([])
    ax[0,0].set_yticks([2.0,3.0]);ax[1,0].set_yticks([3.0,4.0]);ax[2,0].set_yticks([4.0,5.0]);ax[3,0].set_yticks([5.0,6.0]);ax[4,0].set_yticks([8.0,9.0]);ax[5,0].set_yticks([11.0,12.0]);ax[6,0].set_yticks([15.0,16.0]);ax[7,0].set_yticks([17.0,19.0])
    ax[0,1].set_yticks([]);ax[1,1].set_yticks([]);ax[2,1].set_yticks([]);ax[3,1].set_yticks([]);ax[4,1].set_yticks([]);ax[5,1].set_yticks([]);ax[6,1].set_yticks([]);ax[7,1].set_yticks([])
    f.subplots_adjust(wspace=0, hspace=0)
    plt.show()

sys.exit()
################### LINE PLOT #################
line=1
#jets=3:52
#jete=4:03
if(line==1):
    x_space=['03:48:16','03:50:38','03:53:00','03:55:27','03:57:49','04:00:20','04:02:42','04:04:48']
    plt.plot(sxx[0][0:8].mean(axis=0),'o-',label='101 MHz')
    plt.plot(sxx[0][28*4:28*5].mean(axis=0),'o-',label='165 MHz')
    plt.plot(sxx[0][28*9:28*10].mean(axis=0),'o-',label='298 MHz')
    plt.axvline(x=448,color='k',linewidth=5,linestyle='-')
    plt.axvline(x=1768,color='k',linewidth=5,linestyle='--')
    plt.fill_betweenx(np.linspace(0,100,8), 100, 300, facecolor='blue', alpha=0.2)
    plt.fill_betweenx(np.linspace(0,100,8), 400, 700, facecolor='blue', alpha=0.2)
    plt.fill_betweenx(np.linspace(0,100,8), 850, 1100, facecolor='blue', alpha=0.2)
    plt.fill_betweenx(np.linspace(0,100,8), 1350, 1600, facecolor='blue', alpha=0.2)
    plt.legend()
    plt.xticks(np.linspace(0,1984,8),x_space)
    plt.xlabel('Time (HH:MM:SS) UT')
    plt.ylabel('Flux (SFU)')
    plt.show()


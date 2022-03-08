import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from astropy.io import fits
from sunpy.map import Map
from matplotlib.patches import Circle
from scipy.io import readsav
from surya.utils import main as ut
import more_itertools as mit
import itertools

def extract_outliers(Tb_sub,n):
    #Tbmap=d*0;s=np.where(d!=0);s0=abs(s[0][1:]-s[0][:-1]);s1=abs(s[1][1:]-s[1][:-1]);ds=s0-s1
    #Tbmap=d*0;s=np.where(d!=0);s0=abs(s[0][1:]-s[0][:-1]);s1=abs(s[1][1:]-s[1][:-1]);ds=s0-s1
    #dsidx=np.where((ds<-2) | (ds>2))[0]+1;dsidx=np.insert(dsidx,0,0);dsidx=np.insert(dsidx,len(dsidx),len(t[0]))
    d5s=Tb_sub*1.0;n5s=Tb_sub*1.0;d5s_=d5s*0
    for i in range(d5s.shape[0]): 
        s=np.nanstd(d5s[i])
        d5s[i][np.where(d5s[i]<n*s)]=0
        n5s[i][np.where(n5s[i]>-1*n*s)]=0
        dd=d5s[i];dd[0:20]=0;dd[:,0:20]=0;dd[180:200]=0;dd[:,180:200]=0;ss=np.where(dd!=0)
        #s0_=np.sqrt(ss[0]**2 + ss[1]**2);s0=np.sort(s0_);s0_arg=np.argsort(s0_);ds=s0[1:]-s0[:-1]#dsidx=np.where(ds>2)[0]+1
        #dsidx=np.insert(dsidx,0,0);dsidx=np.insert(dsidx,len(dsidx),len(ds)-1)
        if(len(ss[0])>5):
            for j in range(len(ss[0])):
                if((dd[ss[0][j]-1,ss[1][j]-1]==0)& (dd[ss[0][j]-1,ss[1][j]]==0)& (dd[ss[0][j]-1,ss[1][j]+1]==0)& (dd[ss[0][j],ss[1][j]-1]==0)& (dd[ss[0][j],ss[1][j]+1]==0)& (dd[ss[0][j]+1,ss[1][j]-1]==0)&(dd[ss[0][j]+1,ss[1][j]]==0)&(dd[ss[0][j]+1,ss[1][j]+1]==0)):
                   dd[ss[0][j],ss[1][j]]=0
            tt=np.where(dd!=0)   
            for j in range(len(tt[0])):
                ddd=dd[tt[0][j]-2:tt[0][j]+2,tt[1][j]-2:tt[1][j]+2]
                if(np.where(ddd.flatten()!=0)[0].shape[0]<4):
                    dd[tt[0][j]-2:tt[0][j]+2,tt[1][j]-2:tt[1][j]+2]=0
        else:
            dd=dd*0
        d5s_[i]=dd
    return d5s_,n5s
                    
def remove_outliers(Tb,std,l,h):
    idx=np.where((std<h) & (std>l))
    outtb=Tb[idx]
    return outtb,idx[0]



forward_100MHz=readsav('/nas08-data02/rohit/20151203_forward/20151203_100MHz_psimas.sav')
pres=forward_100MHz['modsolstruct'][0][0]
dens=forward_100MHz['modsolstruct'][0][1]
temp=forward_100MHz['modsolstruct'][0][2]
br=forward_100MHz['modsolstruct'][0][3]
bth=forward_100MHz['modsolstruct'][0][4]
bph=forward_100MHz['modsolstruct'][0][5]
vr=forward_100MHz['modsolstruct'][0][6]
vth=forward_100MHz['modsolstruct'][0][7]
vph=forward_100MHz['modsolstruct'][0][8]
babs=np.sqrt(br*br+bth*bth+bph*bph)
vabs=np.sqrt(vr*vr+vth*vth+vph*vph)
rt=readsav('/nas08-data02/rohit/20151203_forward/RT_params_108MHz.sav')
r=rt['r3dall']
tempobs=rt['tempall']
densobs=rt['densall']
brobs=rt['brall']
bthobs=rt['bthall']
bphobs=rt['bphall']
bobs=np.sqrt(brobs*brobs+bthobs*bthobs+bphobs*bphobs)



filepath='/nas08-data02/rohit/20151203_MWA/Tb_new/'
Tb240,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_187-188.p','rb'))
Tb217,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_169-170.p','rb'))
Tb179,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_139-140.p','rb'))
Tb197,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_153-154.p','rb'))
Tb145,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_113-114.p','rb'))
Tb133,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_103-104.p','rb'))
Tb108,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_084-085.p','rb'))
Tb120,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_093-094.p','rb'))
Tb240=np.array(Tb240);Tb217=np.array(Tb217);Tb179=np.array(Tb179);Tb145=np.array(Tb145);Tb197=np.array(Tb197);Tb108=np.array(Tb108[0:4010]);Tb120=np.array(Tb120);Tb133=np.array(Tb133)
Tb240m=np.nanmean(Tb240[1200:],axis=0);Tb217m=np.nanmean(Tb217[0:500],axis=0);Tb179m=np.nanmean(Tb179[0:500],axis=0);Tb145m=np.nanmean(Tb145[0:500],axis=0)
Tb197m=np.nanmean(Tb197[0:500],axis=0);Tb108m=np.nanmean(Tb108[0:500],axis=0);Tb120m=np.nanmean(Tb120[0:500],axis=0);Tb133m=np.nanmean(Tb133[0:500],axis=0)


baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
freq_filelist=['084-085','093-094','103-104','113-114','139-140','153-154','169-170','187-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
j=7
Tb_write_path='/nas08-data02/rohit/20151203_sub_run_mean_subamp/running_median'
ids='20151203_'
#freq=freq[j]*1.e6#240*1.e6
freq=[108.0,120.0,133.0,145.0,179.0,197.0,217.0,240.0]
#ha=[740,600,500,370,220,210,320,225];la=[660,450,370,270,180,170,200,160]
ha=[960,800,650,450,280,280,270,350];la=[700,600,450,330,200,200,180,200]

def get_num(Tb_sub5s):
    num5s=Tb_sub5s*1.0
    num5s[num5s!=0]=1;num5s_mean=np.mean(num5s,axis=0);num=np.zeros((200,200))
    for i in range(200):
        for j in range(200):
            if(num5s_mean[i,j]!=0):
                num5s_mean[i,j]=1
            t=Tb_sub5s[:,i,j];num[i,j]=len([list(group) for group in mit.consecutive_groups(t.nonzero()[0])])
    num5s1=Tb_sub5s*1.0
    num5s1[num5s1!=0]=1;num5s1_mean=np.sum(num5s1,axis=0)
    num=np.array(num)
    return num5s_mean,num#num5s1_mean

Tb_suball=[0]*8;Tb_sub_std_all=[0]*8;num6s=[0]*8;num6s_mean=[0]*8;tidx=[0]*8;Tb_subf=[0]*8;Tb_sub6s=[0]*8;Tb_sub8s=[0]*8;Tb_subn6s=[0]*8;num8s=[0]*8
numr1=[0]*8;numr2=[0]*8;numr3=[0]*8;numr4=[0]*8;numr5=[0]*8;numr6=[0]*8
for j in range(8):
    print str(freq[j])+' MHz'
    f=freq_filelist[j]#'187-188'
    [Tb_sub,Tb_sub_std,xc,yc,time_string,time_sec,bmaj,bmin,bpa,ndata]=pickle.load(open(Tb_write_path+'/Tb_'+str(ids)+f+'_sub_test.p','rb'))
    Tb_suball[j]=Tb_sub;Tb_sub_std_all[j]=Tb_sub_std
    h=ha[j];l=la[j]
    Tb_subf[j],tidx[j]=remove_outliers(Tb_sub,Tb_sub_std,l,h)
    Tb_sub5s,Tb_subn5s=extract_outliers(Tb_subf[j],5)
    Tb_sub6s[j],Tb_subn6s[j]=extract_outliers(Tb_subf[j],6)
    Tb_sub8s[j],Tb_subn8s=extract_outliers(Tb_subf[j],8)
    num5s,num5s_mean=get_num(Tb_sub5s)
    num6s[j],num6s_mean[j]=get_num(Tb_sub6s[j])
    num8s[j],num8s_mean=get_num(Tb_sub8s[j])
    numr1[j]=int(num6s_mean[j][96,111]);numr2[j]=int(num6s_mean[j][107,117]);numr3[j]=int(num6s_mean[j][99,77]);numr4[j]=int(num6s_mean[j][103,92]);numr5[j]=int(num6s_mean[j][88,95]);numr6[j]=int(num6s_mean[j][82,79])
Tb_suball=np.array(Tb_suball);Tb_sub_std=np.array(Tb_sub_std);Tb_sub6s=np.array(Tb_sub6s);Tb_subn6s=np.array(Tb_subn6s);Tb_sub8s=np.array(Tb_sub8s)


h=740;l=660 # 108 MHz
h=740;l=660 # 120 MHz
h=500;l=370 # 133 MHz
#h=800;l=160 # 240 MHz

k=0;Tbstd_arr=[0]*100
for i in np.linspace(0,1500,100):
    #Tbstd_arr[k]=np.std(Tb_subf[0:int(i)])
    k=k+1 



noise=[0]*8
for j in range(8): 
    noise[j]=[0]*1000
    for i in range(1000):
        noise[j][i]=Tb_subf[j][30:30+i,0:50,:].mean(axis=0).std()
noise=np.array(noise)



l1=np.linspace(0,200,51)#;l2=np.linspace(20,80,13)[1:]
xxx,yyy=np.meshgrid(l1, l1)
Tb_region=[0]*8;Tb_region_mean=[0]*8;Tb_region1=[0]*8;Tb_region2=[0]*8;Tb_region3=[0]*8;Tb_region4=[0]*8;Tb_region5=[0]*8;Tb_region6=[0]*8
Tb_cregion=[0]*8;Tb_cregion_mean=[0]*8;Tb_cregion1=[0]*8;Tb_cregion2=[0]*8;Tb_cregion3=[0]*8;Tb_cregion4=[0]*8;Tb_cregion5=[0]*8;Tb_cregion6=[0]*8
ts1=[0]*8;ts2=[0]*8;ts3=[0]*8;ts4=[0]*8;ts5=[0]*8;ts6=[0]*8;ts7=[0]*8;ts8=[0]*8
#Tb_suball=[Tb_sub108,Tb_sub120,Tb_sub133,Tb_sub145,Tb_sub179,Tb_sub197,Tb_sub217,Tb_sub240]
#Tb_suball=[Tb_sub108_cor[:,50:-50,50:-50],Tb_sub120_cor[:,50:-50,50:-50],Tb_sub133_cor[:,50:-50,50:-50],Tb_sub145_cor[:,50:-50,50:-50],Tb_sub179_cor[:,50:-50,50:-50],Tb_sub197_cor[:,50:-50,50:-50],Tb_sub217_cor[:,50:-50,50:-50],Tb_sub240_cor[:,50:-50,50:-50]]
#Tb_suball=[Tb_sub108_cor_sel[:,50:-50,50:-50],Tb_sub120_cor_sel[:,50:-50,50:-50],Tb_sub133_cor_sel[:,50:-50,50:-50],Tb_sub145_cor_sel[:,50:-50,50:-50],Tb_sub179_cor_sel[:,50:-50,50:-50],Tb_sub197_cor_sel[:,50:-50,50:-50],Tb_sub217_cor_sel[:,50:-50,50:-50],Tb_sub240_cor_sel[:,50:-50,50:-50]]
idxr1=[0]*8;idxr2=[0]*8;idxr3=[0]*8;idxr4=[0]*8;idxr5=[0]*8;idxr6=[0]*8
for j in range(8):
    Tb_region[j]=[0]*50
    for i in range(50):
        Tb_region[j][i]=[0]*50
        for k in range(50):
            Tb_region[j][i][k]=np.nanmean(Tb_sub6s[j][:,int(xxx[k+1][i+1]-3):int(xxx[k+1][i+1])+3,int(yyy[k+1][i+1]-3):int(yyy[k+1][i+1]+3)],axis=(1,2))
    Tb_region_mean[j]=np.nanmean(Tb_region[j],axis=2)
    Tb_region2[j]=np.array(Tb_region[j])[26,28];Tb_region3[j]=np.array(Tb_region[j])[24,18];Tb_region4[j]=np.array(Tb_region[j])[26,23]
    Tb_region5[j]=np.array(Tb_region[j])[21,24];Tb_region6[j]=np.array(Tb_region[j])[20,16];Tb_region1[j]=np.array(Tb_region[j])[23,26]
    Tb_cregion1[j]=Tb_region1[j]*1.0;Tb_cregion2[j]=Tb_region2[j]*1.0;Tb_cregion3[j]=Tb_region3[j]*1.0;Tb_cregion4[j]=Tb_region4[j]*1.0;Tb_cregion5[j]=Tb_region5[j]*1.0;Tb_cregion6[j]=Tb_region6[j]*1.0
    #Tb_cregion1[j][Tb_cregion1[j]<Tblim[j]]=0
    #Tb_cregion2[j][Tb_cregion2[j]<Tblim[j]]=0
    #Tb_cregion3[j][Tb_cregion3[j]<Tblim[j]]=0
    #Tb_cregion4[j][Tb_cregion4[j]<Tblim[j]]=0
    #Tb_cregion5[j][Tb_cregion5[j]<Tblim[j]]=0
    #Tb_cregion6[j][Tb_cregion6[j]<Tblim[j]]=0
    idx1=np.where(Tb_cregion1[j]!=0)[0];diffs = np.diff(idx1) != 1;indexes = np.nonzero(diffs)[0] + 1;ts1[j] = np.split(idx1, indexes)
    idx2=np.where(Tb_cregion2[j]!=0)[0];diffs = np.diff(idx2) != 1;indexes = np.nonzero(diffs)[0] + 1;ts2[j] = np.split(idx2, indexes)
    idx3=np.where(Tb_cregion3[j]!=0)[0];diffs = np.diff(idx3) != 1;indexes = np.nonzero(diffs)[0] + 1;ts3[j] = np.split(idx3, indexes)
    idx4=np.where(Tb_cregion4[j]!=0)[0];diffs = np.diff(idx4) != 1;indexes = np.nonzero(diffs)[0] + 1;ts4[j] = np.split(idx4, indexes)
    idx5=np.where(Tb_cregion5[j]!=0)[0];diffs = np.diff(idx5) != 1;indexes = np.nonzero(diffs)[0] + 1;ts5[j] = np.split(idx5, indexes)
    idx6=np.where(Tb_cregion6[j]!=0)[0];diffs = np.diff(idx6) != 1;indexes = np.nonzero(diffs)[0] + 1;ts6[j] = np.split(idx6, indexes)
    idxr1[j]=idx1;idxr2[j]=idx2;idxr3[j]=idx3;idxr4[j]=idx4;idxr5[j]=idx5;idxr6[j]=idx6

nts1=[0]*8;tts1=[0]*8;nhist1=[0]*8
nts2=[0]*8;tts2=[0]*8;nhist2=[0]*8
nts3=[0]*8;tts3=[0]*8;nhist3=[0]*8
nts4=[0]*8;tts4=[0]*8;nhist4=[0]*8
nts5=[0]*8;tts5=[0]*8;nhist5=[0]*8
nts6=[0]*8;tts6=[0]*8;nhist6=[0]*8
bins=np.arange(1,10)*0.5
for i in range(8):
    nts1[i]=len(ts1[i]);tts1[i]=[0]*(len(ts1[i])-1)
    for j in range(len(ts1[i])-1):
        tts1[i][j]=len(ts1[i][j])*0.5
    nhist1[i]=np.histogram(tts1[i],bins=bins)[0]
    nts2[i]=len(ts2[i]);tts2[i]=[0]*(len(ts2[i])-1)
    for j in range(len(ts2[i])-1):
        tts2[i][j]=len(ts2[i][j])*0.5
    nhist2[i]=np.histogram(tts2[i],bins=bins)[0]
    nts3[i]=len(ts3[i]);tts3[i]=[0]*(len(ts3[i])-1)
    for j in range(len(ts3[i])-1):
        tts3[i][j]=len(ts3[i][j])*0.5
    nhist3[i]=np.histogram(tts3[i],bins=bins)[0]
    nts4[i]=len(ts4[i]);tts4[i]=[0]*(len(ts4[i])-1)
    for j in range(len(ts4[i])-1):
        tts4[i][j]=len(ts4[i][j])*0.5
    nhist4[i]=np.histogram(tts4[i],bins=bins)[0]
    nts5[i]=len(ts5[i]);tts5[i]=[0]*(len(ts5[i])-1)
    for j in range(len(ts5[i])-1):
        tts5[i][j]=len(ts5[i][j])*0.5
    nhist5[i]=np.histogram(tts5[i],bins=bins)[0]
    nts6[i]=len(ts6[i]);tts6[i]=[0]*(len(ts6[i])-1)
    for j in range(len(ts6[i])-1):
        tts6[i][j]=len(ts6[i][j])*0.5
    nhist6[i]=np.histogram(tts6[i],bins=bins)[0]
nhist1=np.array(nhist1);nhist2=np.array(nhist2);nhist3=np.array(nhist3);nhist4=np.array(nhist4);nhist5=np.array(nhist5);nhist6=np.array(nhist6)


def get_fhpair(tidx,Tb_region3,l,m):
    list06=list(set(tidx[l]).intersection(tidx[m]));num=0;Tb_region3_06=[l]*len(list06);idx_06=[m]*len(list06)
    for i in list06:
        if((Tb_region3[l][list06.index(i)]!=0) & (Tb_region3[m][list06.index(i)]!=0)):
            Tb_region3_06[num]=[Tb_region3[l][list06.index(i)],Tb_region3[m][list06.index(i)]];idx_06[num]=[tidx[l][list06.index(i)],tidx[m][list06.index(i)]]
            num=num+1
    Tb_region3_06=[i for i in Tb_region3_06 if i !=0];Tb_region3_06=np.array(Tb_region3_06[0:num])
    idx_06=[i for i in idx_06 if i !=0];idx_06=np.array(idx_06[0:num])
    return idx_06,Tb_region3_06

idx1_06,Tb_region1_06=get_fhpair(tidx,Tb_region1,0,6)
idx1_17,Tb_region1_17=get_fhpair(tidx,Tb_region1,1,7)
idx2_06,Tb_region2_06=get_fhpair(tidx,Tb_region2,0,6)
idx2_17,Tb_region2_17=get_fhpair(tidx,Tb_region2,1,7)
idx3_06,Tb_region3_06=get_fhpair(tidx,Tb_region3,0,6)
idx3_17,Tb_region3_17=get_fhpair(tidx,Tb_region3,1,7)
idx4_06,Tb_region4_06=get_fhpair(tidx,Tb_region4,0,6)
idx4_17,Tb_region4_17=get_fhpair(tidx,Tb_region4,1,7)
idx5_06,Tb_region5_06=get_fhpair(tidx,Tb_region5,0,6)
idx5_17,Tb_region5_17=get_fhpair(tidx,Tb_region5,1,7)
idx6_06,Tb_region6_06=get_fhpair(tidx,Tb_region6,0,6)
idx6_17,Tb_region6_17=get_fhpair(tidx,Tb_region6,1,7)

def get_num_flare(d):
    Tbmap=d*0;s=np.where(d!=0);s0=abs(s[0][1:]-s[0][:-1]);s1=abs(s[1][1:]-s[1][:-1]);ds=s0-s1
    dsidx=np.where((ds<-2) | (ds>2))[0]+1;dsidx=np.insert(dsidx,0,0);dsidx=np.insert(dsidx,len(dsidx),len(s[0]))
    s0out=[0]*(len(dsidx)-1);s1out=[0]*(len(dsidx)-1)
    for ss in range(len(dsidx)-1):
        s0out[ss]=s[0][dsidx[ss]:dsidx[ss+1]]
        s1out[ss]=s[1][dsidx[ss]:dsidx[ss+1]]
    num=len(s0out);Tb=[0]*num;area=[0]*num;num_sel=0
    for i in range(num):
        Tb[i]=np.mean(d[s0out[i],s1out[i]])
        area[i]=len(s0out[i])*2500.0/(200)**2 # In PSF Area
        if(area[i]>0.9):
            Tbmap[s0out[i],s1out[i]]=d[s0out[i],s1out[i]]
            num_sel=num_sel+1
    return Tb,Tbmap,num,area,num_sel

Tball=[Tb108m,Tb120m,Tb133m,Tb145m,Tb179m,Tb197m,Tb217m,Tb240m]

for i in range(8):
    #print freq[i],' &',np.round(Tball[i].max()/1000/1000,2),'&', np.round(Tb_region1[i].sum()/nts1[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region1[i])/1.e3,2),'&', nts1[i],'&', np.round(Tb_region2[i].sum()/nts2[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region2[i])/1.e3,2),'&', nts2[i],'&', \
     #   np.round(Tb_region3[i].sum()/nts3[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region3[i])/1.e3,2),'&',nts3[i],'&',np.round(Tb_region4[i].sum()/nts4[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region4[i][0:500])/1.e3,2),'&',nts4[i], '&', \
      #  np.round(Tb_region5[i].sum()/nts5[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region5[i])/1.e3,2),'&',nts5[i],'&', np.round(Tb_region6[i].sum()/nts6[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region6[i])/1.e3,2),'&',nts6[i],'\\\\ \hline'
    print freq[i],' &',np.round(Tball[i].max()/1000/1000,2),'&', np.round(Tb_region1[i].sum()/numr1[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region1[i])/1.e3,2),'&', numr1[i],'&', np.round(Tb_region2[i].sum()/numr2[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region2[i])/1.e3,2),'&', numr2[i],'&', \
        np.round(Tb_region3[i].sum()/numr3[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region3[i])/1.e3,2),'&',numr3[i],'&',np.round(Tb_region4[i].sum()/numr4[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region4[i][0:500])/1.e3,2),'&',numr4[i], '&', \
        np.round(Tb_region5[i].sum()/numr5[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region5[i])/1.e3,2),'&',numr5[i],'&', np.round(Tb_region6[i].sum()/numr6[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region6[i])/1.e3,2),'&',numr6[i],'\\\\ \hline'

def f2r(f):
    r_=r[50,50,:][ut.find_nearest(np.sqrt(densobs[50,50,:])*9000/1.e6,f)[0]]
    dens_=densobs[50,50,:][ut.find_nearest(np.sqrt(densobs[50,50,:])*9000/1.e6,f)[0]]
    temp_=tempobs[50,50,:][ut.find_nearest(np.sqrt(densobs[50,50,:])*9000/1.e6,f)[0]]
    b=bobs[50,50,:][ut.find_nearest(np.sqrt(densobs[50,50,:])*9000/1.e6,f)[0]]
    return r_*6.9e5,dens_,temp_,b

beampix=20
numall=[0]*8;rall_omegap=[0]*8;numall_str=[0]*8;freq_str=[0]*8;numregion=[0]*8;perregion=[0]*8;omega_p=[0]*8;l=[0]*8
for i in range(8):
    numall[i]=num6s_mean[i].sum()/beampix
    omega_p[i]=np.sqrt(densobs[50,50,:])*9000
    rall_omegap[i]=r[50,50,:][ut.find_nearest(omega_p[i]/1.e6,freq[i])[0]]-1;rall_omegap[7]=0.011
    l[i]=((r[50,50,:][ut.find_nearest(omega_p[i]/1.e6+2.5,freq[i])[0]]-1)-(r[50,50,:][ut.find_nearest(omega_p[i]/1.e6 - 2.5,freq[i])[0]]-1))*6.9e10
    numall_str[i]=str(int(numall[i]))
    freq_str[i]=str(int(freq[i]))
    numregion[i]=nts1[i]+nts2[i]+nts3[i]+nts4[i]
    perregion[i]=int(numregion[i]/numall[i]*100)
l=np.array(l);l[np.where(l==0)]=l[3]

A=1.31e19 # in cgs
A=5.0e17 # 10" by 10" 
tau=[0.9,1.3,1.4,1.5,1.6,2.6,4.0,3.1]
nemwa=[0]*8;tempmwa=[0]*8;bmwa=[0]*8;rmwa=[0]*8;Eb=[0]*8;Eth=[0]*8;Er=[0]*8;Smwa=[0]*8;Wmwa=[0]*8;frac=[0]*8;Eth6s=[0]*8;Wall=[0]*8
for i in range(8):
    rmwa[i],nemwa[i],tempmwa[i],bmwa[i]=f2r(freq[i])
    Eb[i]=bmwa[i]**2 /(8*np.pi)*A*l[i];Eth[i]=nemwa[i]*1.38e-16*tempmwa[i]*A*l[i];Eth6s[i]=nemwa[i]*1.38e-16*A*l[i]*np.nanmean(Tb_sub6s[i].sum(axis=0)/num6s_mean[i])
    Er[i]=Tb_sub6s[i].mean(axis=0).sum()*nemwa[i]*1.38e-16*A*l[i]
    Smwa[i]=Tb_sub6s[i].mean(axis=0).sum()*(freq[i]/1000)**2 *200**2 /1224/1.e4 # mSFU
    #Wmwa[i]=4*np.pi*(1.5e13)**2 * 2.e6*0.5*Smwa[i]*1.e-3*1.e-4*1.e-7
    Wmwa[i]=0.15*(1.5e11)**2 * (5.e6*0.5*Smwa[i]*1.e-3)*1.e-4*1.e-7*np.exp(tau[i]) # Jy --> ergs
    Wall_=0.15*(1.5e11)**2 * (5.e6*0.5*(Tb_sub6s[i].flatten()*(freq[i]/1000)**2 *200**2 /1224/1.e4)*1.e-3)*1.e-4*1.e-7*np.exp(tau[i]) # Jy and ergs -> J
    frac[i]=Wmwa[i]/Eth6s[i].mean()
    Wall[i]=Wall_[Wall_!=0]
rmwa=np.array(rmwa)/1.e3;nemwa=np.array(nemwa);tempmwa=np.array(tempmwa);bmwa=np.array(bmwa);Eb=np.array(Eb);Eth=np.array(Eth);Er=np.array(Er);Smwa=np.array(Smwa);Wmwa=np.array(Wmwa);Wall=np.array(Wall);frac=np.array(frac);Eth6s=np.array(Eth6s)
Wall=list(itertools.chain(*Wall))

alp=1.e-5;tau=[0.9,1.3,1.4,1.5,1.6,2.6,4.0,3.1]
frac=[0]*8;areaall=[0]*8;Tbfreqall=[0]*8;Tb_sub6s_burst=[0]*8;Tb_sub6s_burstmap=[0]*8;Tb_sub6s_burstmap_all=[0]*8
num6s_burst=[0]*8;area_burst=[0]*8;tot_num6s_burst=[0]*8;area_burst_all=[0]*8;num6s_burst_area=[0]*8;tot_num6s_burst_area=[0]*8;Tb_sub6s_mean=[0]*8;Eth_sub6s=[0]*8
for i in range(8): # Frequency
    areaall[i]=[0]*len(Tb_sub6s[i]);Tbfreqall[i]=[0]*len(Tb_sub6s[i]);Tb_sub6s_burst[i]=[0]*len(Tb_sub6s[i]);Tb_sub6s_burstmap[i]=[0]*len(Tb_sub6s[i])
    num6s_burst[i]=[0]*len(Tb_sub6s[i]);area_burst[i]=[0]*len(Tb_sub6s[i]);num6s_burst_area[i]=[0]*len(Tb_sub6s[i])
    for j in range(len(Tb_sub6s[i])): # Time
        areaall[i][j]=len(np.where(Tb_sub6s[i][j]!=0)[0])*2500/((200)**2) # In PSF Area
        Tbfreqall[i][j]=np.mean(Tb_sub6s[i][j][np.where(Tb_sub6s[i][j]!=0)])
        Tb_sub6s_burst[i][j],Tb_sub6s_burstmap[i][j],num6s_burst[i][j],area_burst[i][j],num6s_burst_area[i][j]=get_num_flare(Tb_sub6s[i][j])
    areaall[i]=np.array(areaall[i])[np.isfinite(np.array(areaall[i]))]
    Tbfreqall[i]=np.array(Tbfreqall[i])[np.isfinite(np.array(Tbfreqall[i]))]
    tot_num6s_burst[i]=np.array(num6s_burst[i]).sum()
    tot_num6s_burst_area[i]=np.array(num6s_burst_area[i]).sum()
    area_burst_all[i]=[k for k in area_burst[i] if k != [0]];area_burst_all[i]=list(itertools.chain(*area_burst_all[i]))
    Tb_sub6s_burstmap_all[i]=[k for k in Tb_sub6s_burst[i] if len(k) >1]
    Tb_sub6s_burstmap_all[i]=list(itertools.chain(*Tb_sub6s_burstmap_all[i]))
    Tb_sub6s_mean[i]=np.sum(Tb_sub6s[i],axis=0)
    Eth_sub6s[i]=alp*Tb_sub6s_mean[i]*nemwa[i]*1.38e-16*A*l[i]*np.exp(tau[i])
    frac[i]=Wmwa[i]/Eth_sub6s[i].max()
area_burst_freq=list(itertools.chain(*area_burst_all));Tb_sub6s_mean=np.array(Tb_sub6s_mean);Eth_sub6s=np.array(Eth_sub6s);frac=np.array(frac)
Tb_sub6s_burstmap_all=list(itertools.chain(*Tb_sub6s_burstmap_all))
Wall=list(itertools.chain(*Wall))

area_burst_freq=np.array(area_burst_freq)
area_burst_freq[np.where(area_burst_freq<1.0)] = 1.0
plt.hist(np.array(area_burst_freq),bins=8,align='left',normed=True)
plt.ylabel('Normalised Counts');plt.xlabel('PSF Area');plt.yscale('log');plt.xlim(0.7,3.4)
plt.show()

Tbfr=np.array(Eth_sub6s).flatten();Tbfreqall1=Tbfr[Tbfr!=0]
Tbfreqall1_hist=plt.hist(Tbfreqall1,bins=np.logspace(np.log10(1.e20),np.log10(2.e21),10),normed=True)
#Tb_sub6s_burstmap_all_hist=np.histogram(Tb_sub6s_burstmap_all,bins=15)
z1=np.polyfit(np.log10(Tbfreqall1_hist[1][0:9]),np.log10(Tbfreqall1_hist[0][0:9]), 1,cov=True)
#z=np.polyfit(np.log10(Tb_sub6s_burstmap_all_hist[1][1:6]),np.log10(Tb_sub6s_burstmap_all_hist[0][0:5]), 1)
p1=np.poly1d(z1[0]);em=np.sqrt(np.diag(z1[1])[0]);xfit1=np.linspace(1.e20,2.e21,100);xfit_1=np.log10(xfit1)
yfit_1=p1(xfit_1);yfit1=10**(yfit_1+0.2)
plt.hist(Tbfreqall1,bins=np.logspace(np.log10(1.e20),np.log10(2.e21),10),label='Nonthermal Energy',normed=True)
#plt.hist(Tb_sub6s_burstmap_all,bins=15,label='T$_B$')
plt.plot(xfit1,yfit1,label='Linear fit, Slope='+str(np.round(z1[0][0],2))+'$\pm$'+str(np.round(em,2)),color='k');plt.xscale('log');plt.yscale('log')
plt.ylabel('dN/dE$_{nth}$ (ergs$^{-1}$)');plt.xlabel('Nonthermal Flare Energy (ergs)');plt.legend(loc=1)
plt.show()

Tbfreqall1_hist=plt.hist(Wall,bins=np.logspace(np.log10(1.e14),np.log10(2.e16),10),normed=True)
#Tb_sub6s_burstmap_all_hist=np.histogram(Tb_sub6s_burstmap_all,bins=15)
z1=np.polyfit(np.log10(Tbfreqall1_hist[1][0:9]),np.log10(Tbfreqall1_hist[0][0:9]), 1,cov=True)
#z=np.polyfit(np.log10(Tb_sub6s_burstmap_all_hist[1][1:6]),np.log10(Tb_sub6s_burstmap_all_hist[0][0:5]), 1)
p1=np.poly1d(z1[0]);em=np.sqrt(np.diag(z1[1])[0]);xfit1=np.linspace(1.e14,2.e16,100);xfit_1=np.log10(xfit1)
yfit_1=p1(xfit_1);yfit1=10**(yfit_1+0.2)
plt.hist(Wall,bins=np.logspace(np.log10(1.e14),np.log10(2.e16),10),label='Radiated Energy',normed=True)
#plt.hist(Tb_sub6s_burstmap_all,bins=15,label='T$_B$')
plt.plot(xfit1,yfit1,label='Linear fit, Slope='+str(np.round(z1[0][0],2))+'$\pm$'+str(np.round(em,2)),color='k');plt.xscale('log');plt.yscale('log')
plt.ylabel('dN/dW (ergs$^{-1}$)');plt.xlabel('Radiated Energy (ergs)');plt.legend(loc=1)
plt.show()

#Tbfreqall1=[item for sublist in Tbfreqall for item in sublist]
fit=[0]*8;yfit1=[0]*8;xfit1=[0]*8;ul=[1.e4,1.e4,1.e4,1.e4,1.e4,1.e4,1.e4,1.e4];ur=[9.e5,9.e5,9.e5,9.e5,9.e5,9.e5,9.e5,9.e5]
#m1=[1,1,1,1,1,1,4,5];m2=[9,9,9,9,10,8,9,14];em=[0]*8
m1=[1,1,1,1,1,1,1,1];m2=[10,9,9,10,10,9,11,14];em=[0]*8
for i in range(8):
    Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0]
    #Tbfreqall1_hist=np.histogram(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),normed=True)
    Tbfreqall1_hist=plt.hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),label='T$_B$',normed=True)
    #Tb_sub6s_burstmap_all_hist=np.histogram(Tb_sub6s_burstmap_all,bins=15)
    z1=np.polyfit(np.log10(Tbfreqall1_hist[1][m1[i]+1:m2[i]+1]),np.log10(Tbfreqall1_hist[0][m1[i]:m2[i]]), 1,cov=True)
    #z=np.polyfit(np.log10(Tb_sub6s_burstmap_all_hist[1][1:6]),np.log10(Tb_sub6s_burstmap_all_hist[0][0:5]), 1)
    p1=np.poly1d(z1[0]);em[i]=np.sqrt(np.diag(z1[1])[0]);xfit1[i]=np.linspace(ul[i],ur[i],100);xfit_1=np.log10(xfit1[i])
    yfit_1=p1(xfit_1);yfit1[i]=10**yfit_1;fit[i]=z1[0]

f,ax=plt.subplots(2,4,sharex=True,sharey=True)
i=0;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[0,0].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[0,0].plot(xfit1[i][1:40],yfit1[i][1:40],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[0,0].legend()
i=1;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[0,1].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[0,1].plot(xfit1[i][1:35],yfit1[i][1:35],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[0,1].legend()
i=2;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[0,2].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[0,2].plot(xfit1[i][1:35],yfit1[i][1:35],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[0,2].legend()
i=3;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[0,3].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[0,3].plot(xfit1[i][1:35],yfit1[i][1:35],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[0,3].legend()
i=4;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[1,0].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[1,0].plot(xfit1[i][1:38],yfit1[i][1:38],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[1,0].legend()
i=5;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[1,1].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[1,1].plot(xfit1[i][1:32],yfit1[i][1:32],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[1,1].legend()
i=6;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[1,2].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[1,2].plot(xfit1[i][1:40],yfit1[i][1:40],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[1,2].legend()
i=7;Tbfr=np.array(Tb_sub6s_mean[i]).flatten();Tbfreqall1=Tbfr[Tbfr!=0];ax[1,3].hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),histtype='stepfilled',label=str(freq[i])+' MHz',normed=True)
ax[1,3].plot(xfit1[i][1:150],yfit1[i][1:150],'--',label='$\\alpha$='+str(np.round(fit[i][0],2))+'$\pm$'+str(np.round(em[i],2)),color='k');plt.xscale('log');ax[1,3].legend()
ax[1,0].set_xscale('log');ax[1,0].set_yscale('log')#;ax[1,0].set_ylim(0.8,5001);ax[1,0].set_xlim(-0.09,30)
ax[0,0].set_ylabel('dN/dT$_{B}$ (K$^{-1}$)');ax[1,0].set_ylabel('dN/dT$_{B}$ (K$^{-1}$)');ax[1,0].set_xlabel('T$_{B}$ (K)');ax[1,1].set_xlabel('T$_{B}$ (K)');ax[1,2].set_xlabel('T$_{B}$ (K)');ax[1,3].set_xlabel('T$_{B}$ (K)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()


plt.hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),label='T$_B$',normed=True)
#plt.hist(Tb_sub6s_burstmap_all,bins=15,label='T$_B$')
plt.plot(xfit1,yfit1,label='Linear fit, Slope='+str(np.round(z1[0][0],2))+'$\pm$'+str(np.round(em,2)),color='k');plt.xscale('log');plt.yscale('log')
plt.ylabel('dN/dT$_B$');plt.xlabel('T$_B$ (K)');plt.legend()
plt.show()


Tbfr=np.array(num6s_mean).flatten();Tbfreqall2=Tbfr[Tbfr!=0]
Tbfreqall2_hist=np.histogram(Tbfreqall2,bins=15)
z2=np.polyfit(np.log10(Tbfreqall2_hist[1][1:12]),np.log10(Tbfreqall2_hist[0][0:11]), 1,cov=True)
p2=np.poly1d(z2[0]);em=np.sqrt(np.diag(z2[1])[0]);xfit2=np.linspace(10,2.e2,100);xfit_2=np.log10(xfit2)
yfit_2=p2(xfit_2);yfit2=10**yfit_2
plt.hist(Tbfreqall2,bins=15,label='T$_B$')
plt.plot(xfit2,yfit2,label='Linear fit, Slope='+str(np.round(z2[0][0],2))+'$\pm$'+str(np.round(em,2)),color='k');plt.xscale('log');plt.yscale('log')
plt.ylabel('Occurance Frequency');plt.xlabel('Number of bursts');plt.legend()
plt.show()

aiafile='/nas08-data02/rohit/20151203_EUV/aia.lev1.193A_2015-12-03T03_20_17.84Z.image_lev1.fits'
af=fits.open(aiafile)
h=af[0].header;d=af[0].data
aiamap=Map(aiafile)

lascoc2_file='/nas08-data02/rohit/20151203_LASCO/22572578_1.fts'
c2=fits.open(lascoc2_file)
hc=c2[0].header;dc=c2[0].data

    
for i in range(8):
    print freq[i],' &',np.round(nemwa[i]/1.e8,2),' & ',np.round(bmwa[i],2),'&',np.round(tempmwa[i]/1.e6,2),' & ',np.round(Eb[i]/1.e25,2), ' & ',np.round(Eth[i]/1.e25,2),' & ',np.round(Smwa[i],2),' & ' , \
            np.round(Wmwa[i]/1.e15,2),' & ',np.round(frac[i]/1.e-5,2),' & ',np.round(np.nanmean(Tb_sub6s[i].sum(axis=0)/num6s_mean[i])/1000,2),' & ',np.round(Eth_sub6s[i].max()/1.e20,2),'\\\\ \\hline'

for i in range(8):
    print freq[i],' & ',int(np.array(numall_6s)[i]),' & ',int(np.array(numregion_6s)[i]),' & ',int(numregion_6s[i]/numall_6s[i]*100),'\\\\ \\hline'

Tbfr=np.array(Tb_sub6s_mean[3]).flatten();Tbfreqall1=Tbfr[Tbfr!=0]
Tbfreqall1_hist=np.histogram(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15))
#Tb_sub6s_burstmap_all_hist=np.histogram(Tb_sub6s_burstmap_all,bins=15)
z1=np.polyfit(np.log10(Tbfreqall1_hist[1][6:13]),np.log10(Tbfreqall1_hist[0][5:12]), 1,cov=True)
#z=np.polyfit(np.log10(Tb_sub6s_burstmap_all_hist[1][1:6]),np.log10(Tb_sub6s_burstmap_all_hist[0][0:5]), 1)
p1=np.poly1d(z1[0]);em=np.diag(z1[1])[0];xfit1=np.linspace(5.e4,9.e5,100);xfit_1=np.log10(xfit1)
yfit_1=p1(xfit_1);yfit1=10**yfit_1
plt.hist(Tbfreqall1,bins=np.logspace(np.log10(1.e4),np.log10(1.e6),15),label='T$_B$')
#plt.hist(Tb_sub6s_burstmap_all,bins=15,label='T$_B$')
plt.plot(xfit1,yfit1,label='Linear fit, Slope='+str(np.round(z1[0][0],2))+'$\pm$'+str(np.round(em,2)),color='k');plt.xscale('log');plt.yscale('log')
plt.ylabel('Number');plt.xlabel('T$_B$ (K)');plt.legend()
plt.show()


sys.exit()

f,ax=plt.subplots(2,4,sharex=True,sharey=True)
ax[0,0].hist(Tbfreqall[0],bins=50,histtype='step',label=str(freq[0])+' MHz');ax[0,0].legend()
ax[0,1].hist(Tbfreqall[1],bins=50,histtype='step',label=str(freq[1])+' MHz');ax[0,1].legend()
ax[0,2].hist(Tbfreqall[2],bins=50,histtype='step',label=str(freq[2])+' MHz');ax[0,2].legend()
ax[0,3].hist(Tbfreqall[3],bins=50,histtype='step',label=str(freq[3])+' MHz');ax[0,3].legend()
ax[1,0].hist(Tbfreqall[4],bins=50,histtype='step',label=str(freq[4])+' MHz');ax[1,0].legend()
ax[1,1].hist(Tbfreqall[5],bins=50,histtype='step',label=str(freq[5])+' MHz');ax[1,1].legend()
ax[1,2].hist(Tbfreqall[6],bins=50,histtype='step',label=str(freq[6])+' MHz');ax[1,2].legend()
ax[1,3].hist(Tbfreqall[7],bins=50,histtype='step',label=str(freq[7])+' MHz');ax[1,3].legend()
ax[1,0].set_xscale('log');ax[1,0].set_yscale('log');ax[1,0].set_ylim(0.8,5001);ax[1,0].set_xlim(999,5.e4)
ax[0,0].set_ylabel('Numbers');ax[1,0].set_ylabel('Numbers');ax[1,0].set_xlabel('$T_B$ (K)');ax[1,1].set_xlabel('$T_B$ (K)');ax[1,2].set_xlabel('$T_B$ (K)');ax[1,3].set_xlabel('$T_B$ (K)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

f,ax=plt.subplots(2,4,sharex=True,sharey=True)
ax[0,0].hist(areaall[0],bins=10,histtype='step',label=str(freq[0])+' MHz');ax[0,0].legend()
ax[0,1].hist(areaall[1],bins=10,histtype='step',label=str(freq[1])+' MHz');ax[0,1].legend()
ax[0,2].hist(areaall[2],bins=10,histtype='step',label=str(freq[2])+' MHz');ax[0,2].legend()
ax[0,3].hist(areaall[3],bins=10,histtype='step',label=str(freq[3])+' MHz');ax[0,3].legend()
ax[1,0].hist(areaall[4],bins=10,histtype='step',label=str(freq[4])+' MHz');ax[1,0].legend()
ax[1,1].hist(areaall[5],bins=10,histtype='step',label=str(freq[5])+' MHz');ax[1,1].legend()
ax[1,2].hist(areaall[6],bins=10,histtype='step',label=str(freq[6])+' MHz');ax[1,2].legend()
ax[1,3].hist(areaall[7],bins=10,histtype='step',label=str(freq[7])+' MHz');ax[1,3].legend()
ax[1,0].set_xscale('log');ax[1,0].set_yscale('log');ax[1,0].set_ylim(0.8,5001);ax[1,0].set_xlim(-0.09,30)
ax[0,0].set_ylabel('Numbers');ax[1,0].set_ylabel('Numbers');ax[1,0].set_xlabel('Size (PSF)');ax[1,1].set_xlabel('Size (PSF)');ax[1,2].set_xlabel('Size (PSF)');ax[1,3].set_xlabel('Size (PSF)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(r[38,56,:]-1,np.sqrt(densobs[38,56,:])*9000/1.e6,'o-',label='Plasma Frequency (Region 3)')
ax.axvspan(0.027, 0.092, alpha=0.5,color='gray');ax.legend()
ax.set_ylabel('Frequency (MHz)');ax.set_xlabel('Coronal Height ($R_{\odot}$)')
plt.show()

f,ax=plt.subplots(1,1)
f.subplots_adjust(right=0.75)
ax.plot(r[50,50,:]-1,densobs[50,50,:]/1.e8,'o-',color='k');ax.set_ylabel('Electron Density ($\\times 10^{8}$ cm$^{-3}$)')
ax1=ax.twinx();ax1.plot(r[50,50,:]-1,tempobs[50,50,:]/1.e6,'o-',color='red');ax1.set_ylabel('Temperature (MK)')
ax1.yaxis.label.set_color('red');ax1.tick_params(axis='y',colors='red');ax1.spines["right"].set_edgecolor('red')
ax2=ax.twinx();ax2.plot(r[50,50,:]-1,bobs[50,50,:],'o-',color='green');ax2.set_ylabel('Magnetic field (G)')
ax2.spines["right"].set_position(("axes", 1.1));ax2.yaxis.label.set_color('green');ax2.tick_params(axis='y',colors='green');ax2.spines["right"].set_edgecolor('green')
ax.set_xlabel('Coronal Height ($R_{\odot}$)')
plt.show()



def convert_f2r(ax):
    y1, y2 = ax.get_ylim()
    ax1.set_ylim(f2r(y1), f2r(y2))
    ax1.figure.canvas.draw()

numall_6s=np.sum(np.array(num6s_mean),axis=(1,2))/25.
numregion_6s=np.array(numr1)+np.array(numr2)+np.array(numr3)+np.array(numr4)
perregion_6s=np.round(numregion_6s/numall_6s*100,1)
f,ax=plt.subplots(1,1)
#ax.plot(numall,freq,'o',label='Burst Numbers')
ax.barh(freq_str,np.array(numall_6s),color='r',label='Total bursts')
ax.barh(freq_str,np.array(numregion_6s),color='b',label='Bursts (Region 1 to 4)')
ax1=ax.twinx()
#ax.callbacks.connect("ylim_changed", convert_f2r)
#ax1.set_ylim([100,250]);ax.set_ylim([100,250]);ax1.set_yscale('linear')
ax1.set_ylim([0,8]);ax1.set_yscale('linear')
ax.set_ylabel('Frequency (MHz)');ax.set_xlabel('Number of Bursts');ax1.set_ylabel('Coronal Height (Mm)')
ax1.set_yticks(np.arange(8)+0.5);ax1.set_yticklabels(np.round(np.array(rall_omegap)*6.9e5/1.e3,1))
ax.set_xscale('log');ax.set_xlim([9.e0,1.1e5])
for i in range(8):
    ax.text(20, -0.5+i,str(perregion_6s[i])+'%',ha='center', va='bottom',color='white')
ax.legend()
plt.show()

for i in range(8):
    plt.plot(np.arange(1,1001),noise[i],'o-',label=str(freq[i])+' MHz')
plt.plot(np.arange(1,1000),1000/np.sqrt(np.arange(1,1000)),'o-',label='$\sqrt{N}$')
plt.legend()
plt.xlabel('Number of averaged images')
plt.ylabel('Noise (K)')
plt.show()


plt.plot(Tb_sub_std[idx],'o-')
plt.show()



k=2
lTb=-3000;hTb=20000
t=np.arange(len(Tb_region1[k]))*0.5#;Tb_region1[k][(Tb_region1[k]<lTb) | (Tb_region1[k]>hTb)]=0;Tb_region2[k][(Tb_region2[k]<lTb) | (Tb_region2[k]>hTb)]=0
t1=np.arange(1740)*0.5
#Tb_region3[k][(Tb_region3[k]<lTb) | (Tb_region3[k]>hTb)]=0;Tb_region4[k][(Tb_region4[k]<lTb) | (Tb_region4[k]>hTb)]=0;Tb_region5[k][(Tb_region5[k]<lTb) | (Tb_region5[k]>hTb)]=0
#Tb_region6[k][(Tb_region6[k]<lTb) | (Tb_region6[k]>hTb)]=0;Tb_region7[k][(Tb_region7[k]<lTb) | (Tb_region7[k]>hTb)]=0;Tb_region8[k][(Tb_region8[k]<lTb) | (Tb_region8[k]>hTb)]=0
f,ax=plt.subplots(6,1,figsize=(20,10),sharex=True)
ax[0].plot(tidx[7]*0.5,Tb_region1[7],'-',label='240 MHz');ax[0].set_ylim(-2.e3,3.0e3);ax[0].legend(loc=2,fontsize=10)
ax[0].plot(tidx[6]*0.5,Tb_region1[6],'-',label='217 MHz');ax[0].set_ylim(-2.e3,3.0e3);ax[0].legend(loc=2,fontsize=10)
ax[1].plot(tidx[7]*0.5,Tb_region2[7],'-',label='240 MHz');ax[1].set_ylim(-2.e3,2.0e4);ax[1].legend(loc=2,fontsize=10)
ax[1].plot(tidx[6]*0.5,Tb_region2[6],'-',label='217 MHz');ax[1].set_ylim(-2.e3,2.0e4);ax[1].legend(loc=2,fontsize=10)
ax[1].axvspan(1500, 1650, alpha=0.5,color='gray')
ax[2].plot(tidx[7]*0.5,Tb_region3[7],'-',label='240 MHz');ax[2].set_ylim(-2.e3,2.e4);ax[2].legend(loc=2,fontsize=10)
ax[2].plot(tidx[6]*0.5,Tb_region3[6],'-',label='217 MHz');ax[2].set_ylim(-2.e3,1.e4);ax[2].legend(loc=2,fontsize=10)
ax[3].plot(tidx[0]*0.5,Tb_region4[0],'-',label='108 MHz');ax[3].set_ylim(-2.e3,4.e3);ax[3].legend(loc=2,fontsize=10)
ax[3].plot(tidx[4]*0.5,Tb_region4[4],'-',label='179 MHz');ax[3].set_ylim(-2.e3,4.e3);ax[3].legend(loc=2,fontsize=10)
ax[3].plot(tidx[5]*0.5,Tb_region4[5],'-',label='197 MHz');ax[3].set_ylim(-2.e3,4.e3);ax[3].legend(loc=2,fontsize=10)
ax[4].plot(tidx[0]*0.5,Tb_region5[0],'-',label='108 MHz');ax[4].set_ylim(-2.e3,3.e3);ax[4].legend(loc=2,fontsize=10)
ax[4].plot(tidx[1]*0.5,Tb_region5[1],'-',label='120 MHz');ax[4].set_ylim(-1.e3,3.e3);ax[4].legend(loc=2,fontsize=10)
ax[5].plot(tidx[0]*0.5,Tb_region6[0],'-',label='108 MHz');ax[5].set_ylim(-1.e3,3.e3);ax[5].legend(loc=2,fontsize=10)
ax[5].plot(tidx[1]*0.5,Tb_region6[1],'-',label='120 MHz');ax[5].set_ylim(-1.e3,3.e3);ax[5].legend(loc=2,fontsize=10)
ax[3].annotate("", xy=(972/2, 3.e3), xytext=(600, 3.e3), arrowprops=dict(arrowstyle="->"))
ax[0].annotate("Region 1",xy=(1300,2.e3),fontsize=15)
ax[1].annotate("Region 2",xy=(1300,9.e3),fontsize=15)
ax[2].annotate("Region 3",xy=(1300,6.e3),fontsize=15)
ax[3].annotate("Region 4",xy=(1300,2.e3),fontsize=15)
ax[4].annotate("Region 5",xy=(1300,1.e3),fontsize=15)
ax[5].annotate("Region 6",xy=(1300,1.e3),fontsize=15)
ax[3].set_ylabel('$T_B$ (K)');ax[5].set_xlabel('Time (sec)')
plt.show()

omegap=9000*np.sqrt(densobs)/1.e6;omegab=2.8*bobs
f,ax=plt.subplots(1,1)
ax.plot(r[50,50],densobs[50,50],'o-',label='$n_e$')
ax1=ax.twinx()
ax1.plot(r[50,50],bobs[50,50],'o-',label='$B$')
plt.show()

f,ax=plt.subplots(2,1,figsize=(8,8));ax0=ax[0];ax1=ax[1]
ax0.plot(tidx[7]*0.5,Tb_region2[7]/1.e3,'o-',label='240 MHz')
ax0.legend();ax0.set_xlim([1500,1650])
ax0.set_ylabel('$T_B$ (kK)');ax0.set_xlabel('Time (sec)')
im1=ax1.imshow(Tb_suball[:,:,110,95]/1.e3,aspect='auto',origin=0,vmin=0.1,vmax=60);ax1.set_xlim([964,972])
ax1.set_xticks([964,966,968,970,972]);ax1.set_xticklabels(['487','488','489','490','491','492']);ax1.set_xlabel('Time (sec)');ax1.set_ylabel('Frequency (MHz)')
ax1.set_yticks(np.linspace(0,7,8));ax1.set_yticklabels(['108','120','133','145','179','197','217','240'])
ax1.annotate("Region 4", xy=(966, 6), xytext=(966, 6),color='white')
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im1, cax=cax, orientation='vertical',label='$T_B$ (kK)')
plt.show()

f,ax=plt.subplots(1,1)
ax.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax.imshow(Tb240m/np.nanmax(Tb240m),origin=0,extent=[-2500,2500,-2500,2500],alpha=0.7,cmap='jet')#;ax11.text(-1000,1100,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
ax.contour(Tb108m/np.nanmax(Tb108m),[0.001],extent=[-2500,2500,-2500,2500],  colors='tab:red',linewidths=2)#;ax11.text(-1000,1100,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
circ=Circle((800,400),150,color='k',fill=False);ax.annotate("2",xy=(792,402),fontsize=18,color='k');ax.add_patch(circ)
circ=Circle((-1050,150),150,color='k',fill=False);ax.annotate("3",xy=(-1062,152),fontsize=18,color='k');ax.add_patch(circ)
circ=Circle((-250,250),150,color='k',fill=False);ax.annotate("4",xy=(-262,252),fontsize=18,color='k');ax.add_patch(circ)
circ=Circle((-250,-400),150,color='k',fill=False);ax.annotate("5",xy=(-262,-402),fontsize=18,color='k');ax.add_patch(circ)
circ=Circle((-1200,-850),150,color='k',fill=False);ax.annotate("6",xy=(-1212,-852),fontsize=18,color='k');ax.add_patch(circ)
circ=Circle((512,-192),150,color='k',fill=False);ax.annotate("1",xy=(510,-194),fontsize=18,color='k');ax.add_patch(circ)
ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)');ax.set_xlim(-1950,1950);ax.set_ylim(-1950,1950)
plt.show()


f,ax=plt.subplots(2,3)
ax[0,0].imshow(nhist1,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,0].set_yticks(np.linspace(0,7,8));ax[0,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,0].set_xticks(np.linspace(0,8,9));ax[0,0].set_xticklabels(bins)
ax[0,0].text(6, 1, 'Region 1', horizontalalignment='center', verticalalignment='center',fontsize=20)
ax[0,1].imshow(nhist2,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,1].set_yticks(np.linspace(0,7,8));ax[0,1].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,1].set_xticks(np.linspace(0,8,9));ax[0,1].set_xticklabels(bins)
ax[0,1].text(6, 1, 'Region 2', horizontalalignment='center', verticalalignment='center',fontsize=20)
im1=ax[0,2].imshow(nhist3,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,2].set_yticks(np.linspace(0,7,8));ax[0,2].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,2].set_xticks(np.linspace(0,8,9));ax[0,2].set_xticklabels(bins)
ax[0,2].text(6, 1, 'Region 3', horizontalalignment='center', verticalalignment='center',fontsize=20)
ax[1,0].imshow(nhist4,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,0].set_yticks(np.linspace(0,7,8));ax[1,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,0].set_xticks(np.linspace(0,8,9));ax[1,0].set_xticklabels(bins)
ax[1,0].text(6, 1, 'Region 4', horizontalalignment='center', verticalalignment='center',fontsize=20)
ax[1,1].imshow(nhist5,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,1].set_yticks(np.linspace(0,7,8));ax[1,1].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,1].set_xticks(np.linspace(0,8,9));ax[1,1].set_xticklabels(bins)
ax[1,1].text(6, 1, 'Region 5', horizontalalignment='center', verticalalignment='center',fontsize=20)
ax[1,2].imshow(nhist6,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,2].set_yticks(np.linspace(0,7,8));ax[1,2].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,2].set_xticks(np.linspace(0,8,9));ax[1,2].set_xticklabels(bins)
ax[1,2].text(6, 1, 'Region 6', horizontalalignment='center', verticalalignment='center',fontsize=20)
divider = make_axes_locatable(ax[0,2])
cax = divider.append_axes('right', size='5%', pad=0.05)
f.colorbar(im1, cax=cax, orientation='vertical')
ax[0,0].set_ylabel('Frequency (MHz)');ax[1,0].set_ylabel('Frequency (MHz)');ax[1,1].set_xlabel('Duration bins (sec)')
plt.show()


from mpl_toolkits.axes_grid1 import make_axes_locatable
f,ax=plt.subplots(2,2)
im00=ax[0,0].imshow(Tb_subf[0].mean(axis=0),origin=0,aspect='auto');ax[0,0].set_title('Time averaged map')
im01=ax[0,1].imshow(Tb_sub6s[0].mean(axis=0),origin=0,aspect='auto');ax[0,1].set_title('+6-sigma Time averaged map')
im10=ax[1,0].imshow(Tb_subn6s[0].mean(axis=0),origin=0,aspect='auto');ax[1,0].set_title('-6-sigma Time averaged map')
im11=ax[1,1].imshow(num6s[0],origin=0,aspect='auto');ax[1,1].set_title('Areas of +6-sigma emission')
divider = make_axes_locatable(ax[0,0]);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im00, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax[0,1]);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im01, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax[1,0]);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im10, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax[1,1]);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im11, cax=cax, orientation='vertical')
plt.show()


f,ax=plt.subplots(2,1)
ax[0].plot(idx,Tb_sub6s[:,110,95],'o-')
ax[0].plot(idx,Tb_subn6s[:,110,95],'o-')
ax[0].plot(Tb_sub[:,110,95],'o-')
ax[1].plot(idx,Tb_sub6s[:,97,108],'o-')
ax[1].plot(idx,Tb_subn6s[:,97,108],'o-')
ax[1].plot(Tb_sub[:,97,108],'o-')
plt.show()

f,ax=plt.subplots(1,1,figsize=(8, 8))
a1=Tb_sub6s[0][8];a1[a1<1]=np.nan;a2=Tb_sub6s[4][19];a2[a2<1]=np.nan;a3=Tb_sub6s[7][11];a3[a3<1]=np.nan
im=ax.imshow(a1,extent=[-5000,5000,-5000,5000],origin=0,vmin=2500,vmax=4000,cmap='Reds',alpha=0.9)
im=ax.imshow(a2,extent=[-5000,5000,-5000,5000],origin=0,vmin=1500,vmax=3000,cmap='Greens',alpha=0.9)
im=ax.imshow(a3,extent=[-5000,5000,-5000,5000],origin=0,vmin=1500,vmax=2000,cmap='Blues',alpha=0.9)
ax.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='brown',linewidths=3)
ax.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='b',linewidths=3)
ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
ax.set_xlim(-2500,2500);ax.set_ylim(-2500,2500)
c=plt.Circle((0,0),960,fill=False,color='k');ax.add_artist(c)
#divider = make_axes_locatable(ax);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im, cax=cax, label='$T_B$ (K)', orientation='vertical')
ax.annotate('108 MHz',xy=(-50,2000),fontsize=19,color='brown')
ax.annotate('179 MHz',xy=(-50,1800),fontsize=19,color='g')
ax.annotate('240 MHz',xy=(-50,1600),fontsize=19,color='b')
plt.show()


f,ax=plt.subplots(2,4,figsize=(25, 10));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
im00=ax00.imshow(Tb_subf[0].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-200,vmax=200,cmap='coolwarm')
ax00.contour(Tb_subf[0].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[0],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],colors='white',linewidths=3)
ax00.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
im01=ax01.imshow(Tb_subf[1].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-100,vmax=100,cmap='coolwarm')
ax01.contour(Tb_subf[1].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[1],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax01.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
im02=ax02.imshow(Tb_subf[2].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-100,vmax=100,cmap='coolwarm')
ax02.contour(Tb_subf[2].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[2],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax02.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(Tb_subf[3].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-100,vmax=100,cmap='coolwarm')
ax03.contour(Tb_subf[3].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[3],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax03.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
im10=ax10.imshow(Tb_subf[4].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-50,vmax=50,cmap='coolwarm')
ax10.contour(Tb_subf[4].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[4],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax10.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
im11=ax11.imshow(Tb_subf[5].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-50,vmax=50,cmap='coolwarm')
ax11.contour(Tb_subf[5].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[5],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax11.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
im12=ax12.imshow(Tb_subf[6].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-50,vmax=50,cmap='coolwarm')
ax12.contour(Tb_subf[6].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[6],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax12.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(Tb_subf[7].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=-50,vmax=50,cmap='coolwarm')
ax13.contour(Tb_subf[7].mean(axis=0), [np.nanstd(np.nanmean(Tb_subf[7],axis=0)[0:50,:])*5.0],extent=[-5000,5000,-5000,5000],  colors='white',linewidths=3)
ax13.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax00);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im00, cax=cax, label='', orientation='vertical');divider = make_axes_locatable(ax01);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im01, cax=cax, label='', orientation='vertical')
divider = make_axes_locatable(ax02);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im02, cax=cax, label='', orientation='vertical');divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='T$_B$ (K)', orientation='vertical')
divider = make_axes_locatable(ax10);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im10, cax=cax, label='', orientation='vertical');divider = make_axes_locatable(ax11);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im11, cax=cax, label='', orientation='vertical')
divider = make_axes_locatable(ax12);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im12, cax=cax, label='', orientation='vertical');divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='T$_B$ (K)', orientation='vertical')
plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/paper_fig/meanTb.png',dpi=100)
plt.show()

Tbsub6smean=[0]*8
for i in range(8):
    Tbsub6smean[i]=Tb_sub6s[i].sum(axis=0)/num6s_mean[i]/1000;Tbsub6smean[i][0:70,:]=np.nan;Tbsub6smean[i][-70:,:]=np.nan;Tbsub6smean[i][:,0:70]=np.nan;Tbsub6smean[i][:,-70:]=np.nan
f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(Tbsub6smean[0],extent=[-5000,5000,-5000,5000],origin=0,vmin=5,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax00.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
ax01.imshow(Tbsub6smean[1],extent=[-5000,5000,-5000,5000],origin=0,vmin=5,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax01.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
ax02.imshow(Tbsub6smean[2],extent=[-5000,5000,-5000,5000],origin=0,vmin=5,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax02.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(Tbsub6smean[3],extent=[-5000,5000,-5000,5000],origin=0,vmin=5,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax03.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.imshow(Tbsub6smean[4],extent=[-5000,5000,-5000,5000],origin=0,vmin=2,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax10.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
ax11.imshow(Tbsub6smean[5],extent=[-5000,5000,-5000,5000],origin=0,vmin=2,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax11.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
ax12.imshow(Tbsub6smean[6],extent=[-5000,5000,-5000,5000],origin=0,vmin=2,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax12.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(Tbsub6smean[7],extent=[-5000,5000,-5000,5000],origin=0,vmin=2,vmax=9,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax13.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='$T_B$ (kK)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='$T_B$ (kK)', orientation='vertical')
plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/paper_fig/meanTb_6s.png',dpi=100)
plt.show()


f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(Tb_sub6s[0].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax00.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
ax01.imshow(Tb_sub6s[1].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax01.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
ax02.imshow(Tb_sub6s[2].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax02.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(Tb_sub6s[3].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax03.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.imshow(Tb_sub6s[4].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax10.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
ax11.imshow(Tb_sub6s[5].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax11.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
ax12.imshow(Tb_sub6s[6].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax12.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(Tb_sub6s[7].std(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=100,vmax=900,cmap='YlOrRd',norm=matplotlib.colors.LogNorm())
ax13.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='$T_B$ (K)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='$T_B$ (K)', orientation='vertical')
plt.show()

f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(num6s_mean[0],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax00.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
ax01.imshow(num6s_mean[1],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax01.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
ax02.imshow(num6s_mean[2],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax02.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(num6s_mean[3],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax03.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.imshow(num6s_mean[4],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax10.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
ax11.imshow(num6s_mean[5],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax11.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
ax12.imshow(num6s_mean[6],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax12.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(num6s_mean[7],extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='YlOrBr',norm=matplotlib.colors.LogNorm())
ax13.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='Number', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='Number', orientation='vertical')
plt.show()

f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(Tb_sub6s[0].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax00.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
ax01.imshow(Tb_sub6s[1].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax01.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
ax02.imshow(Tb_sub6s[2].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax02.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(Tb_sub6s[3].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax03.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.imshow(Tb_sub6s[4].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax10.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
ax11.imshow(Tb_sub6s[5].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax11.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
ax12.imshow(Tb_sub6s[6].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax12.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(Tb_sub6s[7].mean(axis=0),extent=[-5000,5000,-5000,5000],origin=0,vmin=1,vmax=90,cmap='coolwarm',norm=matplotlib.colors.LogNorm())
ax13.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='g',linewidths=3)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='$T_B$ (K)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='$T_B$ (K)', orientation='vertical')
plt.show()

i=10
for i in range(3550):
    ii="%04d"%i
    f,ax=plt.subplots(2,2,figsize=(10, 10),sharex=True,sharey=True);ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
    ax00.imshow(d,aspect='auto',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax01.imshow(d,aspect='auto',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax10.imshow(d,aspect='auto',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax11.imshow(d,aspect='auto',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax00.contour(Tb_suball[7][i],np.array([8,9,10])*np.nanstd(Tb_suball[7][i]),extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=4,alpha=0.7);ax00.text(-1700,-1700,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue');ax00.set_xlim([-1800,1800]);ax00.set_ylim([-1800,1800])
    ax00.contour(Tb_suball[6][i],np.array([8,9,10])*np.nanstd(Tb_suball[6][i]),extent=[-5000,5000,-5000,5000],  colors='tab:orange',linewidths=4,alpha=0.7);ax00.text(-1700,1500,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange');ax00.set_xlim([-1800,1800]);ax00.set_ylim([-1800,1800])
    ax00.contour(Tb240m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:blue',linewidths=1,linestyles='--');ax00.contour(Tb217m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:orange',linewidths=1,linestyles='--')
    ax01.contour(Tb_suball[5][i],np.array([8,9,10])*np.nanstd(Tb_suball[5][i]),extent=[-5000,5000,-5000,5000],  colors='tab:cyan',linewidths=4,alpha=0.7);ax01.text(-1700,-1700,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan');ax01.set_xlim([-1800,1800]);ax01.set_ylim([-1800,1800])
    ax01.contour(Tb_suball[4][i],np.array([8,9,10])*np.nanstd(Tb_suball[4][i]),extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=4,alpha=0.7);ax01.text(-1700,1500,'179 MHz',fontweight='bold',fontsize=15,color='tab:green');ax01.set_xlim([-1800,1800]);ax01.set_ylim([-1800,1800])
    ax01.contour(Tb197m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:cyan',linewidths=1,linestyles='--');ax01.contour(Tb179m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:green',linewidths=1,linestyles='--')
    ax10.contour(Tb_suball[3][i],np.array([8,9,10])*np.nanstd(Tb_suball[3][i]),extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=4,alpha=0.7);ax10.text(-1700,-1700,'145 MHz',fontweight='bold',fontsize=15,color='tab:red');ax10.set_xlim([-1800,1800]);ax10.set_ylim([-1800,1800])
    ax10.contour(Tb_suball[2][i],np.array([8,9,10])*np.nanstd(Tb_suball[2][i]),extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=4,alpha=0.7);ax10.text(-1700,1500,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple');ax10.set_xlim([-1800,1800]);ax10.set_ylim([-1800,1800])
    ax10.contour(Tb145m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:red',linewidths=1,linestyles='--');ax10.contour(Tb133m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:purple',linewidths=1,linestyles='--')
    ax11.contour(Tb_suball[1][i],np.array([8,9,10])*np.nanstd(Tb_suball[1][i]),extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=4,alpha=0.7);ax11.text(-1700,-1700,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink');ax11.set_xlim([-1800,1800]);ax11.set_ylim([-1800,1800])
    ax11.contour(Tb_suball[0][i],np.array([8,9,10])*np.nanstd(Tb_suball[0][i]),extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=4,alpha=0.7);ax11.text(-1700,1500,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive');ax11.set_xlim([-1800,1800]);ax11.set_ylim([-1800,1800])
    ax11.contour(Tb120m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:pink',linewidths=1,linestyles='--');ax11.contour(Tb108m, [2.e4],extent=[-2500,2500,-2500,2500],  colors='tab:olive',linewidths=1,linestyles='--')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=0.0);ax00.set_title('Time: '+str(i*0.5)+' sec')
    ax00.set_ylabel('Solar Y (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax10.set_xlabel('Solar X (arcsec)')
    plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/pngs_running/i_'+str(ii)+'.png',dpi=70)
    plt.close()



f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(Tb108m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='white')
ax01.imshow(Tb120m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='white')
ax02.imshow(Tb133m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='white')
im03=ax03.imshow(Tb145m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='white')
ax10.imshow(Tb179m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='white')
ax11.imshow(Tb197m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='white')
ax12.imshow(Tb217m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='white')
im13=ax13.imshow(Tb240m/1.e6,extent=[-2500,2500,-2500,2500],origin=0,vmin=0,vmax=0.4,cmap='jet')
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='white')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='T$_B$ (MK)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='T$_B$ (MK)', orientation='vertical')
plt.show()

l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype(np.float16);dc[np.where(dc==0)]=np.nan
levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.75,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
f,ax=plt.subplots(2,4,figsize=(16, 8));ax13=ax[0,0];ax12=ax[0,1];ax11=ax[0,2];ax10=ax[0,3]
ax03=ax[1,0];ax02=ax[1,1];ax01=ax[1,2];ax00=ax[1,3]
ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax00.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax00.contour(Tb_subf[7].mean(axis=0)/np.nanmax(Tb_subf[7].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=4,alpha=0.7);ax00.text(-1000,-1100,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax01.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax01.contour(Tb_subf[6].mean(axis=0)/np.nanmax(Tb_subf[6].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:orange',linewidths=4,alpha=0.7);ax01.text(-1000,-1100,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax02.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax02.contour(Tb_subf[5].mean(axis=0)/np.nanmax(Tb_subf[5].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:cyan',linewidths=4,alpha=0.8);ax02.text(-1000,-1100,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax03.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax03.contour(Tb_subf[4].mean(axis=0)/np.nanmax(Tb_subf[4].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=4);ax03.text(-1000,-1100,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax10.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax10.contour(Tb_subf[3].mean(axis=0)/np.nanmax(Tb_subf[3].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=4,alpha=0.8);ax10.text(-1000,-1100,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax11.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax11.contour(Tb_subf[2].mean(axis=0)/np.nanmax(Tb_subf[2].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=4);ax11.text(-1000,-1100,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax12.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax12.contour(Tb_subf[1].mean(axis=0)/np.nanmax(Tb_subf[1].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=4);ax12.text(-1000,-1100,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax13.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax13.contour(Tb_subf[0].mean(axis=0)/np.nanmax(Tb_subf[0].mean(axis=0)),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=4);ax13.text(-1000,-1100,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
ax00.set_xlim(-1500,1500);ax00.set_ylim(-1500,1500);ax01.set_xlim(-1500,1500);ax01.set_ylim(-1500,1500);ax02.set_xlim(-1500,1500);ax02.set_ylim(-1500,1500);ax03.set_xlim(-1500,1500);ax03.set_ylim(-1500,1500);ax10.set_xlim(-1500,1500);ax10.set_ylim(-1500,1500);ax11.set_xlim(-1500,1500);ax11.set_ylim(-1500,1500);ax12.set_xlim(-1500,1500);ax12.set_ylim(-1500,1500);ax13.set_xlim(-1500,1500);ax13.set_ylim(-1500,1500)
ax00.set_xlabel('Solar X (arcsec)');ax01.set_xlabel('Solar X (arcsec)');ax02.set_xlabel('Solar X (arcsec)');ax03.set_xlabel('Solar X (arcsec)');ax13.set_ylabel('Solar Y (arcsec)');ax03.set_ylabel('Solar Y (arcsec)')
#plt.xlim(-3500,3500);plt.ylim(-3500,3500);ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
plt.show()


l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype(np.float16);dc[np.where(dc==0)]=np.nan
levels_sub=[0.1,0.5,0.9];levels_sub1=[0.5,0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
f,ax=plt.subplots(2,4,figsize=(16, 8));ax13=ax[0,0];ax12=ax[0,1];ax11=ax[0,2];ax10=ax[0,3]
ax03=ax[1,0];ax02=ax[1,1];ax01=ax[1,2];ax00=ax[1,3]
ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax00.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax00.contour(Tb_sub6s[7].mean(axis=0)/np.nanmax(Tb_sub6s[7].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=4,alpha=0.7);ax00.text(-1000,1100,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax01.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax01.contour(Tb_sub6s[6].mean(axis=0)/np.nanmax(Tb_sub6s[6].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:orange',linewidths=4,alpha=0.7);ax01.text(-1000,1100,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax02.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax02.contour(Tb_sub6s[5].mean(axis=0)/np.nanmax(Tb_sub6s[5].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:cyan',linewidths=4,alpha=0.8);ax02.text(-1000,1100,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax03.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax03.contour(Tb_sub6s[4].mean(axis=0)/np.nanmax(Tb_sub6s[4].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=4);ax03.text(-1000,1100,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax10.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax10.contour(Tb_sub6s[3].mean(axis=0)/np.nanmax(Tb_sub6s[3].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=4,alpha=0.8);ax10.text(-1000,1100,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax11.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax11.contour(Tb_sub6s[2].mean(axis=0)/np.nanmax(Tb_sub6s[2].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=4);ax11.text(-1000,1100,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax12.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax12.contour(Tb_sub6s[1].mean(axis=0)/np.nanmax(Tb_sub6s[1].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=4);ax12.text(-1000,1100,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax13.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax13.contour(Tb_sub6s[0].mean(axis=0)/np.nanmax(Tb_sub6s[0].mean(axis=0)),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=4);ax13.text(-1000,1100,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
ax00.set_xlim(-1500,1500);ax00.set_ylim(-1500,1500);ax01.set_xlim(-1500,1500);ax01.set_ylim(-1500,1500);ax02.set_xlim(-1500,1500);ax02.set_ylim(-1500,1500);ax03.set_xlim(-1500,1500);ax03.set_ylim(-1500,1500);ax10.set_xlim(-1500,1500);ax10.set_ylim(-1500,1500);ax11.set_xlim(-1500,1500);ax11.set_ylim(-1500,1500);ax12.set_xlim(-1500,1500);ax12.set_ylim(-1500,1500);ax13.set_xlim(-1500,1500);ax13.set_ylim(-1500,1500)
ax00.set_xlabel('Solar X (arcsec)');ax01.set_xlabel('Solar X (arcsec)');ax02.set_xlabel('Solar X (arcsec)');ax03.set_xlabel('Solar X (arcsec)');ax13.set_ylabel('Solar Y (arcsec)');ax03.set_ylabel('Solar Y (arcsec)')
#plt.xlim(-3500,3500);plt.ylim(-3500,3500);ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
plt.show()






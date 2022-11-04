import numpy as np
import glob
import pickle
import matplotlib.pyplot as plt

b='000-010';fid='187-188'
fl=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
bl=['000-008','000-009','000-010','008-009','008-010','009-010']
freq=np.array([108,120,133,145,161,179,197,217,240])
Tb=[0]*6;nccf=[0]*6;corr=[0]*6
k=0
for b in bl:
    j=0;Tb[k]=[0]*9;nccf[k]=[0]*9;corr[k]=[0]*9
    for fid in fl:
        listp=sorted(glob.glob('/media/rohit/MWA/20140914/pickle/f*'+str(fid)+'*T'+b+'*.p'))
        Tb[k][j]=[0]*len(listp);nccf[k][j]=[0]*len(listp);corr[k][j]=[0]*len(listp)
        for i in range(len(listp)):
            aa=pickle.load(open(listp[i],'rb'),encoding='bytes')
            Tb[k][j][i]=aa[17][3][0][1:3].mean(axis=0);nccf[k][j][i]=aa[5][0][1:3].mean(axis=0);corr[k][j][i]=aa[17][2][0][1:3].mean(axis=0)
        Tb[k][j]=np.array(Tb[k][j]);nccf[k][j]=np.array(nccf[k][j]);corr[k][j]=np.array(corr[k][j])
        j=j+1
    k=k+1
Tb=np.array(Tb);corr=np.array(corr);nccf=np.array(nccf)
Tb=Tb.reshape(6,9,1325).mean(axis=0)
#pickle.dump([central_freq,chan_list,auto_t1,auto_t2,cross,ncross,phase_ncross,u,v,w,azi_pointing,ele_pointing,ph_az,ph_el,start_time,mid_time,end_time,[0,0,corr_factor,S_sun,T_sun,Un_Tbeam_Sun,Temp_beam_sun,fringe_factor,T_baseline,Tsky_integrated]]

Tbimg=[0]*3;Tbimgmax=[0]*3
for i in range(3):
    bb=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_'+fl[i*4]+'.p','rb'),encoding='bytes')
    Tbimg[i]=np.array(bb[0]);Tbimgmax[i]=np.array(bb[0]).max(axis=(1,2))

############## LASCO

lasco=pickle.load(open('/media/rohit/MWA/20140914/lasco/lasco_maps.p','rb'),encoding='bytes')


for i in range(3):
    plt.plot(np.arange(len(Tbimgmax[i]))*12,Tbimgmax[i]/1.e6,label=str(freq[i*4])+' MHz')
plt.xlabel('Time (s)');plt.ylabel('T$_B$ (MK)');plt.legend()
plt.show()


#plt.plot(corr[0],'o-');plt.show()
plt.imshow(Tb,origin='lower',aspect='auto',interpolation='None',vmax=200,vmin=0)
plt.colorbar(label='Flux (SFU)');plt.yticks(np.arange(9),freq);
plt.xticks([0,150,300,450,600,750,900,1050,1200],['01:30','02:00','02:30','03:00','03:30','04:00','04:30','05:00','05:30'])
plt.xlabel('Time (HH:MM UT)');plt.ylabel('Frequency (MHz)')
plt.show()

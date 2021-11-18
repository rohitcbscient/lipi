import os
import numpy as np
import matplotlib.pyplot as plt
import pickle

f=['L190330','L190331','L190401','L190402','L190403','L190404','L190405','L190406','L190407','L190408','L190409','L190410','L190411']
freq=[245,410,610,1415,2695,4995,8800,15400]
time=[0]*len(f);time_s=[0]*len(f);flux_245=[0]*len(f)
flux_410=[0]*len(f);flux_610=[0]*len(f);flux_1415=[0]*len(f)
flux_2695=[0]*len(f);flux_4995=[0]*len(f);flux_8800=[0]*len(f);flux_15400=[0]*len(f)
j=0
for filename in f:
    print 'Reading '+str(filename)+'.SRD ....'
    ff=open('/media/rohit/MWA/RSTN_2019/'+str(filename)+'.SRD',"r")
    l = ff.readlines()
    n=len(l)
    time=[0]*n
    time_s=[0]*n
    flux_245[j]=[0]*n;flux_410[j]=[0]*n;flux_610[j]=[0]*n;flux_1415[j]=[0]*n
    flux_2695[j]=[0]*n;flux_4995[j]=[0]*n;flux_8800[j]=[0]*n;flux_15400[j]=[0]*n
    for i in range(n):
            time[i]=l[i].split(' ')[0]
            time_s[i]=float(l[i].split(' ')[0][0:2])*3600+float(l[i].split(' ')[0][2:4])*60+float(l[i].split(' ')[0][4:6])
            l_=l[i].split(' ')
            if(l_[1][0]!='/'):
                flux_245[j][i]=int(l_[1])
            else:
                flux_245[j][i]=np.nan
            if(l_[2][0]!='/'):
                flux_410[j][i]=int(l_[2])
            else:
                flux_410[j][i]=np.nan
            if(l_[3][0]!='/'):
                flux_610[j][i]=int(l_[3])
            else:
                flux_610[j][i]=np.nan
            if(l_[4][0]!='/'):
                flux_1415[j][i]=int(l_[4])
            else:
                flux_1415[j][i]=np.nan
            if(l_[5][0]!='/'):
                flux_2695[j][i]=int(l_[5])
            else:
                flux_2695[j][i]=np.nan
            if(l_[6][0]!='/'):
                flux_4995[j][i]=int(l_[6])
            else:
                flux_4995[j][i]=np.nan
            if(l_[7][0]!='/'):
                flux_8800[j][i]=int(l_[7])
            else:
                flux_8800[j][i]=np.nan
            if(l_[8][0]!='/'):
                flux_15400[j][i]=int(l_[8])
            else:
                flux_15400[j][i]=np.nan
    j=j+1

flux_245=np.array(flux_245);flux_410=np.array(flux_410);flux_610=np.array(flux_610)
flux_1415=np.array(flux_1415);flux_2695=np.array(flux_2695);flux_4995=np.array(flux_4995);flux_8800=np.array(flux_8800);flux_15400=np.array(flux_15400)

plt.plot(np.arange(len(f)),np.nanmean(flux_245,axis=1),'o',color='k')
plt.errorbar(np.arange(len(f)),np.nanmean(flux_245,axis=1),yerr=np.nanstd(flux_245,axis=1),color='k',label='RSTN Flux at 245 MHz')
#plt.plot(np.arange(len(f)),np.nanmean(flux_410,axis=1),'o',color='r')
#plt.errorbar(np.arange(len(f)),np.nanmean(flux_410,axis=1),yerr=np.nanstd(flux_410,axis=1),color='r',label='410 MHz')
plt.plot([6],[16.5],'s',color='red',markersize=15,label='MWA Flux Density at 240 MHz')
#plt.plot(np.arange(len(f)),np.nanmean(flux_610,axis=1),'o',color='r')
#plt.errorbar(np.arange(len(f)),np.nanmean(flux_610,axis=1),yerr=np.nanstd(flux_610,axis=1),color='r',label='610 MHz')
plt.xlim(-0.5,13.5)
plt.legend(loc=2)
plt.xticks(np.arange(13),['20190330','20190331','20190401','20190402','20190403','20190404','20190405','20190406','20190407','20190408','20190409','20190410','20190411'])
plt.xlabel('Time (dates)')
plt.ylabel('Flux (SFU)')
plt.show()

sys.exit()
mean_flux=np.median(np.array(flux_245)[np.where(np.array(flux_245)>0)])
print 'Median Flux for date '+str(date)+' is '+str(mean_flux)+' SFU at 245 MHz'
filename_=filename.split('L')[1]
#pickle.dump([time,time_s,flux_245],open('/Data/rohit/reserve/RSTN_data/RSTN_20'+str(filename_)+'.p','wb'))
#plt.plot(np.array(time_s),np.array(flux_245),'o-')
#plt.show()

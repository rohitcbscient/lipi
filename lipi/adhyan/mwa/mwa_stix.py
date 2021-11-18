import json
import numpy as np
import matplotlib.pyplot as plt
import glob
import csv
import pickle
from scipy.io import readsav

def read_csv(f):
    with open(f, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            print(f'\t{row["name"]} works in the {row["department"]} department, and was born in {row["birthday month"]}.')
            line_count += 1
        print(f'Processed {line_count} lines.')

with open('/media/rohit/MWA/MWA_STIX_2020/xrays-7-day.json') as json_file:
        data = json.load(json_file)

n=len(data);j1=0;j2=0
time_s=[0]*int(n/2);time_str=[0]*int(n/2);en_18=[0]*int(n/2);en_54=[0]*int(n/2)
for i in range(n):
    day=int(data[i]['time_tag'].split('T')[0].split('-')[-1])-15
    if(data[i]['energy']=='0.05-0.4nm'):
        en_54[j1]=data[i]['flux']
        time_str[j1]=data[i]['time_tag'].split('T')[1]
        time_s[j1]=day*86400+int(time_str[j1].split(':')[0])*3600+int(time_str[j1].split(':')[1])*60
        j1=j1+1
    if(data[i]['energy']=='0.1-0.8nm'):
        en_18[j2]=data[i]['flux']
        j2=j2+1

time_s=np.array(time_s);en_18=np.array(en_18);en_54=np.array(en_54)

csv16=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1612020.txt')
csv17=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1712020.txt')
csv18=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1812020.txt')
csv19=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1912020.txt')
csv20=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_2012020.txt')
csv21=np.loadtxt('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_2112020.txt')
csv21=csv21[0:7200,:]
# 4 - 10 keV - Counts / 4 s 10 - 15 keV - Counts / 4 s  15 - 25 keV - Counts / 4 s  25 - 50 keV - Counts / 4 s  50 - 84 keV - Counts / 4 s
time_s=time_s[1133:9773];en_18=en_18[1133:9773];en_54=en_54[1133:9773]
# 1133 is 00:00 UT for 16th Nov 2020
time_s=time_s-time_s[0]
time_s=time_s.reshape(6,1440)[:,0:60*8]
en_18=en_18.reshape(6,1440)[:,0:60*8];en_54=en_54.reshape(6,1440)[:,0:60*8]

sublist94=sorted(glob.glob('/media/rohit/MWA/20201116_STIX_EUV/AR1/*94*sub.sav'))
sublist171=sorted(glob.glob('/media/rohit/MWA/20201116_STIX_EUV/AR1/*171*sub.sav'))
bdiff94=[0]*len(sublist94);submap94=[0]*len(sublist94)
bdiff171=[0]*len(sublist171);submap171=[0]*len(sublist171)
for i in range(len(sublist94)):
    subdata=readsav(sublist94[i])
    #submap94[i]=subdata['drot_map_'][0][0]
    submap94[i]=subdata['submap'][0][0]
    submap94[i]=submap94[i][0:668,0:667]
    bdiff94[i]=submap94[i]-submap94[0]
for i in range(len(sublist171)):
    subdata=readsav(sublist171[i])
    #submap171[i]=subdata['drot_map_'][0][0]
    submap171[i]=subdata['submap'][0][0]
    submap171[i]=submap171[i][0:668,0:667]
    bdiff171[i]=submap171[i]-submap171[0]

submap94=np.array(submap94);bdiff94=np.array(bdiff94)
max94=np.max(submap94,axis=(1,2))
submap171=np.array(submap171);bdiff171=np.array(bdiff171)
max171=np.max(submap171,axis=(1,2))
mean171=np.mean(submap171[:,250:550,250:550],axis=(1,2))
mean94=np.mean(submap94[:,250:550,250:550],axis=(1,2))
tsaia94=23000+np.arange(330)*12;max171p=max171[115:];max94p=max94[115:];mean171p=mean171[115:];mean94p=mean94[115:]
tsaia171=23000+np.arange(335)*12

freq=np.array([  79.36,   85.76,   93.44,   99.84,  107.52,  113.92,  121.6 ,128.  ,  135.68,  142.08,  149.76,  156.16,  163.84,  170.24,177.92,  184.32,  192.  ,  198.4 ,  206.08,  212.48,  220.16,226.56,  234.24,  240.64])
mwa=pickle.load(open('/media/rohit/MWA/MWA_STIX_2020/0600-0700/ds_median.p','rb'),encoding="bytes")
mwap=mwa[:,41:];tsmwa=23000+np.arange(169)*10

#################################################### Plotting #########################################################33

t8=np.arange(7200)*4
i=0;j=22
plt.plot(t8,csv16[:,0]/np.max(csv16[:,0]),'o-',label='STIX (4-10 keV)')
plt.plot(t8,csv16[:,1]/np.max(csv16[:,1]),'o-',label='STIX (10-15 keV)')
plt.plot(t8,csv16[:,2]/np.max(csv16[:,2]),'o-',label='STIX (15-25 keV)')
plt.plot(tsmwa,mwap[j]/np.nanmax(mwap[j])/3.0+0.5,'o-',label='MWA '+str(freq[j])+' MHz')
plt.plot(tsaia94,mean94p-mean94p[0],'o-',label='AIA 94 $\AA$')
#plt.plot(tsaia171,mean171p-mean171p[0],'o-',label='AIA 171 $\AA$')
plt.plot(time_s[i]+35,en_54[i]/np.max(en_54[i]),'o-',label='GOES (0.5-4.0 $\AA$)')
plt.plot(time_s[i]+35,en_18[i]/np.max(en_18[i]),'o-',label='GOES (1.0-8.0 $\AA$)')
plt.legend(loc=2)
plt.xlabel('Time (sec)')
plt.title('16-11-2020')
plt.xlim(23000,25000) # 23000 --> 06:23:20
plt.xticks(23000+np.array([0,240,480,720,960,1200,1440,1680,1920]),['06:23:20','06:27:20','06:31:20','06:35:20','06:39:20','06:43:20','06:47:20','06:51:20','06:55:20'])
plt.show()


plt.imshow(mwap,aspect='auto',origin=0,interpolation='none',cmap='YlGnBu',vmin=0,vmax=0.03)
plt.yticks(np.arange(24)[::4],freq[::4])
plt.xticks(np.array([0,24,48,72,96,120,144,168,192]),['06:23:20','06:27:20','06:31:20','06:35:20','06:39:20','06:43:20','06:47:20','06:51:20','06:55:20'])
plt.xlabel('Time (HH:MM:SS) UT');plt.ylabel('Freqeuncy (MHz)')
plt.show()

#with open('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1612020.csv', 'rb') as csvfile_16:
#    csv16=csv.reader(csvfile_16)

#with open('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1712020.csv', 'rb') as csvfile_17
#with open('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1812020.csv', 'rb') as csvfile_18
#with open('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_1912020.csv', 'rb') as csvfile_19
#with open('/media/rohit/MWA/MWA_STIX_2020/stix_quick_look_light_curves_2012020.csv', 'rb') as csvfile_20



import numpy as np
import matplotlib.pyplot as plt


filename='/media/rohit/MWA/20140914/rstn_20140914.SRD'
print 'Reading '+str(filename)+'....'
ff=open(str(filename),"r")
l = ff.readlines()
n=30700
time=[0]*n
time_s=[0]*n
flux_245=[0]*n;flux_410=[0]*n;flux_610=[0]*n;flux_1415=[0]*n
flux_2695=[0]*n;flux_4995=[0]*n;flux_8800=[0]*n;flux_15400=[0]*n
for i in range(n):
        time[i]=l[i].split(' ')[0]
        time_s[i]=float(l[i].split(' ')[0][0:2])*3600+float(l[i].split(' ')[0][2:4])*60+float(l[i].split(' ')[0][4:6])
        l_=l[i].split(' ')
        if(l_[1][0]!='/'):
            flux_245[i]=float(l_[1][0]+'.'+l_[1][1]+l_[1][2])*10**int(l_[1][3])
        else:
            flux_245[i]=np.nan
        if(l_[2][0]!='/'):
            flux_410[i]=float(l_[2][0]+'.'+l_[2][1]+l_[2][2])*10**int(l_[2][3])
        else:
            flux_410[i]=np.nan
        if(l_[3][0]!='/'):
            flux_610[i]=float(l_[3][0]+'.'+l_[3][1]+l_[3][2])*10**int(l_[3][3])
        else:
            flux_610[i]=np.nan
        if(l_[4][0]!='/'):
            flux_1415[i]=float(l_[4][0]+'.'+l_[4][1]+l_[4][2])*10**int(l_[4][3])
        else:
            flux_1415[i]=np.nan
        if(l_[5][0]!='/'):
            flux_2695[i]=float(l_[5][0]+'.'+l_[5][1]+l_[5][2])*10**int(l_[5][3])
        else:
            flux_2695[i]=np.nan
        if(l_[6][0]!='/'):
            flux_4995[i]=float(l_[6][0]+'.'+l_[6][1]+l_[6][2])*10**int(l_[6][3])
        else:
            flux_4995[i]=np.nan
        if(l_[7][0]!='/'):
            flux_8800[i]=float(l_[7][0]+'.'+l_[7][1]+l_[7][2])*10**int(l_[7][3])
        else:
            flux_8800[i]=np.nan
        if(l_[8][0]!='/'):
            flux_15400[i]=float(l_[8][0]+'.'+l_[8][1]+l_[8][2])*10**int(l_[8][3])
        else:
            flux_15400[i]=np.nan

flux_245=np.array(flux_245);flux_410=np.array(flux_410);flux_610=np.array(flux_610)
flux_1415=np.array(flux_1415);flux_2695=np.array(flux_2695);flux_4995=np.array(flux_4995);flux_8800=np.array(flux_8800);flux_15400=np.array(flux_15400)
flux_all=np.array([flux_245,flux_410,flux_610,flux_1415,flux_2695,flux_4995,flux_8800,flux_15400])
flux_all_diff=flux_all-np.array(list(flux_all[:,0])*n).reshape(n,8).swapaxes(0,1)
freq_rstn=[245,410,610,1415,2695,4995,8800,15400]

plt.plot(time_s,flux_245,'-',label='245 MHz')
plt.plot(time_s,flux_410,'-',label='410 MHz')
plt.plot(time_s,flux_610,'-',label='610 MHz')
plt.plot(time_s,flux_1415,'-',label='1415 MHz')
plt.plot(time_s,flux_2695,'-',label='2695 MHz')
plt.plot(time_s,flux_4995,'-',label='4995 MHz')
plt.plot(time_s,flux_8800,'-',label='8800 MHz')
plt.plot(time_s,flux_15400,'-',label='15400 MHz')
plt.axvline(x=8130,linestyle='--',label='Maximum at 1.5 GHz',color='k')
plt.legend(),plt.xlabel('Time (s)');plt.ylabel('Solar Flux Value (SFU)');plt.show()


plt.plot(freq_rstn,flux_all_diff[:,8132],'o-')
plt.xlabel('Frequency (MHz)');plt.ylabel('Solar Flux Value (SFU)')
plt.show()


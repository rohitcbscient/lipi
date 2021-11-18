import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from scipy.optimize import curve_fit

data=readsav('/home/i4ds1807205/20151203/20151203_240MHz_psimas.sav')
dens=data['modsolstruct']['dens'][0]
temp=data['modsolstruct']['temp'][0]
br=data['modsolstruct']['br'][0]
bth=data['modsolstruct']['bth'][0]
bph=data['modsolstruct']['bph'][0]
r=data['gridpramsstruct']['rpos'][0]

npr=r[70:,50]-1
spr=r[:29,50]-1
np_dens=dens[70:,50]
sp_dens=dens[:29,50]
np_temp=temp[70:,50]
sp_temp=temp[:29,50]

np_densfit=curve_fit(lambda t,a,b,c: a*np.log(c*t)+b,  npr,  np.log(np_dens),p0=(-6,17,1))
np_tempfit=curve_fit(lambda t,a,b,c: a*np.tanh(b*t)+c, npr, np_temp,p0=(1.e6,1,0))
sp_densfit=curve_fit(lambda t,a,b,c: a*np.log(c*t)+b,  spr,  np.log(sp_dens),p0=(-6,17,1))
sp_tempfit=curve_fit(lambda t,a,b,c: a*np.tanh(b*t)+c, spr, sp_temp,p0=(1.e6,1,0))

xt=np.linspace(0,3,100000)
np_densfity=np.exp(np_densfit[0][0]*np.log(xt*np_densfit[0][2])+np_densfit[0][1])
np_tempfity=np_tempfit[0][0]*np.tanh(xt*np_tempfit[0][1])+np_tempfit[0][2]
sp_densfity=np.exp(sp_densfit[0][0]*np.log(xt*sp_densfit[0][2])+sp_densfit[0][1])
sp_tempfity=sp_tempfit[0][0]*np.tanh(xt*sp_tempfit[0][1])+sp_tempfit[0][2]


tcor=120
tpho=1.0
ztr=12.5
wtr=0.6
dx=0.2
gm=5./3.
x=np.arange(100000)*dx
gx=-1/gm*np.cos(0.5*np.pi*x/x[-1])
te=tpho+0.5*(tcor-tpho)*(np.tanh((x-ztr)/wtr)+1)
ro=[0]*(len(x)-1)
ro[0]=1
for i in range(1,len(x)-1):
    #ro[i]=ro[i-1]*(te[i-1]+0.5*gm*gx[i-1]*dx)/(te[i] - 0.5*gm*gx[i-1]*dx)
    if(x[i]>ztr):
        ro[i]=np.exp(np_densfit[0][0]*np.log(x[i]*(1/3500.)*np_densfit[0][2])+np_densfit[0][1])/1.e17
    else:
        ro[i]=(te[i-1]+0.5*gm*gx[i-1]*dx)/(te[i] - 0.5*gm*gx[i-1]*dx)*ro[i-1]


plot_densfit=1
if(plot_densfit):
    plt.plot(npr,np_dens,'o-')
    plt.plot(spr,sp_dens,'o-')
    plt.plot(xt,np_densfity)
    plt.plot(x[1:]*2.e7/7.e10,np.array(ro)*1.e17,'--')
    plt.show()

plot_tempfit=0
if(plot_tempfit):
    plt.plot(npr,np_temp,'o-')
    plt.plot(spr,sp_temp,'o-')
    plt.plot(x,tempfity)
    plt.show()

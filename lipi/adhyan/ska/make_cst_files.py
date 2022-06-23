# Script to write the cst files
import numpy as np

line1='Theta [deg.]  Phi   [deg.]  Abs(Dir.)[dBi   ]   Abs(Theta)[dBi   ]  Phase(Theta)[deg.]  Abs(Phi  )[dBi   ]  Phase(Phi  )[deg.]  Ax.Ratio[dB    ]    '
line2='------------------------------------------------------------------------------------------------------------------------------------------------------'
aa=np.loadtxt('run5.cst',skiprows=2)
np.savetxt('run6.cst',aa, delimiter=" ", header=line1+"\n"+line2,comments='')
theta=aa[:,0].reshape(360,181);phi=aa[:,1].reshape(360,181);absdir=aa[:,2].reshape(360,181)
abstheta=aa[:,3].reshape(360,181);phasetheta=aa[:,4].reshape(360,181);absphi=aa[:,5].reshape(360,181);phasephi=aa[:,6].reshape(360,181)
ax_ratio=aa[:,7].reshape(360,181)

theta_pb=57.5 # arcmin
beam=np.zeros((360*181,8))
thphi_arr=np.meshgrid(np.linspace(0,181,181),np.linspace(0,360,360))
beam[:,0]=thphi_arr[1].flatten();beam[:,1]=thphi_arr[0].flatten()
ab=(np.cos(1.189*np.pi*((thphi_arr[0].flatten()-90)/(theta_pb/60.)))/(1-4*(1.189*(thphi_arr[0].flatten()-90)/(theta_pb/60.))**2))**2;beam[:,2]=ab
np.savetxt('run6.cst',beam, delimiter=" ", header=line1+"\n"+line2,comments='')

ab=(np.cos(1.189*np.pi*(rho/(theta_pb/60.)))/(1-4*(1.189*rho/(theta_pb/60.))**2))**2
plt.plot(rho,ab,'o-')
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.pcolormesh(phi[:,0]*np.pi/180, theta[0][0:90]*np.pi/180, absdir.swapaxes(0,1)[0:90,:]) #X,Y & data2D must all be same dimensions
plt.show()

f=plt.figure()
ax = f.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.plot(phi[:,0]*np.pi/180,absdir[:,0]*np.pi/180)
plt.show()

ax.plot(aa[:,3],projection='polar')



# Script to write the cst files
import numpy as np

line1='Theta [deg.]  Phi   [deg.]  Abs(Dir.)[dBi   ]   Abs(Theta)[dBi   ]  Phase(Theta)[deg.]  Abs(Phi  )[dBi   ]  Phase(Phi  )[deg.]  Ax.Ratio[dB    ]    '
line2='------------------------------------------------------------------------------------------------------------------------------------------------------'
aa=np.loadtxt('run5.cst',skiprows=2)
np.savetxt('run6.cst',aa, delimiter=" ", header=line1+"\n"+line2,comments='')
theta=aa[:,0].reshape(360,181);phi=aa[:,1].reshape(360,181);absdir=aa[:,2].reshape(360,181)
abstheta=aa[:,3].reshape(360,181);phasetheta=aa[:,4].reshape(360,181);absphi=aa[:,5].reshape(360,181);phasephi=aa[:,6].reshape(360,181)
ax_ratio=aa[:,7].reshape(360,181)

fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.pcolormesh(phi[:,0], theta[0], absdir.swapaxes(0,1)) #X,Y & data2D must all be same dimensions
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(aa[:,3],projection='polar')



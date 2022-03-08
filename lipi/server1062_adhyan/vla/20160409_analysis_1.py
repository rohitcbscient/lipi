import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/nas08-data02/rohit/scripts/general/plt_style.py')
data=np.load('/nas08-data02/rohit/20160409/2200/sun_L_20160409T220400-222600UT.50ms.ms.dspec.npz')
dataLL=data['spec'][1][0]
dataRR=data['spec'][0][0]
freq=np.round(data['freq']/1.e9,3)

plt.imshow(np.log10(dataRR),aspect='auto',origin=0,extent=[0,26400,0,506])
plt.yticks(np.arange(506)[::50],freq[::50])
plt.xticks(np.arange(26400)[::5000],np.arange(26400)[::5000]*0.05)
plt.xlabel('Time (sec) (ST: 22:04:00 UT)');plt.ylabel('Frequency (GHz)')
plt.show()

plt.plot(np.arange(26400)*0.05,dataRR[0:128].mean(axis=0),label='Freq: 0.994 -'+str(freq[128])+' GHz')
plt.xlabel('Time (sec) (ST: 22:04:00 UT)');plt.ylabel('Intensity');plt.legend()
plt.show()



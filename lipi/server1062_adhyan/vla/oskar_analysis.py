import numpy as np
import matplotlib.pyplot as plt

fil='/nas08-data02/rohit/oskar_tests/beam_study/telescope.tm/layout.txt'
layout_all=np.loadtxt(fil,delimiter=',')
layout=[0]*30
for i in range(30):
    ii="%03d" % i
    layout[i]=np.loadtxt('/nas08-data02/rohit/oskar_tests/beam_study/telescope.tm/station'+str(ii)+'/layout.txt',delimiter=',')
layout=np.array(layout);layall=layout.reshape((30*2587,2))
plt.plot(layout_all[:,0],layout_all[:,1],'o')
#plt.plot(layout[1,:,0],layout[1,:,1],'o')
#plt.plot(layall[:,0],layall[:,1],'o')
plt.xlabel('X-Distance (m)');plt.ylabel('Y-Distance (m)');plt.show()


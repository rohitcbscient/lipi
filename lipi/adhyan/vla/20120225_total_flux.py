import numpy as np
import pickle

sumall,freq=pickle.load(open('/media/rohit/VLA/20120225_sub/total_flux_20120225_max.p','rb'))
sumall_diff=sumall-sumall[:,0][:,None]
#sumall_diff[:,200:]=1.e8
sumall_diff[75:82]=sumall_diff[75:82]*sumall_diff[74]/sumall_diff[75]
sumall_diff[82:89]=sumall_diff[82:89]*sumall_diff[79]/sumall_diff[82]
sumall_diff[0:10]=sumall_diff[10:20]
sumall_diff=np.hstack((sumall_diff,np.zeros((96,358))))
sumall_diff[:,806:926]=sumall_diff[:,326:446]
sumall_diff[:,926:934]=-10
sumall_diff[:,934:1054]=sumall_diff[:,326:446]
sumall_diff[:,1054:1062]=-10
sumall_diff[:,1062:1182]=sumall_diff[:,326:446]
sumall_diff[:,1182:]=-10


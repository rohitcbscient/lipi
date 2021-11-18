import numpy as np
from scipy.io import readsav
import glob
import matplotlib.pyplot as plt


sublist=sorted(glob.glob('*94*rot.sav'))
bdiff94=[0]*len(sublist)
submap94=[0]*len(sublist)
for i in range(len(sublist)):
    subdata=readsav(sublist[i])
    submap94[i]=subdata['drot_map_'][0][0]
    submap94[i]=submap94[i][0:668,0:667]
    bdiff94[i]=submap94[i]-submap94[0]

submap94=np.array(submap94);bdiff94=np.array(bdiff94)



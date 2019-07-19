from lmfit import Parameters
import numpy as np



var_list=['A','d','T0','eps','kappa','emin','emax','ebrk','del1','del2',
          'nth','B','theta','fs','delf','ed','nf','paf','lca','bda','delmu',
          'a4','fc','fw','mat','q']
default_values=np.array([1.e18,1.e8,10e6,0.05,4.0,16,0.005,1.0,1.0,6.5,7.5,2.e10,5.e7,200,60,1.e9,0.002,3,400,3,90,0,0.2,1,-1,-1,1,2])


var=['B','nb']
varval=[300,1.e8]
varmin=[100,1.e7]
varmax=[500,1.e9]
delvar=[50,2.e7]

param=Parameters()
VARY=True
for i in range(len(var)):
    param.add_many((var[i],varval[i],VARY,varmin[i],varmax[i],None,delvar[i]))





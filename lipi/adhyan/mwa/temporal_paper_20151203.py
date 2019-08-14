import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
import pickle
from scipy import signal
from surya.utils import main as ut
from surya.utils import model as md
import glob
import matplotlib.cm as cm


def get_forward(list_):
    freq=[0]*len(list_)
    Tb=[0]*len(list_)
    flux=[0]*len(list_)
    ll=0
    for l in list_:
        data=readsav(l)
        freq[ll]=float(l.split('_')[1].split('MHz')[0])
        dx=data['quantmap'][0][4]*16*60
        Tb_=data['stokesstruct'][0][0]
        size=dx*Tb_.shape[0]
        Tb[ll]=np.mean(Tb_)
        flux[ll]=np.mean(ut.Tb2flux(Tb[ll], size, size, freq[ll]/1000.))
        ll=ll+1
    return freq,flux,Tb

dir_='/media/rohit/VLA/20151203_MWA/Tb_new/'

#flist=[108,120,132,145,160,179,196,217,240]
#flabel=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
flist=[108,132,145,160,179,196,217,240]
flabel=['084-085','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
base=['000-008','000-009','000-010','008-009','008-010','009-010']
mwa_flux=[0]*len(base)
i=0
for b in base:
    print b
    mwa_flux[i]=[0]*len(flabel)
    j=0
    for f in flabel:
        d=pickle.load(open('/home/i4ds1807205/20151203/pickle/flux_V1_1133149192-%b'+str(f)+'_T'+str(b)+'.p','rb'))
        d_=d[17][3][0][np.isfinite(d[17][3][0])]
        mwa_flux[i][j]=np.mean(d_)
        #print mwa_flux[i][j]
        j=j+1
    i=i+1

emwa_flux=np.std(np.array(mwa_flux),axis=0)
mwa_flux=np.mean(np.array(mwa_flux),axis=0)
Te=1.e6 # in K
Tb_all=[0]*len(flabel)
flux_mwa=[0]*len(flabel)
tau=[0]*len(flabel)
k=0
for f in flabel:
    print flist[k],' MHz'
    data=pickle.load(open(dir_+'Tb_1133149192-'+f+'.p','r'))
    Tb_all[k]=data[0][0:577]
    bmin[k]=data[3][0]*60
    bmax[k]=data[4][0]*60
    flux_mwa[k]=np.mean(Tb_all[k],axis=0)*mwa_flux[k]/np.sum(Tb_all[k][0])#ut.Tb2flux(Tb_all[k][0], 50, 50, flist[k]/1000.)
    #print flist[k],np.max(Tb_fwd[k]),np.max(Tb_all[k]),np.std(Tb_all[k])#,np.mean(tau[k]),np.std(tau[k])
    k=k+1

Tb_all=np.array(Tb_all)


################### PRINTS #########################
print_tau=1
if(print_tau):
    for f in flist:
        print f 

################### PLOTS ##########################

plot_sigma=1
if(plot_sigma):
    plt.plot(flist,sig)


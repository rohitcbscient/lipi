import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import scipy.io as io
import glob
from scipy.io import readsav
from sys import path as sys_path
from surya.aia.dn2dem_pos import dn2dem_pos
import astropy.time as atime
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
from astropy.io import fits
from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import register, update_pointing
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_dem(dn_in,trin,gains,dn2ph,rdnse,exptime,d1,m1,s1,root2pi,temps):
    # Get rid of the b in the string name (byte vs utf stuff....)
    for i in np.arange(len(trin['channels'])):
        trin['channels'][i]=trin['channels'][i].decode("utf-8")
        #print(trin['channels'])
        # Get the temperature response functions in the correct form for demreg
    tresp_logt=np.array(trin['logt'])
    dem_mod=(d1/(root2pi*s1))*np.exp(-(tresp_logt-m1)**2/(2*s1**2))
    nt=len(tresp_logt)
    nf=len(trin['tr'][:])
    trmatrix=np.zeros((nt,nf))
    for i in range(0,nf):
        trmatrix[:,i]=trin['tr'][i]
    # Now work out the DN/s/px
    # For AIA responses all are dlogt=0.05
    tresp_dlogt=np.full(nt,0.05)
    tc_full=np.zeros([nt,nf])
    for i in range(0,nf):
        tc_full[:,i]=dem_mod*trmatrix[:,i]*10**tresp_logt*np.log(10**tresp_dlogt)
    #dn_in=np.sum(tc_full,0)
    dn_in=np.array(dnpers)
    #print('dn_in: ',dn_in)
    # And the associated uncertainty (no systematics)
    dn0=dn_in*exptime
    shotnoise=(dn2ph*dn0)**0.5/dn2ph/exptime
    # error in DN/s/px
    edn_in=(rdnse**2+shotnoise**2)**0.5 
    #print('edn_in: ',edn_in)
    #  Setup the T binning for DEM solution
    # Temperature bin mid-points for DEM plotting
    mlogt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) for i in np.arange(0,len(temps)-1)])
    # Now work out the DEM - investigate 3 standard ways of running
    # 1. Default - reg runs twice, 1st time to work out weight for constraint matrix, then regs with that
    #         Best option if don't know what doing, hence its the default 
    dem0,edem0,elogt0,chisq0,dn_reg0=dn2dem_pos(dn_in,edn_in,trmatrix,tresp_logt,temps,max_iter=100) #gloci=0 is default behaviour
    # 2. EMloci - reg runs once, works out weight for constraint matrix as min of EM Loci, then regs with that
    #        If some of your filters have a sharper T response (lines or X-ray obs) might be useful to try
    #dem1,edem1,elogt1,chisq1,dn_reg1=dn2dem_pos(dn_in,edn_in,trmatrix,tresp_logt,temps,gloci=1)
    # 3. User weight - reg runs once, user provide weight for constraint matrix L, then regs with that
    #         If have an idea of what DEM might look like could try a rough form of it (though check vs 1, 2 above)
    # As working with synthetic data from DEM model, could try weighting by this model, interp on output DEM T bins
    demwght0=10**np.interp(mlogt,tresp_logt,np.log10(dem_mod))
    #dem2,edem2,elogt2,chisq2,dn_reg2=dn2dem_pos(dn_in,edn_in,trmatrix,tresp_logt,temps,dem_norm0=demwght0/max(demwght0))
    return dem0,edem0,elogt0,chisq0,dn_reg0,mlogt,dem_mod,tresp_logt

# For some DEM model (i.e. a Gaussian) produce the synthetic DN/s/px for each AIA channel
d1=4e22
m1=6.5
s1=0.15
root2pi=(2.*math.pi)**0.5
gains=np.array([18.3,17.6,17.7,18.3,18.3,17.6])
dn2ph=gains*np.array([94,131,171,193,211,335])/3397.
rdnse=np.array([1.14,1.18,1.15,1.20,1.20,1.18])
exptime=np.array([2.901334, 2.000201, 1.999635, 2.901215, 2.900814, 2.900825])
temps=np.logspace(5.7,7.1,num=35)
files=sorted(glob.glob('/home/rohit/DEM/*.fits'));dnpers=[0]*6
files_sav=sorted(glob.glob('/home/rohit/DEM/*.fits_lev1.5.sav'))
dd=readsav('aia_arr_02_23_lev1.5.sav')
#ln=[1808,1809];mn=[485,486]
#ln=np.arange(4096);mn=np.arange(4096)
ln=1200+np.arange(400);mn=3100+np.arange(400)
dem0=[0]*len(ln);edem0=[0]*len(ln) # ln---> Y-coord; mn---> X-coord
for l in range(len(ln)):
    dem0[l]=[0]*len(mn);edem0[l]=[0]*len(mn)
    for m in range(len(mn)):
        #for i in range(6):
        #aa=fits.open(files[i]);h=aa[0].header;d=aa[0].data
        #exptime[i]=h['EXPTIME']
        #bb=readsav(files_sav[i])['fullsunmap'][0][0]
        dnpers=dd['aia_arr'][:,ln[l],mn[m]]/exptime
        if(len(np.where(dnpers==0))!=6):
            trin=io.readsav('/home/rohit/DEM/demreg-master/python/aia_tresp_en.dat')
            dem0[l][m],edem0[l][m],elogt0,chisq0,dn_reg0,mlogt,dem_mod,tresp_logt=get_dem(dnpers,trin,gains,dn2ph,rdnse,exptime,d1,m1,s1,root2pi,temps)
        else:
            dem0[l][m],edem0[l][m]=0,0
        print(l,m,np.max(dem0[l][m]))
dem0=np.array(dem0);edem0=np.array(edem0)
pickle.dump([dem0,edem0],open('/home/rohit/DEM/20140914_02_23_aiadem_py.p','wb'))
sys.exit()

dem1=pickle.load(open('20140914_aiadem_py.p','rb'));dem2=pickle.load(open('20140914_02_23_aiadem_py.p','rb'))
f,ax=plt.subplots(3,1,sharex=True);ax0=ax[0];ax1=ax[1];ax2=ax[2]
im0=ax0.imshow(dem1[0][:,:,:].sum(axis=2),aspect='auto',origin='lower',vmin=1.e21,vmax=1.e22,cmap='YlGn')
divider = make_axes_locatable(ax0);cax = divider.append_axes("right", size="2%", pad=0.0);f.colorbar(im0,cax=cax,label='EM (cm$^{-5}$)')
im1=ax1.imshow(dem2[0][:,:,:].sum(axis=2),aspect='auto',origin='lower',vmin=1.e21,vmax=1.e22,cmap='YlGn')
divider = make_axes_locatable(ax1);cax = divider.append_axes("right", size="2%", pad=0.0);f.colorbar(im1,cax=cax,label='EM (cm$^{-5}$)')
im2=ax2.imshow(dem2[0][:,:,:].sum(axis=2)/dem1[0][:,:,:].sum(axis=2),aspect='auto',origin='lower',vmin=0.5,vmax=1.5,cmap='coolwarm')
divider = make_axes_locatable(ax2);cax = divider.append_axes("right", size="2%", pad=0.0);f.colorbar(im2,cax=cax,label='Ratio (2/1)')
plt.show()


dem_idl=readsav("/home/rohit/DEM/emcube_2014.sav")['emcube']
status_idl=readsav("/home/rohit/DEM/status_2014.sav")['status']
coeff_idl=readsav("/home/rohit/DEM/coeff_2014.sav")['coeff']

sys.exit()
matplotlib.rcParams['font.size'] = 16
clrs=['darkgreen','darkcyan','gold','sienna','indianred','darkslateblue']
# Do the plot
fig = plt.figure(figsize=(8, 7))
for i in np.arange(6):
    plt.semilogy(tresp_logt,trmatrix[:,i],label=trin['channels'][i],color=clrs[i],lw=4)
plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
plt.ylabel('$\mathrm{AIA\;Response\;[DN\;s^{-1}\;px^{-1}\;cm^5]}$')
plt.ylim([2e-29,5e-24])
plt.xlim([5.0,8.0])
plt.legend(ncol=2,prop={'size': 16})
plt.rcParams.update({'font.size': 16})
plt.grid(True,which='both',lw=0.5,color='gainsboro')
plt.show()


# # Check what the DEM model looks like
# fig = plt.figure(figsize=(8, 4.5))
# plt.plot(tresp_logt,dem_mod)
# plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
# plt.ylabel('$\mathrm{DEM\;[cm^{-5}\;K^{-1}]}$')
# plt.ylim([2e20,4e23])
# plt.xlim([5.7,7.2])
# plt.rcParams.update({'font.size': 16})
# plt.yscale('log')
# plt.show()




yr=[2e19,4e23]
xr=[5.7,7.2]#np.log10([min(temps),max(temps)])
fig = plt.figure(figsize=(8, 4.5))
plt.errorbar(mlogt,dem0[5,5],xerr=elogt0,yerr=edem0[5,5],fmt='or',\
    ecolor='lightcoral', elinewidth=3, capsize=0,label='Def Self LWght')
plt.plot(tresp_logt,dem_mod,'--',color='grey')
plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
plt.ylabel('$\mathrm{DEM\;[cm^{-5}\;K^{-1}]}$')
plt.ylim(yr)
plt.xlim(xr)
plt.rcParams.update({'font.size': 16})
plt.yscale('log')
plt.legend()
# plt.savefig('demregpy_aiasyn_slw.png',bbox_inches='tight')
plt.show()


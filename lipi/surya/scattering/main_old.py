import numpy as np
import math
from scipy.io import readsav

def lnlike(theta, x, y, yerr):
    B,nb,nth,emin,d,lnf = theta
    model=residual_model(theta)
    inv_sigma2 = 1.0/(yerr**2 + model**2)
    print(np.sum((y-model)**2*inv_sigma2),np.sum(yerr),np.sum(model))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    B, nb, nth, emin, d, lnf = theta
    if 100 < B < 300 and 1.e7 < nb < 1.e9 and 1.e9 <nth<1.e11 and 0.001 <emin< 0.01 and 1.e7 <d< 1.e9 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

def lnprob(theta,x,y,yerr):
    lp=lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp+lnlike(theta,x,y,yerr)



#init_val=np.array([B_true,nb_true,nth_true,emin_true,d_true,f_true])
#nll = lambda *args: -lnlike(*args)
#bnds = ((0, None), (0, None),(0,None),(0,None),(0,None),(None,None))
#result = op.minimize(nll,[B_true, nb_true, nth_true, emin_true, d_true, np.log(f_true)], args=(freq_data, Tb_data, Tb_noise),method='SLSQP',bounds=bnds)

#ndim, nwalkers = 6, 50
#pos = [init_val + (init_val/1)*np.random.randn(ndim) for i in range(nwalkers)]

#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(freq, Tb_data, Tb_noise))
#print('#### OPTIMIZING...')
#sampler.run_mcmc(pos, 500)
#samples = sampler.chain[:, :, :].reshape((-1, ndim))
#pickle.dump(sampler.chain,open('mcmc_20120225_a.p','wb'))
#fig = corner.corner(samples, labels=["B", "nb", "nth","emin","depth","$\ln,f$"],
#fig.savefig("mcmc_20120225_a.png")




data=readsav('/home/i4ds1807205/20151203/20151203_240MHz_psimas.sav')

#POS: Plane of Sky
thpos=aa['gridpramsstruct'][0][3]
phipos=aa['gridpramsstruct'][0][4]
phpos=aa['gridpramsstruct'][0][5]

mdtor=np.pi/180.
dx=aa['gridpramsstruct'][0][8]
dy=aa['gridpramsstruct'][0][9]
xmin,xmax=aa['gridpramsstruct'][0][10]
ymin,ymax=aa['gridpramsstruct'][0][11]
ngrid=256
x=np.linspace(xmin,xmax,ngrid)
y=np.linspace(ymin,ymax,ngrid)
r=np.sqrt(x*x+y*y)
th=acos(y/x)
cmer=aa['gridpramsstruct'][0][6]
phio=aa['gridpramsstruct'][0][7]
ph=(cmer-phio)*mdtor

#MODEL

pres=aa['modsolstruct'][0][0]
dens=aa['modsolstruct'][0][1]
temp=aa['modsolstruct'][0][2]
br=aa['modsolstruct'][0][3]
bth=aa['modsolstruct'][0][4]
bph=aa['modsolstruct'][0][5]
vr=aa['modsolstruct'][0][6]
vth=aa['modsolstruct'][0][7]
vph=aa['modsolstruct'][0][8]

#LOS
losmin=aa['lospramsstruct'][1]
losint=aa['lospramsstruct'][2]
nlos=aa['lospramsstruct'][3]


btot=np.sqrt(br*br+bth*bth+bph*bph)
blos=sin()


np.ones(nlos)




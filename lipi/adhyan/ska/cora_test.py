import numpy as np
from cora.signal import corr21cm
from cora.foreground import galaxy
import healpy
import matplotlib.pyplot as plt
import healpy as hp
from mpl_toolkits.axes_grid1 import make_axes_locatable

cr = galaxy.FullSkySynchrotron()
aps1 = cr.angular_powerspectrum(np.arange(1000), 800.0, 800.0)
assert len(aps1) == 1000
assert np.allclose(aps1.sum(), 75.47681191093129, rtol=1e-7)  # Calculated for commit 02f4d1cd3f402d
fa = np.linspace(400.0, 800.0, 64)
aps2 = cr.angular_powerspectrum(np.arange(1000)[:, None, None], fa[None, :, None], fa[None, None, :])
assert aps2.shape == (1000, 64, 64)

gpolsky=cr.getpolsky()
gsky=cr.getsky()
healpy.visufunc.mollview(gsky[0])
plt.show()

ipix=np.arange(100)
theta, phi = hp.pix2ang(nside, ipix)
ra = np.rad2deg(phi);dec = np.rad2deg(0.5 * np.pi - theta)
pix=hp.ang2pix(nside,theta,phi)


##### Select a circle

ra = 45
dec = 45
radius = 10.1 # in degrees
theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)
radius = np.deg2rad(radius)
xyz = hp.ang2vec(theta, phi) # Cartisian Coordinates
ipix_disc = hp.query_disc(nside, xyz, radius)

emptyhealpix=np.zeros(gsky.shape[1])
emptyhealpix[ipix_disc]=1.0
healpy.visufunc.mollview(emptyhealpix)
plt.show()

##################

from cora.util import hputil
alms=hputil.sphtrans_real(gsky[0],lside=12,lmax=50)
cl = hp.anafast(gsky[0], lmax=50)
ell = np.arange(len(cl));gsky_inv=[0]*50
for l in range(50):
    print(l);gsky_inv[l]=[0]*50
    for m in range(50):
        alms1=alms*0;alms1=alms1+alms[l,m]
        alm_packed=hputil.pack_alm(alms1,lmax=50)
        gsky_inv[l][m] = hp.alm2map(alm_packed,nside=128)
gsky_inv=np.array(gsky_inv)
alm_packed_direct=hputil.pack_alm(alms,lmax=50)
gsky_inv_direct=hp.alm2map(alm_packed_direct,nside=128)
ylm = healpy.alm2map(alm_packed_direct, nside=128, pol=False, sigma=0.0, lmax=50, mmax=50)

alms_sq = alms.real**2 + alms.imag**2

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
plt.axes(ax1)
hp.visufunc.mollview(gsky_inv.mean(axis=(0,1)),hold=True,title='a$_{lm}$(l=0,m=0)')
plt.axes(ax2)
hp.visufunc.mollview(gsky_inv_direct,hold=True,title='Full Synchroton Sky [CORA]')
im3=ax3.imshow(alms_sq,origin='lower',vmin=0,vmax=0.1);ax3.set_ylabel('m');ax3.set_xlabel('l');ax3.set_title('a$_{lm}$')
ax4.plot(ell,cl,'o-');ax4.set_xlabel('l');ax4.set_ylabel('$C_{l}$');ax4.set_title('Power Spectrum')
divider = make_axes_locatable(ax3);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im3, cax=cax, orientation='vertical')
plt.show()




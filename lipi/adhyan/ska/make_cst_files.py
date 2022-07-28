# Script to write the cst files
import numpy as np
from katbeam import JimBeam
import matplotlib.pylab as plt
from astropy.io import fits


def get_meerkat_uhfbeam(f,pol,beamextent):
    beam=JimBeam('MKAT-AA-UHF-JIM-2020');freqlist=uhfbeam.freqMHzlist
    margin=np.linspace(-beamextent/2.,beamextent/2.,int(beamextent*2))
    x,y=np.meshgrid(margin,margin)
    freqMHz_idx=np.where(freqlist==freqlist.flat[np.abs(freqlist - f).argmin()])[0][0]
    freqMHz=freqlist[freqMHz_idx]
    if pol=='H':
        beampixels=beam.HH(x,y,freqMHz)
    elif pol=='V':
        beampixels=beam.VV(x,y,freqMHz)
    else:
        beampixels=beam.I(x,y,freqMHz)
        pol='I'
    return beampixels

def show_beam(beampixels,beamextent,freq,pol):
    plt.imshow(beampixels,extent=[-beamextent/2,beamextent/2,-beamextent/2,beamextent/2])
    plt.title('%s pol beam\nfor %s at %dMHz'%(pol,'',freq))
    plt.xlabel('deg');plt.ylabel('deg');plt.colorbar()
    plt.show()

import eidos
from eidos.create_beam import zernike_parameters,save_fits
from eidos.parallelize import parmap
from eidos.spatial import recon_par,jones_to_mueller_all
meerkat_beam_coeff_ah='/home/rohit/eidos/eidos/data/meerkat_beam_coeffs_ah_zp_dct.npy'
meerkat_beam_coeff_em='/home/rohit/eidos/eidos/data/meerkat_beam_coeffs_em_zp_dct.npy'
npix=500;dia=10;thres=0
params_em, freqs = zernike_parameters(meerkat_beam_coeff_em, npix, dia, thres) # Params 1st index is channel; 2nd and 3rd indices are zernike coeefient (shape is (2,2,2,20)); 4th index is Npix and 5th threshold
params_ah, freqs = zernike_parameters(meerkat_beam_coeff_ah, npix, dia, thres) # Params 1st index is channel; 2nd and 3rd indices are zernike coeefient (shape is (2,2,2,20)); 4th index is Npix and 5th threshold
ch=0
B_em = recon_par(params_em[ch,:])
B_ah = recon_par(params_ah[ch,:])
thphi_arr_eidos = np.meshgrid(np.linspace(-5,5,npix),np.linspace(-5,5,npix))
B_em[:,:,0:2,:]=0+0j;B_em[:,:,-2:,:]=0+0j;B_em[:,:,:,-2:]=0+0j;B_em[:,:,:,0:2]=0+0j
B_ah[:,:,0:2,:]=0+0j;B_ah[:,:,-2:,:]=0+0j;B_ah[:,:,:,-2:]=0+0j;B_ah[:,:,:,0:2]=0+0j


###################################
f,ax=plt.subplots(2,2);ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
ax00.imshow(10*np.log10(np.abs(B_ah[0,0])),aspect='auto',origin='lower',extent=[-5,5,-5,5]);ax00.set_title('E$_{00}^{h}$')
ax01.imshow(10*np.log10(np.abs(B_ah[0,1])),aspect='auto',origin='lower',extent=[-5,5,-5,5]);ax01.set_title('E$_{01}^{h}$')
ax10.imshow(10*np.log10(np.abs(B_ah[1,0])),aspect='auto',origin='lower',extent=[-5,5,-5,5]);ax10.set_title('E$_{10}^{h}$')
im=ax11.imshow(10*np.log10(np.abs(B_ah[1,1])),aspect='auto',origin='lower',extent=[-5,5,-5,5]);ax11.set_title('E$_{11}^{h}$')
ax10.set_xlabel('Deg');ax00.set_ylabel('Deg')
ax11.set_xlabel('Deg');ax10.set_ylabel('Deg')
plt.colorbar(im)
plt.show()

f,ax = plt.subplots(2,1);ax0=ax[0];ax1=ax[1]
ax0.plot(np.linspace(-5,5,npix),10*np.log10(np.abs(B_ah[0,0]))[250],'o-',label='AH')
ax0.plot(np.linspace(-5,5,npix),10*np.log10(np.abs(B_em[0,0]))[250],'o-',label='EM')
ax1.plot(np.linspace(-5,5,npix),10*np.log10(np.abs(B_em[0,0]))[250]-10*np.log10(np.abs(B_ah[0,0]))[250],'o-',label='Residual')
ax1.set_xlabel('Distance from center (deg)');ax0.set_ylabel('Power (dB)');ax0.legend()
plt.show()


B=B_ah
if len(B.shape)==4: B = np.expand_dims(B, axis=0)
data_M = jones_to_mueller_all(B)
data = np.zeros((B.shape[0],1,1,B.shape[3],B.shape[4]), dtype=np.complex)
eidos.create_beam.main(['-p', '500','-d','10','-f', '1.e9'])
re_=fits.open('primary_beam_mh_1000000000MHz_10deg_re.fits');re_beam=re_[0].data[0]
im_=fits.open('primary_beam_mh_1000000000MHz_10deg_im.fits');im_beam=im_[0].data[0]
abs_beam=np.sqrt(re_beam*re_beam+im_beam*im_beam)

bb=np.load('/home/rohit/eidos/eidos/data/meerkat_beam_coeffs_ah_zp_dct.npy',encoding='latin1',allow_pickle=True).item()
#dict_keys(['nu', 'dctc', 'dcti', 'zi'])

######################################################
#f,ax=plt.subplots(4,4)
#ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
#ax10=ax[0,0];ax11=ax[0,1];ax12=ax[1,2];ax13=ax[1,3]
#plt.imshow(abs_beam,aspect='auto',origin='lower')
#plt.show()


line1='Theta [deg.]  Phi   [deg.]  Abs(Dir.)[dBi   ]   Abs(Theta)[dBi   ]  Phase(Theta)[deg.]  Abs(Phi  )[dBi   ]  Phase(Phi  )[deg.]  Ax.Ratio[dB    ]    '
line2='------------------------------------------------------------------------------------------------------------------------------------------------------'
line3='Theta [deg.]  Phi   [deg.]  Abs(Dir.)[dBi   ]   Abs(Horiz)[dBi   ]  Phase(Horiz)[deg.]  Abs(Vert  )[dBi   ]  Phase(Vert  )[deg.]  Ax.Ratio[dB    ]  '
line3='Theta [deg.]  Phi   [deg.]  Abs(Dir.)[dBi   ]   Horiz(Abs)[dBi   ]  Horiz(Phase)[deg.]  Vert(Abs)[dBi   ]  Vert(Phase )[deg. ]  Ax.Ratio[dB    ]  '
aa=np.loadtxt('run5.cst',skiprows=2)
np.savetxt('run6.cst',aa, delimiter=" ", header=line1+"\n"+line2,comments='')
theta=aa[:,0].reshape(360,181);phi=aa[:,1].reshape(360,181);absdir=aa[:,2].reshape(360,181)
abstheta=aa[:,3].reshape(360,181);phasetheta=aa[:,4].reshape(360,181);absphi=aa[:,5].reshape(360,181);phasephi=aa[:,6].reshape(360,181)
ax_ratio=aa[:,7].reshape(360,181)

theta_pb=57.5 # arcmin
beam=np.zeros((360*181,8))
thphi_arr=np.meshgrid(np.linspace(-90,90,181),np.linspace(-180,179,360))

th_abs=[0]*360;phi_abs=[0]*360;th_ph=[0]*360;phi_ph=[0]*360;ax_ratio_meerkat=[0]*360
for i in range(360):
    th_abs[i]=[0]*181;phi_abs[i]=[0]*181;th_ph[i]=[0]*181;phi_ph[i]=[0]*181;ax_ratio_meerkat[i]=[0]*181
    for j in range(181):
        t0=np.abs(thphi_arr_eidos[0][0] - thphi_arr[0][i][j]).argmin();t1=np.abs(thphi_arr_eidos[1][:,0] - thphi_arr[1][i][j]).argmin()
        th_abs[i][j] = np.abs(B_ah[0][0][t0][t1])
        th_ph[i][j]= np.rad2deg(np.angle(B_ah[0][0][t0][t1]))
        phi_abs[i][j]=np.abs(B_ah[0][1][t0][t1])
        phi_ph[i][j]= np.rad2deg(np.angle(B_ah[0][1][t0][t1]))
        ax_ratio_meerkat[i][j]= 0.0
th_abs=np.array(th_abs)
phi_abs=np.array(phi_abs)
th_ph=np.array(th_ph)
phi_ph=np.array(phi_ph)
ax_ratio_meerkat=np.array(ax_ratio_meerkat)
beam=np.zeros((360*181,8))
beam[:,1]=thphi_arr[1].flatten()+180;beam[:,0]=thphi_arr[0].flatten()+90;beam[:,3]=th_abs.flatten();beam[:,4]=th_ph.flatten();beam[:,5]=phi_abs.flatten();beam[:,6]=phi_ph.flatten()
beam[:,7]=ax_ratio_meerkat.flatten()
beam[np.where(beam==0)]=1.e-7
np.savetxt('run9.cst',beam, delimiter=" ", header=line1+"\n"+line2,comments='')

write_eidos_beam=1
if(write_eidos_beam):
    B=B_ah;cc=1;B[np.where(B==0)]=1.e-10
    beam=np.zeros((250000,8))
    r_arr_eidos=np.sqrt(thphi_arr_eidos[1].flatten()**2 + thphi_arr_eidos[0].flatten()**2)
    ang_arr_eidos=np.rad2deg(np.arctan2(thphi_arr_eidos[1].flatten(),thphi_arr_eidos[0].flatten()))
    beam[:,1]=ang_arr_eidos;beam[:,0]=r_arr_eidos;beam[:,3]=np.abs(B[0][0]).flatten();beam[:,4]=np.rad2deg(np.angle(B[0][0])).flatten();beam[:,5]=np.abs(B[0][1]).flatten();beam[:,6]=np.rad2deg(np.angle(B[0][1])).flatten()
    beam[:,7]=np.zeros(len(thphi_arr_eidos[1].flatten()))
    #beam[np.where(beam==0)]=1.e-8
    np.savetxt('run8.cst',beam, delimiter=" ", header=line3+"\n"+line2,comments='')

sys.exit()
beam=np.zeros((360*181,8))
beam[:,0]=thphi_arr[1].flatten();beam[:,1]=thphi_arr[0].flatten()
#ab=(np.cos(1.189*np.pi*((thphi_arr[0].flatten()-90)/(theta_pb/60.)))/(1-4*(1.189*(thphi_arr[0].flatten()-90)/(theta_pb/60.))**2))**2;beam[:,2]=ab
ab_gauss=np.exp(-0.5*((thphi_arr[0].flatten()-90)**2)/((theta_pb/60.)**2))*(1/(np.sqrt(2*np.pi)*(theta_pb/60.)));ab_gauss=ab_gauss/np.max(ab_gauss);beam[:,2]=ab_gauss
np.savetxt('run7.cst',beam, delimiter=" ", header=line1+"\n"+line2,comments='')

H_beampixels=get_meerkat_uhfbeam(720,'H',180)
V_beampixels=get_meerkat_uhfbeam(720,'H',180)
show_beam(H_beampixels,180,720,'H')


ab=(np.cos(1.189*np.pi*(rho/(theta_pb/60.)))/(1-4*(1.189*rho/(theta_pb/60.))**2))**2
plt.plot(rho,ab,'o-')
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.pcolormesh(phi[:,0]*np.pi/180, theta[0][0:90]*np.pi/180, absdir.swapaxes(0,1)[0:90,:]) #X,Y & data2D must all be same dimensions
plt.show()

f=plt.figure()
ax = f.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.plot(phi[:,0]*np.pi/180,absdir[:,0]*np.pi/180)
plt.show()

ax.plot(aa[:,3],projection='polar')



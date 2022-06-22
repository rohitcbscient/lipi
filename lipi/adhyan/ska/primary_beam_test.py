import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable


aa=fits.open('/home/rohit/simulations/primary_beam/beam_test_S0000_TIME_SEP_CHAN_SEP_AMP_XX.fits')
data=aa[0].data[0][0]
print(np.nanmax(data),np.where(data==np.nanmax(data)))
plt.imshow(data);plt.show()


pb=np.loadtxt("/home/rohit/simulations/primary_beam/azel0.beam")
pb[:,2]=np.arange(5,10,len(pb[:,2]));pb[:,3]=np.linspace(0,180,len(pb[:,3]))

gr=np.meshgrid(np.arange(0,90,1),np.arange(0,360,1))
pb1=np.vstack((gr[0].flatten(),gr[1].flatten(),np.linspace(5,10,32400),np.linspace(0,180,32400),np.linspace(0,180,32400),np.linspace(0,180,32400),np.linspace(0,180,32400),np.linspace(0,180,32400)))
pb2=np.vstack((gr[0].flatten(),gr[1].flatten()))
np.savetxt('/home/rohit/simulations/primary_beam/cst.beam',pb2,delimiter=' ')
np.savetxt('/home/rohit/simulations/primary_beam/scalar.beam',pb1,delimiter=' ')

mf=fits.open('/home/rohit/simulations/primary_beam/rascil_feko/MID_FEKO_VP_B1_45_0765_real.fits')
wcs = WCS(mf[0].header)
mid_feko_real=mf[0].data[0][0]
mf=fits.open('/home/rohit/simulations/primary_beam/rascil_feko/MID_FEKO_VP_B1_45_0765_imag.fits')
mid_feko_imag=mf[0].data[0][0]
f,ax=plt.subplots(1,3);ax0=ax[0];ax1=ax[1];ax2=ax[2]
im0=ax0.imshow(mid_feko_real,origin='lower');ax0.set_title('Real')
im1=ax1.imshow(mid_feko_imag,origin='lower');ax1.set_title('Imaginary')
im2=ax2.imshow(np.sqrt(mid_feko_imag**2 +mid_feko_real**2),origin='lower');ax2.set_title('Amplitude')
divider = make_axes_locatable(ax0);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im0, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax1);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im1, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax2);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im2, cax=cax, orientation='vertical')
plt.show()

plt.subplot(projection=wcs,slices=('y', 'x',0,0))
plt.imshow(np.sqrt(mid_feko_imag**2 +mid_feko_real**2),origin='lower')
plt.colorbar()
plt.show()

img=fits.open('/home/rohit/simulations/primary_beam/test_image_I.fits')
image=img[0].data[0]
print(np.nanmax(image),np.where(image==np.nanmax(image)),np.nanstd(image))
plt.imshow(image);plt.show()

imageall=[0]*5
for i in range(5):
    img=fits.open(f"/home/rohit/simulations/primary_beam/run{i+1}/run{i+1}_I.fits")
    imageall[i]=img[0].data[0]

fig,ax=plt.subplots(2,3)
ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2]
im00=ax00.imshow(imageall[0],origin='lower');ax00.set_title('Run 1')
im01=ax01.imshow(imageall[1],origin='lower');ax01.set_title('Run 2')
im10=ax10.imshow(imageall[2],origin='lower');ax10.set_title('Run 3')
im11=ax11.imshow(imageall[3],origin='lower');ax11.set_title('Run 4')
im12=ax12.imshow(imageall[4],origin='lower');ax12.set_title('Run 5')
divider = make_axes_locatable(ax00);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im00, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax01);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im01, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax10);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im10, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax11);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im11, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax12);cax = divider.append_axes('right', size='5%', pad=0.05);fig.colorbar(im12, cax=cax, orientation='vertical')
plt.show()

img_spl=fits.open('/home/rohit/simulations/primary_beam/test_image_I_with_spline.fits')
image_spl=img_spl[0].data[0]
print(np.nanmax(image_spl),np.where(image_spl==np.nanmax(image_spl)),np.nanstd(image_spl))
plt.imshow(image_spl);plt.show()

with open('telescope.tm/element_pattern_fit_x_0_1000.bin', mode='rb') as file:
    fileContent = file.read();val=[0]*len(fileContent)
    for i in range(len(fileContent)):
        val[i]=fileContent[i]



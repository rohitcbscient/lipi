import pickle 
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import numpy as np
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import pickle

aa=pickle.load(open('/media/rohit/MWA/20151203_sub_amp/Tb/Tb_20151203_187-188_sub.p'))
a=fits.open('/media/rohit/MWA/20151203_EUV/aia.lev1.193A_2015-12-03T03_10_17.84Z.image_lev1.fits');ad=a[0].data
d=aa[0]
d[d<6000]=0

nd=[0]*d.shape[0]
k=0
for i in range(d.shape[0]):
    if(d[i].mean()!=0):
        nd[k]=d[i]
        k=k+1
nd=np.array(nd[0:209])

for i in range(209):
    ii='%02d' %i
    lev=np.array([0.1,0.3,0.5,0.7,0.9,0.99])
    f,ax=plt.subplots(1,1,figsize=(15,15))
    ax.imshow(ad,origin=0,extent=[-1228,1228,-1228,1228],cmap='gist_stern',vmin=0,vmax=1200)
    ax.contourf(nd[i]/np.nanmax(nd[i]),origin='lower',extent=[-5000,5000,-5000,5000],levels=lev,alpha=0.8)
    ax.set_xlim(-1228,1228);ax.set_ylim(-1228,1228)
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelbottom=False)
    arr_lena = mpimg.imread('/home/i4ds1807205/Pictures/02.jpg')
    #arr_lena1=arr_lena[17:-17,25:-25,:]
    imagebox = OffsetImage(arr_lena,zoom=0.42)
    ab = AnnotationBbox(imagebox, (-900, -1000))
    ax.add_artist(ab)
    arr_lena1 = mpimg.imread('/home/i4ds1807205/Pictures/03.jpg')
    arr_lena1=arr_lena1[:,::-1,:]
    imagebox1 = OffsetImage(arr_lena1,zoom=0.42)
    ab1 = AnnotationBbox(imagebox1, (950, -1000))
    ax.add_artist(ab1)
    ax.text(-1150, 1000, '$\mathcal{HAPPY}$', style='normal',size=50, color='red',weight='bold',
                    bbox={'facecolor': 'yellow', 'alpha': 0.8, 'pad': 10})
    ax.text(-1150, 800, '$\mathcal{DIWALI}$', style='oblique',size=50, color='red',
                    bbox={'facecolor': 'yellow', 'alpha': 0.8, 'pad': 10})
    ax.text(900,1100,'AIA 193 $\AA$',color='white',size=15,style='italic')
    ax.text(900,1050,'MWA 240 MHz',color='white',size=15,style='italic')
    plt.axis('off')
    plt.savefig('/media/rohit/MWA/20151203_sub_amp/for_diwali/'+str(ii)+'.png',dpi=100)
    plt.close()

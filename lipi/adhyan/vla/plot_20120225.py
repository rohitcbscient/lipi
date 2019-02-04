import matplotlib.pyplot as plt
import numpy as np
import setting_20120225

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')


def euv_vla(cmap,vlasubmap,xl,xr,yl,yr):
    lev=np.linspace(0.35,0.95,10)*7
    plt.imshow(cmap[40],origin=True,extent=[xl,xr,yl,yr],cmap='hot')
    plt.contour(vlasubmap[::-1,:],extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='white')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    #plt.title(t)
    plt.grid(True)
    plt.show()


if __name__=='__main__':
    main();
else:
    print 'plotting module for 20120225....'



from astropy.io import fits
import glob
from astropy.io.fits import getdata
from sunpy.map import Map
import pickle

listlasco=sorted(glob.glob('/media/rohit/MWA/20140914/lasco/*.fts'))

d=[0]*len(listlasco);smap=[0]*len(listlasco);time=[0]*len(listlasco)
for i in range(len(listlasco)):
    aa=fits.open(listlasco[i])
    d[i], h = getdata(listlasco[i], 0, header=True)
    time[i]=h['TIME-OBS']
    dd=d[i]-d[0]
    smap[i]=Map(dd,h)
pickle.dump([smap,time],open('/media/rohit/MWA/20140914/lasco/lasco_maps.p','wb'))




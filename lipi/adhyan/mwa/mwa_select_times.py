from astropy.table import Table
from astropy.time import Time
import numpy as np
tab = Table.read("/media/rohit/MWA/MWA-data-selection/phase1.fits")
times = Time(tab['obsid'], format='gps').utc
year = '2014'
month = '10'

tt=list(times.isot)
obsid=[]
time_obsid = []
for t in tt:
    if((t.split('-')[0]==year) and (t.split('-')[1]==month) ):
        idx = tt.index(t) 
        obsid.append(tab['obsid'][idx])
        time_obsid.append(times[idx].isot)
np.savetxt(year+"_"+month+'.txt',np.array(obsid).T,fmt='%d')

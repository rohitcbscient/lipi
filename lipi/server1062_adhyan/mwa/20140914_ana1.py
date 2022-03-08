import numpy as np
from matplotlib import pyplot as plt
import glob
import sunpy.map
from sunpy.instr.aia import aiaprep
from sunpy.net import Fido, attrs as a
import pickle
from astropy.coordinates import SkyCoord
from astropy import units as u

import warnings

spwlist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
idlist=["1094693472","1094693776","1094694072","1094694376","1094694672","1094694968","1094695272","1094695576","1094695872","1094696168","1094696472","1094696768","1094697072","1094697376",
        "1094697672","1094697968","1094698272","1094698568","1094698872","1094699176","1094699472","1094699768","1094700072","1094700368","1094700672","1094700976","1094701272","1094701568","1094701872","1094702168","1094702472","1094702768","1094703072", 
        "1094703368","1094703672","1094703976","1094704272","1094704568","1094704872","1094705168","1094705472","1094705776","1094706072","1094706368","1094706672","1094706968","1094707272","1094707576","1094707872","1094708168","1094708472","1094708768","1094709072"]

result = Fido.search(a.Time('2014-09-14T02:45:00', '2014-09-14T06:01:00'),a.Instrument("hmi"), a.Wavelength(193*u.angstrom),a.vso.Sample(12*u.second))
result = Fido.search(a.Time('2014/09/14 02:45:00', '2014/09/14 06:01:00'),a.Instrument.hmi, a.Physobs.los_magnetic_field)
file_download = Fido.fetch(result, site='ROB',path='.')


list171=sorted(glob.glob('aia*171*lev1.fits'));n=len(list171)
aia=[0]*n;maxaia=[0]*n;aiab=[0]*n;aiar=[0]*n
aia[0]=aiaprep(sunpy.map.Map(list171[0]))
for i in range(20,n):
    out='map_'+list171[i].split('1.')[0]+'1.5.p'
    outbd='bmap_'+list171[i].split('1.')[0]+'1.5.p'
    outrd='rmap_'+list171[i].split('1.')[0]+'1.5.p'
    aia1 = sunpy.map.Map(list171[i])
    aia[i] = aiaprep(aia1);maxaia[i]=aia[i].data.max();aiar1=aiaprep(sunpy.map.Map(list171[i-20]))
    aiab[i]=sunpy.map.Map(aia[i].data-aia[0].data,aia[i].meta);aiar[i]=sunpy.map.Map(aia[i].data-aiar1.data,aia[i].meta)
    pickle.dump([aia[i],aiab[i],aiar[i]],open('/nas08-data02/rohit/20140914/sunpy/'+out,'wb'))

from matplotlib.animation import FuncAnimation
fig=plt.figure()


def animate(i):
    aiab[i].plot()
    return line
    


import numpy as np
import os

os.system('rm -rf Gauss_point.cl')
os.system('rm -rf complist_only')
cl.done()
cl.addcomponent(dir="J2000 10h12m00.00s -27d00m00.0s", shape='point',flux=10000.0, fluxunit='Jy', freq='0.24GHz')
cl.addcomponent(dir="J2000 10h10m00.00s -28d00m00.0s", shape='point',flux=10000.0, fluxunit='Jy', freq='0.24GHz')
cl.addcomponent(dir="J2000 10h05m00.00s -29d00m00.0s", shape='point',flux=10000.0, fluxunit='Jy', freq='0.24GHz')
#cl.addcomponent(dir="J2000 10h03m00.00s -30d00m00.0s", shape='point',flux=10000.0, fluxunit='Jy', freq='0.24GHz')
cl.addcomponent(dir="J2000 10h00m00.00s -30d00m00.0s", flux=20000.0, fluxunit='Jy', freq='0.240GHz', shape="disk", majoraxis="50arcmin", minoraxis='50arcmin', positionangle='45.0deg')
#cl.setshape(1,'disk','16.0arcmin')
#
#cl.addcomponent(dir="J2000 10h00m00.08s -30d00m02.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
#cl.addcomponent(dir="J2000 09h59m59.92s -29d59m58.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
#cl.addcomponent(dir="J2000 10h00m00.40s -29d59m55.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
#cl.addcomponent(dir="J2000 09h59m59.60s -30d00m05.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
cl.rename('Gauss_point.cl')
cl.done()

default("simobserve")
project = "complist_only"
complist = 'Gauss_point.cl'
compwidth = '0.002GHz'
direction = "J2000 10h00m00.0s -30d00m00.0s" 
obsmode = "int"
antennalist = 'mwa_phase-I.cfg'
totaltime = "0.5s"
mapsize = "400arcmin"
simobserve()

default("simanalyze")
project = "complist_only"
vis="complist_only.mwa_phase-I.ms"
#imsize = [2048]
imdirection = "J2000 10h00m00.0s -30d00m00.0s" 
cell = '50arcsec'
niter = 1
threshold = ''
analyze = True
simanalyze()

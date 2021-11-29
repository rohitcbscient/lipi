import numpy as np
import os
import pyproj
import sys


wgs84 = pyproj.Proj(projparams = 'epsg:4326')
InputGrid = pyproj.Proj(projparams = 'epsg:3857')
x1, y1 = -11705274.6374,4826473.6922
pyproj.transform(InputGrid, wgs84, x1, y1)

isn2004=pyproj.CRS("+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=1700000 +y_0=300000 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1")   
# Define some common projections using EPSG codes 
wgs84=pyproj.CRS("EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
pyproj.transform(wgs84, isn2004, 63.983, -19.700)


sys.exit()

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

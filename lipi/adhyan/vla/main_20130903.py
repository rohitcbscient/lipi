from sunpy import coordinates as co
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame as bf
from sunpy.coordinates import frames
import matplotlib.pyplot as plt
#sc = SkyCoord(349.4*u.deg, -5.64*u.deg, 1.54*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')rdinates import frames

#co.transformations.hgs_to_hcc()

# hgs: Heliographic Stonyhurst
# hgc: Heliographic Carrington
# hcc: Heliocentric Cartisean
# hpc: Helioprojective Cartisean
# hcrs: Heliographic Stonyhurst with z-axis alligned with the rotation axis

# hgc_to_hgs

# hgs_to_hcc

# hcc_to_hpc
#sc = SkyCoord(349.4*u.deg, -5.64*u.deg, 1.54*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
#obstime="2013/09/03T03:50:00"

#hpsc=sc.transform_to(frames.Helioprojective)
#co.transformations.hgc_to_hgs(sc,fram)

B_lcl=np.genfromtxt('/home/i4ds1807205/Dropbox/20130903/large_closed_field_line_data.txt', skip_header=2)
dist_lcl=B_lcl[:,1]
lat_lcl=B_lcl[:,4]
lon_lcl=B_lcl[:,5]
sc_lcl = SkyCoord(lon_lcl*u.deg, lat_lcl*u.deg, dist_lcl*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
hpsc_lcl=sc_lcl.transform_to(frames.Helioprojective)
x_hpsc_lcl=hpsc_lcl.Tx
y_hpsc_lcl=hpsc_lcl.Ty


B_mcl=np.genfromtxt('/home/i4ds1807205/Dropbox/20130903/mid-size_closed_field_line_data.txt', skip_header=2)
dist_mcl=B_mcl[:,1]
lat_mcl=B_mcl[:,4]
lon_mcl=B_mcl[:,5]
sc_mcl = SkyCoord(lon_mcl*u.deg, lat_mcl*u.deg, dist_mcl*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
hpsc_mcl=sc_mcl.transform_to(frames.Helioprojective)
x_hpsc_mcl=hpsc_mcl.Tx
y_hpsc_mcl=hpsc_mcl.Ty

B_open=np.genfromtxt('/home/i4ds1807205/Dropbox/20130903/open_field_line_data.txt', skip_header=2)
dist_open=B_open[:,1]
lat_open=B_open[:,4]
lon_open=B_open[:,5]
sc_open = SkyCoord(lon_open*u.deg, lat_open*u.deg, dist_open*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
hpsc_open=sc_open.transform_to(frames.Helioprojective)
x_hpsc_open=hpsc_open.Tx
y_hpsc_open=hpsc_open.Ty

B_cl1=np.genfromtxt('/home/i4ds1807205/Dropbox/20130903/closed-line-1_data.txt', skip_header=2)
dist_cl1=B_cl1[:,1]
lat_cl1=B_cl1[:,4]
lon_cl1=B_cl1[:,5]
sc_cl1 = SkyCoord(lon_cl1*u.deg, lat_cl1*u.deg, dist_cl1*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
hpsc_cl1=sc_cl1.transform_to(frames.Helioprojective)
x_hpsc_cl1=hpsc_cl1.Tx
y_hpsc_cl1=hpsc_cl1.Ty

B_cl2=np.genfromtxt('/home/i4ds1807205/Dropbox/20130903/closed-line-2_data.txt', skip_header=2)
dist_cl2=B_cl2[:,1]
lat_cl2=B_cl2[:,4]
lon_cl2=B_cl2[:,5]
sc_cl2 = SkyCoord(lon_cl2*u.deg, lat_cl2*u.deg, dist_cl2*695508.0*u.km, obstime="2013/09/03T03:50:00",frame='heliographic_carrington')
hpsc_cl2=sc_cl2.transform_to(frames.Helioprojective)
x_hpsc_cl2=hpsc_cl2.Tx
y_hpsc_cl2=hpsc_cl2.Ty

np.savetxt('/home/i4ds1807205/Dropbox/20130903/20130903_hcs_B_cl1.txt',np.array([x_hpsc_cl1,y_hpsc_cl1]).T)
np.savetxt('/home/i4ds1807205/Dropbox/20130903/20130903_hcs_B_cl2.txt',np.array([x_hpsc_cl2,y_hpsc_cl2]).T)
np.savetxt('/home/i4ds1807205/Dropbox/20130903/20130903_hcs_B_lcl.txt',np.array([x_hpsc_lcl,y_hpsc_lcl]).T)
np.savetxt('/home/i4ds1807205/Dropbox/20130903/20130903_hcs_B_mcl.txt',np.array([x_hpsc_mcl,y_hpsc_mcl]).T)
np.savetxt('/home/i4ds1807205/Dropbox/20130903/20130903_hcs_B_open.txt',np.array([x_hpsc_open,y_hpsc_open]).T)
plt.plot(x_hpsc_lcl,y_hpsc_lcl,'o-')
plt.plot(x_hpsc_mcl,y_hpsc_mcl,'o-')
plt.plot(x_hpsc_open,y_hpsc_open,'o-')
plt.show()


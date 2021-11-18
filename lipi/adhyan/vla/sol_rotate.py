import numpy as np
import sunpy
from sunpy.map import Map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates import frames
from astropy.wcs import WCS
from reproject import reproject_interp
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord




f='2021_03_21__23_59_15_006__STEREO-A_SECCHI_EUVI_171.jp2'
g='2021_03_21__18_59_57_350__SDO_AIA_AIA_171.jp2'
map_euvi=Map(f)
map_aia=Map(g)
out_shape=(4096,4096)

out_header=sunpy.map.make_fitswcs_header(out_shape,map_aia.reference_coordinate,scale=u.Quantity(map_euvi.scale),rotation_matrix=map_euvi.rotation_matrix,wavelength=map_euvi.wavelength)
#mars_header = sunpy.map.make_fitswcs_header(out_shape,mars_ref_coord,scale=u.Quantity(map_aia.scale),rotation_matrix=map_aia.rotation_matrix,instrument="AIA",wavelength=map_aia.wavelength)
out_wcs = WCS(out_header)
output, footprint = reproject_interp(map_euvi, out_wcs, out_shape)
outmap = sunpy.map.Map(output, out_header)
outmap.plot_settings = map_euvi.plot_settings
#######################
mars = get_body_heliographic_stonyhurst('mars', '2001-02-03')
mars_ref_coord = SkyCoord(map_aia.reference_coordinate.Tx,map_aia.reference_coordinate.Ty,obstime=map_aia.reference_coordinate.obstime,observer=mars,frame="helioprojective")

mars_header = sunpy.map.make_fitswcs_header(out_shape,mars_ref_coord,scale=u.Quantity(map_aia.scale),rotation_matrix=map_aia.rotation_matrix,instrument="AIA",wavelength=map_aia.wavelength)
mars_wcs = WCS(mars_header)
mars_output, footprint = reproject_interp(map_aia, mars_wcs, out_shape)
marsmap = sunpy.map.Map((mars_output, mars_header))

########################
so_cor = SkyCoord(-107.042, 1.406, 0.699954, unit=('deg','deg','AU'), obstime='2021-03-21T18:59:57.350', frame='heliographic_stonyhurst')
so_ref_coord=SkyCoord(map_euvi.reference_coordinate.Tx,map_euvi.reference_coordinate.Ty,obstime=map_euvi.reference_coordinate.obstime,observer=so_cor,frame="helioprojective")
so_header=sunpy.map.make_fitswcs_header(out_shape,so_ref_coord,
                    scale=u.Quantity(map_euvi.scale),instrument="STIX",observatory="SO",wavelength=map_euvi.wavelength)
so_wcs = WCS(so_header)
so, footprint = reproject_interp(map_euvi, so_wcs, out_shape)
somap = sunpy.map.Map(so, so_header)
########################

so_ref_coorda=SkyCoord(map_aia.reference_coordinate.Tx,map_aia.reference_coordinate.Ty,obstime=map_aia.reference_coordinate.obstime,observer=so_cor,frame="helioprojective")
so_headera=sunpy.map.make_fitswcs_header(out_shape,so_ref_coorda,
                    scale=u.Quantity(map_aia.scale),instrument="AIA",observatory="SDO",wavelength=map_aia.wavelength)
so_wcsa = WCS(so_headera)
soa, footprint = reproject_interp(map_aia, so_wcsa, out_shape)
somap_a = sunpy.map.Map(soa, so_headera)

fig = plt.figure()
ax1 = fig.add_subplot(1, 4, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 4, 2, projection=map_euvi)
map_euvi.plot(axes=ax2, title='EUVI image')
ax3 = fig.add_subplot(1, 4, 3, projection=outmap)
outmap.plot(axes=ax3, title='EUVI image as seen from SDO')
ax4 = fig.add_subplot(1, 4, 4, projection=somap)
somap.plot(axes=ax4,cmap='euvi171', title='EUVI image as seen from SO')
#ax5 = fig.add_subplot(1, 5, 5, projection=marsmap)
#marsmap.plot(axes=ax5,cmap='jet', title='AIA image as seen from Mars')
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(1, 3, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 3, 2, projection=map_euvi)
map_euvi.plot(axes=ax2, title='EUVI image')
ax3 = fig.add_subplot(1, 3, 3, projection=somap)
somap.plot(axes=ax3, cmap='euvi171',title='EUVI image as seen from SO')
plt.show()

sys.exit()










h=m.reference_coordinate.transform_to(frames.HeliographicCarrington())

out_shape = (4096, 4096)
out_header = sunpy.map.make_fitswcs_header(out_shape,m.reference_coordinate,scale=u.Quantity(m.scale),instrument='EUVI',wavelength=171*u.Angstrom)
out_wcs = WCS(out_header)
#earth = get_body_heliographic_stonyhurst('earth', m.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(m, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))




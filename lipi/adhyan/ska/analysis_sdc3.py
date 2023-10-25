import oskar
from oskar.measurement_set import MeasurementSet
from karabo.imaging.imager import Imager
from karabo.simulation.visibility import Visibility
import numpy as np

res_filename='/media/rohit/Seagate_Expansion_Drive/SDC3/sdc3point_ch750_4h1d.MS'

res_vis=Visibility.read_from_file(res_filename)

imager = oskar.Imager()
imager.set(fov_deg=4, image_size=2048)
imager.set(input_file=)
output = imager.run(return_images=1)
image = output["images"][0]

imager = Imager(res_vis,imaging_npixel=2048,imaging_cellsize=np.deg2rad(50./3600),imaging_weighting='natural')
dirty = imager.get_dirty_image()
dirty.write_to_file('/media/rohit/Seagate_Expansion_Drive/SDC3/residual_sdc3point_ch750_4h1d.fits')


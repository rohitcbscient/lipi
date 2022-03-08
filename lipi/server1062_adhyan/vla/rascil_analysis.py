import numpy as np
#%env RASCIL_DATA='/nas08-data02/rohit/rascil'
import rascil
from rascil.processing_components.simulation import create_named_configuration
from astropy.coordinates import SkyCoord
from astropy import units as u
from rascil.processing_components.skycomponent import create_skycomponent
from rascil.processing_components import create_test_image
from rascil.data_models.polarisation import PolarisationFrame
from rascil.processing_components.visibility import create_blockvisibility
from astropy.time import Time

#lowcore = create_named_configuration('LOWBD2-CORE')
lowcore = create_named_configuration('LOWBD2')
#utc_time=Time("2021-09-01T05:00:00",format='isot', scale='utc')
times = np.linspace(-3, +3, 13) * (np.pi / 12.0)

frequency = np.array([1e8])
channel_bandwidth = np.array([1e7])

# Define the component and give it some polarisation and spectral behaviour
f = np.array([100.0])
flux = np.array([f])

phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=35.0 * u.deg, frame='icrs', equinox='J2000')
compabsdirection = SkyCoord(ra=17.0 * u.deg, dec=36.5 * u.deg, frame='icrs', equinox='J2000')

comp = create_skycomponent(flux=flux, frequency=frequency, direction=compabsdirection,
                                polarisation_frame=PolarisationFrame('stokesI'))
#image = create_test_image(frequency=frequency, phasecentre=phasecentre,cellsize=0.001,polarisation_frame=PolarisationFrame('stokesI'))

vis = create_blockvisibility(lowcore, times=times, frequency=frequency,
                             channel_bandwidth=channel_bandwidth,
                             phasecentre=phasecentre, weight=1,
                             polarisation_frame=PolarisationFrame('stokesI'),
                             integration_time=1.0)



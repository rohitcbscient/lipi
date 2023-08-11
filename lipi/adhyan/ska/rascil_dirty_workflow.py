import os
import karabo
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import numpy as np
from ska_sdp_datamodels.science_data_model.polarisation_model import PolarisationFrame
from ska_sdp_datamodels.visibility import create_visibility
from ska_sdp_datamodels.configuration import (
    create_named_configuration,
)
import pandas
from ska_sdp_func_python.image.deconvolution import (
    deconvolve_cube,
)
from ska_sdp_datamodels.configuration.config_coordinate_support import (
    hadec_to_azel,
    xyz_at_latitude,
    xyz_to_uvw,
)

from ska_sdp_datamodels.sky_model.sky_model import (
    SkyComponent,
    SkyModel,
)
from rascil.workflows.rsexecute.execution_support.rsexecute import rsexecute
from ska_sdp_func_python.imaging import (
    invert_visibility,
    invert_ng,
    predict_visibility,
    invert_visibility,
    create_image_from_visibility,
)  # IMPORTANT

from ska_sdp_func_python.sky_component import (
    insert_skycomponent,
    apply_beam_to_skycomponent,
    filter_skycomponents_by_flux,
)
from rascil.processing_components.imaging.primary_beams import create_pb

from rascil.processing_components import create_test_image
from rascil.processing_components.simulation import (
    ingest_unittest_visibility,
    create_unittest_model,
    create_unittest_components,
)

from ska_sdp_datamodels.visibility.vis_utils import (
    calculate_transit_time,
    generate_baselines,
)

results_dir = "/home/rohit/simulations/rascil_results/"
# Construct LOW core configuration
lowr3 = create_named_configuration("MID", rmax=750.0)
# We create the visibility. This just makes the uvw, time, antenna1, antenna2,
# weight columns in a table. We subsequently fill the visibility value in by
# a predict step.
times = np.arange(10)
frequency = np.array([1e8])
channel_bandwidth = np.array([1e6])
phasecentre = SkyCoord(
    ra=+20.0 * u.deg, dec=-30.0 * u.deg, frame="icrs", equinox="J2000"
)
utc_time_obs = Time("2024-01-01T10:00:00", format="isot", scale="utc")
stime = calculate_transit_time(lowr3.location, utc_time_obs, phasecentre)
vt = create_visibility(
    lowr3,
    times,
    frequency,
    channel_bandwidth=channel_bandwidth,
    weight=1.0,
    phasecentre=phasecentre,
    polarisation_frame=PolarisationFrame("stokesI"),
    utc_time=utc_time_obs,
)

# advice = advise_wide_field(
#    vt, guard_band_image=3.0, delA=0.1, oversampling_synthesised_beam=4.0
# )
# cellsize = advice["cellsize"]
cellsize = 0.000552  # In radians
# Read the venerable test image, constructing a RASCIL Image
m31image = create_test_image(
    cellsize=cellsize, frequency=frequency, phasecentre=vt.phasecentre
)
vt = predict_visibility(vt, m31image, context="2d")

model = create_image_from_visibility(vt, cellsize=cellsize, npixel=512)
print(model.pixels.data.max())
dirty, sumwt = invert_visibility(vt, model, context="2d")
print(dirty.pixels.data.max())
psf, sumwt = invert_ng(vt, model, context="2d", dopsf=True)
comp, residual = deconvolve_cube(
    dirty,
    psf,
    niter=10000,
    threshold=0.001,
    fractional_threshold=0.001,
    window_shape="quarter",
    gain=0.7,
    scales=[0, 3, 10, 30],
)


imagename = "imaging_dirty.fits"
psffilename = "psf.fits"
input_model="input_model.fits"
print(os.path.join(results_dir, imagename))
dirty.image_acc.export_to_fits(fits_file=os.path.join(results_dir, imagename))
psf.image_acc.export_to_fits(fits_file=os.path.join(results_dir, psffilename))
m31image.image_acc.export_to_fits(fits_file=os.path.join(results_dir, input_model))



###-----------------------------------------------------------------------------

results_dir = "/home/rohit/simulations/rascil_results/"
# Construct LOW core configuration
lowr3 = create_named_configuration("MID", rmax=750.0)
# We create the visibility. This just makes the uvw, time, antenna1, antenna2,
# weight columns in a table. We subsequently fill the visibility value in by
# a predict step.
times = np.arange(10)
frequency = np.array([1e8])
channel_bandwidth = np.array([1e6])
phasecentre = SkyCoord(
    ra=+20.0 * u.deg, dec=-30.0 * u.deg, frame="icrs", equinox="J2000"
)
utc_time_obs = Time("2024-01-01T10:00:00", format="isot", scale="utc")
stime = calculate_transit_time(lowr3.location, utc_time_obs, phasecentre)
vt = create_visibility(
    lowr3,
    times,
    frequency,
    channel_bandwidth=channel_bandwidth,
    weight=1.0,
    phasecentre=phasecentre,
    polarisation_frame=PolarisationFrame("stokesI"),
    utc_time=utc_time_obs,
)


""" latitude = lowr3.location.geodetic[1].to("rad").value
ants_xyz = lowr3["xyz"].data
ants_xyz = xyz_at_latitude(ants_xyz, latitude)
nants = len(lowr3["names"].data)
ant_pos = xyz_to_uvw(ants_xyz, 0, phasecentre.dec.rad)
baselines = pandas.MultiIndex.from_tuples(
    generate_baselines(nants), names=("antenna1", "antenna2")
) """

skycomps = list()
compabsdirection = SkyCoord(ra=20.0 * u.deg, dec=-30.0 * u.deg, frame='icrs', equinox='J2000')
flux_point=np.array([np.array([100.0])])
skycomps.append(SkyComponent(direction=compabsdirection,
                        flux=flux_point,
                        frequency=frequency,
                        shape="Point",
                        polarisation_frame=PolarisationFrame('stokesI'),
                    ))

beam = create_pb(model, telescope="MID", use_local=False)
beam_data=beam.data_vars.variables['pixels']
skycomps = apply_beam_to_skycomponent(skycomps, beam)
i=0;nfreqwin = 1
channelwidth=np.array([10])
bvis_list = [
    rsexecute.execute(ingest_unittest_visibility, nout=1)(
        lowr3,
        [frequency[i]],
        [channelwidth[i]],
        times,
        PolarisationFrame('stokesI'),
        phasecentre,
        zerow=False,
    )
    for i in range(nfreqwin)
]

rsexecute.execute(create_unittest_components, nout=1)(
            lowr3,
            [frequency[i]],
            [channelwidth[i]],
            times,
            PolarisationFrame('stokesI'),
            phasecentre,
            zerow=False,
        )

components_list = [
        rsexecute.execute(create_unittest_components)(
            model_imagelist[freqwin],
            flux[freqwin, :][numpy.newaxis, :],
            offset=(offset, 0.0),
        )
        for freqwin, m in enumerate(model_imagelist)
    ]
    components_list = rsexecute.persist(components_list)
    
################-----------------
sc = create_low_test_skycomponents_from_gleam(flux_limit=1.0,
                                            polarisation_frame=PolarisationFrame("stokesIQUV"),
                                            frequency=frequency, kind='cubic',
                                            phasecentre=phasecentre,
                                            radius=0.1)
model = create_image_from_visibility(vis, cellsize=0.001, npixel=512, frequency=frequency,
                                    polarisation_frame=PolarisationFrame('stokesIQUV'))

bm = create_low_test_beam(model=model)
sc = apply_beam_to_skycomponent(sc, bm)
vis = dft_skycomponent_visibility(vis, sc)


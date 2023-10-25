import numpy as np
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from datetime import datetime, timedelta


sky = SkyModel()
sky_data=np.zeros((1,13))
sky_data[:,0]=20;sky_data[:,1]=-30;sky_data[:,2]=1
sky.add_point_sources(sky_data)
ms_path='test_ms.ms'
vis_path='test_vis.vis'

telescope = Telescope.get_SKA1_MID_Telescope()
path_telescope = 'telescope.tm'
telescope = Telescope.read_from_file(path_telescope)
simulation = InterferometerSimulation(
            ms_file_path=ms_path,
            vis_path=vis_path,
            channel_bandwidth_hz=1e6,
            time_average_sec=1,
            noise_enable=True,
            noise_seed="time",
            noise_freq="Range",
            noise_rms="Range",
            noise_start_freq=1.0e9,
            noise_inc_freq=1.0e8,
            noise_number_freq=24,
            noise_rms_start=5000,
            noise_rms_end=10000,
        )
observation = Observation(
            phase_centre_ra_deg=20.0,
            start_date_and_time=datetime(2022, 9, 1, 23, 00, 00, 521489),
            length=timedelta(hours=0, minutes=0, seconds=1, milliseconds=0),
            phase_centre_dec_deg=-30.5,
            number_of_time_steps=1,
            start_frequency_hz=1.0e9,
            frequency_increment_hz=1e6,
            number_of_channels=1,
        )

visibility = simulation.run_simulation(telescope, sky, observation)


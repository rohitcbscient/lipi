import casacore.tables as ct
from oskar.measurement_set import MeasurementSet
import oskar
import os
import numpy as np
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pytest
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.sourcedetection.evaluation import SourceDetectionEvaluation
from karabo.sourcedetection.result import (
    PyBDSFSourceDetectionResult,
    SourceDetectionResult,
)
from karabo.test.conftest import NNImageDiffCallable, TFiles
from astropy.io import fits
import bdsf
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import Angle
import astropy.units as u
from datetime import datetime, timedelta
import tools21cm as t2c
import sys
from astropy.wcs import WCS

path='/scratch/snx3000/rsharma/galactic_fg_test/'
filename0=path+'residual_sdc3point_ch0000_04.MS'
dc=ct.table(filename0)
uvw =dc.getcol("UVW")
ms_data=	dc_uvw =dc.getcol("DATA")
for k in range(1,5):
	kk="%04d" % k
	filename_low=path+'residual_sdc3point_ch'+kk+'_04.MS'
	#filename_low='/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_'+kk+'.MS'
	dc=ct.table(filename_low)
	dc_data=dc.getcol('DATA')
	ms_data=np.nanmean(np.vstack((ms_data,dc_data[:,0,0])),axis=0)

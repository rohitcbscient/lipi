import casacore.tables as ct
import numpy as np
import sys
import os

i = int(os.environ['SLURM_ARRAY_TASK_ID'])+790

ii="%04d" % i
skadc=ct.table('/store/ska/sk01/sdc3-new/MS/ZW3_IFRQ_'+ii+'.ms')
skadc_u = skadc.getcol("UVW")[:,0].reshape(1440,130816)
skadc_v = skadc.getcol("UVW")[:,1].reshape(1440,130816)
skadc_w = skadc.getcol("UVW")[:,2].reshape(1440,130816)
np.save('/scratch/snx3000/rsharma/ska_data_reduced/SKADC_U_'+ii+'.npy',skadc_u)
np.save('/scratch/snx3000/rsharma/ska_data_reduced/SKADC_V_'+ii+'.npy',skadc_v)
np.save('/scratch/snx3000/rsharma/ska_data_reduced/SKADC_W_'+ii+'.npy',skadc_w)

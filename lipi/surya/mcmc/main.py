import numpy as np


ParmLocal=29*[0] # the array of the volume element parameters
ParmLocal[0] = 3.e18  # Area, cm^2
ParmLocal[1] =depth#1e8  # Depth, cm (will be changed later!)
ParmLocal[2] =10e6   # T_0, K
ParmLocal[3] =0.05  # \eps (not used in this example)
ParmLocal[4] =4.0   # \kappa (not used in this example)
ParmLocal[5] =16    # number of integration nodes in energy (continuous code)
ParmLocal[6] = emin#0.005   # E_min, MeV
ParmLocal[7] =1.0  # E_max, MeV
ParmLocal[8] =1.0   # E_break, MeV (not used in this example)
ParmLocal[9] = 5.5   # \delta_1
ParmLocal[10]= 7.0   # \delta_2 (not used in this example)
ParmLocal[11]=nth#2e10   # n_0 - thermal electron density, cm^{-3}
ParmLocal[12]=nb#5e7   # n_b - nonthermal electron density, cm^{-3}
ParmLocal[13]=B#300   # B - magnetic field, G
ParmLocal[14]=60    # theta - the viewing angle, degrees (will be changed later!)
ParmLocal[15]=fs#1e9   # starting frequency to calculate spectrum, Hz
ParmLocal[16]=delf  # logarithmic step in frequency
ParmLocal[17]=3     # distribution over energy (PLW is chosen)
ParmLocal[18]=Nf    # number of frequencies (specified above)
ParmLocal[19]=3     # distribution over pitch-angle (GLC is chosen)
ParmLocal[20]=90    # loss-cone boundary, degrees
ParmLocal[21]=0     # beam direction (degrees) in GAU and SGA (not used in this example)
ParmLocal[22]=0.2   # \Delta\mu
ParmLocal[23]=1     # a_4 in SGA (not used in this example)
ParmLocal[25]=-1    # f^C_cr
ParmLocal[26]=ParmLocal[25]    # f^WH_cr
ParmLocal[27]=1     # matching on
ParmLocal[28]=2     # Q-optimization on

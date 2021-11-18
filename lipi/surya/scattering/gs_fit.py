# This is an example of calculating the left- and right-polarized intensity spectra from an inhomogeneous source.
# In this example, the magnetic field direction (relative to the line-of-sight) varies linearly from 80 to 110 degrees
# along the line-of-sight. All other parameters are assumed to be constant, but they are allowed to be variable as well.
# Function GET_MW uses only open intervals, messages are off, and the parameter \Delta\eta is set to the default value. 

import MWTransfer
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav

def mwt(depth,emin,nb,nth,B,fs,delf,nn):
#def mwt(depth,emin,B,fs,delf,nn):
    Nf=nn    # number of frequencies
    NSteps=30 # number of nodes along the line-of-sight
    # Firstly, we create one-dimensional array describing the properties of a single volume element.
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
    # Now, we create two-dimensional array which contains information about all volume elements along the line-of-sight.
    Parms=NSteps*[None] # the array of input parameters
    for i in range(NSteps):
        Parms[i]=ParmLocal[:]               # The parameters of all volume elements are initially the same...
        Parms[i][1]=Parms[i][1]/NSteps      # But the length of an elementary interval should equal the source depth
                                            # divided by the number of intervals (all intervals are assumed to have the
                                            # same length)... 
        Parms[i][14]=80.0+30.0*i/(NSteps-1) # And, finally, the viewing angle varies from 80 to 110 degrees.

    res=MWTransfer.GET_MW(Parms) # calling the main code

    # extracting the data from the output array
    if res:
        f=res[0]   # emission frequency (GHz)
        I_L=res[1] # left-hand polarized emission intensity, sfu (as observed at the Earth)
        I_R=res[2] # right-hand polarized emission intensity, sfu (as observed at the Earth)
    else:
        print("Calculation error!")

    # Now, we print the results.
    res_print=0
    if res_print:
        print("f=", f)
        print("I_L=", I_L)
        print("I_R=", I_R)
    S=np.array(I_L)+np.array(I_R)
    sfu2cgs=1.0e-19
    freq=np.array(f)*1.e9
    vc=2.998e10
    kb=1.38e-16
    rad2asec=180.0*3600/(np.pi)
    sr=ParmLocal[0]/(7.25e7)**2/rad2asec**2
    Tb=S*sfu2cgs*vc**2/(2.0*kb*freq*freq*sr)
    return Tb,freq



def mwt_test(B):
    print(B)
    Nf=400    # number of frequencies
    NSteps=30 # number of nodes along the line-of-sight
    # Firstly, we create one-dimensional array describing the properties of a single volume element.
    ParmLocal=29*[0] # the array of the volume element parameters
    ParmLocal[0] = 1e18  # Area, cm^2
    ParmLocal[1] =1e8  # Depth, cm (will be changed later!)
    ParmLocal[2] =10e6   # T_0, K
    ParmLocal[3] =0.05  # \eps (not used in this example)
    ParmLocal[4] =4.0   # \kappa (not used in this example)
    ParmLocal[5] =16    # number of integration nodes in energy (continuous code)
    ParmLocal[6] =0.005   # E_min, MeV
    ParmLocal[7] =1.0  # E_max, MeV
    ParmLocal[8] =1.0   # E_break, MeV (not used in this example)
    ParmLocal[9] =6.5   # \delta_1
    ParmLocal[10]=7.5   # \delta_2 (not used in this example)
    ParmLocal[11]=2e10   # n_0 - thermal electron density, cm^{-3}
    ParmLocal[12]=5e7   # n_b - nonthermal electron density, cm^{-3}
    ParmLocal[13]=B#300   # B - magnetic field, G
    ParmLocal[14]=60    # theta - the viewing angle, degrees (will be changed later!)
    ParmLocal[15]=1e9   # starting frequency to calculate spectrum, Hz
    ParmLocal[16]=0.002  # logarithmic step in frequency
    ParmLocal[17]=3     # distribution over energy (PLW is chosen)
    ParmLocal[18]=400    # number of frequencies (specified above)
    ParmLocal[19]=3     # distribution over pitch-angle (GLC is chosen)
    ParmLocal[20]=90    # loss-cone boundary, degrees
    ParmLocal[21]=0     # beam direction (degrees) in GAU and SGA (not used in this example)
    ParmLocal[22]=0.2   # \Delta\mu
    ParmLocal[23]=1     # a_4 in SGA (not used in this example)
    ParmLocal[25]=-1    # f^C_cr
    ParmLocal[26]=ParmLocal[25]    # f^WH_cr
    ParmLocal[27]=1     # matching on
    ParmLocal[28]=2     # Q-optimization on
    # Now, we create two-dimensional array which contains information about all volume elements along the line-of-sight.
    Parms=NSteps*[None] # the array of input parameters
    for i in range(NSteps):
        Parms[i]=ParmLocal[:]               # The parameters of all volume elements are initially the same...
        Parms[i][1]=Parms[i][1]/NSteps      # But the length of an elementary interval should equal the source depth
                                            # divided by the number of intervals (all intervals are assumed to have the
                                            # same length)... 
        Parms[i][14]=80.0+30.0*i/(NSteps-1) # And, finally, the viewing angle varies from 80 to 110 degrees.

    res=MWTransfer.GET_MW(Parms) # calling the main code
    # extracting the data from the output array
    if res:
        f=res[0]   # emission frequency (GHz)
        I_L=res[1] # left-hand polarized emission intensity, sfu (as observed at the Earth)
        I_R=res[2] # right-hand polarized emission intensity, sfu (as observed at the Earth)
    else:
        print("Calculation error!")

    # Now, we print the results.
    res_print=0
    if res_print:
        print("f=", f)
        print("I_L=", I_L)
        print("I_R=", I_R)
    S=np.array(I_L)+np.array(I_R)
    sfu2cgs=1.0e-19
    freq=np.array(f)*1.e9
    vc=2.998e10
    kb=1.38e-16
    rad2asec=180.0*3600/(np.pi)
    sr=ParmLocal[0]/(7.25e7)**2/rad2asec**2
    Tb=S*sfu2cgs*vc**2/(2.0*kb*freq*freq*sr)
    return Tb,freq

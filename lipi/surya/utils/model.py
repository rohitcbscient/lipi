# NAME  quiet_sun
#
# PURPOSE
#    1) calculates newkirk
#    2) calculates saito
#    3) 
#
# OUTPUTS
#    1) nk_ne,nk_r,nk_fp
#    2) st_ne,st_r,st_fp
#
# INPUTS
#
# HISTORY
#    written 2019-04-15 R. Sharma


import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
from scipy import interpolate
from scipy.integrate import quad


def plot_freq_vs_distance():
	plt.plot(freq_mwa,r_ne_mwa)
	plt.ylabel('Distance from Sun\'s centre (in $R_{\odot}$)')
	plt.xlabel('Frequency (MHz)')
	plt.show()
def plot_all_models():	
	plt.plot(r,ne_nk,'-',label='Newkirk')
	plt.plot(r,ne_nk4,'-',label='Newkirk*4')
	plt.plot(r,ne_st,'-',label='Saito')
	plt.plot(r,ne_st10,'-',label='Saito*10')
	plt.xlabel('Distance from Sun\'s centre (in $R_{\odot}$)')
	plt.ylabel('Electron Density (cm$^{-3}$)')
	plt.legend(loc=1)
	plt.show()

	
def freq2dens(f):
	'''
	name: Newkirk model
	input: plasma frequency (MHz)
	onput: electron density (cm^-3)
	'''
	ne=((f*1.e6)/9000)**2
	return ne

def nk_freq2r(f,n):
	'''
	name: Newkirk model
	input: f:plasma frequency (MHz),n:number to mutilply to density
	output: r:height (Rsun from centre), r_array, ne_array
	'''
	ne=freq2dens(f)
	r_=np.linspace(1,3,1000) # Units: Rsun
	ne_=n*(4.2*10**(4+(4.32/r_))) # Unit: cm^{-3}
	f = interpolate.interp1d(ne_,r_)
	r=f(ne)
	return r,ne,r_,ne_


def st_freq2r(f,n):
	'''
	name: Newkirk model
	input: f:plasma frequency (MHz),n:number to mutilply to density
	output: r:height (Rsun from centre), r_array, ne_array
	'''
	ne=freq2dens(f)
	r_=np.linspace(1,3,1000) # Units: Rsun
	ne_=n*(1.36e12*(1/r_**2.14)+1.68e14*(1/r_**6.13))*1.0e-6 # Unit: cm^{-3}
	f = interpolate.interp1d(ne_,r_)
	r=f(ne)
	return r,ne,r_,ne_

def refractive_index(ne,freq):
	'''
	name: computes refractive index
	input:ne(cm^-3),freq(MHz)
	output:n
	'''
	n=1-80.6*1.e6*ne/(freq*1.e6)**2
	return n	

def kappa(freq,ne,Te):
	'''
	name: compute optical depth
	input:freq(MHz),ne(cm^-3),Te(in K)
	output:kappa
	'''
	n=refractive_index(ne,freq)
	kappa_=(1.e-11)*ne*ne*Te**(-1.5)/(freq*freq*1.e12*n)
	return kappa_

def integrand(f,r):
    ne_nk=(4.2*10**(4+(4.32/r)))
    fpe=9000*np.sqrt(ne_nk)
    return fpe #fpe**4 /(f*f-fpe*fpe)**2




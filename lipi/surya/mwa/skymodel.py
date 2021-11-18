import numpy as np
import pickle
import ephem
import time
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import griddata as grd
##############################################################################################################################################
def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):	# Find xy index of haslam map that correspons to a (Alt,Azimuth) pair
	distance = np.sqrt((y_array-y_point)**2 + (x_array-x_point)**2)
	idy,idx = np.where(distance==distance.min())
	return idy[0],idx[0],distance 

def gal2local(ra,dec,obs):
	body = ephem.FixedBody()
	body._ra = ra
	body._dec = dec
	body.compute(obs)
	return body.alt*180/np.pi,body.az*180/np.pi

def gal2azel(haslam_ra,haslam_dec,obs):			# Convert Haslam map coordinates into Az-El coordinates
	alt_haslam=np.zeros((2047,4096))
	azi_haslam=np.zeros((2047,4096))
	print 'Calculating az-el..'
	for i in range(2047):
		for j in range(4096):
			alt_haslam[i,j], azi_haslam[i,j] =gal2local(haslam_ra[i,j],haslam_dec[i,j],obs)
	return (alt_haslam,azi_haslam)

def Tsky_408(haslam_gal,haslam_ra,haslam_dec,array_lon,array_lat,obs_date,azimuth_bins,elevation_bins,spec_index,theta_resol,phi_resol,file_):	# Convert Tsky values and spectral index values from galactic coordinates to Az-El coordinates
	obs = ephem.Observer()	
	obs.lon = str(array_lon)
	obs.lat = str(array_lat)
	obs.date =  str(obs_date)
	alt_has,azi_has=gal2azel(haslam_ra,haslam_dec,obs)
	Tsky =np.zeros((azimuth_bins,elevation_bins))
	#spectr_index = np.ones((azimuth_bins,elevation_bins))*2.55
	alt_flat = alt_has.flatten()
	azi_flat = azi_has.flatten()
	grid_x,grid_y = np.meshgrid(np.arange(0,360,phi_resol),np.arange(0,90,theta_resol))
	points = np.array((azi_flat,alt_flat)).T
	Tsky = grd(points,haslam_gal[0:2047].flatten(),(grid_x,grid_y),method='cubic').T
	Tsky[0] = Tsky[-1]		# To make azimuth a circular coordinate
	spectr_index = grd(points,spec_index[0:2047].flatten(),(grid_x,grid_y),method='cubic').T
	spectr_index[0] = spectr_index[-1]
	spectr_index[np.isnan(spectr_index)]=2.55
# 	spectr_index[0] = spectr_index[-1]
	pickle.dump([alt_has,azi_has,Tsky,spectr_index],open(str(file_)+'.p','wb'))
	return alt_has,azi_has,Tsky,spectr_index

def sky_model(T_sk, spectr_index,freq):		# Power law scaling to get Tsky at freq
	Tsk = T_sk*((408.0/freq)**spectr_index)
	return Tsk
###############################################################################################################################################

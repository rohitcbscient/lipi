import pickle
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import glob
import itertools
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import RectSphereBivariateSpline
import thread
import time
from surya.mwa import skymodel
#import skymodel
import ephem
import datetime
import astropy
from surya.mwa import primarybeam
from surya.mwa import fringefactor
import os
import math
import sys
#import plot
###############################################################################################################################################
def midtime_calculate(start_datetime,end_datetime):		#Calculate midtime of observation
	starttime = start_datetime.split(' ')[1]
	date = start_datetime.split(' ')[0]		# Assuming date remains the same throughout the duration of the observation
	endtime = end_datetime.split(' ')[1]
	start_seconds = int(starttime.split(':')[0])*3600 + int(starttime.split(':')[1])*60 +  int(starttime.split(':')[2])
	end_seconds = int(endtime.split(':')[0])*3600 + int(endtime.split(':')[1])*60 +  int(endtime.split(':')[2])
	mid_temp = (start_seconds+end_seconds)/2
	mid_hour = mid_temp/3600
	mid_temp -= 3600*mid_hour
	mid_minutes = mid_temp/60
	mid_temp -= mid_minutes*60
	mid_seconds = mid_temp
	mid_time = date + ' ' + str(mid_hour) + ':' + str(mid_minutes) + ':' + str(mid_seconds)
	return mid_time

### STEP 2
def Trec_interpolate(flist,trec):		#Spline interpolation for Trec at intermediate frequencies
	trec_interp = interp1d(trec[0], trec[1])
	Trec=trec_interp(flist)
	return Trec

def Tpickup_interpolate(flist,tpickup):		#Spline interpolation for Tpickup at intermediate frequencies
	tpickup_interp = interp1d(tpickup[0], tpickup[1])
	Tpickup=tpickup_interp(flist)
	return Tpickup

###STEP 3

def solar_coords(array_lon,array_lat,array_elevation,start_time):
	rad2deg = float(180.0/math.pi)
	obs = ephem.Observer()
	obs.lon = str(array_lon)
	obs.lat = str(array_lat)
	obs.elevation = array_elevation
	obs.date = start_time
	sun = ephem.Sun()
	sun.compute(obs)
	return sun.alt*rad2deg,sun.az*rad2deg,sun.ra*rad2deg,sun.dec*rad2deg

def solid_angle(azimuth_bins,elevation_bins):		# Solid angle area element
	area=np.zeros((azimuth_bins,elevation_bins))
	grid_x, grid_y = np.mgrid[0:azimuth_bins+1:361j,0:elevation_bins+1:91j]
	for i in range(azimuth_bins):
		for j in range(elevation_bins):
			area[i,j]=(grid_x[i+1,j]-grid_x[i,j])*(grid_y[i,j+1]-grid_y[i,j])*np.cos(grid_y[i,j]*np.pi/180)/((180/np.pi)**2) # radian**2
	return area	


#### STEP 4

class ds:    
    def __init__(self,T_sun,S_sun,Temp_beam_sun,Un_Tbeam_Sun,corr_factor,fringe_factor,T_baseline,Tsky_integrated):
        self.Tsun = T_sun
	self.Ssun = S_sun
	self.Temp_beam_sun=Temp_beam_sun
	self.Un_Tbeam_Sun=Un_Tbeam_Sun
	self.corr_factor=corr_factor
	self.T_baseline=T_baseline
	self.Tsky_integrated=Tsky_integrated
	self.fringe_factor=fringe_factor
        self.I = T_sun[0]**2+T_sun[3]**2
	self.Q = T_sun[0]**2-T_sun[3]**2

#####################################################################################################################
def star_az_el(StartYear,StartMonth,StartDate,StartHour,StartMin,StartSec,ra,dec,lon,lat,ele):
	rad2deg = float(180.0/math.pi)
	StartTime = datetime.datetime(StartYear, StartMonth, StartDate, StartHour, StartMin, StartSec)
	StartTimeStr = StartTime.strftime("%Y-%m-%d,%H:%M:%S")
	array_128 = ephem.Observer()
	# Coordinates of Tile011MWA
	array_128.long, array_128.lat = str(lon), str(lat)
	array_128.elevation = ele
	array_128.date = StartTime.strftime("%Y/%m/%d %H:%M:%S")
        star = ephem.FixedBody()
	eq = ephem.Equatorial(ra*np.pi/180, dec*np.pi/180, epoch='2000')
        star._ra = eq.ra
        star._dec = eq.dec
	star.compute(array_128)
	star_al=float(star.alt)*rad2deg
	star_az=float(star.az)*rad2deg
	return star_az,star_al,StartTimeStr

def solar_az_el(StartYear,StartMonth,StartDate,StartHour,StartMin,StartSec):
	rad2deg = float(180.0/math.pi)
	StartTime = datetime.datetime(StartYear, StartMonth, StartDate, StartHour, StartMin, StartSec)
	StartTimeStr = StartTime.strftime("%Y-%m-%d,%H:%M:%S")
	array_128 = ephem.Observer()
	# Coordinates of Tile011MWA
	array_128.long, array_128.lat = '116.66931', '-26.54675'
	array_128.elevation = 377.83
	array_128.date = StartTime.strftime("%Y/%m/%d %H:%M:%S")
	Sun = ephem.Sun(array_128.date)
	Sun.compute(array_128)
	sun_al=float(Sun.alt)*rad2deg
	sun_az=float(Sun.az)*rad2deg
	return sun_az,sun_al,StartTimeStr

def radec2azel(ra,dec,date):
	obs = ephem.Observer()
	obs.lon = '116.6708'
	obs.lat = '-26.7033'
	obs.elevation = 377.83
	obs.date =  str(date)
	eq = ephem.Equatorial(ra*np.pi/180, dec*np.pi/180, epoch='2000')
	body = ephem.FixedBody()
	body._ra = eq.ra
	body._dec = eq.dec
	body.compute(obs)
	return body.alt*180/np.pi,body.az*180/np.pi

###############################################################################################################################################
def tsun_computation(DS,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec):
	if (os.path.isfile(str(WORKING_DIR)+'flux_V1_'+DS.split('.')[0] + '.p')==False):
		phi_resolution = 1.0		# Resolution in azimuth (in degrees)
		theta_resolution = 1.0		# Resolution in elevation (in degrees)
		azimuth_bins = int(360/phi_resolution)
		elevation_bins = int(90/theta_resolution)
		solar_size = 22.5		#Radius of Sun (in arcmin)
		solidangle_sun = (np.pi)*((solar_size*np.pi/(180*60))**2)	#Solid angle subtended by Sun for solar radius = 22.5'
		kb=1.38e-23		# Boltzmann constant
		c=3.e8			# Speed of light in vacuum
		time1 = time.time()
		Trec_array=pickle.load(open(rec_path,'r'))		# Min. freq. = 69.59MHz, Max freq. = 300.26MHz
		Tpickup_array = pickle.load(open(grd_path,'r'))	# Min. freq.= 77.5MHz, Max. freq. = 299.116MHz Interpolation works only for intermediate frequencies
		#haslam_gal,spec_index,gal_lat,gal_lon,haslam_ra,haslam_dec = pickle.load(open('/home/rohit/scripts/flux_calibration/haslam_spec_gal.p','r'))
		haslam_gal,spec_index,gal_lat,gal_lon,haslam_ra,haslam_dec = pickle.load(open(spidx_path,'r'))
		solidangle = solid_angle(azimuth_bins,elevation_bins)
		print '######## Working with ',str(DS_DIR),str(DS),' ###########'
		os.chdir(DS_DIR)
		#central_freq,chan_list,auto_t1,auto_t2,cross,ncross,phase_ncross,u,v,w,azi_pointing,ele_pointing,ph_ra,ph_dec,start_time,end_time = pickle.load(open(DS,'r'))
		central_freq,chan_list,auto_t1,auto_t2,cross,ncross,phase_ncross,u,v,w,azi_pointing,ele_pointing,start_time,end_time = pickle.load(open(DS,'r'))
		if(mwa_phase==2):
			ncross=np.mean(ncross,axis=2)
		ncross=np.mean(ncross,axis=2)
		ph_ra,ph_dec=azi_pointing,ele_pointing
		#azi_pointing,ele_pointing,ph_ra,ph_dec = 	308.66, 40.61,	158.79, 8.72
		mid_time = midtime_calculate(start_time,end_time)
		print 'Mid Time: ',mid_time
		ph_el,ph_az= ele_pointing,azi_pointing#radec2azel(ph_ra*180/np.pi,ph_dec*180/np.pi,mid_time)
		print 'Phase centre elevation: ',ph_el,'Phase centre azimuth: ',ph_az
		os.chdir(WORKING_DIR)
		no_pol = len(ncross)		#No. of polarisations
		no_chan = len(central_freq)		#No. of channels
		Trec = np.zeros(no_chan)
		Tpickup = np.zeros(no_chan)
		T_sun = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		S_sun = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		Un_Tbeam_Sun = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		corr_factor = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		fringe_factor = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		Tsky_integrated = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		T_baseline = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		Temp_beam_sun = np.zeros((no_pol,len(ncross[0]),len(ncross[0][0])))
		alt_sun,az_sun,ra_sun, dec_sun = solar_coords(array_lon,array_lat,array_elevation,mid_time)	# All quantities in degrees
		Trec = Trec_interpolate(central_freq,Trec_array)
		Tpickup = Tpickup_interpolate(central_freq,Tpickup_array)
		#alt_haslam_s,azi_haslam_s,T_sky_s,spectral_index_array = skymodel.Tsky_408(haslam_gal,haslam_ra,haslam_dec,array_lon,array_lat,start_time,azimuth_bins,elevation_bins,spec_index,theta_resolution,phi_resolution)
		#alt_haslam_e,azi_haslam_e,T_sky_e,spectral_index_array = skymodel.Tsky_408(haslam_gal,haslam_ra,haslam_dec,array_lon,array_lat,end_time,azimuth_bins,elevation_bins,spec_index,theta_resolution,phi_resolution)

		#print T_sky
		# Sky map at 408 MHz and spectral index map in (Az,Alt) coordinates 
		for ch in range(no_chan):
			freq = central_freq[ch]
			primary_beam = primarybeam.beam(np.round(azi_pointing,2),np.round(ele_pointing,2),freq,azimuth_bins,elevation_bins,BEAM_DIR,WORKING_DIR)
			norm_primary_beam  = (10**(primary_beam[8]/10))/np.max(10**(primary_beam[8]/10))
			#norm_primary_beam = primarybeam.normalise(primary_beam,no_pol,azimuth_bins,elevation_bins)#Normalising primary beam as per polarisation
			print '\nch=',ch
			pola=np.array([0,3])
			nn=ncross.shape[2]-1#74
			for p in range(2):
				p=pola[p]
				print '\np=',p,' FILE: ',str(DS), len(ncross[p][ch]), nn
				for t in range(0,len(ncross[p][ch]),nn):			# Working with mean values of u,v,w per minute
					t_index=range(0,len(ncross[p][ch]),nn).index(t)
					date=start_time.split(' ')[0].split('/')
					hhmmss=start_time.split(' ')[1].split(':')
					sec=int(hhmmss[2])+37*(t_index+1)
					hour_=int(hhmmss[0])
					print 'sec:',sec, 'Phase centre elevation: ',ph_el,'Phase centre azimuth: ',ph_az
					for j in range(len(range(0,len(ncross[p][ch]),nn))):
						if(sec==60):
							sec=sec-1
						if(60*j<sec<60*(j+1)):
							sec_=sec-60*j
							min_=int(hhmmss[1])+j
							if(min_>=60):
								hour_=	int(hhmmss[0])+1
								min_= min_-60
                                        if(ifsun):
					    az_sun,alt_sun,StartTimeStr = solar_az_el(int(date[0]),int(date[1]),int(date[2]),hour_,min_,sec_)
                                        if(ifsun!=1):    
                                            az_sun,alt_sun,StartTimeStr = star_az_el(int(date[0]),int(date[1]),int(date[2]),hour_,min_,sec_,star_ra,star_dec,array_lon,array_lat,array_elevation)
                                        print 'Time: ',StartTimeStr,' Source Azimuth:',az_sun,' Source Altitude:',alt_sun
					#beam_factor = primarybeam.factor(norm_primary_beam[p],az_sun,alt_sun,azimuth_bins,elevation_bins)
					beam_factor = primarybeam.factor(norm_primary_beam,az_sun,alt_sun,azimuth_bins,elevation_bins)
					#plot_polar_dot(norm_primary_beam,'',az_sun,alt_sun)
					print 'Beam Factor :',beam_factor
					this_time=StartTimeStr.split(',')[0].split('-')[0]+'/'+StartTimeStr.split(',')[0].split('-')[1]+'/'+StartTimeStr.split(',')[0].split('-')[2]+' '+StartTimeStr.split(',')[1]
					if (os.path.isfile(str(HASLAM_DIR)+str(StartTimeStr)+'.p')== False):
						print 'Creating haslam array....'
						alt_haslam_m,azi_haslam_m,T_sky_m,spectral_index_array = skymodel.Tsky_408(haslam_gal,haslam_ra,haslam_dec,array_lon,array_lat,this_time,azimuth_bins,elevation_bins,spec_index,theta_resolution,phi_resolution,str(HASLAM_DIR)+str(StartTimeStr)) # np.ones((2048,4096)),np.ones((2048,4096)),np.ones((360,90)),np.ones((360,90))#
					alt_haslam_m,azi_haslam_m,T_sky_m,spectral_index_array = pickle.load(open(str(HASLAM_DIR)+str(StartTimeStr)+'.p','r'))
					T_sky_m[np.isnan(T_sky_m)] = 0
					Tsky_m = skymodel.sky_model(T_sky_m, spectral_index_array,freq)
                                        if(ifsun!=1):
                                            Tsky_m[int(az_sun)-2:int(az_sun)+2,int(alt_sun)-2:int(alt_sun)+2]=0
					beam_area,sum_beamarea = primarybeam.beamarea(norm_primary_beam,solidangle)
					#beam_area,sum_beamarea = primarybeam.beamarea(norm_primary_beam[p],solidangle)
					tgb_m = Tsky_m*beam_area/sum_beamarea
					T_sky_avg_m = np.sum(tgb_m)/1000		# Galactic background contribution (in K)
					print 'Beam averaged Galactic temperature: ',T_sky_avg_m,' K'
					Tsky_integrated[p][ch][t:(t+nn)] = T_sky_avg_m
					#print '********** Fringe factor **********************'
					fringe_factor_,fringe_cos, fringe_sin = fringefactor.fraction(ph_az,ph_el,az_sun,alt_sun,freq,u[t:(t+nn)],v[t:(t+nn)],w[t:(t+nn)],solar_size,azimuth_bins,elevation_bins,t)
					fringe_factor[p][ch][t:(t+nn)]=fringe_factor_
					tgb_cosine=tgb_m*fringe_cos[0:360,0:90][:,::-1]				# COSINE FRINGE ON SKY
					tgb_sin=tgb_m*fringe_sin[0:360,0:90][:,::-1]				# SINE FRINGE ON SKY
					tgb_abs=np.sqrt(np.sum(tgb_sin)**2 + np.sum(tgb_cosine)**2) 		# ABSOLUTE FRINGE ON THE SKY
					T_sky_baseline=np.sum(tgb_abs)/1000 					# Converting from mK to K
					print 'Fringe ::',T_sky_baseline#,tgb_abs,ph_az,ph_el,az_sun,alt_sun,freq,np.mean(u[t:(t+nn)]),np.mean(v[t:(t+nn)]),np.mean(w[t:(t+nn)])		
					#plot.plot_polar(tgb_sin,'')
					T_baseline[p][ch][t:(t+nn)] = T_sky_baseline
					corr_factor[p][ch][t:(t+nn)] = beam_factor#*fringe_factor[p][ch][t:(t+nn)]
					#ncross[p][ch][t:(t+nn)] = 0.5 ### PUT YOUR NCCF ###
					Un_Tbeam_Sun[p][ch][t:(t+nn)]=(ncross[p][ch][t:(t+nn)]*(Tsky_integrated[p][ch][t:(t+nn)]+Trec[ch]+Tpickup[ch])-T_sky_baseline)/(1-ncross[p][ch][t:(t+nn)])
			    	        Tbeam_sun = corr_factor[p][ch][t:(t+nn)]*Un_Tbeam_Sun[p][ch][t:(t+nn)]
					Temp_beam_sun[p][ch][t:(t+nn)] = Tbeam_sun
					T_sun[p][ch][t:(t+nn)] = Tbeam_sun*sum_beamarea/solidangle_sun	# in K
					S_sun[p][ch][t:(t+nn)] = (2*kb*Tbeam_sun*sum_beamarea/(c/(freq*1.e6))**2)/1.e-22 	# in SFU
					print 'FLUX: ',np.nanmean(S_sun[p][ch][t:(t+nn)]),' NCCF: ',np.nanmean(ncross[p][ch][t:(t+nn)]),' Sky+rec+grd: ',np.nanmean(Tsky_integrated[p][ch][t:(t+nn)]+Trec[ch]+Tpickup[ch]),'Sky contribution: ',np.nanmean(Tsky_integrated[p][ch][t:(t+nn)]),'Baseline: ',np.mean(T_baseline[p][ch][t:(t+nn)])
					#print Tbeam_sun
					#plt.plot(Un_Tbeam_Sun[p][ch][t:(t+62)])
					#plt.plot(Tbeam_sun)
					#plt.show()
				#plt.plot(Un_Tbeam_Sun[p][ch])
				#plt.plot(Temp_beam_sun[p][ch])
				#plt.show()
			#sys.exit()
			dynamic_spectrum = ds(T_sun,S_sun,Temp_beam_sun,Un_Tbeam_Sun,corr_factor,fringe_factor,T_baseline,Tsky_integrated)
			os.chdir(DS_DIR)
			pickle.dump([central_freq,chan_list,auto_t1,auto_t2,cross,ncross,phase_ncross,u,v,w,azi_pointing,ele_pointing,ph_az,ph_el,start_time,mid_time,end_time,[0,0,corr_factor,S_sun,T_sun,Un_Tbeam_Sun,Temp_beam_sun,fringe_factor,T_baseline,Tsky_integrated]],open('flux_V1_'+DS.split('.')[0] + '.p','wb'))
			#pickle.dump([Tsky_integrated,T_baseline,Tbeam_sun,Temp_beam_sun,T_sun,S_sun],open('check_flux.dat.p','wb'))
		    	os.chdir(WORKING_DIR)
		#time2=time.time()
		#print (time2-time1)
		###############################################################################################################################################

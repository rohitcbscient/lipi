from astropy.io import ascii
import os
import numpy as np
import pickle
import thread
from scipy.interpolate import interp1d
import subprocess
from scipy.interpolate import interp2d
###############################################################################################################################################
def create_read_beam_old(azimuth_bins,elevation_bins,freq,azimuth,elevation):		# Reading from the beam files created using Alan's code	
	#Creating .beam file 
	os.system('pwd')
	filename = 'azel'+str(freq)+'_az'+str(azimuth)+'_el'+str(elevation)+'.beam'
	if (os.path.isfile(filename)== False):
		sp = subprocess.Popen(' ./azelbeam azel'+str(freq)+'_final.txt -az '+str('%0.2f'% (azimuth))+' -el '+str('%0.2f'% (elevation))+' -outf '+str(freq),shell=True)
		sp.wait()
		os.system(' mv '+'azel'+str(freq)+'.beam '+str(filename))
		os.system(' mv '+'azel.pos '+str(filename)+'.ps')
		os.system('print '+str(filename))
	'''
	creates 360,90 array with 360 rows
	'''
	ele = elevation_bins + 1
	print filename
	beam = ascii.read(filename)
	amp_x0=np.array(beam['col5']).reshape(azimuth_bins,ele)[:,1:ele]
	ph_x0=np.array(beam['col6']).reshape(azimuth_bins,ele)[:,1:ele]
	amp_x1=np.array(beam['col7']).reshape(azimuth_bins,ele)[:,1:ele]
	ph_x1=np.array(beam['col8']).reshape(azimuth_bins,ele)[:,1:ele]
	amp_y0=np.array(beam['col9']).reshape(azimuth_bins,ele)[:,1:ele]
	ph_y0=np.array(beam['col10']).reshape(azimuth_bins,ele)[:,1:ele]
	amp_y1=np.array(beam['col11']).reshape(azimuth_bins,ele)[:,1:ele]
	ph_y1=np.array(beam['col12']).reshape(azimuth_bins,ele)[:,1:ele]
	beamdb=np.array(beam['col13']).reshape(azimuth_bins,ele)[:,1:ele]	
	return amp_x0,ph_x0,amp_x1,ph_x1,amp_y0,ph_y0,amp_y1,ph_y1,beamdb


def create_read_beam(azimuth_bins,elevation_bins,freq,azimuth,elevation):		# Reading from the beam files created using Alan's code	
	#Creating .beam file 
	os.system('pwd')
	filename = 'azel'+str(freq)+'_az'+str(azimuth)+'_el'+str(elevation)+'.beam'
	if (os.path.isfile(filename)== False):
		sp = subprocess.Popen(' ./azelbeam azel'+str(freq)+'_final.txt -az '+str('%0.2f'% (azimuth))+' -el '+str('%0.2f'% (elevation))+' -outf '+str(freq),shell=True)
		sp.wait()
		os.system(' mv '+'azel'+str(freq)+'.beam '+str(filename))
		os.system(' mv '+'azel.pos '+str(filename)+'.ps')
		os.system('print '+str(filename))
	'''
	creates 360,90 array with 360 rows
	'''
	ele = elevation_bins + 1
	print filename
	beam = np.loadtxt(filename)
	amp_x0=beam[:,2].reshape(azimuth_bins,ele)[:,1:ele]
	ph_x0=beam[:,3].reshape(azimuth_bins,ele)[:,1:ele]
	amp_x1=beam[:,4].reshape(azimuth_bins,ele)[:,1:ele]
	ph_x1=beam[:,5].reshape(azimuth_bins,ele)[:,1:ele]
	amp_y0=beam[:,6].reshape(azimuth_bins,ele)[:,1:ele]
	ph_y0=beam[:,7].reshape(azimuth_bins,ele)[:,1:ele]
	amp_y1=beam[:,8].reshape(azimuth_bins,ele)[:,1:ele]
	ph_y1=beam[:,9].reshape(azimuth_bins,ele)[:,1:ele]
	beamdb=beam[:,10].reshape(azimuth_bins,ele)[:,1:ele]
	return amp_x0,ph_x0,amp_x1,ph_x1,amp_y0,ph_y0,amp_y1,ph_y1,beamdb

def interpolate(yi,yi1,freq,freqi,freqi1,azimuth_bins,elevation_bins):
	ynew=np.zeros((azimuth_bins,elevation_bins))
	for i in range(azimuth_bins):
		for j in range(elevation_bins):
			alpha=(yi1[i,j]-yi[i,j])/(freqi1-freqi)
			ynew[i,j]=yi[i,j]+(freq-freqi)*alpha
	return ynew

def beam(azi_pointing,ele_pointing,freq,azimuth_bins,elevation_bins,BEAM_DIR,WORKING_DIR):
	flist=100+200*np.arange(21)/20
	a=freq-flist
	fli=flist[np.where((a>=0)&(a<10))][0]
	if (freq>300):
		fli=300
		fli1=310
	else:
		fli1=flist[np.where((a>=0)&(a<10))[0][0]+1]
	print fli,fli1,freq
	os.chdir(BEAM_DIR)
	f1=create_read_beam(azimuth_bins,elevation_bins,fli,azi_pointing,ele_pointing)
	f2=create_read_beam(azimuth_bins,elevation_bins,fli1,azi_pointing,ele_pointing)		# Primary beam paramters at nearest frequencies fli and fli1
	os.chdir(WORKING_DIR)
	amp_x0n=interpolate(f1[0],f2[0],freq,fli,fli1,azimuth_bins,elevation_bins)
	ph_x0n=interpolate(f1[1],f2[1],freq,fli,fli1,azimuth_bins,elevation_bins)
	amp_x1n=interpolate(f1[2],f2[2],freq,fli,fli1,azimuth_bins,elevation_bins)
	ph_x1n=interpolate(f1[3],f2[3],freq,fli,fli1,azimuth_bins,elevation_bins)
	amp_y0n=interpolate(f1[4],f2[4],freq,fli,fli1,azimuth_bins,elevation_bins)
	ph_y0n=interpolate(f1[5],f2[5],freq,fli,fli1,azimuth_bins,elevation_bins)
	amp_y1n=interpolate(f1[6],f2[6],freq,fli,fli1,azimuth_bins,elevation_bins)
	ph_y1n=interpolate(f1[7],f2[7],freq,fli,fli1,azimuth_bins,elevation_bins)
	beamdb=interpolate(f1[8],f2[8],freq,fli,fli1,azimuth_bins,elevation_bins)
	primary_beam = [amp_x0n,ph_x0n,amp_x1n,ph_x1n,amp_y0n,ph_y0n,amp_y1n,ph_y1n,beamdb]  #Computing primary beam paramters at required frequency by interpolation
	return primary_beam	

def normalise(primary_beam,no_pol,azimuth_bins,elevation_bins):
	nprimarybeam = np.zeros((no_pol,azimuth_bins,elevation_bins))
	nprimarybeam[0] = (primary_beam[0]*primary_beam[0]+primary_beam[2]*primary_beam[2])/max((primary_beam[0]*primary_beam[0]+primary_beam[2]*primary_beam[2]).flatten())  	# XX polarisation
	nprimarybeam[1] = (primary_beam[0]*primary_beam[0]+primary_beam[6]*primary_beam[6])/max((primary_beam[0]*primary_beam[0]+primary_beam[6]*primary_beam[6]).flatten())  	# XY polarisation
        nprimarybeam[2] = (primary_beam[4]*primary_beam[4]+primary_beam[2]*primary_beam[2])/max((primary_beam[4]*primary_beam[4]+primary_beam[2]*primary_beam[2]).flatten())  	# YX polarisation
        nprimarybeam[3] = (primary_beam[4]*primary_beam[4]+primary_beam[6]*primary_beam[6])/max((primary_beam[4]*primary_beam[4]+primary_beam[6]*primary_beam[6]).flatten())  # YY polarisation
	return nprimarybeam

def beamarea(nprim_beam,solidangle):
	beam_area = nprim_beam*solidangle
	sum_beamarea = np.sum(beam_area)
	return beam_area,sum_beamarea

def factor(norm_prim_beam,az_sun,alt_sun,azimuth_bins,elevation_bins):
    f = interp2d(np.arange(elevation_bins),np.arange(azimuth_bins), norm_prim_beam, kind='cubic')
    factor=1.0/f(alt_sun,az_sun)[0]
    return factor

###############################################################################################################################################

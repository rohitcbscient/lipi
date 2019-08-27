import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import ephem
################################################################################################################################################
def radec2azel(ra,dec,array_lon,array_lat,starttime):		# Coordinate conversion from (RA,Dec) to (Alt, Az)
	obs = ephem.Observer()
	obs.lon = array_lon
	obs.lat = array_lat
	obs.date =  starttime
	eq = ephem.Equatorial(ra*np.pi/180, dec*np.pi/180, epoch='2000')
	body = ephem.FixedBody()
	body._ra = eq.ra
	body._dec = eq.dec
	body.compute(obs)
	return body.alt*180/np.pi,body.az*180/np.pi

def spherical2cartesian(ro,th,phi):		# Conversion from spherical polar coordinates to cartesian ccoordinates
	x=ro*np.sin(th)*np.cos(phi)
	y=ro*np.sin(th)*np.sin(phi)
	z=ro*np.cos(th)
	return x,y,z

def cartesian2spherical(x,y,z):			# Conversion from cartesian to spherical polar coordinates
	ro=np.sqrt(x*x+y*y+z*z)
	th=np.arctan(y/x)
	phi=np.arccos(z/np.sqrt(x*x+y*y+z*z))
	return ro,th,phi

def vectorrotationcartisean_x(x,y,z,angle):	# Vector rotation by angle theta in clockwise direction about x-axis
	xo=x
	yo=y*np.cos(angle)-z*np.sin(angle)
	zo=y*np.sin(angle)+z*np.cos(angle)
	return xo,yo,zo

def vectorrotationcartisean_y(x,y,z,angle):	# Vector rotation by angle theta in clockwise direction about y-axis
	xo=x*np.cos(angle)+z*np.sin(angle)
	yo=y
	zo=-x*np.sin(angle)+z*np.cos(angle)
	return xo,yo,zo

def vectorrotationcartisean_z(x,y,z,angle):	# Vector rotation by angle theta in clockwise direction about z-axis
	xo=x*np.cos(angle)-y*np.sin(angle)
	yo=x*np.sin(angle)+y*np.cos(angle)
	zo=z
	return xo,yo,zo
	
def fringe_fraction(um,vm,wm,freq,s,az,el,az_sun,el_sun,dl,dm):		# Fraction captured by baseline
	'''
	fraction(u,v,w,freq,s,dx,dy) captured
	s: Half size of Radio Sun
	dx:  (along the azimuth direction) in deg
	dy: (along the elevation direction) in deg
	'''
	print 'dl: ',dl,'dm: ',dm ,'Sun el-az:',az_sun,el_sun
	el_sun=90-el_sun
	el=90-el
	arcmin22=s*np.pi/(60*180.)
	lnew,mnew= np.meshgrid(-1*arcmin22*(-1+np.arange(200)/100.),-1*arcmin22*(-1+np.arange(200)/100.))
	cosine=np.zeros((200,200))
	sine=np.zeros((200,200))
	one_circle=np.zeros((200,200))
	l=0
	for i in range(200):
		for j in range(200):
			if((lnew[i,j]*lnew[i,j]+mnew[i,j]*mnew[i,j])<arcmin22*arcmin22):
				lnn=lnew[i,j]+dl*(np.pi/180)
				mnn=mnew[i,j]+dm*(np.pi/180)
				nnn=np.sqrt(1-lnn*lnn-mnn*mnn)
				cosine[i,j]=np.cos(2*np.pi*(um*lnn+vm*mnn+wm*(nnn-1)))
				sine[i,j]=np.sin(2*np.pi*(um*lnn+vm*mnn+wm*(nnn-1)))
				one_circle[i,j]=1
				l=l+1
	abs_sun=np.sqrt(np.sum(cosine)*np.sum(cosine)+np.sum(sine)*np.sum(sine))
	frac_=(abs_sun/l) 
	print 'FRINGE CONTRIBUTION (in percentage):: ',frac_*100.0
	return frac_
	
def lmn(az,el,azimuth_bins,elevation_bins):			# Computing l,m,n
	za=90-el
	l=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	m=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	n=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	cart_x=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	cart_y=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	cart_z=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	ncart_x=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	ncart_y=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	ncart_z=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	nx=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	ny=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	nz=np.zeros((azimuth_bins+1,2*elevation_bins+1))
	for i in range(azimuth_bins+1):
		for j in range(2*elevation_bins+1):
			cart_x[i,j]=spherical2cartesian(1,j*np.pi/180,np.pi*i/180)[0]#np.sin(j*np.pi/180)*np.cos(i*np.pi/180)
			cart_y[i,j]=spherical2cartesian(1,j*np.pi/180,i*np.pi/180)[1]#np.sin(j*np.pi/180)*np.sin(i*np.pi/180)
			cart_z[i,j]=spherical2cartesian(1,j*np.pi/180,i*np.pi/180)[2]#np.cos(j*np.pi/180)
			nx[i,j],ny[i,j],nz[i,j]=vectorrotationcartisean_y(cart_x[i,j],cart_y[i,j],cart_z[i,j],-1*za*(np.pi/180)) # ROTATION IN ALTITUDE (ROTATION FROM ZENITH TO GIVEN ALT-AZ)
			ncart_x[i,j],ncart_y[i,j],ncart_z[i,j]=vectorrotationcartisean_x(nx[i,j],ny[i,j],nz[i,j],az*(np.pi/180)) # ROTATION IN AZIMUTH 
			m[i,j]=ncart_y[i,j]
			l[i,j]=ncart_x[i,j]
			n[i,j]=ncart_z[i,j]
	return	l,m,n	

def fringe_frac_baseline(um,vm,wm,freq,s,l,m,az,el,az_sun,el_sun):
	print str(freq)+' MHz'
	lamd=300.0/freq	
	u,v,w=um/lamd,vm/lamd,wm/lamd
	dl=(l[int(az),int(el)]+l[int(az)+1,int(el)]+l[int(az),int(el)+1]+l[int(az)+1,int(el)+1])/4
	dm=(m[int(az),int(el)]+m[int(az)+1,int(el)]+m[int(az),int(el)+1]+m[int(az)+1,int(el)+1])/4
	frac_=fringe.fraction(u,v,w,freq,s,m,az,el,az_sun,el_sun,dl,dm)[0]
	return (frac_,u,v,w)

def fraction(azi_pointing,ele_pointing,az_sun,el_sun,freq,u,v,w,solar_size,azimuth_bins,elevation_bins,t):
	l,m,n = lmn(azi_pointing,ele_pointing,azimuth_bins,elevation_bins)
	lambd = 300.0/freq		# Frequency in MHz		# using mean value of u,v,w per minute
	u_mean = np.mean(u)
	v_mean = np.mean(v)
	w_mean = np.mean(w)
	u_mean,v_mean,w_mean=u_mean/lambd,v_mean/lambd,w_mean/lambd
	dl=(l[int(azi_pointing),int(ele_pointing)]+l[int(azi_pointing)+1,int(ele_pointing)]+l[int(azi_pointing),int(ele_pointing)+1]+l[int(azi_pointing)+1,int(ele_pointing)+1])/4
	dm=(m[int(azi_pointing),int(ele_pointing)]+m[int(azi_pointing)+1,int(ele_pointing)]+m[int(azi_pointing),int(ele_pointing)+1]+m[int(azi_pointing)+1,int(ele_pointing)+1])/4
	print 'Time index: ',t,'dm: ',dm
	frac_ = fringe_fraction(u_mean,v_mean,w_mean,freq,solar_size,azi_pointing,ele_pointing,az_sun,el_sun,dl,dm)
	fringe_cos=np.cos(2*np.pi*(l*u_mean+m*v_mean+w_mean*(n-1)))				# FRINGE COSINE
	fringe_sin=np.sin(2*np.pi*(l*u_mean+m*v_mean+w_mean*(n-1)))				# FRINGE SINE	
	return (1/frac_,fringe_cos,fringe_sin)
	
################################################################################################################################################
'''
phase_el,phase_az = radec2azel(ra_sun,dec_sun,array_lon,array_lat,starttime)
l,m,n=fringe.lmn(phase_az,phase_el,azimuth_bins,elevation_bins)        			# COMPUTING THE L,M,N FOR GIVEN AZ-EL DIRECTION
fraction_baseline[nb][p],u,v,w=fringe_frac_baseline(um,vm,wm,freq,solar_size,l,m,phase_az,phase_el,az_sun,el_sun)
		fringe_factor=(100.0/fraction_baseline[nb][p])	
		fringe_cos=np.cos(2*np.pi*(l*u+m*v+w*(n-1)))				# FRINGE COSINE
		fringe_sin=np.sin(2*np.pi*(l*u+m*v+w*(n-1)))				# FRINGE SINE		
fcos=fringe_cos[0:360,0:90][:,::-1]						
'''

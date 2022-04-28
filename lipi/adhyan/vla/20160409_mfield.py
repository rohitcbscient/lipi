import numpy as np
from astropy.io import fits
from surya.utils import Bextrap
from astropy import units as u
import matplotlib.pyplot as plt
import pickle
from astropy.coordinates import SkyCoord

file_name_fr="/media/rohit/VLA/paraview/radio_loop_new.csv"
file_name_ls="/media/rohit/VLA/paraview/large_scale.csv"
file_name_euv="/media/rohit/VLA/paraview/EUV_loop2.csv"
file_name_radio='/media/rohit/VLA/paraview/EUV_loop_radio.csv'
x,y,z,bx,by,bz=Bextrap.get_fieldlines(file_name_fr)
ff='/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav'
bsmap,bcarrmap=Bextrap.get_gxs_sav2hpp(ff,'2016-04-09T18:34:17.00')
x,y,z,bx,by,bz,b_hp,b_hp_pix,b_carr,b_carr_pix=Bextrap.transform_fieldlines(x,y,z,bx,by,bz,'2016/04/09T18:45:00',bsmap[0].wcs,bcarrmap[0])
babs=np.sqrt(bx**2 + by**2 + bz**2)
####
s=2
B1000=1000./2.8/s;B1500=1500./2.8/s # in G
idx1000=np.where((babs<B1000+1) & (babs>B1000-1))
idx1500=np.where((babs<B1500+1) & (babs>B1500-1))
z_1000=z[idx1000].max()*1.4;z_1500=z[idx1500].max()*1.4

####
fermi_energy=15 # keV
omegapbygyro=0.3
R=z_1000-z_1500
v=np.sqrt(fermi_energy/(0.5*9.1e-31)*(1.e3*1.6e-19))/1.e6 # Mm/s
B_low=200;B_top=40;pitch_angle=np.logspace(0,2,10000);htop=80; s_cl=z_1000 # Mm/s
loop_length = np.pi*htop; B_scale_height=loop_length*0.5/np.e # Mm
alpha0=np.arcsin(np.sqrt(B_top/B_low))*180/np.pi
alphac=np.arcsin(np.sqrt(B_top/B1000))*180/np.pi
mirror_point = B_top-B_scale_height/np.tan(pitch_angle*3.14159/180);mirror_point[mirror_point<0] = 0
mirror_ratio=B_low/B_top
bounce_time = 2*np.pi*B_scale_height/(v*np.sin(pitch_angle*np.pi/180))
tcc=bounce_time/np.pi*np.arccos((htop-s_cl)/(htop-mirror_point))
epsilonD = 1/((np.pi*loop_length*1000)/(2*726*18)) # 18" radius
alpha_c=np.arcsin(np.sqrt(B_top/B1000))*180/np.pi
lambda_ei=1.e4*5.e6**2/5.e8/1.e8 # in Mm/sec
t_collision=lambda_ei/v # in sec
w=tcc/t_collision
t_loss=0.5*bounce_time
t_loss = t_collision*(0.5*np.pi/np.arccos((htop-s_cl)/(htop-mirror_point)))


taup=5.1;tau_growth=6.3e-4
tau_diff_eff=(taup/(2*np.pi))**2 /tau_growth
tau_diff=tau_diff_eff*epsilonD
I_ratio=np.sqrt(tau_growth/tau_diff_eff) # Ampltitude of perturbations in the number distribution to density


plt.plot(pitch_angle,t_loss,color='b',label='$\\tau_{loss}$')
plt.axvline(x=alpha_c,color='r',label='$\\alpha_C$')
plt.axhline(y=t_collision,color='green',linestyle='-',label='$\\tau_{coll}$')
plt.ylabel('Time (sec)');plt.xlabel('$\\alpha_0$ (degrees)');plt.legend()#;plt.ylim(0.01,1)
plt.yscale('log');plt.show()

cm = plt.cm.get_cmap('RdYlBu');sc = plt.scatter(x*1.4-200,z*1.4, c=babs, vmin=0, vmax=400, s=35, cmap=cm);plt.colorbar(sc,label='|B| (G)')
plt.plot(np.linspace(20,55,10),np.ones(10)*z_1000,'-',color='k',label='Gyrofrequency layer (s=2)')
plt.plot(np.linspace(20,55,10),np.ones(10)*z_1500,'-',color='k')
plt.plot(np.linspace(20,55,10),np.ones(10)*mirror_point.min(),'--',color='b',label = 'Mirror point')
#plt.plot(np.linspace(220,255,10),np.ones(10)*mirror_point.max(),'--',color='b')
plt.xlabel('Distance (Mm)');plt.ylabel('Height (Mm)');plt.legend(loc=2)
plt.show()

#####
xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r.p','rb'),encoding = 'latin1')
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs.p','rb'),encoding = 'latin1')
qsmaxTbr=np.array(qsmaxTbr)
x0,y0,sigx0,sigy0,rot0,tb0,x1,y1,sigx1,sigy1,rot1,tb1=np.load('/media/rohit/VLA/20160409/blob/all_params.npz.npy')
corr_t=np.load('/media/rohit/VLA/20160409/correlation_t.npz.npy')
xlaia=-948;xraia=-648;ylaia=70;yraia=370


plot_fields=1
if(plot_fields):
    plt.close()
    fig = plt.figure(figsize=(16, 8))
    ax0 = fig.add_subplot(121,projection=bsmap[0]);ax1 = fig.add_subplot(122,projection=bcarrmap[0])
    p0=bsmap[0].plot(axes=ax0,aspect='auto')
    p1=bcarrmap[0].plot(axes=ax1,aspect='auto')
    hp_lon = b_hp.Tx.value * u.arcsec
    hp_lat = b_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
    ax0.plot_coord(seeds0, color='tab:olive', marker='o', markersize=10)
    #ca_lon = b_proj[0] * u.deg
    #ca_lat = b_proj[1] * u.deg
    seeds1 = SkyCoord(b_carr.lon.ravel(), b_carr.lat.ravel(),frame=bcarrmap[0].coordinate_frame)
    ax1.plot_coord(seeds1, color='tab:olive', marker='o', markersize=10)
    #plt.plot(xcr90[15][1000:1400],ycr90[15][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[10][1000:1400],ycr90[10][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[5][1000:1400],ycr90[5][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[1][1000:1400],ycr90[1][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[12][500:600],ycr90[12][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[10][500:600],ycr90[10][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[5][500:600],ycr90[5][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[1][500:600],ycr90[1][500:600],'o',alpha=0.2)
    #ax0.set_xlim(-860,-680);ax0.set_ylim(180,360)
    plt.show()

hp_lon = b_hp.Tx.value * u.arcsec
hp_lat = b_hp.Ty.value* u.arcsec
fig = plt.figure()
ax = plt.subplot(projection=bsmap[0])
bsmap[0].plot(axes=ax)
seeds = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
ax.plot_coord(seeds, color='white', marker='o', linewidth=0)


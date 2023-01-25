import numpy as np
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from surya.utils import Bextrap
from surya.utils import main as ut
import pickle
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm

h=np.arange(200)*1400/1.e3 # in Mm
h1=np.arange(200)*500/1.e3 # in Mm
sh=np.linspace(10,500,100)
ux=np.linspace(0,50,50)
dem_arr=[0]*100;ne=[0]*100;ne_dem=[0]*100;wpe_dem=[0]*100;omegap_dem=[0]*100
tem=2.6e27

for i in range(100):
        dem_arr[i]=[0]*200;ne[i]=[0]*200
        ne[i]=1.16e17*(np.exp(-1*h1*1000/sh[i]))+9.5e8*np.exp(-1*h1*1000/31.e3)#+4.2e4 *10**(4.32/((695700.+h*1000)/695700.))


#####################333
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
B1000=1000./2.8/s;B1500=1500./2.8/s;B1200 = 1200./2.8/s # in G
idx1000=np.where((babs<B1000+1) & (babs>B1000-1))
idx1500=np.where((babs<B1500+1) & (babs>B1500-1))
idx1200=np.where((babs<B1200+1) & (babs>B1200-1))
z_1000=z[idx1000].max()*1.4;z_1500=z[idx1500].max()*1.4;z_1200=z[idx1200].max()*1.4

####
penergy=np.linspace(1,20,100)
#fermi_energy=1 # keV
pitch_angle=np.logspace(0,np.log10(90),1000)
omegapbygyro=0.3
R=z_1000-z_1500
B_low = 200;
B_top = 40
htop = 80;
s_cl1000 = loop_length * 0.5-z_1000  # Mm/s
s_cl1200 = loop_length * 0.5-z_1200  # Mm/s
s_cl1500 = loop_length * 0.5-z_1500  # Mm/s
loop_length = np.pi * htop;
B_scale_height = loop_length * 0.5 / np.e  # Mm
alpha0 = np.arcsin(np.sqrt(B_top / B_low)) * 180 / np.pi
alphac = np.arcsin(np.sqrt(B_top / B1000)) * 180 / np.pi
mirror_point = B_scale_height / np.tan(pitch_angle * 3.14159 / 180);
mirror_point[mirror_point >loop_length*0.5] = loop_length*0.5
mirror_ratio = B_low / B_top

f,ax=plt.subplots(1,1)
ax.plot(pitch_angle,mirror_point,'o-',label='Mirror Point')
ax.axhline(y=s_cl1000,label='1.0 GHz',color='r')
ax.axhline(y=s_cl1500,label='1.5 GHz',color='g')
ax.axvline(x=70,label='$\\alpha=70^o$',color='k')
ax.axvline(x=85,label='$\\alpha=85^o$',color='brown')
ax.legend()
ax.set_xlabel('Pitch Angle (deg)');ax.set_ylabel('$s$ (Mm)')
plt.show()


ne_1d=ne[28];ne_1p5GHz=ne_1d[75];ne_1d_low=1.5e8
v=300*np.sqrt(1-(1/(0.002*penergy+1)**2))
lambda_turb=lambda_ei*np.e
tcc=[0]*len(penergy);bounce_time=[0]*len(penergy);w=[0]*len(penergy);t_loss=[0]*len(penergy)
slength=[0]*len(penergy);slength_turb=[0]*len(penergy);
i=0
for fermi_energy in penergy:
    #v[i]=np.sqrt(fermi_energy/(0.5*9.1e-31)*(1.e3*1.6e-19))/1.e6 # Mm/s
    T=fermi_energy*1.16e7
    #idx=np.where(mirror_point>20)
    bounce_time[i] = 2*np.pi*B_scale_height/(v[i]*np.sin(pitch_angle*np.pi/180))
    tcc[i]=bounce_time[i]/np.pi*np.arccos((s_cl)/(mirror_point))
    epsilonD = 1/((np.pi*loop_length*1000)/(2*726*18)) # 18" radius
    Tt=2.e6;lambda_ei=5210*Tt**2/(ne_1p5GHz)/1.e8 #(CGS) in Mm ne[28][72] is ne for the 1.5 GHz
    w[i]=tcc[i]/(lambda_ei/v[i])
    t_loss[i] = (lambda_ei/v[i]) * (0.5 * np.pi / np.arccos((s_cl)/(mirror_point)))
    slength[i]=10*(fermi_energy**2/400)*(10**11/ne_1d_low)*(np.cos(pitch_angle*np.pi/180))
    slength_turb[i]=np.sqrt(lambda_turb*1.e8/9.3e-36/ne_1d_low)*penergy
    i=i+1
tcc=np.array(tcc);bounce_time=np.array(bounce_time)
v=np.array(v);w=np.array(w);slength=np.array(slength)
slength[slength>60] = np.nan;slength[slength<20] = np.nan
slength_p=slength*1.0
t_loss=np.array(t_loss)
slength_p[0:9]=np.nan; slength_p[19:]=np.nan
slength_p[np.isfinite(slength_p)]=1.0
slength_turb=np.array(slength_turb)

K=1.e-36

slength_diff=np.sqrt(lambda_ei/K/ne_1d_low)*penergy/1.e8

#--- Bounce Time
f,ax=plt.subplots(1,1)
im=ax.imshow(bounce_time,aspect='auto',origin='lower',cmap='jet',norm=LogNorm(vmin=5, vmax=200))
ax.contour(bounce_time,levels=[40],color='k')
ax.set_xscale('linear');ax.set_xlim([500,980])
tickxidx=[200,500,600,700,800,980];ax.set_xticks([])
ax.set_xticks(tickxidx)
ax.set_xticklabels(np.round(pitch_angle[tickxidx],1))
ax.set_yticks(np.arange(len(penergy))[::10])
ax.set_yticklabels(np.round(penergy[::10],1))
ax.set_ylabel('Electron Energy (keV)');ax.set_xlabel('Pitch Angle (deg)')
f.colorbar(im,label='Bounce Time (sec)')
plt.show()

#--- tcc
f,ax=plt.subplots(1,1)
im=ax.imshow(tcc,aspect='auto',origin='lower',cmap='jet',norm=LogNorm(vmin=0.5, vmax=3))
ax.set_xscale('linear');ax.set_xlim([500,980])
tickxidx=[200,500,600,850,900,950,980];ax.set_xticks([])
ax.set_xticks(tickxidx)
ax.set_xticklabels(np.round(pitch_angle[tickxidx],1))
ax.set_yticks(np.arange(len(penergy))[::10])
ax.set_yticklabels(np.round(penergy[::10],1))
ax.set_ylabel('Electron Energy (keV)');ax.set_xlabel('Pitch Angle (deg)')
f.colorbar(im,label='$t_{cc}$ (sec)')
plt.show()

#--- t_loss
f,ax=plt.subplots(1,1)
im=ax.imshow(t_loss,aspect='auto',origin='lower',cmap='jet',norm=LogNorm(vmin=0.01, vmax=0.5))
ax.set_xscale('linear');ax.set_xlim([500,980])
tickxidx=[200,300,400,450,500,500,600,650,700,750,800];ax.set_xticks([])
ax.set_xticks(tickxidx)
ax.set_xticklabels(np.round(pitch_angle[tickxidx],1))
ax.set_yticks(np.arange(len(penergy))[::10])
ax.set_yticklabels(np.round(penergy[::10],1))
ax.set_ylabel('Electron Energy (keV)');ax.set_xlabel('Pitch Angle (deg)')
f.colorbar(im,label='$t_{loss}$ (sec)')
plt.show()


#--- Stopping Length
f,ax=plt.subplots(1,1)
im=ax.imshow(slength,aspect='auto',origin='lower',cmap='spring',norm=LogNorm(vmin=5, vmax=100))
ax.set_xscale('linear');ax.set_xlim([800,990])
tickxidx=[200,500,600,700,750,850,900,950,980];ax.set_xticks([])
ax.set_xticks(tickxidx)
ax.set_xticklabels(np.round(pitch_angle[tickxidx],1))
ax.set_yticks(np.arange(len(penergy))[::10])
ax.set_yticklabels(np.round(penergy[::10],1))
ax.set_ylabel('Electron Energy (keV)');ax.set_xlabel('Pitch Angle (deg)')
f.colorbar(im,label='$L_{stop}$ (Mm)')
ax.imshow(slength_p,aspect='auto',origin='lower',alpha=0.4,cmap='winter')
#ax.axhline(y=9,color='brown')
ax.axhline(y=18,color='brown')
ax.axhline(y=9,color='red',label='$E_{min}$=2.4 keV',linewidth=2,linestyle='--')
ax.axvline(x=925,color='blue',label='$\\alpha=65^o$')
ax.axvline(x=990,color='brown',label='$\\alpha=86^o$')
ax.axvline(x=975,color='black',label='$\\alpha=81^o$',linewidth=2,linestyle='--')
ax.legend(loc=2)
plt.show()




epsilonD = 1/((np.pi*loop_length*1000)/(2*726*18)) # 18" radius
Tt=2.e6;lambda_ei=5210*Tt**2/(ne[28][72])/1.e8 #(CGS) in Mm ne[28][72] is ne for the 1.5 GHz 
t_collision=lambda_ei/v # in sec
w=tcc/t_collision
t_loss=0.5*bounce_time
t_loss = t_collision*(0.5*np.pi/np.arccos((htop-s_cl)/(htop-mirror_point)))

tcc=(Talpha0*np.arccos(scl/s1))

# In CGS
fraction_flare=A/(4*np.pi*L**2)
nnth=fraction_flare*3.e9
v0=3.e9;S=120;m0=9.1e-28;n0=3.e8
L=1.2e10;Gamma=1.e17;alp0=0.4;masing_length=6.e9
Gamma=(2*np.pi*0.005*1.e9)
facc=(v0/(L*Gamma))
W = facc*n0*m0*v0*v0*alp0**3
A=np.pi*(18*726)**2*1.e10;V=A*masing_length
rad_energy=S*1.e-19*A*1.e9*2.5
a=8;eta1=0.5*(a-3)*np.pi/alp0*np.sin(alp0)*np.sin(alp0)*np.cos(alp0)*np.cos(alp0)*alp0*alp0
print('W (ergs/cm^3):',W,'facc:',facc,'W (ergs):',W*V,'Rad Ener:',rad_energy)
print('Number of masing cells:',A/(rad_energy/W/masing_length),'Gamma:',Gamma/1.e5)
print('Ratio:',W*V/rad_energy)
del_alph=(rad_energy/V/facc/n0/m0/v0/v0)**(1./3)

taup=5.1;tau_growth=6.3e-4
tau_diff_eff=(taup/(2*np.pi))**2 /tau_growth
tau_diff=tau_diff_eff*epsilonD
I_ratio=np.sqrt(tau_growth/tau_diff_eff) # Ampltitude of perturbations in the number distribution to density

tof0=100 # sec
alpha_trap0 = np.rad2deg(np.arcsin(2.e9/tof0/v0))

plt.plot(pitch_angle[idx],mirror_point[idx],'o-')
plt.axvline(x=alpha0,color='k',label='$\\alpha_0$')
plt.axhline(y=z_1500,color='r',label='1.5 GHz')
plt.axhline(y=z_1000,color='b',label='1.0 GHz')
plt.legend();plt.xlabel('$\\alpha$ (degrees)');plt.ylabel('Coronal Height (Mm)')
plt.show()

#CGS
me=9.1e-28;kb=1.23e-16;v0=3.e9;n0=1.e9;omega=1.e9;alpha0=85;Tb=2.e7;c=3.e10
incoh_Tb=me/kb*v0*v0/2/np.pi
del_alpha=alpha0*np.sqrt(incoh_Tb*n0*2*np.pi*c*c/omega/v0/Tb)

plt.plot(pitch_angle,t_loss,color='b',label='$\\tau_{loss}$')
plt.axvline(x=alpha0,color='r',linestyle='-',label='$\\alpha_0$')
plt.axvline(x=alphac,color='r',linestyle='--',label='$\\alpha_C$')
plt.axhline(y=t_collision,color='green',linestyle='-',label='$\\tau_{coll}$')
plt.plot(pitch_angle,tcc,color='black',linestyle='-',label='$t_{cc}$')
plt.ylabel('Time (sec)');plt.xlabel('$\\alpha_0$ (degrees)');plt.legend()#;plt.ylim(0.01,1)
plt.yscale('log');plt.show()

cm = plt.cm.get_cmap('nipy_spectral');sc = plt.scatter(x*1.4-200,z*1.4, c=babs, vmin=0, vmax=500, s=35, cmap=cm,alpha=0.2)
cb=plt.colorbar(sc,label='|B| (G)');cb.set_alpha(1);cb.draw_all()
#idx1=np.where((babs<201)&(babs>199));plt.scatter(x[idx1]*1.4-200,z[idx1]*1.4, c='k')
plt.plot(np.linspace(20,55,10),np.ones(10)*z_1000,'-',color='r',label='1 GHz Gyrofrequency layer (s=2)')
plt.plot(np.linspace(20,55,10),np.ones(10)*z_1200,'-',color='brown',label='1.2 GHz Gyrofrequency layer (s=2)')
plt.plot(np.linspace(20,55,10),np.ones(10)*z_1500,'-',color='orange',label='1.5 GHz Gyrofrequency layer (s=2)')
plt.plot(np.linspace(20,55,10),np.ones(10)*35,'--',color='k',label = 'Mirroring layers')
plt.plot(np.linspace(-95,-50,10),np.ones(10)*20,'--',color='k')
plt.xlabel('Distance (Mm)');plt.ylabel('Height (Mm)');plt.legend(loc=2);plt.ylim(0,200)
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


#####

penergy=np.linspace(1,20,100)
pitch_angle_array=np.linspace(70,89,100)
Mm2cm=1.e8
volume=np.pi*(18/2.)**2*60
stopping_depth=0.92e17*penergy**2
ne=6.e8
stopping_length=2.e17*penergy**2/ne/1.e8
#stopping_length=stopping_depth*volume*Mm2cm**3
ve=300*np.sqrt(1-(1/(0.002*penergy+1)**2))
travel_time=np.zeros((len(ve),len(pitch_angle_array)))
for i in range(len(ve)):
    for j in range(len(pitch_angle_array)):
        travel_time[i][j]=115./(ve[i]*np.cos(pitch_angle_array[j]*np.pi/180))

energy_stopping_length=np.sqrt(60/2.e17*ne*1.e8) # stopping width = 60 Mm
energy_stopping_idx=ut.find_nearest(penergy,energy_stopping_length)[0]

ne_low=1.e8
energy_stopping_length_low=np.sqrt(115/2.e17*ne_low*1.e8) # stopping width = 60 Mm
energy_stopping_idx_low=ut.find_nearest(penergy,energy_stopping_length_low)[0]
travel_time1=1.0*travel_time;travel_time1[travel_time1<40] = np.nan;travel_time1[travel_time1>50] = np.nan
travel_time1[0:energy_stopping_idx_low]=np.nan; travel_time1[energy_stopping_idx:]=np.nan
travel_time1[np.isfinite(travel_time1)]=1.0

f,ax=plt.subplots(1,1)
im=ax.imshow(np.log10(travel_time),origin='lower',aspect='auto',cmap='coolwarm')
ax.contour(np.log10(travel_time),[np.log10(40)],color='k',linewidth=4)
ax.contour(np.log10(travel_time),[np.log10(50)],color='k',linewidth=4)
ax.imshow(travel_time1,origin='lower',alpha=0.5,cmap='Wistia')
ax.axhline(y=energy_stopping_idx_low,color='brown')
ax.axhline(y=energy_stopping_idx,color='brown')
f.colorbar(im,label='$log_{10}$(Travel Time) (s)')
#plt.plot([70.5],[17],'o',markersize=10,color='green')
ax.set_xticks(np.arange(len(pitch_angle_array))[::10]);ax.set_xticklabels(np.round(pitch_angle_array,1)[::10])
ax.set_yticks(np.arange(len(penergy))[::10]);ax.set_yticklabels(np.round(penergy,1)[::10])
ax.set_ylabel('Energy (keV)');ax.set_xlabel('Pitch Angle (deg)')
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(penergy,stopping_length,'o-')
ax.axvline(x=energy_stopping_length)
plt.show()

Talpha0=2*np.pi*30/(ve*np.sin(28*np.pi/180))
s1=90;scl=80
tcc=(Talpha0*np.arccos(scl/s1))
t_collision_arr=lambda_ei/ve
w=tcc/t_collision_arr
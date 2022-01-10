import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import sys
import glob
from astropy.io import fits
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
import pickle
from sunpy import sun
from dateutil import parser
from reproject import reproject_interp
from sunpy.coordinates import frames
from astropy.wcs import WCS
from astropy.wcs import utils
import csv
from sunpy.coordinates import get_body_heliographic_stonyhurst
from surya.utils import model
from sunpy.map import make_fitswcs_header, Map
from scipy.io import readsav
import os
rename=1
if(rename):
    lst=glob.glob('hmi.m_45s.*.fits')
    for l in lst:
        out=l.split('.')[0]+'_'+l.split('.')[2]+l.split('.')[3]+l.split('.')[4].split('_')[0]+'_'+l.split('.')[4].split('_')[1]+l.split('.')[4].split('_')[2]+l.split('.')[4].split('_')[3]+'_magnetogram.fits'
        os.system('mv '+l+' '+out)
        

hmifile='/sdata/20140914_hmi/hmi.m_45s.2014.09.14_01_55_30_TAI.magnetogram.fits'
hmifile='/sdata/20140914_hmi/hmi.m_45s.2014.09.14_02_45_45_TAI.magnetogram.fits'
hmimap=Map(hmifile)
hmid=hmimap.data[::-1,::-1]
hmid[np.where(hmid<-5000)]=0
hmimap=Map(hmid,hmimap.meta)

#expolB_=readsav('/sdata/20140914_hmi/20140914_magnetic_extrapolation_1arcsec.sav');expolB=expolB_['box']
#expolB_1=readsav('/sdata/20140914_hmi/20140914_magnetic_extrapolation_024450.sav');expolB1=expolB_1['box']
listB=sorted(glob.glob('/sdata/20140914_hmi/hmi/*cube.sav'))
for i in range(len(listB)):
    ii="%04d" % i
    expolB_=readsav(listB[i]);expolB=expolB_['bout']
    bx,by,bz=expolB[0],expolB[1],expolB[2]
    babs=np.sqrt(bx*bx+by*by+bz*bz)
    # Generate the grid
    dim=bx.shape
    zz,yy,xx=np.mgrid[0:dim[0],0:dim[1],0:dim[2]]
    pts = np.empty(bx.shape + (3,), dtype=np.int32)
    pts[..., 0] = xx
    pts[..., 1] = yy
    pts[..., 2] = zz
    vectors = np.empty(bx.shape + (3,), dtype=np.float64)
    vectors[..., 0] = bx
    vectors[..., 1] = by
    vectors[..., 2] = bz
    # We reorder the points and vectors so this is as per VTK's
    # requirement of x first, y next and z last.
    pts = pts.transpose(2, 1, 0, 3).copy()
    pts.shape = pts.size // 3, 3
    vectors = vectors.transpose(2, 1, 0, 3).copy()
    vectors.shape = vectors.size // 3, 3
    sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)
    sg.point_data.vectors = vectors
    sg.point_data.vectors.name = 'Magnetic Field'
    absv=np.sqrt(vectors[:,0]**2+vectors[:,1]**2+vectors[:,2]**2)
    sg.point_data.scalars = absv
    sg.point_data.scalars.name = 'Magnitude'
    write_data(sg, listB[i]+'_'+str(ii)+'.vtk')


sys.exit()

k=0 # Layer
index=expolB['index'][0]
crval1=index['crval1']
crval2=index['crval2']
crpix1=index['crpix1']
crpix2=index['crpix2']
ctype1=index['ctype1']
ctype2=index['ctype2']
cdelt1=index['cdelt1']
cdelt2=index['cdelt2']
cunit1=index['cunit1']
cunit2=index['cunit2']
hdu = fits.PrimaryHDU(babs[k])
list_all=list(index.dtype.names)
list_all.remove('COMMENT')
list_all.remove('HISTORY')
list_all.remove('BITPIX')
list_all.remove('NAXIS')
list_all.remove('DATE_D$OBS')
index['WCSNAME'],index['CTYPE1'],index['CUNIT1'],index['CTYPE2'],index['CUNIT2']=['Carrington-Heliographic'],['CRLN-CEA'],['deg'],['CRLT-CEA'],['deg']
index['DATE_OBS']=['2014-09-14T01:55:30']
ii=0
for idx in list_all:
    #print idx
    #hdu.header.append((idx,index[list_all[ii]][0],[]))
    hdu.header.update({str(idx):index[list_all[ii]][0]})
    ii=ii+1

#source_height=1400;hdu.header.update({'NAXIS':3});hdu.header.append(('CRPIX3',0));hdu.header.append(('CRVAL3',695700.0));hdu.header.append(('CTYPE3','HECH'));hdu.header.append(('CUNIT3','km'));hdu.header.append(('CDELT1',1400))
#hdu.header.append(('RSUN_REF',source_height))
hdu.data=babs[0]
hhdu=hdu.header
hdul = fits.HDUList([hdu])
mymap=Map(babs[k],hhdu)
hp_coord=mymap.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
hp_hcc=mymap.reference_coordinate.transform_to(frames.Heliocentric(observer="earth"))

sys.exit()
############################

out_shape = (300, 300)
out_header = sunpy.map.make_fitswcs_header(out_shape,hp_coord)
out_wcs = WCS(out_header)
earth = get_body_heliographic_stonyhurst('earth', mymap.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(mymap, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))
outmap.plot_settings = mymap.plot_settings
aiamap=Map('/media/rohit/VLA/20160409_EUV/171/ssw_cutout_20160409_184434_AIA_171_.fts')

#file_name_fr="/media/rohit/VLA/paraview/field_radio.csv"
file_name_fr="/media/rohit/VLA/paraview/radio_loop_new.csv"
file_name_ls="/media/rohit/VLA/paraview/large_scale.csv"
file_name_euv="/media/rohit/VLA/paraview/EUV_loop2.csv"
file_name_radio='/media/rohit/VLA/paraview/EUV_loop_radio.csv'

def get_fieldlines(file_name,out_wcs):
    file_ = open(file_name)
    with open(file_name) as f:
        lines=f.readlines()
    numline = len(lines)-1
    bx=[0]*numline;by=[0]*numline;bz=[0]*numline
    x=[0]*numline;y=[0]*numline;z=[0]*numline
    vtkidx=[0]*numline;i=0
    for row in lines[1:]:
        row=row.split(',')
        vtkidx[i]=float(row[1])
        bx[i]=float(row[1]);by[i]=float(row[2]);bz[i]=float(row[3])
        x[i]=float(row[5]);y[i]=float(row[6]);z[i]=float(row[7])
        i=i+1
    vtkidx=np.array(vtkidx)
    x=np.array(x);y=np.array(y);z=np.array(z)
    bx=np.array(bx);by=np.array(by);bz=np.array(bz)
    #with open(file_name, 'r') as csvfile:
    #    #reader = csv.reader(csvfile, skipinitialspace=True)
    #    reader = csv.reader(csvfile)
    #    print(reader.line_num)
    #    bx=[0]*numline;by=[0]*numline;bz=[0]*numline
    #    x=[0]*numline;y=[0]*numline;z=[0]*numline
    #    vtkidx=[0]*numline
    #    i=0
    #    for row in reader:
    #        if(i!=0):
    #            vtkidx[i]=float(row[1])
    #            bx[i]=float(row[3]);by[i]=float(row[4]);bz[i]=float(row[5])
    #            x[i]=float(row[10]);y[i]=float(row[11]);z[i]=float(row[12])
    #        i=i+1
    #    vtkidx=np.array(vtkidx)
    #    x=np.array(x);y=np.array(y);z=np.array(z)
    #    bx=np.array(bx);by=np.array(by);bz=np.array(bz)
    #xkm=((x[1:]-hhdu['CRPIX1'])*(1400/724.*np.cos(45.58*np.pi/180))+out_header['crval1'])*724
    #ykm=((y[1:]-hhdu['CRPIX2'])*(1400/724.*np.cos(45.58*np.pi/180))+out_header['crval2'])*724
    #r=np.sqrt(x[1:]**2 + y[1:]**2 + z[1:]**2);phix=np.arccos(x[1:]/r);phiy=np.arccos(y[1:]/r)
    #xkm=(x[1:]-hhdu['CRPIX1'])*1400*np.cos(phix)+out_header['crval1']*725;ykm=(y[1:]-hhdu['CRPIX2'])*1400*np.cos(phiy)+out_header['crval2']*725
    #zkm=695700.0+z[1:]*1400.0#;r=np.sqrt(xkm*xkm+ykm*ykm+zkm*zkm);theta=np.arccos(xkm/r);phi=np.arccos(ykm/(r*np.sin(theta)))
    #frame_out = SkyCoord(x[1:] * u.pix, y[1:] * u.pix, obstime='2016/04/09T18:45:00', observer="earth",frame="heliographic_carrington")
    dd1=mymap.pixel_to_world(x[1:]*u.pix,y[1:]*u.pix);dd=dd1.transform_to("heliocentric")
    xkm=dd.cartesian.x.value;ykm=dd.cartesian.y.value;zkm=dd.cartesian.z.value
    sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime="2016/04/09T18:45:00", observer="earth", frame="heliocentric")
    b_hp=sc.transform_to(frames.Helioprojective(observer='earth'))
    b_proj = utils.skycoord_to_pixel(b_hp, out_wcs)
    return x,y,z,bx,by,bz,b_hp,b_proj

#idxsort=sorted(range(len(vtkidx)), key=lambda k: vtkidx[k])
#x,y,z,bx,by,bz,b_hp_fr,b_proj_fr=get_fieldlines(file_name_radio,out_wcs)
babs_=np.sqrt(bx*bx+by*by+bz*bz)
x,y,z,bx,by,bz,b_hp_fr,b_proj_fr=get_fieldlines(file_name_fr,out_wcs)
#x,y,z,bx,by,bz,b_hp_ls,b_proj_ls=get_fieldlines(file_name_ls,out_wcs)
#x,y,z,bx,by,bz,b_hp_euv,b_proj_euv=get_fieldlines(file_name_euv,out_wcs)
tb0_pk1=tb0[:,840:910];tb0_pk2=tb0[:,930:1000];tb0_pk3=tb0[:,1040:1100];tb0_pk4=tb0[:,1130:1200];tb0_pk5=tb0[:,1270:1350];tb0_pk6=tb0[:,1400:1460]
tb1_pk1=tb1[:,840:910];tb1_pk2=tb1[:,930:1000];tb1_pk3=tb1[:,1040:1100];tb1_pk4=tb1[:,1130:1200];tb1_pk5=tb1[:,1270:1350];tb1_pk6=tb1[:,1400:1460]
y0_pk1=y0[:,840:910];y0_pk2=y0[:,930:1000];y0_pk3=y0[:,1040:1100];y0_pk4=y0[:,1130:1200];y0_pk5=y0[:,1270:1350];y0_pk6=y0[:,1400:1460]
y1_pk1=y1[:,840:910];y1_pk2=y1[:,930:1000];y1_pk3=y1[:,1040:1100];y1_pk4=y1[:,1130:1200];y1_pk5=y1[:,1270:1350];y1_pk6=y1[:,1400:1460]
x0_pk1=x0[:,840:910];x0_pk2=x0[:,930:1000];x0_pk3=x0[:,1040:1100];x0_pk4=x0[:,1130:1200];x0_pk5=x0[:,1270:1350];x0_pk6=x0[:,1400:1460]
x1_pk1=x1[:,840:910];x1_pk2=x1[:,930:1000];x1_pk3=x1[:,1040:1100];x1_pk4=x1[:,1130:1200];x1_pk5=x1[:,1270:1350];x1_pk6=x1[:,1400:1460]
listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_3/*spw.*0-15*.FITS'))[0:2000];mapp=Map(listvla_r[800]);dcp_data1=mapp.data;dcp_data1[np.isnan(dcp_data1)]=0;dcp_data1[np.where(dcp_data1<0)]=0
listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_5/*spw.*0-15*.FITS'))[0:2000];mapp=Map(listvla_r[800]);dcp_data2=mapp.data;dcp_data2[np.isnan(dcp_data2)]=0;dcp_data2[np.where(dcp_data2<0)]=0
listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.*0-15*.FITS'))[0:2000];mapp=Map(listvla_r[800]);dcp_data0=mapp.data;dcp_data0[np.isnan(dcp_data0)]=0;dcp_data0[np.where(dcp_data0<0)]=0
dd0=Map(dcp_data0,mapp.meta);dd1=Map(dcp_data1,mapp.meta);dd2=Map(dcp_data2,mapp.meta);qmapp=Map(listvla_r[350]);qmapp_data=qmapp.data;qmapp_data[np.isnan(qmapp_data)]=0;qmapp0=Map(qmapp_data,qmapp.meta)
xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
xlaia1=aiamap.center.Tx.value-2.0*int(aiamap.data.shape[0]/2);xraia1=aiamap.center.Tx.value+2.0*int(aiamap.data.shape[0]/2);ylaia1=aiamap.center.Ty.value-2.0*int(aiamap.data.shape[1]/2);yraia1=aiamap.center.Ty.value+2.0*int(aiamap.data.shape[0]/2)
b_proj_rad = utils.skycoord_to_pixel(b_hp_fr, mapp.wcs)

idx1=np.where((babs_>100) & (babs_<400))[0]-1


xlaia=-948;xraia=-648;ylaia=70;yraia=370;j=0
plot_fields=1
if(plot_fields):
    f,ax0=plt.subplots(1,1)
    p=outmap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #ax0.scatter(b_hp_fr.Tx.value[idx1],b_hp_fr.Ty.value[idx1],s=5,color='tab:olive')
    ax0.scatter(b_hp_fr.Tx.value,b_hp_fr.Ty.value,s=5,color='tab:olive')
    ax0.contour(dd0.data/np.nanmax(dd0.data),levels=[0.6,0.8,0.9],colors='r',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    ax0.contour(dd1.data/np.nanmax(dd1.data),levels=[0.6,0.8,0.9],colors='g',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    ax0.contour(dd2.data/np.nanmax(dd2.data),levels=[0.6,0.8,0.9],colors='b',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    #plt.plot(x0[j],y0[j],'o')
    #plt.plot(x1[j],y1[j],'o')
    ax0.set_xlim(-860,-680);ax0.set_ylim(180,360)
    plt.show()

plot_fields=1
if(plot_fields):
    cm = plt.cm.get_cmap('YlGnBu')
    f,ax0=plt.subplots(1,1)
    p=outmap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #sc=ax0.scatter(b_hp_fr.Tx.value[idx1],b_hp_fr.Ty.value[idx1],s=5,c=babs_[idx1],vmin=100,vmax=400,cmap=cm)
    sc=ax0.scatter(b_hp_fr.Tx.value[idx1],b_hp_fr.Ty.value[idx1],s=5,c=z[idx1],vmin=100,vmax=400,cmap=cm)
    ax0.contour(dd0.data/np.nanmax(dd0.data),levels=[0.6,0.8,0.9],colors='r',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    ax0.contour(dd1.data/np.nanmax(dd1.data),levels=[0.6,0.8,0.9],colors='g',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    ax0.contour(dd2.data/np.nanmax(dd2.data),levels=[0.6,0.8,0.9],colors='b',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    #plt.plot(x0[j],y0[j],'o')
    #plt.plot(x1[j],y1[j],'o')
    ax0.set_xlim(-860,-680);ax0.set_ylim(180,360);plt.colorbar(sc,label='|B| (G)')
    plt.show()

plot_fields=1
if(plot_fields):
    f,ax0=plt.subplots(1,1)
    p=outmap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #ax0.scatter(b_hp_fr.Tx.value[idx1],b_hp_fr.Ty.value[idx1],s=5,color='tab:olive')
    ax0.scatter(b_hp_fr.Tx.value,b_hp_fr.Ty.value,s=5,color='tab:olive')
    #plt.plot(xcr90[15][1000:1400],ycr90[15][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[10][1000:1400],ycr90[10][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[5][1000:1400],ycr90[5][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[1][1000:1400],ycr90[1][1000:1400],'o',alpha=0.2)
    plt.plot(xcr90[12][500:600],ycr90[12][500:600],'o',alpha=0.2)
    plt.plot(xcr90[10][500:600],ycr90[10][500:600],'o',alpha=0.2)
    plt.plot(xcr90[5][500:600],ycr90[5][500:600],'o',alpha=0.2)
    plt.plot(xcr90[1][500:600],ycr90[1][500:600],'o',alpha=0.2)
    ax0.set_xlim(-860,-680);ax0.set_ylim(180,360)
    plt.show()

plot_fields=1
if(plot_fields):
    f,ax0=plt.subplots(1,1)
    p=aiamap.plot(axes=ax0,aspect='auto')
    ax0.contour(qmapp0.data/np.nanmax(qmapp0.data),levels=[0.1,0.2,0.3,0.5,0.7,0.9],colors='r',linewidths=4,origin='lower',extent=[xlvla,xrvla,ylvla,yrvla])
    #plt.plot(x0[j],y0[j],'o')
    #plt.plot(x1[j],y1[j],'o')
    ax0.set_xlim(-1000,-600);ax0.set_ylim(0,400)
    plt.show()

freq=np.round(np.linspace(0.994,2.006,32),3)
plot_paper_timeseries=1
Tball0=np.hstack((qsmaxTbr,maxTbr[0]));Tball0[2132:2170]=np.nan;Tball0[4531:4570]=np.nan;qsmaxTbr[2132:2170]=np.nan;qsmaxTbr[4531:4570]=np.nan
if(plot_paper_timeseries):
    f,ax=plt.subplots(1,1)
    ax.plot(Tball0/1.e6,'o-',label=str(freq[0])+' GHz')
    ax.set_xticks(np.arange(len(Tball0))[::1000]);ax.set_xticklabels(['18:44:00','18:44:50','18:45:40','18:46:30','18:47:20','18:48:10','18:49:00','18:49:50'])
    ax.legend();ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS UT)')
    plt.show()
if(plot_paper_timeseries):
    f,ax=plt.subplots(1,1)
    ax.plot(np.array(qsmaxTbr)/1.e6,'o-',label=str(freq[0])+' GHz')
    ax.set_xticks(np.arange(len(qsmaxTbr))[::1000]);ax.set_xticklabels(['18:44:00','18:44:50','18:45:40','18:46:30','18:47:20'])
    ax.legend();ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS UT)')
    plt.show()
if(plot_paper_timeseries):
    f,ax=plt.subplots(1,1)
    ax.plot(np.array(maxTbr[0])/1.e6,'o-',label=str(freq[0])+' GHz')
    ax.set_xticks(np.arange(len(maxTbr[0]))[::500]);ax.set_xticklabels(['18:48:00','18:48:25','18:48:50','18:49:15','18:49:45'])
    ax.legend();ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS UT)')
    plt.show()

#np.savetxt([b_proj_fr,b_proj_ls,b_proj_euv],'/media/rohit/VLA/paraview/bproj.txt')
#pickle.dump([[b_hp_fr.Tx.value[0],b_hp_fr.Ty.value[1]],[b_hp_ls.Tx.value[0],b_hp_ls.Ty.value[1]],[b_hp_euv.Tx.value[0],b_hp_euv.Ty.value[1]]],open('/media/rohit/VLA/paraview/bproj.p','wb'),protocol=2)

b1=280;b2=170
idxb1=np.where((babs<b1+1) & (babs>b1-1))
idxb2=np.where((babs<b2+1) & (babs>b2-1))

hh2_0=np.array([43.4,40.5,40.1,39.2,38.5,37.7,37.1,36.3,36,35.6,35,34.1,33.5,33.1,32.6,32.1,31.6,31.1,30.8,30.3,29.8,29.1,28.8,28.3,27.9,27.5,27.2,27,26.5,26.2,25.8])
babs_h2=freq/0.0028/2;hh2=[0]*len(babs_h2)
for i in range(len(babs_h2)):
    hh2[i]=ut.find_nearest(babs[:,160,172],babs_h2[i])[0]*1400/1000.


sys.exit()

###########################

omegap=2.8*babs[:,160,172]/1000 # in GHz
h=np.arange(200)*1400/1.e3 # in Mm
hscale=5000;k=10
ne=1.16e17*np.exp(-1000*h/hscale)+1.e9
fp=9000*np.sqrt(ne)/1.e9
ne5=1.16e17*(1+(h*1000/(10*500)))**(-10+1)
fp5=9000*np.sqrt(ne5)/1.e9
ne2=1.16e17*(1+(h*1000/(10*200)))**(-10+1)
fp2=9000*np.sqrt(ne2)/1.e9
ne3=1.16e17*(1+(h*1000/(10*300)))**(-10+1)
fp3=9000*np.sqrt(ne3)/1.e9
plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
plt.plot(h,fp,'-',label='Plasma Freqeuncy ($1.e9+1.16e17*exp(-h/h0)$)')
plt.plot(h,fp5,'-',label='Plasma Freqeuncy ($1.16e17*(1+h/(k*h0))^{-k+1}$) h0=500 km')
plt.plot(h,fp2,'-',label='h0=200 km')
plt.plot(h,fp3,'-',label='h0=300 km')
plt.axhline(y=1,color='k',linewidth=2)
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Freqeuncy (GHz)'),plt.legend()
plt.yscale('log');plt.xlim([0,80]);plt.ylim([0.01,100])
plt.show()

plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
plt.axhline(y=1.4);plt.axhline(y=0.99)
plt.show()
### CANS
d=readsav('/media/rohit/VLA/20160409_sim/code/output/test_param.sav')
x=d['x']*2.e7/1.e8;ro=d['ro']*1.e17;temp=d['te']*1.e4
fpro=9000*np.sqrt(ro)/1.e9
f,ax=plt.subplots(1,1)
ax.plot(x,fpro[0],'o-',label='Plasma Frequency')
ax.plot(h,omegap,'o-',label='Gyro-resonance Freqeuncy')
ax.set_xlabel('Coronal Height (Mm)');ax.set_ylabel('Freqeuncy (GHz)'),ax.legend()
plt.show()

### FORWARD

ff=readsav('/media/rohit/VLA/20160409_EUV/forward/RT_params.sav')
dens=ff['densall'];r=ff['r3dall'];temp=ff['tempall']
plt.plot((r[81,154]-1)*700,dens[81,154],'o-')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Density ($cm^{-3}$)')
plt.show()
plt.plot((r[81,154]-1)*700,9000*np.sqrt(dens[81,154])/1.e6,'o-')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Plasma Freqeuncy (MHz)')
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(x,ro[0],'o-',label='Density')
ax1=ax.twinx()
ax1.plot(x,temp[0]/1.e6,'o-',label='Temperature',color='k')
ax.set_ylabel('Density (cm^-3)');ax1.set_ylabel('Temperature (MK)');ax1.legend()
ax.set_xlabel('Coronal Height (Mm)');ax.set_yscale('log')
plt.show()


dem=readsav('/media/rohit/VLA/20160409/demmap.sav')
xc=dem['demmap']['xc'][0];yc=dem['demmap']['yc'][0];densdata=dem['demmap']['data']
demptemp=[0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0,11.0,14.0]

dem_region=[0]*10
dem_region1=[0]*10
for i in range(10):
    dem_region[i]=densdata[i][80:130,20:35].mean()
    dem_region1[i]=densdata[i][120:160,40:55].mean()

#nk=model.nk_freq2r(np.arange(100)+1000,1)
R=np.linspace(2000,280.e3,1000)
nkne=4.2e4 *10**(4.32/((695700.+R)/695700.))
th=10*np.pi/180.
saitone=(3.09e8*(1-0.5*np.sin(th))/((695700.+R)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+R)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+R)/695700.)**2.5)
nkfp=9000*np.sqrt(nkne)
saitofp=9000*np.sqrt(saitone)

#plt.plot(R,nkne)
plt.plot(R,nkfp/1.e6,label='Newkirk')
plt.plot(R,saitofp/1.e6,label='SAITO')
plt.legend()
#plt.xlabel('Coronal Height (km)');plt.ylabel('Electron Density (ccm)')
plt.xlabel('Coronal Height (km)');plt.ylabel('Plasma Frequency (MHz)')
plt.show()

ne25=1.e9+1.16e17*np.exp(-1*(h*1000/500.))
fp25=9000*np.sqrt(ne25)/1.e9
plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
#plt.plot(h,fpu,'-',label='Plasma Freqeuncy ($1.16e17*(1+h/(k*h0))^{-k+1}$) h0=400 km')
plt.plot(h,fp25,'-',label='Plasma Freqeuncy ($1.e6+1.16e17*e^{(-h/h0)}$) h0=110 km')
#plt.axvline(x=25,color='k',linewidth=4,label='25 Mm')
#plt.axhline(y=1.0,color='r',linewidth=4, label='1 GHz')
#plt.axhline(y=0.4,color='g',linewidth=4, label='400 MHz')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Freqeuncy (GHz)'),plt.legend()
plt.yscale('log');plt.xlim([0,50]);plt.ylim([0.01,100])
plt.show()

plt.imshow(densdata[3],aspect='auto',origin=0,extent=[xc-0.6*83,xc+83,yc-0.6*83,yc+0.6*83])
#plt.plot(xcrmax[0],ycrmax[0],'o',color='k')
plt.title('2.0-3.0 MK')
plt.colorbar(label='EM (cm$^{-5}$)')
plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
plt.show()


h=np.arange(200)*1400/1.e3 # in Mm
h1=np.arange(200)*500/1.e3 # in Mm
sh=np.linspace(10,500,100)
ux=np.linspace(0,50,50)
dem_arr=[0]*100;ne=[0]*100;ne_dem=[0]*100;wpe_dem=[0]*100;omegap_dem=[0]*100
tem=2.6e27

for i in range(100):
    dem_arr[i]=[0]*200;ne[i]=[0]*200
    #ne[i]=1.16e17*(1+(h*1000/(10*sh[i])))**(-10+1)
    ne[i]=1.16e17*(np.exp(-1*h1*1000/sh[i]))+9.5e8*np.exp(-1*h1*1000/31.e3)#+4.2e4 *10**(4.32/((695700.+h*1000)/695700.))
    #ne[i]=4.2e4 *10**(4.32/((695700.+h*1000)/695700.))+1.16e17*(np.exp(-1*h*1000/sh[i]))
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+1.e6
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+(3.09e8*(1-0.5*np.sin(th))/((695700.+h*1000)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+h*1000)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+h*1000)/695700.)**2.5)
    for j in range(200):
        dem_arr[i][j]=np.sum(ne[i][j:]*ne[i][j:])*(h[1]-h[0])*1.e8
    idx=ut.find_nearest(dem_arr[i],tem)[0]
    ne_dem[i]=ne[i][idx];wpe_dem[i]=9000*np.sqrt(ne_dem[i])/1.e9
    omegap_dem[i]=omegap[idx]

dem_arr=np.array(dem_arr);ne=np.array(ne);ne_dem=np.array(ne_dem);wpe_dem=np.array(wpe_dem);omegap_dem=np.array(omegap_dem)
#dem_arr[np.where((dem_arr<2.e27) & (dem_arr>1.e27))]=np.nan

s=2;sf=np.math.factorial(s);e=1.e-19;me=9.1e-31;c=3.e8;nu=1.e9;mu=1;th=45
lb=0.10e6;ratio=wpe_dem/omegap_dem
def optical_depth(s,ratio):
    Ftho=(np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*(np.sin(th*np.pi/180)**2+4*np.cos(th*np.pi/180)**2 +np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))**2)/(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2+np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))
    Fthx=(np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*(np.sin(th*np.pi/180)**2+4*np.cos(th*np.pi/180)**2 -np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))**2)/(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2-np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))
    tauo=np.pi*e*e*ne[28]*1.e9*lb*s*s*s**(2*(s-1))*np.sin(th*np.pi/180)*Ftho/(2*me*c*nu*sf*(2*mu)**(s-1))
    taux=np.pi*e*e*ne[28]*1.e9*lb*s*s*s**(2*(s-1))*np.sin(th*np.pi/180)*Fthx/(2*me*c*nu*sf*(2*mu)**(s-1))
    return tauo,taux

tauo2,taux2=optical_depth(2,ratio)
tauo3,taux3=optical_depth(3,ratio)
tauo4,taux4=optical_depth(5,ratio)


plt.imshow(np.log10(dem_arr),origin=0,aspect='auto',extent=[h[0],h[-1],sh[0],sh[-1]])
plt.colorbar(label='Log(EM(cm$^{-5}$))')
plt.ylabel('Scale Height (km)');plt.xlabel('Source Height (Mm)')
plt.show()

plt.plot(sh,wpe_dem/omegap_dem,'o-')
plt.ylabel('$\omega_p/\Omega_e$');plt.xlabel('Scale Height (km)')
plt.show()

plt.plot(sh,dem_arr[:,10],'o-')
plt.axhline(y=tem,color='k')
plt.xscale('log');plt.yscale('log')
plt.xlabel('Scale Height (km)');plt.ylabel('TEM ($cm^{-5}$)')
plt.show()


f,ax=plt.subplots(2,2)
ax[0,0].plot(h,babs[:,160,172],'o-',label='Magnetic Field')
ax[0,0].set_ylabel('|B| (Gauss)')
ax[0,1].plot(h,ne[28],'o-',label='Electron density')
ax[0,1].set_ylabel(' n$_e$ ($cm^{-3}$)')
ax[1,0].plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax[1,0].set_ylabel('Frequency (GHz)')
ax[1,1].plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax[1,1].plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax[1,1].set_ylabel('Frequency (GHz)')
ax[1,0].set_xlabel('Coronal Height (Mm)');ax[1,1].set_xlabel('Coronal Height (Mm)')
ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax.plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax1=ax.twinx();ax.legend()
ax1.plot(h,tauo3,'o--',color='k',label='O-mode, s=3')
ax1.plot(h,taux3,'s--',color='k',label='X-mode, s=3')
ax1.axhline(y=1.0,color='k')
ax1.legend(loc=2);ax.set_xscale('log');ax.set_yscale('log');ax1.set_xscale('log');ax1.set_yscale('log')
ax.set_xlabel('Coronal Heights (Mm)');ax.set_ylabel('Freqeuncy (GHz)');ax1.set_ylabel('$\\tau$')
plt.show()

vA=babs[:,160,172]/np.sqrt(4*np.pi*ne[28]*9.1e-28)/1.e5
f,ax=plt.subplots(1,1)
ax.plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax.plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax1=ax.twinx();ax.legend()
ax1.plot(h,vA,'-',color='k')
ax.set_xlabel('Coronal Heights (Mm)');ax.set_ylabel('Freqeuncy (GHz)');ax1.set_ylabel('V${_A}$ km/s')
ax1.legend(loc=2);ax.set_xscale('log');ax.set_yscale('log');ax1.set_xscale('log');ax1.set_yscale('log')
plt.show()

plt.plot(h,omegap,'o-',label='Gyroresonce Frequency')
plt.plot(h1,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Frequency h0=150 km')
#plt.axhline(y=1.4,color='k');plt.axhline(y=0.99,color='k')
#plt.axvline(x=20.0,color='k');plt.axvline(x=26.0,color='k')
plt.axvline(x=hh2_0[0],color='k');plt.axvline(x=hh2_0[16],color='k');plt.axvline(x=hh2_0[30],color='k')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Frequency (GHz)'),plt.legend()
plt.yscale('log');plt.xscale('log');plt.xlim([0,100]);plt.ylim([0.01,100])
plt.show()


f,ax=plt.subplots(1,1)
ax.plot(h,ne[28],'o-',label='Inside Density')
ax.plot(h,ne[28]*(2/2.6),'o-',label='Outside Density')
ax.legend();ax.set_ylim(1.e8,1.e9);ax.set_xscale('log')
ax.set_xlim(1.e1,1.e2)
ax.set_ylabel('Density ($cm^{-3}$)');ax.set_xlabel('Coronal Height (Mm)')
plt.show()

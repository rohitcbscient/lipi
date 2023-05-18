import numpy as np
import os
import sys
import pyproj
import scipy.spatial.transform as te
import glob

# Path to Configuration: ~/simulations/configuration


def get_xyz(filename):
    filecontent=open(filename, 'r')
    lines=filecontent.read().split('\n');lines = list(filter(None, lines))
    i=0;x=[0]*len(lines);y=[0]*len(lines);z=[0]*len(lines)
    for l in lines:
        if(l[0]!='#'):
            line=l.replace("\t"," ").split(' ')
            line=[i for i in line if i != '']
            x[i]=np.float(line[0]);y[i]=np.float(line[1]);z[i]=np.float(line[2])
        i=i+1
    x = [i for i in x if i != 0];y = [i for i in y if i != 0];z = [i for i in z if i != 0]
    x=np.array(x);y=np.array(y);z=np.array(z)
    return x,y,z

def convert_ecef2latlong(x,y,z):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
    return lon, lat, alt

def convert_ecef2enu(x,y,z,lat0, lon0, alt0):
    transformer = pyproj.Transformer.from_crs({"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
                                {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},)
    x_org, y_org, z_org = transformer.transform(lon0,lat0,alt0,radians=False)
    vec=np.array([[ x-x_org, y-y_org, z-z_org]]).T
    rot1 =  te.Rotation.from_euler('x', -(90-lat0), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  te.Rotation.from_euler('z', -(90+lon0), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rotMatrix = rot1.dot(rot3)
    enu = rotMatrix.dot(vec).T.ravel()
    return enu.T.reshape(x.shape[0],3)


tel_param={'alma':[-23.0234,-67.7538,5050],'vla':[34.078749167,-107.6177275,2120],'askap':[-26.6961,116.6369,125],
        'meerkat':[-30.72111,21.4111,1000],'mkatplus':[-30.72111,21.4111,1000],'ska1low':[-50,0,1],'WSRT':[52.9145,6.6031,16],
        'sma':[19.8242,-155.4781,4080],'GMRT':[19.0919,74.0506,656],'pdbi':[44.6330,5.8965,2550],'vlba':[19.8016,155.4556,10],
        'aca':[-23.0234,-67.7538,5050],'atca':[-30.31278,149.56278,273],'carma':[37.2385,-118.3041,2196],'lofar':[52.9089,6.8689,15],
        'mwa':[-26.7033,116.6711,125],'ngvla':[34,-107.6,1000],'ska1mid':[-30.72111,21.4111,1000],'ska1low':[-26.82472208,116.7644482,120]} # Lat, Lon, alt

cfg_list=sorted(glob.glob('*.cfg'))
for filename in [cfg_list[173]]:
    tel_list=list(tel_param.keys())
    tel=filename.split('.')[0].split('_')[0].split('-')[0]
    if((tel!='alma') & (tel!='aca') & (tel!='sma') & (tel!='meerkat') & (tel!='mwa')):
        x,y,z=get_xyz(filename)
        lon, lat, alt = convert_ecef2latlong(x,y,z)
        lat0,lon0,alt0=tel_param[tel]
        enu_coord=convert_ecef2enu(x,y,z,lat0, lon0, alt0)
        enu_x,enu_y,enu_z=enu_coord[:,0],enu_coord[:,1],enu_coord[:,2]
    else:
        lat0,lon0,alt0=tel_param[tel]
        enu_x,enu_y,enu_z=get_xyz(filename)
    os.system('rm -rf '+filename.split('.cfg')[0]+'.tm')
    os.system('mkdir '+filename.split('.cfg')[0]+'.tm')
    station_num=len(enu_x)
    #np.savetxt(filename.split('.cfg')[0]+'.tm/layout.txt',np.array((enu_x,enu_y)).T)
    np.savetxt(filename.split('.cfg')[0]+'.tm/layout.txt',np.array((enu_x,enu_y,enu_z)).T)
    np.savetxt(filename.split('.cfg')[0]+'.tm/position.txt',np.c_[lon0,lat0])
    for i in range(station_num):
        ii="%03d"%i
        os.system('mkdir '+filename.split('.cfg')[0]+'.tm/station'+ii)
        os.system('cp /home/rohit/simulations/layout_station.txt '+filename.split('.cfg')[0]+'.tm/station'+ii+'/layout.txt')




sys.exit()

lay=np.loadtxt('meerkat_config.txt')
station_num=lay.shape[0]
name='meerkat'
os.system('mkdir '+name+'.tm')
os.system('cp meerkat_config.txt '+name+'.tm')
for i in range(station_num):
    ii="%03d" %i
    os.system('mkdir '+name+'.tm/station'+ii)
    os.system('cp layout_station.txt '+name+'.tm/station'+ii+'/layout.txt')
    #np.savetxt(name+'.tm/station'+ii+'/layout.txt',np.array([0,0]).T)



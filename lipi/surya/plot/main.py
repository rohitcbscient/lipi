import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np


def make_circle(axx,x,y,res,xshape,s):
    '''
    axx: axis
    x: x-coordinate in 0,1 values
    s: radius in arcsec
    '''
    r=(s*60.0/res)*(1.0/xshape)
    #axx=fig.add_subplot(1,1,1)
    circ=plt.Circle((x,y), radius=r, color='red', linewidth=3,fill=False,transform=axx.transAxes)
    axx.add_patch(circ)

def add_beam(ax,xcenter,ycenter,width, height,angle):
    theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
    x = 0.5 * width * np.cos(theta)
    y = 0.5 * height * np.sin(theta)
    rtheta = np.radians(angle)
    R = np.array([[np.cos(rtheta), -np.sin(rtheta)],[np.sin(rtheta),np.cos(rtheta)],])
    x, y = np.dot(R, np.array([x, y]))
    x += xcenter
    y += ycenter
    ax.fill(x, y, alpha=0.0, facecolor='yellow', edgecolor='grey', linewidth=0.2, zorder=1)

    e1 = patches.Ellipse((xcenter, ycenter), width, height,
                 angle=angle, linewidth=2, fill=False, zorder=2)
    ax.add_patch(e1)

def plot_polar_dot(data,colorbar,az_sun,el_sun):#,title):
        '''
        Takes in an array in AZ, EL (and not ZA) and produces a plot
        in AZ, ZA coordinates.
        '''
        za=90-el_sun
        azimuths = np.radians(np.linspace(0, 360, 360))
        zeniths = np.arange(90, 0, -1)
        r, theta = np.meshgrid(zeniths, azimuths)
        plt.rcParams.update({'font.weight':'bold' })
        plt.rc('grid', c='0.1', ls='-', lw=0.5)
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.tick_params(axis='y', colors='black', which='both')
        ax1=ax.contourf(theta, r, 10*np.log10(data), 50, cmap=plt.cm.jet, interpolation='none')
        #ax1=ax.contourf(theta, r, data, 50, norm=LogNorm(vmin=1e-7,vmax=1), cmap=plt.cm.jet, interpolation='none')
        #ax1=ax.contourf(theta, r, data, 50, cmap=plt.cm.jet, interpolation='none')
        ax.scatter([np.radians(az_sun)], [za], marker='*',s=140, c='orange')
        #ax.set_title(str(title))
        ax.set_theta_offset(np.pi/2)
        ax.set_rmax(90)
        ax.set_rmin(1)
        #tick=[-1,-0.25,-0.5,-0.75, 0,0.25,0.5,0.75, 1]
        tick=[0,-10,-20,-30,-40,-50,-60,-70]
        cbar = plt.colorbar(ax1,orientation='vertical',ticks=tick)
        cbar.ax.set_ylabel(str(colorbar),fontweight='bold')
        plt.show()



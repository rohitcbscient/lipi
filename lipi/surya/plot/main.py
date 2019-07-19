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
	ax.fill(x, y, alpha=0.4, facecolor='yellow', edgecolor='yellow', linewidth=0.5, zorder=1)

        e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)
        ax.add_patch(e1)


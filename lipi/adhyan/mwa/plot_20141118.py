import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
from matplotlib import patches


def ds_plot(Tsun_plot,ytick_list,xtick_list,ylabel_list,xlabel_list,vmn,vmx,cmtitle,tick):
        flux_space=put_spaces(Tsun_plot)
        print flux_space.shape
        f = plt.gcf()
        plt.rcParams['axes.linewidth'] = 1
        f, ax2=plt.subplots(1, 1, sharey=True)
        im2=ax2.imshow(flux_space,interpolation='none',aspect='auto',cmap='jet',norm=LogNorm(vmin=vmn,vmax=vmx))
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        cbar2 = plt.colorbar(im2, cax=cax2, format="%.2f",ticks=tick)
        cbar2.set_label(cmtitle,fontweight='bold')
        #ax2.set_xlim([d1,d2])
        xa = ax2.get_xaxis()
        ax2.set_xlabel("Time (HH:MM:SS UT)",fontweight='bold')
        ax2.set_ylabel("Frequency (MHz)",fontweight='bold')
        ax2.tick_params('both', length=5, width=2, which='major')
        ax2.tick_params('both', length=5, width=2, which='minor')
        ax2.set_yticks(ytick_list)
        ax2.set_xticks(xtick_list)
        ax2.set_yticklabels(ylabel_list)
        ax2.set_xticklabels(xlabel_list)
        plt.rc('font',weight='bold')
        #plt.savefig(filename)
        #xa.set_major_formatter(FuncFormatter(lambda x, pos:str(start + datetime.timedelta(seconds=x))))
        f.show()


def put_spaces(flux):
        n=flux.shape[0]
        flux_space=flux.reshape(flux.shape[0]*flux.shape[1],flux.shape[2])
        for i in range(n-1):
                flux_space = np.insert(flux_space, 64*(i+1)+24*i, np.zeros((24,flux.shape[2])), 0)
        return  flux_space

def plot_map(d,title,bmin,bmax,levels):
        fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=False,figsize=(9,8))
        im1=plt.imshow(d,cmap='hot',origin='lower',extent=[-2.75,2.75,-2.75,2.75],norm=LogNorm(vmin=1.e3,vmax=1.e6))
        beam_centre=[-2,-2]
        e1 = patches.Ellipse((beam_centre[0], beam_centre[1]), bmin,bmax, angle=0, linewidth=2, fill=False, zorder=2,color='magenta')
        ax1.add_patch(e1)
        x2=np.arange(d.shape[0])
        y2=np.arange(d.shape[1])
        X2,Y2=np.meshgrid(x2,y2)
        X21=(X2-d.shape[0]*0.5)*99/3600
        Y21=(Y2-d.shape[1]*0.5)*99/3600
        #levels=np.array([0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])*np.max(mwa)
        #levels=np.array([0.009,0.09,0.9])*np.max(d[np.isfinite(d)])
        plt.contour(X21,Y21,d, levels, hold='on', colors='cyan',linewidths=2)
        plt.xlim([-2.75,2.75])
        plt.ylim([-2.75,2.75])
        plt.grid(True)
        plt.xlabel('X (arcsecs)')
        plt.ylabel('Y (arcsecs)')
        plt.title(title)
        #plt.annotate(str(freq[i])+' MHz', xy=(-300, 1400),  xycoords='data',xytext=(-300, 1400))
        plt.show()


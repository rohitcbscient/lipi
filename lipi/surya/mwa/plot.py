import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter

def plot_polar(data,colorbar):
	'''
	Takes in an array in AZ, EL (and not ZA) and produces a plot
	in AZ, ZA coordinates.
	'''
	azimuths = np.radians(np.linspace(0, 360, 360))
	zeniths = np.arange(90, 0, -1)
	r, theta = np.meshgrid(zeniths, azimuths)
	plt.rcParams.update({'font.weight':'bold' })
	fig, ax = plt.subplots(figsize=(20,5),subplot_kw=dict(projection='polar'))
	#ax1=ax.contourf(theta, r, data, 50, locator=ticker.LogLocator(), cmap=plt.cm.jet, interpolation='none')
	ax1=ax.contourf(theta, r, data, 50, cmap=plt.cm.jet, interpolation='none')
	ax.set_theta_offset(np.pi/2)
	cbar = plt.colorbar(ax1,orientation='vertical')
	cbar.ax.set_ylabel(str(colorbar),fontweight='bold',fontsize=20)
	plt.show()
	
def plot_polar_dot_contour(data,data1,colorbar,az_sun,el_sun):
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
	ax.tick_params(axis='y', colors='yellow', which='both')
	#ax1=ax.contourf(theta, r, data, 50, locator=ticker.LogLocator(), cmap=plt.cm.jet, interpolation='none')
	ax1=ax.contourf(theta, r, data,np.linspace(0, 600, 50), cmap=plt.cm.jet, interpolation='none')
	#ax1=ax.contourf(theta, r, data,np.linspace(-1, 1, 50), cmap=plt.cm.winter, interpolation='none')
	#ax1=ax.contourf(theta, r, data, 50, norm=LogNorm(), cmap=plt.cm.jet, interpolation='none')
	ax2=ax.contour(theta, r,data1, levels=[0.9,0.09, 0.009],colors = 'red',linewidths=3)
	ax.scatter([np.radians(az_sun)], [za], marker='*',s=200, c='orange')
	ax.set_theta_offset(np.pi/2)
	ax.set_rmax(90)
	ax.set_rmin(1)
	#tick=[-1,-0.25,-0.5,-0.75, 0,0.25,0.5,0.75, 1]
	tick=np.arange(7)*100
	cbar = plt.colorbar(ax1,orientation='vertical',ticks=tick)
	cbar.ax.set_ylabel(str(colorbar),fontweight='bold',fontsize=20)
	plt.show()	
	
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
	cbar.ax.set_ylabel(str(colorbar),fontweight='bold',fontsize=20)
	plt.show()	
	
def plot_polar_log(data,colorbar,az,el):
	'''
	Takes in an array in AZ, EL (and not ZA) and produces a plot
	in AZ, ZA coordinates.
	'''
	azimuths = np.radians(np.linspace(0, 360, 360))
	zeniths = np.arange(90, 0, -1)
	r, theta = np.meshgrid(zeniths, azimuths)
	fig = plt.gcf()
	fig.set_size_inches(8.5, 6.5)
	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
#	ax1=ax.contourf(theta, r, data, 50, locator=ticker.LogLocator(), cmap=plt.cm.jet, interpolation='none')
	#levels = np.power(10, -7.0*np.arange(20)[::-1]/20) # nprimary beam
	levels = -70.0*np.arange(50)[::-1]/50
	ax1=ax.contourf(theta, r, data, levels,cmap=plt.cm.jet, interpolation='none')
	#ax.scatter(np.radians(az), 90-np.array(el), marker='*',s=140, c='orange')
	plt.plot(np.radians(az), 90-np.array(el), c='cyan',linewidth=4)
	tick=[0,-10,-20,-30, -40,-50,-60]	
	cbar = plt.colorbar(ax1, ticks=tick,format="%.0f")
	cbar.set_label(colorbar,fontweight='bold')
	ax.set_theta_offset(np.pi/2)
	ax.set_rmax(90)
	ax.set_rmin(1)
	plt.show()
	
def plot_fraction(freq,frac_flux11_22,frac_flux21_22,frac_flux21_23,xx,yy):
	plt.plot(freq,frac_flux21_23,'o-',c='red',label='Tile021-Tile023',linewidth=4,markersize=10)
	plt.plot(freq,frac_flux21_22,'o-',c='green',label='Tile021-Tile022',linewidth=4,markersize=10)
	plt.plot(freq,frac_flux11_22,'o-',c='blue',label='Tile011-Tile022',linewidth=4,markersize=10)
	plt.plot(freq,xx,'o-',c='black',label='X',linewidth=4,markersize=10)
	plt.plot(freq,yy,'o--',c='black',label='Y',linewidth=4,markersize=10)
	plt.rcParams.update({'font.weight':'bold' })
	plt.rc('bold')
	plt.xlim(80,300)
	plt.xlabel('Frequency (MHz)',fontweight='bold')
	plt.ylabel('Flux Captured ($\%$)',fontweight='bold')
	plt.legend(loc=3)
	plt.show()	
	
	
def ds_plot(Tsun_plot,Tn,filename,cmtitle):
	ytick_list=[0,28,42,70,84,112,126,154,168,196,210,238,252,280,294,322,336,364,378,406]
	#flist_space=[102.25,103.50,115.00,117.50,129.00,131.50,145.00,147.50,165.00,167.50,186.00,188.50,211.00,213.50,239.25,240.50,271.25,272.50,296.00,298.50]
	flist_space=[101.12,102.40,115.20,117.76, 129.28, 131.84, 145.92,148.48,165.12,167.68,186.88,189.44,211.20,213.76,238.08,239.36,270.08,271.36,296.96, 298.24]
	xtick_list=[0,35,70,105,140,175,210,245]
	#x_space=['03:40','03:50','04:00','04:10','04:20','04:30','04:40']
	x_space=['03:56:32','03:57:07','03:57:42','03:58:17','03:58:52','03:59:27','04:00:02','04:00:40']
	#x_space=['03:56:32','04:00:40','04:04:48','04:08:56','04:13:04','04:17:21','04:21:20']
	f = plt.gcf()
	plt.rcParams['axes.linewidth'] = 1
	f, ax2=plt.subplots(1, 1, sharey=True)
	im2=ax2.imshow(Tsun_plot[::-1,:]/Tn,interpolation='none',aspect='auto',cmap='jet',extent=[0,490/2,0,406],norm=LogNorm(vmin=0.1, vmax=35))
	divider2 = make_axes_locatable(ax2)
	#tick=[10,50,100,150,200,250,300,400,500]
	tick=[0.1,0.3,1,3,10,30]
	#tick=[0.01,0.03,0.05,0.10,0.30,0.50]
	cax2 = divider2.append_axes("right", size="5%", pad=0.05)
	cbar2 = plt.colorbar(im2, cax=cax2, format="%.2f",ticks=tick)
	cbar2.set_label(cmtitle,fontweight='bold')
	#ax2.set_xlim([d1,d2])
	xa = ax2.get_xaxis()
	ax2.set_xlabel("03-September-2013 (HH:MM UT)",fontweight='bold')
	ax2.set_ylabel("Frequency (MHz)",fontweight='bold')
	ax2.tick_params('both', length=5, width=2, which='major')
	ax2.tick_params('both', length=5, width=2, which='minor')
	ax2.set_yticks(ytick_list)
	ax2.set_xticks(xtick_list)
	ax2.set_yticklabels(flist_space)
	ax2.set_xticklabels(x_space)
	plt.rc('font',weight='bold')
	#plt.savefig(filename)
	#xa.set_major_formatter(FuncFormatter(lambda x, pos:str(start + datetime.timedelta(seconds=x))))
	f.show()
	
def ds_plot_conti(Tsun_plot,Tn,filename,cmtitle):
	#ytick_list=[0,28,42,70,84,112,126,154,168,196,210,238,252,280,294,322,336,364,378,406]
	ytick_list=406*np.arange(20)/20
	flist_space=103.8+np.arange(20)*30/20.0
	f = plt.gcf()
	plt.rcParams['axes.linewidth'] = 1
	f, ax2=plt.subplots(1, 1, sharey=True)
	im2=ax2.imshow(Tsun_plot[::-1,:]/Tn,interpolation='none',aspect='auto',cmap='jet',extent=[0,480/2,0,406],norm=LogNorm(vmin=0.1, vmax=200))
	divider2 = make_axes_locatable(ax2)
	#tick=[10,50,100,150,200,250,300,400,500]
	tick=[0.1,0.3,1,3,10,30,100,300]	
	#tick=[0.01,0.03,0.05,0.10,0.30,0.50]
	cax2 = divider2.append_axes("right", size="5%", pad=0.05)
	cbar2 = plt.colorbar(im2, cax=cax2, format="%.2f",ticks=tick)
	cbar2.set_label(cmtitle,fontweight='bold')
	#ax2.set_xlim([d1,d2])
	xa = ax2.get_xaxis()
	ax2.set_xlabel("Time (in sec)",fontweight='bold')
	ax2.set_ylabel("Frequency (MHz)",fontweight='bold')
	ax2.set_yticks(ytick_list)
	ax2.set_yticklabels(flist_space)
	ax2.tick_params('both', length=5, width=2, which='major')
	ax2.tick_params('both', length=5, width=2, which='minor')
	plt.rc('font',weight='bold')
	#plt.savefig(filename)
	#xa.set_major_formatter(FuncFormatter(lambda x, pos:str(start + datetime.timedelta(seconds=x))))
	f.show()	
	
def temp_scatter(beam_ave11_22xx,beam_ave11_22yy,flist):
	color=['b','g','r','c','m','y','k','b','g','r']
	marker=['o','o','o','o','o','o','o','v','v','v']
	plt.rcParams.update({'font.weight':'bold' })
	plt.rc('grid', c='0.1', ls='-', lw=0.5)
	for i in range(10):
		plt.plot(beam_ave11_22xx[i*42:42*(i+1)-14,:].flatten(),beam_ave11_22yy[i*42:42*(i+1)-14,:].flatten(),c=color[i],alpha=0.7,linestyle='',marker=marker[i],markersize=10,label=str(flist[i])+' MHz')
	plt.xlabel('XX (K)',fontweight='bold')
	plt.ylabel('YY (K)',fontweight='bold')
	plt.legend(loc=4)
	plt.rc('bold')
	plt.show()			
	
def ratio_histograms1(ratio_11_22by21_22xx,ratio_11_22by21_23xx,ratio_21_22by21_23xx,ratio_11_22by22_23xx,ratio_21_22by22_23xx,ratio_21_23by22_23xx,b,ylim):
	b1=int((max(ratio_11_22by21_22xx)-min(ratio_11_22by21_22xx))/0.004)
	b2=int((max(ratio_11_22by21_23xx)-min(ratio_11_22by21_23xx))/0.004)
	b3=int((max(ratio_21_22by21_23xx)-min(ratio_21_22by21_23xx))/0.004)
	b4=int((max(ratio_11_22by22_23xx)-min(ratio_11_22by22_23xx))/0.004)
	b5=int((max(ratio_21_22by22_23xx)-min(ratio_21_22by22_23xx))/0.004)
	b6=int((max(ratio_21_23by22_23xx)-min(ratio_21_23by22_23xx))/0.004)
	plt.hist(ratio_11_22by21_22xx,bins=b1,histtype='step',color='blue',label='Tile011-Tile022 / Tile021-Tile022',lw=3)
	plt.hist(ratio_11_22by21_23xx,bins=b2,histtype='step',color='red',label='Tile011-Tile022 / Tile021-Tile023',lw=3)
	plt.hist(ratio_21_22by21_23xx,bins=b3,histtype='step',color='green',label='Tile021-Tile022 / Tile021-Tile023',lw=3)
	plt.hist(ratio_11_22by22_23xx,bins=b4,histtype='step',color='black',label='Tile011-Tile022 / Tile022-Tile023',lw=3)
	plt.hist(ratio_21_22by22_23xx,bins=b5,histtype='step',color='magenta',label='Tile021-Tile022 / Tile022-Tile023',lw=3)
	plt.hist(ratio_21_23by22_23xx,bins=b6,histtype='step',color='orange',label='Tile021-Tile023 / Tile022-Tile023',lw=3)
	plt.rcParams.update({'font.weight':'bold' })
	plt.legend(loc=1)
	plt.xlabel('Ratio',fontweight='bold')
	plt.ylabel('Number of data points',fontweight='bold')
	plt.ylim(0,2400)
	#plt.savefig('four_baseline_240MHz.eps')
	#plt.yscale('log')
	#plt.xticks(0.25*np.arange(40))
	plt.show()
	
def ratio_histograms(ratio_11_22by21_22xx,ratio_11_22by21_23xx,ratio_21_22by21_23xx,ratio_11_22by22_23xx,ratio_21_22by22_23xx,ratio_21_23by22_23xx,b,ylim):
	plt.hist(ratio_11_22by21_22xx,bins=b,histtype='step',label='Tile011-Tile022 / Tile021-Tile022',lw=3)
	plt.hist(ratio_11_22by21_23xx,bins=b,histtype='step',label='Tile011-Tile022 / Tile021-Tile023',lw=3)
	plt.hist(ratio_21_22by21_23xx,bins=b,histtype='step',label='Tile021-Tile022 / Tile021-Tile023',lw=3)
	plt.hist(ratio_11_22by22_23xx,bins=b,histtype='step',label='Tile011-Tile022 / Tile022-Tile023',lw=3)
	plt.hist(ratio_21_22by22_23xx,bins=b,histtype='step',label='Tile021-Tile022 / Tile022-Tile023',lw=3)
	plt.hist(ratio_21_23by22_23xx,bins=b,histtype='step',label='Tile021-Tile023 / Tile022-Tile023',lw=3)
	plt.rcParams.update({'font.weight':'bold' })
	plt.legend(loc=1)
	plt.xlabel('Ratio',fontweight='bold')
	plt.ylabel('Number of data points',fontweight='bold')
	plt.ylim(0,ylim)
	#plt.savefig('four_baseline_240MHz.eps')
	#plt.yscale('log')
	#plt.xticks(0.25*np.arange(40))
	plt.show()	
	
def plot_median(freq_list,median_11_22by21_22xx,median_11_22by21_23xx,median_21_22by21_23xx,median_11_22by22_23xx,median_21_22by22_23xx,median_21_23by22_23xx):
	plt.plot(freq_list,median_11_22by21_22xx,'o-',c='blue',label='Tile011-Tile022 / Tile021-Tile022',lw=3)
	plt.plot(freq_list,median_11_22by21_23xx,'o-',c='red',label='Tile011-Tile022 / Tile021-Tile023',lw=3)
	plt.plot(freq_list,median_21_22by21_23xx,'o-',c='green',label='Tile021-Tile022 / Tile021-Tile023',lw=3)
	plt.plot(freq_list,median_11_22by22_23xx,'o-',c='black',label='Tile011-Tile022 / Tile022-Tile023',lw=3)	
	plt.plot(freq_list,median_21_22by22_23xx,'o-',c='magenta',label='Tile021-Tile022 / Tile022-Tile023',lw=3)	
	plt.plot(freq_list,median_21_23by22_23xx,'o-',c='orange',label='Tile021-Tile023 / Tile022-Tile023',lw=3)	
	plt.rcParams.update({'font.weight':'bold' })
	plt.legend(loc=1)
	plt.xlabel('Frequency (MHz)',fontweight='bold')
	#plt.ylabel('Median',fontweight='bold')
	plt.ylabel('FWHM',fontweight='bold')
	#plt.ylim(0.5,1.6)
	plt.show()
	
def plot_median_errors(freq_list,median_11_22by21_22xx,median_11_22by21_23xx,median_21_22by21_23xx,median_11_22by22_23xx,median_21_22by22_23xx,median_21_23by22_23xx,e_11_22by21_22xx,e_11_22by21_23xx,e_21_22by21_23xx,e_11_22by22_23xx,e_21_22by22_23xx,e_21_23by22_23xx):
	plt.figure()
	plt.rcParams.update({'font.weight':'bold' })
	plt.errorbar(freq_list,median_11_22by21_22xx,yerr=e_11_22by21_22xx,fmt='o-',c='blue',label='Tile011-Tile022 / Tile021-Tile022',lw=3)
	plt.errorbar(freq_list,median_11_22by21_23xx,yerr=e_11_22by21_23xx,fmt='o-',c='red',label='Tile011-Tile022 / Tile021-Tile023',lw=3)
	plt.errorbar(freq_list,median_21_22by21_23xx,yerr=e_21_22by21_23xx,fmt='o-',c='green',label='Tile021-Tile022 / Tile021-Tile023',lw=3)
	plt.errorbar(freq_list,median_11_22by22_23xx,yerr=e_11_22by22_23xx,fmt='o-',c='black',label='Tile011-Tile022 / Tile022-Tile023',lw=3)	
	plt.errorbar(freq_list,median_21_22by22_23xx,yerr=e_21_22by22_23xx,fmt='o-',c='magenta',label='Tile021-Tile022 / Tile022-Tile023',lw=3)	
	plt.errorbar(freq_list,median_21_23by22_23xx,yerr=e_21_23by22_23xx,fmt='o-',c='orange',label='Tile021-Tile023 / Tile022-Tile023',lw=3)	
	plt.legend(loc=1)
	plt.xlim(90,300)
	plt.xlabel('Frequency (MHz)',fontweight='bold')
	plt.ylabel('Median',fontweight='bold')
	plt.show()	
	
def plot_uvw(u,v,w,filename):
	plt.plot(u,label='u',marker='o',color='blue')
	plt.plot(v,label='v',marker='o',color='red')
	plt.plot(w,label='w',marker='o',color='green')
	plt.legend(loc=1)
	plt.savefig(str(filename))
	plt.close()
	
def plot_azelsun(az,el,filename):
	plt.plot(az,label='AZ',marker='o',color='blue')
	plt.plot(el,label='EL',marker='o',color='red')
	plt.legend(loc=1)
	plt.savefig(str(filename))
	plt.close()
	
def plot_azelphases(az,el,filename):
	plt.plot(az,label='AZ',marker='o',color='blue')
	plt.plot(el,label='EL',marker='o',color='red')
	plt.legend(loc=1)
	plt.savefig(str(filename))
	plt.close()
	
def plot_correctionfact(data,data1,filename):
	plt.plot(data,label='Fringe factor',marker='o',color='red')
	plt.plot(data1,label='NCCF factor',marker='o',color='blue')
	plt.legend(loc=1)
	plt.savefig(str(filename))
	plt.close()	
		

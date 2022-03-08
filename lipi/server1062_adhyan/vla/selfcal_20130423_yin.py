## Script for self-calibration--spw7--L band
selfpath="/nas08-data02/vladata/20130423/L-Band/selfcal/yin_test/"
selfcal0=selfpath+'spw7self0'+'.ms'
selfcal1=selfpath+'spw7self1'+'.ms'
selfcal2=selfpath+'spw7self2'+'.ms'
selfcal3=selfpath+'spw7self3'+'.ms'
selfcal4=selfpath+'spw7self4'+'.ms'

trange='2013/04/23/20:45:17~2013/04/23/20:45:23'
refant='ea07'
antennas='ea01,ea03,ea06,ea07,ea08,ea09,ea17,ea18,ea19,ea21,ea23,ea25,ea26,ea27'
slfcal_table0 = selfpath+'selfcal.G0'
slfcal_table1 = selfpath+'selfcal.G1'
slfcal_table2 = selfpath+'selfcal.G2'
slfcal_table3 = selfpath+'selfcal.G3'
slfcal_table4 = selfpath+'selfcal.G4'


split(vis="/nas08-data02/vladata/20130423/L-Band/selfcal/yin_test/20130423T2030-2050.L.50ms.selfcal.ms",outputvis=selfcal0,spw='7',datacolumn='data',timerange='20:40:00~20:46:00')
clearcal(vis=selfcal0)
#clean(vis=selfcal0,imagename=selfpath+'selinit0',spw='0:3~60',timerange=trange,interactive=False,usescratch=False,npercycle=50\
#		,imsize= [512, 512],cell= ['5.0arcsec', '5.0arcsec'],phasecenter= 'J2000 0.548505138896rad 0.222291421955rad',stokes='RRLL",niter=50)
clean(vis=selfcal0,imagename=selfpath+'ori0',spw='0:3~60',timerange=trange,interactive=True,usescratch=False,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

clean(vis=selfcal0,imagename=selfpath+'selinit0',spw='0:3~60',timerange=trange,interactive=False,usescratch=False,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)
##first round
clean(vis=selfcal0,imagename=selfpath+'selinit',spw='0:3~60',timerange=trange,interactive=True,usescratch=False,npercycle=50\
		,imsize= [512, 512],cell= ['5.0arcsec', '5.0arcsec'],phasecenter= 'J2000 0.548505138896rad 0.222291421955rad',stokes="RRLL",niter=50)
		

clean(vis=selfcal0,imagename=selfpath+'sel0',spw='0:3~60',timerange=trange,interactive=True,usescratch=True,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

gaincal(vis=selfcal0, refant=refant,caltable=slfcal_table0,spw='0:28~32',\
		selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
		minblperant=3,minsnr=2,append=True)

##ant ea19 has no solution
#plotcal(caltable=slfcal_table0,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
#iteration='spw',subplot=111,figfile=slfcal_table0+'_phase.png')

plotcal(caltable=slfcal_table0,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
		iteration='',subplot=111,showgui=False,figfile=slfcal_table0+'_phase.png',showflags=False,overplot=False,clearpanel="Auto")

##calonly here
applycal(vis=selfcal0,gaintable=slfcal_table0,spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal1,datacolumn='corrected')



##2nd

clean(vis=selfcal1,imagename=selfpath+'sel1',spw='0:3~60',timerange=trange,interactive=True,usescratch=True,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

gaincal(vis=selfcal1, refant=refant,caltable=slfcal_table1,spw='0:28~32',\
		selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
		minblperant=3,minsnr=3,append=True)

#plotcal(caltable=slfcal_table1,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
	#	iteration='spw',subplot=111)

plotcal(caltable=slfcal_table1,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
		iteration='',subplot=111,showgui=False,figfile=slfcal_table1+'_phase.png',showflags=False,overplot=False,clearpanel="Auto")


clearcal(vis=selfcal0)
applycal(vis=selfcal0,gaintable=[slfcal_table0,slfcal_table1],spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal2,datacolumn='corrected')

##3rd
clean(vis=selfcal2,imagename=selfpath+'sel2',spw= '0:3~60',timerange=trange,interactive=True,usescratch=True,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

gaincal(vis=selfcal2, refant=refant,caltable=slfcal_table2,spw='0:28~32',\
		selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
		minblperant=3,minsnr=3,append=True)

#plotcal(caltable=slfcal_table2,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
#		iteration='spw',subplot=111)

plotcal(caltable=slfcal_table2,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
		iteration='',subplot=111,showgui=False,figfile=slfcal_table2+'_phase.png',showflags=False,overplot=False,clearpanel="Auto")

clearcal(vis=selfcal0)
applycal(vis=selfcal0,gaintable=[slfcal_table0,slfcal_table1,slfcal_table2],spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal3,datacolumn='corrected')


##4th
clean(vis=selfcal3,imagename=selfpath+'sel3',spw= '0:3~60',timerange=trange,interactive=True,usescratch=True,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

gaincal(vis=selfcal3, refant=refant,caltable=slfcal_table3,spw='0:28~32',\
		selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
		minblperant=3,minsnr=3,append=True)

#plotcal(caltable=slfcal_table2,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
#		iteration='spw',subplot=111)

plotcal(caltable=slfcal_table3,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
		iteration='',subplot=111,showgui=False,figfile=slfcal_table3+'_phase.png',showflags=False,overplot=False,clearpanel="Auto")

clearcal(vis=selfcal0)
applycal(vis=selfcal0,gaintable=[slfcal_table0,slfcal_table1,slfcal_table2,slfcal_table3],spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal4,datacolumn='corrected')




##5th
clean(vis=selfcal4,imagename=selfpath+'sel4',spw= '0:3~60',timerange=trange,interactive=True,usescratch=True,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL",niter=50)

gaincal(vis=selfcal4, refant=refant,caltable=slfcal_table4,spw='0:28~32',\
		selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
		minblperant=3,minsnr=3,append=True)

#plotcal(caltable=slfcal_table2,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
#		iteration='spw',subplot=111)

plotcal(caltable=slfcal_table4,antenna=antennas,xaxis='antenna',yaxis='phase',plotrange=[-1,-1,-180,180],
		iteration='',subplot=111,showgui=False,figfile=slfcal_table4+'_phase.png',showflags=False,overplot=False,clearpanel="Auto")

clearcal(vis=selfcal0)
applycal(vis=selfcal0,gaintable=[slfcal_table0,slfcal_table1,slfcal_table2,slfcal_table3],spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal4,datacolumn='corrected')





##all antennas phase are flat, no need for the deeper round

clearcal(vis=selfcal0)
applycal(vis=selfcal0,gaintable=[slfcal_table0,slfcal_table1,slfcal_table2],spw='0',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis=selfcal0,outputvis=selfcal3,datacolumn='corrected')
## do the final
clean(vis=selfcal3,imagename=selfpath+'self',spw= '0:3~60',timerange=trange,interactive=False,usescratch=False,npercycle=50\
		,imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL")
		
clean(vis=selfcal3,imagename=selfpath+'selfte',spw= '0:3~60',timerange=trange,interactive=False,usescratch=False,npercycle=50\
		,antenna='ea01,ea03,ea06,ea07,ea08,ea09,ea17,ea18,ea21,ea23,ea25,ea26,ea27',imsize= [256, 256],cell= ['4.0arcsec', '4.0arcsec'],phasecenter= 'J2000 02h04m51.211 12d43m12.47',stokes="RRLL")

split(vis="/srg/yjluo/xixi/20130423/20130423T2030-2050.L.50ms.cal.ms/",outputvis='20130423T2040-2050.L.50ms.cal.ms',datacolumn='data',timerange='20:40:00~20:50:00')
applycal(vis='20130423T2040-2050.L.50ms.cal.ms',gaintable=[slfcal_table0,slfcal_table1,slfcal_table2],spw='7',selectdata=True,\
		antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
split(vis='20130423T2040-2050.L.50ms.cal.ms',outputvis='20130423T2040-2050.L.50ms.selfcal.ms',datacolumn='corrected')

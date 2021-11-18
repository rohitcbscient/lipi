

# For 20130423 L-band data
# RA-DEC
# RA: 01 19 09.77 
# DEC: +08 20 56.8

#### Telescope flag #####

vis='20130423_L_T2036-2046_50ms.ms'

scan0='51' # Flag system configuration scans and scans of the start of the field change. 
flagdata(vis=vis,scan=scan0)

quackinterval_= 6.0 # Flag the time at the start of calibrator scans
flagdata(vis=vis,mode='quack',quackmode='beg',quackinterval=quackinterval_)

#### Priori Cal Step ####
hanningsmooth=(vis=vis,datacolumn='data',field='',outputvis='20130423_L_T2036-2046_50ms_hs.ms')
vishs='20130423_L_T2036-2046_50ms_hs.ms'

gencal(vis=vishs,caltable=visant,caltype='antpos')
gencal(vis=vishs,caltable=visgc,caltype='gc')
gencal(vis=vishs,caltable=visrq,caltype='rq')

# Phase vs time; for 0...7 channel 20~25 60~65, etc; time ave: 10 sec; iteration over baselines; fix phase axis -180 to 180 
spw1='0:20~25,1:20~25,2:20~25,3:20~25,4:55~60,5:20~25,6:20~25,7:55~60'
refant='ea06'
gaincal(vis=vishs,caltable='20130423_L_T2036-2046_50ms.ph0',field=0,spw=spw1,gaintype='G',gaintable=['20130423_L_T2036-2046_50ms_hs.rq'],refant=refant,calmode='p',solint='10s',minsnr=3)

# phase vs channel; iteration over baselines; field 0; spw 0...7, fix phase axis -180, 180
# delay calibration is done to calculate changes from the ph0 
# delay reults do not change much
#spw2=0:17~19;29~31;46~63,1:0~10;25~30;41~43;60~63,2:0~13,4:10~63,5:10~63,7:20~55 , # for the delay calibration
spw2='0:4~16;20~21;25~28;32~45;51~60,1:3~7;10~61,2:3~7;12~63,3,4:0~14;35~39;43~47;57~61,5:7~11;15~31;35~63,6,7:0~25;37~47;56~63'
# YL: Gaincal output 1/1/1: prepared/tried/success delay calibrations...
# plotcal range: -1,-1,10,10 points displayed on terminal

# ea07: 15 ea23: 12 ea27: 12 ea21: 15 ea25: 15

# Plot delay vs channel for different spw and frequency and see the abrupt jumps.
# Flag:ea23, LL, spw: 4~7
# When we apply same spw as delay to the bandpass, we see less number of solutions in the bandpass, i.e. we revise the spw selection
# 




import sys
import os

def define_TIMES(MSNAME,TIME_STEP):
#returns a list of time slices, taking steps defined in TIME_STEPS 
        ms.open(MSNAME)
        ms.selectinit(datadescid=0)
        tt = ms.getdata('time')
        ms.close()
        times = tt['time']#[::TIME_STEP]
        dt = float(tt['time'][1])-float(tt['time'][0])
        print "#### First time stamp %s" % (tt['time'][0])
        print "#### Integration time %f" % (dt)
        print "#### Number of time stamps %d" % (len(tt['time']))
        print "#### Last time stamp %s" % (tt['time'][len(tt['time'])-1])
        print "#### No. of time stamps to image %d" % (len(times))
        return times
#03:24:33.0 - 03:29:29.0 
#03:29:37.0 - 03:34:33.0 
#03:34:33.0 - 03:39:29.0
#03:39:37.0 - 03:44:33.0
#03:44:33.0 - 03:49:29.0
#03:49:37.0 - 03:54:33.0
#03:54:33.0 - 03:59:29.0 

MSNAME='merged_all_145MHz_chan.ms'
TIME_LIST=define_TIMES(MSNAME,1)
TIME_LIST=TIME_LIST[0]+np.arange(2400)*1.0
#frac=['.0','.5']*2400
i=0
time_list=[0]*(TIME_LIST.shape[0])
time_pair=[0]*(TIME_LIST.shape[0]-1)
imagename=[0]*(TIME_LIST.shape[0]-1)
for i in range(TIME_LIST.shape[0]):
        timed=TIME_LIST[i]-TIME_LIST[0]
        reftime=qa.time(qa.quantity(TIME_LIST[0],'s'),form="ymd")
        reftime_sec=int(reftime[0].split('/')[3].split(':')[0])*3600+int(reftime[0].split('/')[3].split(':')[1])*60+int(reftime[0].split('/')[3].split(':')[2])
        time_list[i]=time.strftime('%H:%M:%S', time.gmtime(np.round(timed+reftime_sec,4)))#+frac[i]

for i in range(TIME_LIST.shape[0]-1):
    time_pair[i]=time_list[i]+'~'+time_list[i+1]

time_pair_mean=[0]*(4560/2)
for i in range(4560/2):
    time_pair_mean[i]=time_list[i]+'~'+time_list[120/2+i]

#time_pair=time_pair[120:592-120]+time_pair[608+120:608+592-120]+time_pair[1200+120:1200+592-120]+time_pair[1808+120:1808+592-120]+time_pair[2400+120:2400+592-120]+time_pair[3008+120:3008+592-120]+time_pair[3600+120:3600+592-120]
#time_pair_mean=time_pair_mean[0:592-240]+time_pair_mean[608:608+592-240]+time_pair_mean[1200:1200+592-240]+time_pair_mean[1808:1808+592-240]+time_pair_mean[2400:2400+592-240]+time_pair_mean[3008:3008+592-240]+time_pair_mean[3600:3600+592-240]
time_pair=time_pair[60:296-60]+time_pair[304+60:304+296-60]+time_pair[600+60:600+296-60]+time_pair[904+60:904+296-60]+time_pair[1200+60:1200+296-60]+time_pair[1504+60:1504+296-60]+time_pair[1800+60:1800+296-60]
time_pair_mean=time_pair_mean[0:296-120]+time_pair_mean[304:304+296-120]+time_pair_mean[600:600+296-120]+time_pair_mean[904:904+296-120]+time_pair_mean[1200:1200+296-120]+time_pair_mean[1504:1504+296-120]+time_pair_mean[1800:1800+296-120]
outputname=[0]*len(time_pair)
for i in range(len(time_pair)):
    ii="%04d"%i
    print 'Time: '+str(ii), time_pair[i], time_pair_mean[i]
    outputname[i]='ms_files/ms_'+str(ii)+'.ms'
    subvs2(vis='merged_all_145MHz_chan.ms',outputvis='ms_files/ms_'+str(ii)+'.ms',timerange=time_pair[i],mode='linear',subtime1=time_pair_mean[i])
    #os.system('rm -rf ee1.ms');os.system('rm -rf ee2.ms');os.system('rm -rf ee3.ms');os.system('rm -rf ee4.ms')
    #split(vis='merged_all.ms',outputvis='ee1.ms',timerange=time_pair[i],datacolumn='data')
    #split(vis='merged_all.ms',outputvis='ee2.ms',timerange=time_pair[i],datacolumn='data')
    #concat(vis=['ee1.ms','ee2.ms'],concatvis='ee3.ms')
    #subvs2(vis='ee3.ms',outputvis='ee4.ms',timerange='',mode='linear',subtime1=time_pair_mean[i])
    #split(vis='ee4.ms',outputvis='ms_files/ms_'+str(ii)+'.ms',timerange=time_pair[i],datacolumn='data')



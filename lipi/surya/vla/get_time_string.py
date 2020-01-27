
def str2sec(s):
    '''
    Input must be in HH:MM:SS format
    '''
    sec=int(s.split(':')[0])*3600+int(s.split(':')[1])*60+float(s.split(':')[2])
    return sec

def get_time_string(start_time,end_time,res):
    ss=str2sec(start_time)
    es=str2sec(end_time)
    n=int((es-ss)/res)
    ts=[0]*n
    ds=[0]*n
    timestr=[0]*n
    time_string=[0]*n
    for i in range(n):
        ts[i]=ss+i*res
        timestr[i]=time.strftime('%H:%M:%S', time.gmtime(ts[i]))
        hour=int(timestr[i].split(':')[0])
        minu=int(timestr[i].split(':')[1])
        sec=int(timestr[i].split(':')[2])
        ds[i]=round(ts[i]-int(ts[i]),4)
        time_string[i]=datetime.datetime(year,mon,day,hour,minu,sec,int(ds[i]*1.e6)).strftime('%H:%M:%S.%f')
    return time_string


####################################################

MSNAME='20130423T2040-2050.L.50ms.selfcal.ms'

year=2013
mon=4
day=23
start_time='20:40:00'
end_time='20:50:00'
res=0.005 # Desired time interval
time_string=get_time_string(start_time,end_time,res)

time_pair=[0]*(len(time_string)-1)
for i in range(len(time_string)-1):
    time_pair[i]=time_string[i]+'~'+time_string[i+1]



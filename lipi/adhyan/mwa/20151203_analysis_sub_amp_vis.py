
import glob
import analysisUtils as au
import pickle
import matplotlib.pyplot as plt

def get_TileID_phase1(TileName):
        """
        From MWA PHASE-I. 
        A dictionary to return the TileID given the TileName as specified 
        in the MS ('Tile???MWA') Use this function rather than using the IDs 
        directly simply to reduce the possibility of human error.
        
        Divya 02Dec2015
        """
        TileDict = {'Tile011MWA' : 0, 'Tile012MWA' : 1, 'Tile013MWA' : 2, 'Tile014MWA' : 3, 'Tile015MWA' : 4, 'Tile016MWA' : 5, 'Tile017MWA' : 6, 'Tile018MWA' : 7, 'Tile021MWA' : 8, 'Tile022MWA' : 9, 'Tile023MWA' : 10, 'Tile024MWA' : 11, 'Tile025MWA' : 12, 'Tile026MWA' : 13, 'Tile027MWA' : 14, 'Tile028MWA' : 15, 'Tile031MWA' : 16,'Tile032MWA':17, 'Tile033MWA' : 18, 'Tile034MWA' : 19, 'Tile035MWA' : 20, 'Tile036MWA' : 21, 'Tile037MWA' : 22, 'Tile038MWA' : 23, 'Tile041MWA' : 24, 'Tile042MWA' : 25, 'Tile043MWA' : 26, 'Tile044MWA' : 27, 'Tile045MWA' : 28, 'Tile046MWA' : 29, 'Tile047MWA' : 30, 'Tile048MWA' : 31, 'Tile051MWA' : 32, 'Tile052MWA' : 33, 'Tile053MWA' : 34, 'Tile054MWA' : 35, 'Tile055MWA' : 36, 'Tile056MWA' : 37, 'Tile057MWA' : 38, 'Tile058MWA' : 39, 'Tile061MWA' : 40, 'Tile062MWA' : 41, 'Tile063MWA' : 42, 'Tile064MWA' : 43, 'Tile065MWA' : 44, 'Tile066MWA' : 45, 'Tile067MWA' : 46, 'Tile068MWA' : 47, 'Tile071MWA' : 48, 'Tile072MWA' : 49, 'Tile073MWA' : 50, 'Tile074MWA' : 51, 'Tile075MWA' : 52, 'Tile076MWA' : 53, 'Tile077MWA' : 54, 'Tile078MWA' : 55, 'Tile081MWA' : 56, 'Tile082MWA' : 57, 'Tile083MWA' : 58, 'Tile084MWA' : 59, 'Tile085MWA' : 60, 'Tile086MWA' : 61, 'Tile087MWA' : 62, 'Tile088MWA' : 63, 'Tile091MWA' : 64, 'Tile092MWA' : 65, 'Tile093MWA' : 66, 'Tile094MWA' : 67, 'Tile095MWA' : 68, 'Tile096MWA' : 69, 'Tile097MWA' : 70, 'Tile098MWA' : 71, 'Tile101MWA' : 72, 'Tile102MWA' : 73, 'Tile103MWA' : 74, 'Tile104MWA' : 75, 'Tile105MWA' : 76, 'Tile106MWA' : 77, 'Tile107MWA' : 78, 'Tile108MWA' : 79, 'Tile111MWA' : 80, 'Tile112MWA' : 81, 'Tile113MWA' : 82, 'Tile114MWA' : 83, 'Tile115MWA' : 84, 'Tile116MWA' : 85, 'Tile117MWA' : 86, 'Tile118MWA' : 87, 'Tile121MWA' : 88, 'Tile122MWA' : 89, 'Tile123MWA' : 90, 'Tile124MWA' : 91, 'Tile125MWA' : 92, 'Tile126MWA' : 93, 'Tile127MWA' : 94, 'Tile128MWA' : 95, 'Tile131MWA' : 96, 'Tile132MWA' : 97, 'Tile133MWA' : 98, 'Tile134MWA' : 99, 'Tile135MWA' : 100, 'Tile136MWA' : 101, 'Tile137MWA' : 102, 'Tile138MWA' : 103, 'Tile141MWA' : 104, 'Tile142MWA' : 105, 'Tile143MWA' : 106, 'Tile144MWA' : 107, 'Tile145MWA' : 108, 'Tile146MWA' : 109, 'Tile147MWA' : 110, 'Tile148MWA' : 111, 'Tile151MWA' : 112, 'Tile152MWA' : 113, 'Tile153MWA' : 114, 'Tile154MWA' : 115, 'Tile155MWA' : 116, 'Tile156MWA' : 117, 'Tile157MWA' : 118, 'Tile158MWA' : 119, 'Tile161MWA' : 120, 'Tile162MWA' : 121, 'Tile163MWA' : 122, 'Tile164MWA' : 123, 'Tile165MWA' : 124, 'Tile166MWA' : 125, 'Tile167MWA' : 126, 'Tile168MWA' : 127}
        return TileDict[TileName]

def get_amp(MSNAME,tile1,tile2,POL_LIST):               # Extraction of amplitude of crosscorrelations
        ms.open(MSNAME)
        casalog.post("#### Tile1 %03d; Tile2 %03d" % (tile1, tile2));
        print '#### Tile1 %03d; Tile2 %03d; POL_LIST %s' % (tile1, tile2, POL_LIST)
        ms.selectinit(datadescid=0) # Reset any earlier selections
        ms.select({'antenna1':[tile1],'antenna2':[tile2]})
        ms.selectpolarization(POL_LIST)
        #amp=ms.getdata(['amplitude'])
        amp=ms.getdata(["data","axis_info"],ifraxis=True)
        ms.close()
        return amp['data']


MSNAME='/nas08-data02/rohit/20151203_sub_run_mean_subamp/running_median/179MHz/merged_all_179MHz_chan.ms'
#datams = mstool()
#ms.open(MSNAME,nomodify=False)
#ms.selectinit(datadescid=0)  
#ms.select({'antenna1':[0,1]}) 
#rec=ms.getdata(["amplitude","phase"],ifraxis=True)
#rec=ms.getdata(["data"],ifraxis=True)
#d=pickle.load(open('data_merged_all_240MHz_chan.p','rb'))
#rec={"data":[]}
#rec['data']=d
#ampl=rec["amplitude"];phas=rec["phase"]

dump_data=0
if(dump_data):
    get_baselines=1
    if(get_baselines):
        #amp=get_amp(MSNAME,0,0,'I')
        bs=au.getBaselineLengths(MSNAME, sort=True)
        bl=[0]*len(bs);amp=[0]*len(bs)
        for i in range(len(bs)):
            base=bs[i][0]
            t1=base.split('-');b1=get_TileID_phase1(t1[0]+'MWA');b2=get_TileID_phase1(t1[1]+'MWA')
            amp[i]=get_amp(MSNAME,b1,b2,'I')
            bl[i]=bs[i][1]
    pickle.dump([amp,np.array(bl),bs],open(MSNAME+'_baseline.p','wb'))
    [amp,bl,bs]=pickle.load(open(MSNAME+'_baseline.p','rb'))

#amps,bl,bs=pickle.load(open('merged_all_108MHz_chan.1.ms_baseline.p','rb'))
#amp,bl,bs=pickle.load(open('merged_all_108MHz_chan.ms_baseline.p','rb'))

amps,bl,bs=pickle.load(open('merged_all_179MHz_chan.1.ms_baseline.p','rb'))
amp,bl,bs=pickle.load(open('merged_all_179MHz_chan.ms_baseline.p','rb'))
re_meds=[0]*len(amps);im_meds=[0]*len(amps);re_medper=[0]*len(amps);im_medper=[0]*len(amps)
re_med=[0]*len(amps);im_med=[0]*len(amps)
for i in range(len(amps)):
    re_med[i]=np.median(amp[i].flatten().real)
    im_med[i]=np.median(amp[i].flatten().imag)
    re_meds[i]=np.median(amps[i].flatten().real)
    im_meds[i]=np.median(amps[i].flatten().imag)
    re_medper[i]=re_med[i]/np.abs(amp[i].flatten().mean().real)
    im_medper[i]=im_med[i]/np.abs(amp[i].flatten().mean().imag)
re_meds=np.array(re_meds)
im_meds=np.array(im_meds)
re_med=np.array(re_med)
im_med=np.array(im_med)
re_medper=np.array(re_medper)
im_medper=np.array(im_medper)
re_hist=np.histogram(re_med,bins=30)
im_hist=np.histogram(im_med,bins=30)
re_hists=np.histogram(re_meds,bins=30)
im_hists=np.histogram(im_meds,bins=30)
reper_hist=np.histogram(re_medper,bins=90000,normed=1)
imper_hist=np.histogram(im_medper,bins=90000,normed=1)
amp_med=np.sqrt(re_med**2 + im_med**2)
amp_meds=np.sqrt(re_meds**2 + im_meds**2)

f,ax=plt.subplots(5,2,figsize=(10,12),sharex=True)
i=27
x0=amp[i].flatten().real[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().real;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[0,0].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax0=ax[0,0].twinx();ax[0,0].legend(loc=3)
ax0.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[0,0].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
x0=amp[i].flatten().imag[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().imag;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[0,1].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax1=ax[0,1].twinx();ax[0,1].legend(loc=3)
ax1.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[0,1].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow');ax[0,0].set_title('Real');ax[0,1].set_title('Imaginary');i=400
x0=amp[i].flatten().real[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().real;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[1,0].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax0=ax[1,0].twinx();ax[1,0].legend(loc=3)
ax0.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[1,0].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
x0=amp[i].flatten().imag[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().imag;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[1,1].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax1=ax[1,1].twinx();ax[1,1].legend(loc=3)
ax1.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[1,1].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow');i=600
x0=amp[i].flatten().real[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().real;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[2,0].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax0=ax[2,0].twinx();ax[2,0].legend(loc=3)
ax0.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[2,0].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
x0=amp[i].flatten().imag[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().imag;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[2,1].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax1=ax[2,1].twinx();ax[2,1].legend(loc=3)
ax1.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[2,1].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow');i=4000
x0=amp[i].flatten().real[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().real;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[3,0].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax0=ax[3,0].twinx();ax[3,0].legend(loc=3)
ax0.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[3,0].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
x0=amp[i].flatten().imag[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().imag;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[3,1].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax1=ax[3,1].twinx();ax[3,1].legend(loc=3)
ax1.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[3,1].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow');i=8000
x0=amp[i].flatten().real[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().real;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[4,0].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax0=ax[4,0].twinx();ax[4,0].legend(loc=3)
ax0.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[4,0].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
x0=amp[i].flatten().imag[0:1776];x0[570:600]=np.nan;x0[1150:1210]=np.nan;x0[1740:]=np.nan
x1=amps[i].flatten().imag;x1[570:600]=np.nan;x1[1150:1210]=np.nan;x1[1740:]=np.nan
ax[4,1].plot(np.arange(1776)*0.5,x0,'o-',markersize=1,color='blue',label=str(np.round(bl[i],2))+' m');ax1=ax[4,1].twinx();ax[4,1].legend(loc=3)
ax1.plot(np.arange(1776)*0.5,x1,'o-',markersize=1,color='red',label=str(np.round(bl[i],2))+' m')
xf=np.arange(len(x0))[np.isfinite(x0)];yf=x0[np.isfinite(x0)];p=np.polyfit(xf,yf,6);nyf=np.polyval(p,np.arange(len(x0)))
ax[4,1].plot(np.arange(len(x0))*0.5,nyf,'-',color='yellow')
ax[0,0].set_ylabel('Amplitude');ax[1,0].set_ylabel('Amplitude');ax[2,0].set_ylabel('Amplitude');ax[3,0].set_ylabel('Amplitude');ax[4,0].set_ylabel('Amplitude');ax[4,0].set_xlabel('Time (sec)');ax[4,1].set_xlabel('Time (sec)')
plt.show()

plt.plot(bl[:8001],abs(re_med),label='Before Subtraction (Real)')
plt.plot(bl[:8001],abs(re_meds),label='After Subtraction (Real)')
plt.plot(bl[:8001],abs(amp_med),label='Before Subtraction')
plt.plot(bl[:8001],abs(amp_meds),label='After Subtraction')
plt.xlabel('Baseline Length (m)');plt.ylabel('Amplitude');plt.legend()
plt.show()

f,ax=plt.subplots(1,2,figsize=(12,6))
h0=np.histogram(amps[27].flatten().real,bins=300)
h1=np.histogram(amps[400].flatten().real,bins=300)
h2=np.histogram(amps[1200].flatten().real,bins=300)
h3=np.histogram(amps[4000].flatten().real,bins=300)
h4=np.histogram(amps[8000].flatten().real,bins=300)
ax[0].plot(h0[1][1:],h0[0],'o-',markersize=4,color='blue',label=str(np.round(bl[27],2))+' m')
ax[0].plot(h1[1][1:],h1[0],'o-',markersize=4,color='cyan',label=str(np.round(bl[400],2))+' m')
ax[0].plot(h2[1][1:],h2[0],'o-',markersize=4,color='green',label=str(np.round(bl[1200],2))+' m')
ax[0].plot(h3[1][1:],h3[0],'o-',markersize=4,color='yellow',label=str(np.round(bl[4000],2))+' m')
ax[0].plot(h4[1][1:],h4[0],'o-',markersize=4,color='red',label=str(np.round(bl[8000],2))+' m')
ax[0].legend();ax[0].set_yscale('log');ax[0].set_xlim(-20,20)
h0=np.histogram(amps[27].flatten().imag,bins=300)
h1=np.histogram(amps[400].flatten().imag,bins=300)
h2=np.histogram(amps[1200].flatten().imag,bins=300)
h3=np.histogram(amps[4000].flatten().imag,bins=300)
h4=np.histogram(amps[8000].flatten().imag,bins=300)
ax[1].plot(h0[1][1:],h0[0],'o-',markersize=4,color='blue',label=str(np.round(bl[27],2))+' m')
ax[1].plot(h1[1][1:],h1[0],'o-',markersize=4,color='cyan',label=str(np.round(bl[400],2))+' m')
ax[1].plot(h2[1][1:],h2[0],'o-',markersize=4,color='green',label=str(np.round(bl[1200],2))+' m')
ax[1].plot(h3[1][1:],h3[0],'o-',markersize=4,color='yellow',label=str(np.round(bl[4000],2))+' m')
ax[1].plot(h4[1][1:],h4[0],'o-',markersize=4,color='red',label=str(np.round(bl[8000],2))+' m')
ax[1].legend();ax[1].set_yscale('log');ax[1].set_xlim(-20,20)
ax[0].set_xlabel('Real Amplitude');ax[1].set_xlabel('Imaginary Amplitude')
ax[0].set_ylabel('Occupancy');ax[1].set_ylabel('Occupancy')
plt.show()

plt.plot(re_hists[1][1:],re_hists[0],'o-',markersize=4,color='red',label='Real')
plt.plot(im_hists[1][1:],im_hists[0],'o-',markersize=4,color='blue',label='Imaginary')
plt.legend();plt.xlabel('')
plt.xlabel('Magnitude');plt.ylabel('Occupancy')
plt.show()
plt.plot(re_hist[1][1:],re_hist[0],'o-',markersize=4,color='red',label='Real')
plt.plot(im_hist[1][1:],im_hist[0],'o-',markersize=4,color='blue',label='Imaginary')
plt.legend();plt.xlabel('')
plt.xlabel('Magnitude');plt.ylabel('Occupancy')
plt.show()

plt.plot(reper_hist[1][1:],reper_hist[0],'o-',markersize=4,color='red',label='Real')
plt.plot(imper_hist[1][1:],imper_hist[0],'o-',markersize=4,color='blue',label='Imaginary')
plt.legend();plt.xlabel('')
plt.xlim(-0.05,0.05);plt.xlabel('Fraction deviation of medians (%)');plt.ylabel('Occupancy')
plt.show()

plt.hist(np.array(im_med),bins=1000,histtype='step')
plt.hist(np.array(im_med),bins=1000,histtype='step')


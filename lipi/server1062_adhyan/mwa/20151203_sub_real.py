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
        amp=ms.getdata(["amplitude","axis_info"],ifraxis=True)
        ms.close()
        return amp['amplitude']

def RunningMedian(x,N):
    idx = np.arange(N) + np.arange(len(x)-N+1)[:,None]
    b = [row for row in x[idx]]
    return np.array(map(np.median,b))

#MSNAME='merged_all_108MHz_chan.1.sub.10.ms'
MSNAME='merged_all_240MHz_chan.2.ms'
get_baselines=0
if(get_baselines):
    #amp=get_amp(MSNAME,0,0,'I')
    bs=au.getBaselineLengths(MSNAME, sort=True)
    bl=[0]*len(bs);amp=[0]*len(bs)
    for i in range(len(bs)):
        base=bs[i][0]
        t1=base.split('-');b1=get_TileID_phase1(t1[0]+'MWA');b2=get_TileID_phase1(t1[1]+'MWA')
        amp[i]=get_amp(MSNAME,b1,b2,'I')
        bl[i]=bs[i][1]

list_base=[57,59,158,159,160,169,170,171,1347,1348,1349,1350,1351,1352,1353,1354,1462,1463,1465,1466,1467,1468,1469]
datams = mstool()
ms.open(MSNAME,nomodify=False)  
#ms.selectinit(datadescid=0)  
#ms.select({'antenna1':[0,1]}) 
rec=ms.getdata(['data'],ifraxis=True)
#d=pickle.load(open('data_merged_all_240MHz_chan.p','rb'))
#rec={"data":[]}
#rec['data']=d
real=rec["data"].real;imag=rec["data"].imag
#pickle.dump(ampl,open('amp_before.p','wb'))
nb=8128;nt=592;n=3;error_=0
NN=31;sig=2.5
#nb=7875;nt=592;n=1
for i in range(real.shape[0]):
    for j in range(real.shape[1]):
        for k in range(0,real.shape[2]):
            yr=real[i][j][k];yi=imag[i][j][k]
            #y_=y.real.reshape((2,592));yamp=[0]*2
            yr_=yr.reshape((n,nt));yi_=yi.reshape((n,nt));yrnew=[0]*n;yinew=[0]*n
            for l in range(n):
                xr=np.arange(len(yr_[l]))
                x0=xr;y0=yr_[l]
                if(error_):
                    y0[np.where((y0<(np.nanmean(y0)-1.5*np.nanstd(y0))) | (y0>(np.nanmean(y0)+1.5*np.nanstd(y0))))]=y0[np.where((y0<(np.nanmean(y0)-1.5*np.nanstd(y0))) | (y0>(np.nanmean(y0)+1.5*np.nanstd(y0))))]*0
                yr0=y0[np.where((y0>(np.nanmean(y0)-sig*np.nanstd(y0))) & (y0<(np.nanmean(y0)+sig*np.nanstd(y0))))]
                xr0=x0[np.where((y0>(np.nanmean(y0)-sig*np.nanstd(y0))) & (y0<(np.nanmean(y0)+sig*np.nanstd(y0))))] # 0 --> Throw away the Outliers
                yr1=y0[np.where((y0<(np.nanmean(y0)-sig*np.nanstd(y0))) | (y0>(np.nanmean(y0)+sig*np.nanstd(y0))))] # 1 --> The Outliers
                xr1=x0[np.where((y0<(np.nanmean(y0)-sig*np.nanstd(y0))) | (y0>(np.nanmean(y0)+sig*np.nanstd(y0))))]
                #xr0=xr[yr_[l]!=0];yr0=yr_[l][yr_[l]!=0]
                if(len(xr0)>3):
                    print i,j,k,n#p=np.polyfit(xr0,yr0,1)#yrnew[l]=np.polyval(p,x0)
                    yrnew[l]=RunningMedian(yr0,NN) # Fit for the curve without outliers
                    for ii in range(len(yrnew[l])):
                        yr_[l][xr0[ii]+int(NN/2)]=yr0[int(NN/2):-int(NN/2)][ii]-yrnew[l][ii] # subtract real from yr0, i.e. without outliers
                    for ii in range(int(NN/2)):
                        yr_[l][xr0[ii]]=0;yr_[l][xr0[-int(NN/2)+ii]]=0
                else:
                    yr_[l]=yr_[l]*0
                if(len(xr1)>0):
                    yr_[l][xr1]=yr_[l][xr1]-np.median(yr0)
                xi=np.arange(len(yi_[l]))
                x0=xi;y0=yi_[l]
                if(error_):
                    y0[np.where((y0<(np.nanmean(y0)-1.5*np.nanstd(y0))) | (y0>(np.nanmean(y0)+1.5*np.nanstd(y0))))]=0*y0[np.where((y0<(np.nanmean(y0)-1.5*np.nanstd(y0))) | (y0>(np.nanmean(y0)+1.5*np.nanstd(y0))))]
                #xi0=xi[yi_[l]!=0];yi0=yi_[l][yi_[l]!=0]
                yi0=y0[np.where((y0>(np.nanmean(y0)-sig*np.nanstd(y0))) & (y0<(np.nanmean(y0)+sig*np.nanstd(y0))))]
                xi0=x0[np.where((y0>(np.nanmean(y0)-sig*np.nanstd(y0))) & (y0<(np.nanmean(y0)+sig*np.nanstd(y0))))]
                yi1=y0[np.where((y0<(np.nanmean(y0)-sig*np.nanstd(y0))) | (y0>(np.nanmean(y0)+sig*np.nanstd(y0))))]
                xi1=x0[np.where((y0<(np.nanmean(y0)-sig*np.nanstd(y0))) | (y0>(np.nanmean(y0)+sig*np.nanstd(y0))))]
                if(len(xi0)>3):
                    #p=np.polyfit(xi0,yi0,1)
                    #yinew[l]=np.polyval(p,x0)
                    yinew[l]=RunningMedian(yi0,NN)
                    for ii in range(len(yinew[l])):
                        yi_[l][xi0[ii]+int(NN/2)]=yi0[int(NN/2):-int(NN/2)][ii]-yinew[l][ii]
                    for ii in range(int(NN/2)):
                        yi_[l][xi0[ii]]=0;yi_[l][xi0[-int(NN/2)+ii]]=0
                else:
                    yi_[l]=yi_[l]*0
                if(len(xi1)>0):
                    yi_[l][xi1]=yi_[l][xi1]-np.median(yi0)
            #yinew=np.hstack(np.array(yinew));yrnew=np.hstack(np.array(yrnew))
            yi=yi_.flatten()
            yr=yr_.flatten()
            real[i][j][k]=yr
            imag[i][j][k]=yi
            if(0):
                if(k in list_base):
                    plt.plot(yr,'o',markersize=5,color='b')
                    plt.plot(yrnew,'-',color='b')
                    plt.plot(yi,'o',markersize=5,color='r')
                    plt.plot(yinew,'-',color='r')
                    plt.plot(yr,'s',markersize=5,color='b')
                    plt.plot(yi,'s',markersize=5,color='r')
                    plt.title('i='+str(i)+' j='+str(j)+' k='+str(k))
                    plt.savefig('pngs/sub_'"%01d"%i+'_'+"%04d"%k+'.png')
                    plt.close()
#real=real.swapaxes(2,3).reshape(4,1,nb*nt*n)
#imag=imag.swapaxes(2,3).reshape(4,1,nb*nt*n)
#real=real.reshape(4,1,nb*nt*n) # sub2
#imag=imag.reshape(4,1,nb*nt*n)
real1=real.swapaxes(2,3).reshape(4,1,nb*nt*n)
imag1=imag.swapaxes(2,3).reshape(4,1,nb*nt*n)
ndata=real1*0j
#ndata=ndata.reshape(4,1,8128*592*2)
for i in range(4):
    for j in range(nb*nt*n):
        ndata[i][0][j]=np.complex(real1[i][0][j],imag1[i][0][j])
nrec={"data":None};nrec["data"]=ndata
ms.open(MSNAME,nomodify=False)
ms.putdata(nrec)  

sys.exit()
pickle.dump([amp,np.array(bl),bs],open('baselines.p','wb'))
[amp,bl,bs]=pickle.load(open('baselines.p','rb'))


amp_=amp[0][0][0][0].reshape((7,592))
ydiff=[0]*7
for i in range(7):
    y=amp_[i]
    x=np.arange(len(y))
    x0=x[y>0]
    y0=y[y>0]
    p=np.polyfit(x0,y0,3)
    ynew=np.polyval(p,x0)
    ydiff[i]=y0-ynew
    #plt.plot(x0,y0,'-')
    #plt.plot(x0,ynew,'-')
    plt.plot(x0,ydiff[i],'-')
    plt.show()


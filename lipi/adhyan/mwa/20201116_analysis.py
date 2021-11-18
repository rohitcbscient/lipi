import numpy as np
import pickle
import glob
import matplotlib.pyplot as plt
import sys

make_pickle=0
if(make_pickle):
    mslist=sorted(glob.glob('*.ms'))
    for MSNAME in mslist:
        print MSNAME
        ms.open(MSNAME)
        ms.selectinit(datadescid=0)
        amp1=ms.getdata(["amplitude","axis_info"],ifraxis=True)
        amp=amp1['amplitude']
        pickle.dump(amp,open(MSNAME+'.amp.p','wb'))


#sys.exit()
ids=['062','067','073','078','084','089','095','100','106','111','117','122','128','133','139','144','150','155','161','166','172','177','183','188']
freq=np.array([  79.36,   85.76,   93.44,   99.84,  107.52,  113.92,  121.6 ,128.  ,  135.68,  142.08,  149.76,  156.16,  163.84,  170.24,177.92,  184.32,  192.  ,  198.4 ,  206.08,  212.48,  220.16,226.56,  234.24,  240.64])
bmedian=[0]*len(ids);bmed=[0]*len(ids)
j=0
for ii in ids:
    filelist=sorted(glob.glob('/media/rohit/MWA/MWA_STIX_2020/0600-0700/'+ii+'/*_'+str(ii)+'-'+str(ii)+'.ms.amp.p'))
    print ii,len(filelist)
    i=0;bmedian[j]=[0]*len(filelist);bmed[j]=[0]*len(filelist)
    for f in filelist:
        amp=pickle.load(open(f,'rb'))
        bmedian[j][i]=np.median(amp[:,:,1000:,:],axis=2)
        #bmed[j][i]=bmedian[j][i].mean(axis=(0,1))
        bmed[j][i]=bmedian[j][i].mean(axis=0)[0]
        i=i+1
    bmedian[j]=np.array(bmedian[j]);bmed[j]=np.array(bmed[j])
    j=j+1
bmedian=np.array(bmedian)
bmed=np.array(bmed)
bmed[:,:,-1]=np.nan;bmed[:,:,0]=np.nan
bmed2d=bmed.reshape(24,13*30)
for i in range(24):
    bmed2d[i]=bmed2d[i]-bmed2d[i][122]

bmed2d_sub=bmed2d[:,90:300]
plt.imshow(bmed2d_sub,aspect='auto',origin=0,interpolation='none',vmin=0,vmax=0.03,cmap='YlGnBu')
plt.show()

filelist_auto=sorted(glob.glob('*T002-007.DS.dat.p'))
i=0;ac=[0]*len(filelist_auto)
for f in filelist_auto:
    aa=pickle.load(open(f,'rb'))
    ac[i]=aa[2][0][0][0]
    i=i+1
ac=np.array(ac)
#nccf=np.array(nccf).reshape(8,3,2,30)
#nccf_bmean=np.mean(nccf,axis=1)
ts=nccf[:,0,:]
plt.plot(ts.flatten(),'o-')
plt.show()



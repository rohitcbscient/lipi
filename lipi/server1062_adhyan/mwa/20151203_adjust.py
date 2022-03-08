import pickle
import numpy as np

base=['000-008','000-009','000-010','008-009','008-010','009-010']
flist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
n=7

for b in base:
    for f in flist:
        filename='20151203_'+f+'_T'+b+'.p'
        data=pickle.load(open(filename,'r'))
        Tb=data[17][3][0]
        nTb=[0]*64
        for i in range(64):
            nTb[i]=[0]*n
            nTb[i][0]=Tb[i,0:592]
            for j in range(n-1):
                fact=np.nanmean(nTb[i][j][560:592])
                #print fact
                nTb[i][j+1]=Tb[i,(j+1)*592:(j+2)*592]+fact-np.nanmean(Tb[i,(j+1)*592:(j+1)*612])
        nTb=np.array(nTb).reshape(64,7*592)
        data[17][3][0]=nTb
        pickle.dump(data,open('new_pickle/20151203_'+f+'_T'+b+'.p','w'))


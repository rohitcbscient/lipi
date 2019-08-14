import pickle
import numpy as np
import glob

def read_gaussian(n,datafile,p,start_chan,end_chan,data_type):
        bic=pickle.load(open(str(datafile)+'_bic_'+str(n),'r'))
        resp=pickle.load(open(str(datafile)+'_res_'+str(n),'r'))
        param=pickle.load(open(str(datafile)+'_param_'+str(n),'r'))
        return bic, resp, param

def read_n_comp(datafile):
        bic_list=sorted(glob.glob(str(datafile)+'_bic*'))
        bn=np.zeros(len(bic_list))
        for i in range(len(bic_list)):
                bn[i]=pickle.load(open(str(bic_list[i]),'r'))
        n_compp=np.where(bn==min(bn))[0][0]+1
        return n_compp, len(bic_list)


import numpy as np
import matplotlib.pyplot as plt
import os
import ms


def get_listobs(vis_):
    listfile_=vis_.split('.ms')[0]+'.listobs'
    if(os.path.isfile(listfile_)==False):
        listobs(vis=vis_,listfile=listfile_)
        print listfile_+" created"

def get_plotants(vis_):
    figfile_=vis.split('.ms')[0]+'.plotants.png'
    if(os.path.isfile(figfile_)==False):
        plotants(vis=vis_,antindex=True,logpos=True,checkbaselines=True) 
        print figfile_+' created'
        plt.close()

def get_amp_phase(vis_,n0,n1,nspw):
    ms.open(vis_)
    print 'Getting amp-phase for '+str(n0)+' & '+str(n1)
    ms.select({'antenna1':[n0],'antenna2':[n1]})
    amp=ms.getdata(['amplitude'])['amplitude']
    phase=ms.getdata(['phase'])['phase']
    ms.close()
    amp=amp.reshape(amp.shape[0],amp.shape[1],nspw,amp.shape[2]/nspw)
    phase=phase.reshape(phase.shape[0],phase.shape[1],nspw,phase.shape[2]/nspw)
    return amp,phase

#vis='2050.1s.cal.ms'
#ap01=get_amp_phase(vis,0,1,8)        
#ap02=get_amp_phase(vis,0,2,8)        
#ap03=get_amp_phase(vis,1,2,8)        
#ap03=get_amp_phase(vis,1,2,8)        
#ap03=get_amp_phase(vis,1,2,8)        




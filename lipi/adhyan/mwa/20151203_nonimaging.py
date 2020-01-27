import sys
from surya.radio import main as rd

MS='1133149192-%b187-188.MS'
ms.open(MS)
nspw=64
n0=0
n1=8
n2=9
n3=10

ms.open(MS)
ms.select({'antenna1':[n0],'antenna2':[n1]})
amp08=ms.getdata(['amplitude'])['amplitude']
phase08=ms.getdata(['phase'])['phase']
ms.close()
ms.open(MS)
ms.select({'antenna1':[n0],'antenna2':[n2]})
amp09=ms.getdata(['amplitude'])['amplitude']
phase09=ms.getdata(['phase'])['phase']
ms.close()
ms.open(MS)
ms.select({'antenna1':[n1],'antenna2':[n2]})
amp89=ms.getdata(['amplitude'])['amplitude']
phase89=ms.getdata(['phase'])['phase']
ms.close()
ms.open(MS)
ms.select({'antenna1':[n0],'antenna2':[n3]})
amp010=ms.getdata(['amplitude'])['amplitude']
phase010=ms.getdata(['phase'])['phase']
ms.close()
ms.open(MS)
ms.select({'antenna1':[n1],'antenna2':[n3]})
amp810=ms.getdata(['amplitude'])['amplitude']
phase810=ms.getdata(['phase'])['phase']
ms.close()
ms.open(MS)
ms.select({'antenna1':[n2],'antenna2':[n3]})
amp910=ms.getdata(['amplitude'])['amplitude']
phase910=ms.getdata(['phase'])['phase']
ms.close()

phase_closure=phase89+phase08-phase09
amp_closure=amp08*amp910/(amp09*amp810)
ampc=amp_closure[0].flatten()[np.isfinite(amp_closure[0].flatten())]


print 'AMP STD: ',np.std(amp08[0].flatten()),np.std(amp09[0].flatten()),np.std(amp89[0].flatten())
print 'PHASE STD: ',np.std(phase08[0].flatten()),np.std(phase09[0].flatten()),np.std(phase89[0].flatten())
print 'PHASE MEAN: ',np.mean(phase08[0].flatten()),np.mean(phase09[0].flatten()),np.mean(phase89[0].flatten())
print 'PHASE CLOSURE: ',np.mean(phase_closure[0].flatten()),np.std(phase_closure[0].flatten())
#amp=amp.reshape(amp.shape[0],amp.shape[1],nspw,amp.shape[2]/nspw)
#phase=phase.reshape(phase.shape[0],phase.shape[1],nspw,phase.shape[2]/nspw)
#amp,phase=rd.get_amp_phase(MS,n0,n1,nspw)


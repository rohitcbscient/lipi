#from sklearn.mixture import BayesianGaussianMixture as GMM
from sklearn.mixture import GaussianMixture as GMM
import pickle
import matplotlib.pyplot as plt
import numpy as np

def gm(data):
	k=10
	bici=1.e8+np.zeros(k)
	X = data
	bic = []
	n_components_range = range(1, k+1)
	for n_components in n_components_range:    
		#g = GMM(covariance_type='full', init_params='random', min_covar=min_coverr,n_components=n_components, n_init=1, n_iter=5000, params='wmc', random_state=1,thresh=0.0001,data_file=datafile)
                g = GMM(covariance_type='full', init_params='random',n_components=n_components, n_init=1, max_iter=5000, random_state=1)
		g.fit(X)
		bici[n_components]=g.bic(X)
		bic.append(g.bic(X))
                print n_components,bici[n_components]
		print 'Convergence state: ',g.converged_
		if (bici[n_components-1] <= bici[n_components]):
			print 'The minimum bic reached='+str(min(bici))+'_at ncomp='+str(np.where(bici==min(bici))[0][0])
			break
		else:	
			continue	
	print 'Convergence state: ',g.converged_
        gg = GMM(covariance_type='full', init_params='random',n_components=np.where(bici==min(bici))[0][0], n_init=1, max_iter=5000, random_state=1)
        gg.fit(X)
        m=gg.means_;c=np.sqrt(gg.covariances_);w=gg.weights_
	return np.where(bici==min(bici))[0][0],bic,m,c,w

data=pickle.load(open('/nas08-data02/rohit/20151203_sub_run/Tb_20151203_125-126_sub.p','rb'))
Tb=data[0]
#datain=np.vstack((Tb[0].flatten(),Tb[0].flatten()))
datain=Tb[22].flatten().reshape(-1,1)
print 'Fitting...'
b=gm(datain)



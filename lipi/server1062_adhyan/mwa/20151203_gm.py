import numpy as np
from sklearn.mixture import BayesianGaussianMixture as GMM
import os

def gm(data,min_val):
	k=10
	bici=1.e8+np.zeros(k)
	X = data
	X=X[np.where(np.greater(X,min_val))]
	bic = []
	n_components_range = range(1, k+1)
	for n_components in n_components_range:    
		#g = GMM(covariance_type='full', init_params='random', min_covar=min_coverr,n_components=n_components, n_init=1, n_iter=5000, params='wmc', random_state=1,thresh=0.0001,data_file=datafile)
                g = GMM(covariance_type='full', init_params='random',n_components=n_components, n_init=1, max_iter=5000, random_state=1)
		g.fit(X)
		bici[n_components]=g.bic(X)
		bic.append(g.bic(X))
		print 'Convergence state: ',g.converged_
		if (bici[n_components-1] <= bici[n_components]):
			print 'The minimum bic reached='+str(min(bici))+'_at ncomp='+str(np.where(bici==min(bici))[0][0])
			break
		else:	
			continue	
	#print 'Convergence state: ',g.converged_	
	allbic=(np.where(bici==min(bici))[0][0]+1)*[0]
	return np.where(bici==min(bici))[0][0],allbic

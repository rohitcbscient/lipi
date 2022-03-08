import pickle

vlacube='/media/rohit/VLA/20120225_cube/'
fcube=open(vlacube+'/centroid_20120225.p','rb')
cube_cent=pickle.load(fcube)
fcube.close()

pickle.load(open('tp.p','r'))
vlasubvs='/media/rohit/VLA/20120225_sub/fits/'
fsub=open(vlasubvs+'/centroid_20120225.p','rb')
subvs_cent=pickle.load(fsub)
fsub.close()
plt.plot(subvs_cent[1],cube_cent[1],'o')
plt.show()


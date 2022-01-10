import sunpy
import matplotlib.pyplot as plt
from sunpy.map import Map
import numpy as np
import plotly.express as px
import pickle
import glob

list_=sorted(glob.glob('map*.p'))
aia=[0]*60;aiab=[0]*60;aiar=[0]*60
for i in range(60):
    aia_,aiab_,aiar_=pickle.load(open(list_[i*5],'rb'))
    aia[i],aiab[i],aiar[i]=aia_.data,aiab_.data,aiar_.data
aia=np.array(aia)[:,:2048,2048:];aiab=np.array(aiab)[:,:2048,2048:];aiar=np.array(aiar)[:,:2048,2048:]

fig = px.imshow(aiar, animation_frame=0, color_continuous_scale='RdBu_r', origin='lower',zmin=-30,zmax=30,labels=dict(animation_frame="Time"))
fig.show()




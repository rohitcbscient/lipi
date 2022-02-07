import numpy as np
import os

lay=np.loadtxt('alma.tm/layout.txt')
station_num=lay.shape[0]
name='alma1'
os.system('mkdir '+name+'.tm')
os.system('cp alma.tm/layout.txt '+name+'.tm')
for i in range(station_num):
    ii="%03d" %i
    os.system('mkdir '+name+'.tm/station'+ii)
    os.system('cp layout_station.txt '+name+'.tm/station'+ii+'/layout.txt')
    #np.savetxt(name+'.tm/station'+ii+'/layout.txt',np.array([0,0]).T)



import numpy as np

def hms2sec_c(time):
    sec=int(time.split(' ')[1].split(':')[0])*3600+int(time.split(' ')[1].split(':')[1])*60+float(time.split(' ')[1].split(':')[2])
    return sec

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]





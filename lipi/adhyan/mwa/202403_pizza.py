import numpy as np

path='/data/rohit/pizza-project/201503/10/240MHz/'

msname = path+'1109993144_ch187-188.ms'
listfile = path+'1109993144_ch187-188.listobs'
#--- Pic A model
calname = path+'1110007544_ch187-188_cal.ms'
calfile = path+'1110007544_ch187-188_cal.listobs'
listobs(msname,listfile=listfile,overwrite=True)
listobs(calname,listfile=calfile,overwrite=True)


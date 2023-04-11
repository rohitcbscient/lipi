import numpy as np
import tools21cm as t2c


x_file = t2c.XfracFile('/sdata/ip6_project_2023/xfrac3d_6.000.bin')
d_file = t2c.DensityFile('/sdata/ip6_project_2023/7.059n_all.dat')

xfrac = x_file.xi
dens  = d_file.cgs_density



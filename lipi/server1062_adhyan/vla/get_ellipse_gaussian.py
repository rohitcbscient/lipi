import numpy as np
import matplotlib.pyplot as plt
import json, glob
import sys, time, os
from scipy.stats import multivariate_normal
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from scipy.constants import *

def integrate_simps (mesh, func):
    nx, ny = func.shape
    px, py = mesh[0][int(nx/2), :], mesh[1][:, int(ny/2)]
    val = simps( simps(func, px), py )
    return val

def normalize_integrate (mesh, func):
    return func / integrate_simps (mesh, func)

def moment (mesh, func, index):
    ix, iy = index[0], index[1]
    g_func = normalize_integrate (mesh, func)
    fxy = g_func * mesh[0]**ix * mesh[1]**iy
    val = integrate_simps (mesh, fxy)
    return val

def moment_seq (mesh, func, num):
    seq = np.empty ([num, num])
    for ix in range (num):
        for iy in range (num):
            seq[ix, iy] = moment (mesh, func, [ix, iy])
    return seq

def get_centroid (mesh, func):
    dx = moment (mesh, func, (1, 0))
    dy = moment (mesh, func, (0, 1))
    return dx, dy

def get_covariance (mesh, func, dxy):
    g_mesh = [mesh[0]-dxy[0], mesh[1]-dxy[1]]
    Mxx = moment (g_mesh, func, (2, 0))
    Myy = moment (g_mesh, func, (0, 2))
    Mxy = moment (g_mesh, func, (1, 1))
    return np.array([[Mxx, Mxy], [Mxy, Myy]])

def plot_contour_sub (mesh, func, loc=[0, 0], title="name", pngfile="./name"):
    sx, sy = loc
    nx, ny = func.shape
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    mx, my = int ( (sy-ys)/dy ), int ( (sx-xs)/dx )
    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    ax.set_aspect('equal')
    ax_x = divider.append_axes("bottom", 1.0, pad=0.5, sharex=ax)
    ax_x.plot (mesh[0][mx, :], func[mx, :])
    ax_x.set_title ("y = {:.2f}".format(sy))
    ax_y = divider.append_axes("right" , 1.0, pad=0.5, sharey=ax)
    ax_y.plot (func[:, my], mesh[1][:, my])
    ax_y.set_title ("x = {:.2f}".format(sx))
    im = ax.contourf(mesh, func, cmap="jet")
    ax.set_title (title)
    plt.colorbar (im, ax=ax, shrink=0.9)
    plt.savefig(pngfile + ".png")

def make_gauss (mesh, sxy, rxy, rot):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = np.exp (-0.5 * (px/rxy[0])**2)
    fy = np.exp (-0.5 * (py/rxy[1])**2)
    return fx * fy

if __name__ == "__main__":
    argvs = sys.argv  
    argc = len(argvs)
    print (argvs)

    nx, ny = 500, 500
    lx, ly = 500, 500
    rx, ry = 40, 25 # Gaussian width
    sx, sy = 50, 10 # Gaussian center
    rot    = 30
    H=1.e7

    px = np.linspace (-1, 1, nx) * lx
    py = np.linspace (-1, 1, ny) * ly
    mesh = np.meshgrid (px, py)
    fxy0 = make_gauss (mesh, [sx, sy], [rx, ry], np.deg2rad(rot))*H 
    s1xy = get_centroid (mesh, fxy0)
    plot_contour_sub (mesh, fxy0, loc=s1xy, title="Original", pngfile="./fxy0")


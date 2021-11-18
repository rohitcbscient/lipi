import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import argrelmax

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle

# 2D Gaussian model
def func(xy, [x0, y0, sigma0, H0],[x1,y1,sigma1,H1]):
    x, y = xy
    A0 = 1.0 / (2 * sigma0**2);A1 = 1.0 / (2 * sigma1**2)
    I = H0 * np.exp(-A0 * ( (x - x0)**2 + (y - y0)**2)) + H1 * np.exp(-A1 * ( (x - x1)**2 + (y - y1)**2))
    return I

# Generate 2D gaussian
def generate(x0, y0, sigma, H):
    x = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    y = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    xx, yy = np.meshgrid(x, y)
    I = func((xx, yy), x0=x0, y0=y0, sigma=sigma, H=H)
    return xx, yy, I

def fit(image,[x0,y0,sigma0,H0],[x1,y1,sigma1,H1]):

    # Prepare fitting
    x = np.arange(0, image.shape[1], 1)
    y = np.arange(0, image.shape[0], 1)
    xx, yy = np.meshgrid(x, y)
    # Guess intial parameters
    initial_guess0 = [x0, y0, sigma0, H0]
    initial_guess1 = [x1, y1, sigma1, H1]
    lower = [0, 0, 0, 0]
    upper = [image.shape[0], image.shape[1], max(*image.shape), image.max() * 2]
    bounds = [lower, upper]
    # FIT
    pred_params, uncert_cov = curve_fit(func, (xx.ravel(), yy.ravel()), image.ravel(),
                                        p0=initial_guess, bounds=bounds)

    # Get residual
    predictions = func((xx, yy), *pred_params)
    rms = np.sqrt(np.mean((image.ravel() - predictions.ravel())**2))

    print("True params : ", true_parameters)
    print("Predicted params : ", pred_params)
    print("Residual : ", rms)

    return pred_params,rms,uncert_cov

def plot(image, params):

    fig, ax = plt.subplots()
    ax.imshow(image, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower')

    ax.scatter(params[0], params[1], s=100, c="red", marker="x")

    circle = Circle((params[0], params[1]), params[2], facecolor='none',
            edgecolor="red", linewidth=1, alpha=0.8)
    ax.add_patch(circle)


xcimax,ycimax,xci90,yci90,maxTbi,Tbi_r1,Tbi_r2,areai50,eTbi=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_i.p','rb'))
Tbi_r1=np.array(Tbi_r1).reshape(32,2000);Tbi_r2=np.array(Tbi_r2).reshape(32,2000);max_Tbi=np.array(maxTbi).reshape(32,2000)
xcimax=np.array(xcimax).reshape(32,2000);ycimax=np.array(ycimax).reshape(32,2000);areai50=np.array(areai50).reshape(32,2000)
xcvmax,ycvmax,xcv90,ycv90,maxTbv,Tbv_r1,Tbv_r2,eTbv=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_v.p','rb'))
xclmax,yclmax,xcl90,ycl90,maxTbl,Tbl_r1,Tbl_r2,eTbl=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_l.p','rb'))
xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r.p','rb'))
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs.p','rb'))
qsx,qsy,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50=np.array(qsx),np.array(qsy),np.array(qsmaxTbr),np.array(qsTbr_r1),np.array(qsTbr_r2),np.array(qsarear50)
ds_=aa=np.load('/media/rohit/VLA/20160409/sun_L_20160409.1s.ms.dspec.npz')
#timevla=allmaps['vla']['timevla']
#time94=allmaps['aia94']['time94'];time131=allmaps['aia131']['time131'];time335=allmaps['aia335']['time335']
#time1600=allmaps['aia1600']['time1600'];time1700=allmaps['aia1700']['time1700']
freq=np.round(np.linspace(0.997,2.005,32),3)


# The fit performs well without bounds
params,rms,uncert_cov = fit(image, with_bounds=True)
plot(image, params)


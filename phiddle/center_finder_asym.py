import glob
import time

from scipy.interpolate import CubicSpline
import numpy as np
import copy as cp
from astropy.modeling import models, fitting


def get_center_asym(data, left_fwhm = 1., right_fwhm = None, window = 21, p=1/4):
    im1 = np.array(data[:,:])
    im = cp.deepcopy(im1)
    if right_fwhm is not None:
        scale = right_fwhm/left_fwhm # < 0
    else:
        scale = 1.
    x = np.array(range(im1.shape[1]))
    correlation = np.zeros(im1.shape[1])
    smoothed = im1
    smoothed = smoothed - np.mean(smoothed, axis=-1).reshape(im.shape[0], 1) 
    smoothed = smoothed/np.max(smoothed)

    for i_center in range(max(1, int(im1.shape[1]*p)), int(im1.shape[1]*(1-p))):
        left = smoothed[:, :i_center]
        right = smoothed[:, i_center:]
        cs = CubicSpline(x[i_center:], right, axis = -1) # interpolate rows
        min_x = np.min(x[i_center:])

        maxval_left = np.max(np.abs(left))

        right_scaled = cs(scale * (x[i_center:]-min_x) + min_x) 
        maxval_right = np.max(np.abs(right_scaled))
        right_scaled = right_scaled / maxval_right * maxval_left

        maxlength = min(left.shape[-1], right.shape[-1])

        left = left[:, i_center-maxlength:i_center]
        right = right_scaled[:, maxlength-1::-1]

        correlation[i_center] += np.sum(left * right) / (maxlength**2)

    wmin = int(round((correlation.shape[0] - window)*0.5))
    wmin = int(len(correlation) * p)
    wmax = int(len(correlation) * (1-p))
    filt = np.zeros(correlation.shape)
    # filt[wmin:correlation.shape[0]-wmin] = 1.
    filt[wmin:wmax] = 1.
    imaximum_conv = np.argmax(np.multiply(correlation, filt)) - 1

    return imaximum_conv  #imaximum excluded because it was throwing error

#@jit(nopython = False)
def fft_smoothing(d,param):
#d is a 1D vector, param is the cutoff frequency
    rft = np.fft.rfft(d)
    rft[:, int(param):] = 0.
    d = np.fft.irfft(rft, d.shape[1])
    return d

# 
# if __name__ == "__main__":
# 
#     for fn in data_files:
#         tau = float(fn.split("_")[3]) 
#         Temp = float(fn.split("_")[4]) 
#         left_fwhm, right_fwhm = left_right_width(tau, Temp)
#         print(tau, Temp, left_fwhm, right_fwhm)
#         print(fn)
#         data = np.load(fn)
#         end_jl = time.time()
#         _, center = get_center_asym(data.T, plotting=True, ROI_search=False, left_fwhm = left_fwhm, right_fwhm =  right_fwhm)
#         end_py = time.time()
#         print(center)
# 

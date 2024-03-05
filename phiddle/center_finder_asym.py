import glob
import time

from scipy.interpolate import CubicSpline
import numpy as np
import copy as cp
from astropy.modeling import models, fitting


def get_center_asym(data, left_fwhm = 1., right_fwhm = 1., window = 51):
    # Using the correlation method
    im1 = np.array(data[:,:])
    im = cp.deepcopy(im1)
    center = []
    weights = []
    weights_max = []
    corr = []
    mid = im.shape[1]
    scale = right_fwhm/left_fwhm # < 0
    x = np.array(range(im1.shape[1]) )
    correlation = np.zeros(im1.shape[1])
    # smoothed = fft_smoothing(im1, 15)
    smoothed = im1
    smoothed = smoothed - np.mean(smoothed, axis = -1).reshape(im.shape[0], 1) 
    #New
    smoothed = smoothed/np.max(smoothed)
    eps = 0.01
    #for i_center in range(int(im1.shape[1]*0.25), int(im1.shape[1]*0.75)):
    for i_center in range(int(im1.shape[1]*0.01), int(im1.shape[1]*0.99)):
        left = smoothed[:, :i_center]
        right = smoothed[:, i_center:]
        cs = CubicSpline(x[i_center:], right, axis = -1) # interpolate rows
        min_x = np.min(x[i_center:])

        #right_scaled = cs(scale * (x[i_center:] -min_x) + min_x)
        maxval_left = np.max(np.abs(left))
        #left = left/maxval_left

        right_scaled = cs(scale * (x[i_center:] -min_x) + min_x) 
        maxval_right = np.max(np.abs(right_scaled))
        right_scaled = right_scaled / maxval_right * maxval_left

        #asds
        #right_scaled = cs(scale * x[i_center:])
        maxlength = min(left.shape[-1], right.shape[-1])
        #print(scale, np.amin(right_scaled), np.amin(left), maxlength)

        left = left[:, i_center-maxlength:i_center]
        #left = (left - np.mean(left)) #/ (np.std(left) + 1.e-8)

        right = right_scaled[:, maxlength-1::-1]
        #right = (right - np.mean(right))#/(np.std(right) + 1.e-8)
        #print(left.shape, right.shape, np.std(left), np.std(right))


        #correlation[i_center] += np.sum(left * right) /maxlength 
        #New
        correlation[i_center] += np.sum(left * right) / (maxlength**2) #* max(maxval_left, maxval_right)

    wmin = int(round((correlation.shape[0] - window)*0.5))
    filt = np.zeros(correlation.shape)
    filt[wmin:correlation.shape[0]-wmin] = 1.
    imaximum_conv = np.argmax(np.multiply(correlation, filt)) - 1

    return imaximum_conv  #imaximum excluded because it was throwing error

#@jit(nopython = False)
def fft_smoothing(d,param):
#d is a 1D vector, param is the cutoff frequency
    rft = np.fft.rfft(d)
    rft[:, int(param):] = 0.
    d = np.fft.irfft(rft, d.shape[1])
    return d

# def reilluminate_image(data):
#     # Fit the data using astropy.modeling
#     p_init = models.Polynomial2D(degree=2)
#     fit_p = fitting.LevMarLSQFitter()
#     x = []
#     y = []
#     z = []
#     size = 20
#     #left bottom
#     for ix in range(0, size, 1): 
#         for iy in range(0, size, 1):
#             x.append(ix) 
#             y.append(iy) 
#             z.append(data[ix, iy])
#     #right bottom
#     for ix in range(data.shape[0]-1, data.shape[0]-size-1, -1): 
#         for iy in range(0, size, 1):
#             x.append(ix) 
#             y.append(iy) 
#             z.append(data[ix, iy])
#     #right top
#     for ix in range(data.shape[0]-1, data.shape[0]-size-1, -1): 
#         for iy in range(data.shape[1]-1, data.shape[1]-size-1, -1):
#             x.append(ix) 
#             y.append(iy) 
#             z.append(data[ix, iy])
#     #left top
#     for ix in range(0, size, 1): 
#         for iy in range(data.shape[1]-1, data.shape[1]-size-1, -1):
#             x.append(ix) 
#             y.append(iy) 
#             z.append(data[ix, iy])
#     
#     with warnings.catch_warnings():
#         # Ignore model linearity warning from the fitter
#         warnings.filterwarnings('ignore', message='Model is linear in parameters',
#                                 category=AstropyUserWarning)
#         p = fit_p(p_init, x, y, z)
#     # Plot the data with the best-fit model
#     x, y = np.mgrid[:data.shape[0], :data.shape[1]]
#     vmin = 0
#     vmax = np.max(data)
#     plt.figure(figsize=(8, 2.5))
#     plt.subplot(1, 3, 1)
#     plt.imshow(data, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
#     plt.title("Data")
#     plt.subplot(1, 3, 2)
#     plt.imshow(p(x, y), origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
#     plt.title("Model")
#     plt.subplot(1, 3, 3)
#     vmax = np.max(data - p(x, y)) 
#     vmin = np.min(data - p(x, y)) 
#     residual = data - p(x, y)
#     plt.imshow(residual, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
#     plt.title("Residual")
#     plt.show()
#     return residual
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

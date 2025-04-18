import math
from types import SimpleNamespace

import numpy as np

# Helper functions
# Should moslty be self explanable

__version__ = "0.2.0"
__date__ = "02/17/2025"
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#dd8dff', '#48971b']


# 1d distribution functions
def oned_gaussian(x, height, x_0, width_x):
    return height * np.exp( -((x-x_0)/width_x)**2 /2) 


def oned_gaussian_func(height, x_0, width_x):
    return lambda x: height * np.exp( -((x-x_0)/width_x)**2 /2)


def lorentz(x, height, x_0, width_x):
    return height / (1+((x-x_0)/width_x)**2)


def two_lorentz(height, x_0, sigma_1, sigma_2):
    return lambda x: ( lorentz(x, height, x_0, sigma_1)*(x<=x_0).astype(int)
                     + lorentz(x, height, x_0, sigma_2)*(x>x_0).astype(int) )
# End of 1d distribution functions


def get_continue_patches(indices):
    """ return [(start_0, width_0), (start_1, width_1), ...] """
    start_width_ls = []
    val = indices[0]
    start_idx = 0
    end_idx = 0
    while end_idx < len(indices)-1:
        end_idx += 1
        val += 1
        if indices[end_idx] != val:
            start_width_ls.append((indices[start_idx], end_idx-start_idx))
            start_idx = end_idx
            val = indices[start_idx]
    # Ending case
    start_width_ls.append((indices[start_idx], end_idx-start_idx+1))
    return start_width_ls


def replace_nan_w_None(ls):
    a = []
    for i in ls:
        if math.isnan(i):
            a.append(None)
        else:
            a.append(i)

     

def collect_data_and_q(h5, conds):
    q = []
    data = []

    for cond in conds:
        q.append(h5[cond]['0']['integrated_1d'][0])

        dim2 = len(h5[cond])
        dim1 = len(h5[cond]['0']['integrated_1d'][0])

        arr = np.zeros((dim1, dim2))
        for i in range(dim2):
            arr[:, i] = h5[cond][str(i)]['integrated_1d'][1]
        data.append(arr)

    return data, q 


def collect_conditions(conds):

    tau_arr = np.zeros(len(conds))
    tpeak_arr = np.zeros(len(conds))

    for idx, cond in enumerate(conds):
        tau_arr[idx], tpeak_arr[idx] = get_condition(cond)

    return tau_arr, tpeak_arr


def get_condition(condition_str):
    _, tau, _, tpeak = condition_str.split('_')
    return float(tau), float(tpeak)


def collect_positions(h5, conds):  # pass conds to make sure the order is correct

    x_arr = np.zeros(len(conds))
    y_arr = np.zeros(len(conds))
    for idx, cond in enumerate(conds):
        x_arr[idx] = h5[cond].attrs['x_center']
        y_arr[idx] = h5[cond].attrs['y_center']

    return x_arr, y_arr


def minmax_norm(data):
    _min = np.min(data)
    _data = data - _min 
    _max = np.max(_data)
    _data /= _max
    return _data, _min, _max

def minmax_denorm(data, _min, _max):
    return data*_max + _min



def index_phase(phase_names, phase_ls):
    return [idx for idx, name in enumerate(phase_ls)
            if name in phase_names]


def find_first_larger(arr, value):
    for idx, val in enumerate(arr):
        if val > value:
            return idx
    return 0


def find_first_smaller(arr, value):
    for idx, val in enumerate(arr[::-1]):
        if val < value:
            return arr.shape[0] - idx
    return 0


def remove_back_slash(s):
    return s.replace('/' , '')


def dict_to_namespace(d):
    # Recursively convert a dictionary to SimpleNamespace
    if isinstance(d, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in d.items()})
    elif isinstance(d, list):
        return [dict_to_namespace(item) for item in d]
    else:
        return d

import numpy as np


COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


def collect_data_and_q(h5, conds):

    data_dict = {}

    for cond in conds:
        cond_dict = {}
        q = h5[cond]['0']['integrated_1d'][0]
        dim2 = len(h5[cond])
        dim1 = len(h5[cond]['0']['integrated_1d'][0])

        arr = np.zeros((dim1, dim2))
        for i in range(dim2):
            arr[:, i] = h5[cond][str(i)]['integrated_1d'][1]
        cond_dict['q'] = q
        cond_dict['data'] = arr
        cond_dict['cond'] = cond
        cond_dict['x'] = h5[cond].attrs['x_center']
        cond_dict['y'] = h5[cond].attrs['y_center']
        if 'fracs' in list(h5[cond].attrs):
            cond_dict['fracs'] = h5[cond].attrs['fracs']
            cond_dict['cations'] = h5[cond].attrs['cations']
        data_dict[cond] = cond_dict

    return data_dict


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
    _data = data - np.min(data)
    _data /= np.max(_data)
    return _data


def index_phase(phase_names, phase_ls):
    return [idx for idx, name in enumerate(phase_ls)
            if name in phase_names]


def find_first_larger(arr, value):
    for idx, val in enumerate(arr):
        if val > value:
            return idx


def find_first_smaller(arr, value):
    for idx, val in enumerate(arr[::-1]):
        if val < value:
            return arr.shape[0] - idx

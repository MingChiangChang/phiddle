import numpy as np


def collect_data_and_q(h5, conds):

    data_dict = {}

    for cond in conds:
        cond_dict = {}
        q = h5[cond]['0']['integrated_1d'][0]
        dim2 = len(h5[cond])
        dim1 = len(h5[cond]['0']['integrated_1d'][0])

        arr = np.zeros((dim1, dim2))
        for i in range(dim2):
            arr[:,i] = h5[cond][str(i)]['integrated_1d'][1]
        cond_dict['q'] = q
        cond_dict['data'] = arr
        cond_dict['cond'] = cond
        cond_dict['x'] = h5[cond].attrs['x_center']
        cond_dict['y'] = h5[cond].attrs['y_center']
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

def collect_positions(h5, conds): # pass conds to make sure the order is correct

    x_arr = np.zeros(len(conds))
    y_arr = np.zeros(len(conds))
    for idx, cond in enumerate(conds):
        x_arr[idx] = h5[cond].attrs['x_center']
        y_arr[idx] = h5[cond].attrs['y_center']

    return x_arr, y_arr

def minmax_norm(data):
    data -= np.min(data)
    data /= np.max(data)
    return data

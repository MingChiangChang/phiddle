import numpy as np


COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#dd8dff', '#48971b']

def lorentz(x, height, x_0, width_x):
    return height / (1+((x-x_0)/width_x)**2)


def two_lorentz(height, x_0, sigma_1, sigma_2):
    return lambda x: ( lorentz(x, height, x_0, sigma_1)*(x<=x_0).astype(int)
                     + lorentz(x, height, x_0, sigma_2)*(x>x_0).astype(int) )


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


# 
def _get_strain(a, a0):
    return (a-a0)/a0


def _regularize(params, mean_θ, std_θ):
    return (params-np.log(mean_θ)) / (np.sqrt(2)*std_θ)


def extend_priors(mean_θ, std_θ, phases):

    total_params = np.sum([c.free_param_num for c in phases])
    full_mean_θ = np.zeros(total_params)
    full_std_θ = np.zeros(total_params)
    start = 0

    for i, cs in enumerate(phases):
        n = cs.cl.free_param_num
        p = 1 if len(mean_θ)/3. == len(phases) else 1
        full_mean_θ[start:start+n] = mean_θ[3*(p-1)] * np.array(cs.cl.get_free_params())
        full_std_θ[start:start+n] = std_θ[3*(p-1)] * np.array(cs.cl.get_free_params())
        full_mean_θ[start+n:start+n+2] = mean_θ[1+3*(p-1):3*p]
        full_std_θ[start+n:start+n+2] = std_θ[1+3*(p-1):3*p]

        if cs.peakprofile.free_param_num > 0:
            full_mean_θ[start + n + 2] = 0.5
            full_std_θ[start + n + 2] = 10.

        start += cs.free_param_num

    return full_mean_θ, full_std_θ

import os
import json

import numpy as np
import pandas as pd
import h5py

from temp_profile import (
        left_right_width_2024, left_right_width_2023,
        temperature_profile_func_dict, width_func_dict
    )
from center_finder_asym import get_center_asym
from util import collect_data_and_q, collect_conditions, collect_positions, get_condition, two_lorentz
from datastructure import LabelData


# TODO: Storing refined lp
class datamodel():

    def __init__(self):
        # Instantiate empty DataFrame so that manuvering through the empty UI
        # Does not cause errors
        self.df_data = {}
        self.size = 1
        self.df_data['Tpeak'] = [0]
        self.df_data['Dwell'] = [0]
        self.df_data['data'] = [np.zeros((1, 1))]
        self.df_data['q'] = [np.zeros((1, 1))]
        self.df_data['x'] = [0]
        self.df_data['y'] = [0]
        self.df_data['fracs'] = [0]
        self.df_data['cation'] = [""]
        self.df_data["phases"] = [[] for _ in range(self.size)]
        self.df_data["refined_lps"] = [[] for _ in range(self.size)]
        self.df_data["refined_lps_uncer"] = [[] for _ in range(self.size)] 
        self.df_data["act"] = [[] for _ in range(self.size)]
        self.df_data["act_uncer"] = [[] for _ in range(self.size)]
        self.df_data["width"] = [[] for _ in range(self.size)]
        self.df_data["width_uncer"] = [[] for _ in range(self.size)]
        self.df_data["is_refined"] = [False for _ in range(self.size)]

        self.df = pd.DataFrame(self.df_data)
        self.left_right_width = left_right_width_2024
        self.labeldata = LabelData()
        self.temperature_profile_year = "2024"
        self.temperature_profile_func_dict = temperature_profile_func_dict
        self.width_func_dict = width_func_dict


    def read_h5(self, file_path, q_path=None):
        self.df_data = {}
        self.h5_path = file_path  # keep for now, may be redundent
        if not os.path.isfile(file_path):
            print("h5 path not found")
            return
        self.h5 = h5py.File(file_path, 'r')['exp']

        self._ind = 0 # State holder.. Should model has state..?
        self.conds = sorted(list(self.h5))
        self.df_data['data'], self.df_data['q'] = collect_data_and_q(self.h5, self.conds)
        self.size = len(self.conds)

        self.df_data['Dwell'], self.df_data['Tpeak'] = collect_conditions(list(self.h5))
        
        self.df_data['x'], self.df_data['y'] = collect_positions(self.h5, self.conds)
        if 'fracs' in list(self.h5[self.conds[0]].attrs):
            self.cations = self.h5[self.conds[0]].attrs['cations'].tolist()

            fracs = []
            for i, cation in enumerate(self.cations):
                _frac = []
                for j, cond in enumerate(self.conds):
                    _frac.append(self.h5[cond].attrs['fracs'][i][0])
                self.df_data[cation] = _frac
                fracs.append(_frac)
            self.df_data['cations'] = [self.cations for _ in range(self.size)]
            self.df_data['fracs'] = (np.array(fracs).T).tolist()
        else:
            self.cations = []


        if 'xx' in list(self.h5[self.conds[0]].attrs):
            xx = []
            for _, cond in enumerate(self.conds):
                xx.append(self.h5[cond].attrs['xx'])
            self.df_data['xx'] = xx

        self.df_data["phases"] = [[] for i in range(self.size)]
        self.df_data["refined_lps"] = [[] for _ in range(self.size)]
        self.df_data["refined_lps_uncer"] = [[] for _ in range(self.size)] 
        self.df_data["act"] = [[] for _ in range(self.size)]
        self.df_data["act_uncer"] = [[] for _ in range(self.size)]
        self.df_data["width"] = [[] for _ in range(self.size)]
        self.df_data["width_uncer"] = [[] for _ in range(self.size)]
        self.df_data["is_refined"] = [False for _ in range(self.size)]
        self.df_data["center_idx"] = [None for _ in range(self.size)]

        self.current_dwell, self.current_tpeak = get_condition(self.conds[self._ind])
        self.df = pd.DataFrame(self.df_data)
        self.labeldata = LabelData()

    def clear_label_data(self):
        self.labeldata = LabelData()


    def get_xaxis_size(self):
        return self.df['data'][self.ind].shape[1]

    def get_current_labeled_indices(self):
        return self.labeldata.get_labeled_indices(self.ind)

    def set_current_center(self, center):
        self.df.at[self.ind, "center_idx"] = center


    def get_cations(self):
        if hasattr(self, 'cations'):
            return self.cations
        return []


    def update(self, ind):
        self.ind = ind

    def set_temp_profile_params_by_year(self, year):
        if year not in temperature_profile_func_dict:
            print("Year specified does not exist. This shouldn't happen.",
                  "Please contact Ming!")
            return
        
        self.temperature_profile_year = year


    def get_lps_update_dict(self):
        return {'refined_lps': [], 
                'refined_lps_uncer': [], 
                'act': [], 
                'act_uncer': [], 
                'width': [], 
                'width_uncer': [] }


    def update_ind(self, ind, data_dict):
        for key, value in data_dict.items():
            self.df.at[ind, key] = value


    def update_refined_lp(self, ind, refined_lp):
        self.df.at[ind, 'refined_lps'] = refined_lp
        self.df.at[ind, 'is_refined'] = True


    def update_refined_lp_uncer(self, ind, refined_lp_uncer):
        self.df.at[ind, 'refined_lps_uncer'] = refined_lp_uncer


    def update_refined_act(self, ind, act):
        self.df.at[ind, 'act'] = act


    def update_refined_act_uncer(self, ind, act_uncer):
        self.df.at[ind, 'act_uncer'] = act_uncer


    def update_refined_width(self, ind, width):
        self.df.at[ind, 'width'] = width 


    def update_refined_width_uncer(self, ind, width_uncer):
        self.df.at[ind, 'width_uncer'] = width_uncer


    def add_to_phase_diagram(self, phase_names):
        if self.df['phases'][self._ind] != phase_names:
            self.df.at[self._ind, 'phases'] = phase_names
            self.df.at[self._ind, 'is_refined'] = False 


    def remove_from_phase_diagram(self):
        self.df.at[self._ind, 'phases'] = []
        self.df.at[self._ind, 'is_refined'] = False 


    def update_phases(self, phases): 
        self.df['phases'] = phases


    def is_refined(self, phase_name):
        indices = self.get_index_with_phase(phase_name)
        return np.all([self.df['is_refined'][ind] for ind in indices])


    def __getitem__(self, ind):
        return self.df.iloc[ind]


    def get_phases(self):
        return self.df['phases'].to_list()


    def get_all_phases(self):
        all_phases = set()
        for phase_ls in self.df['phases']:
            for phase in phase_ls:
                all_phases.add(phase)
        return all_phases


    def get_dict_for_phase_diagram(self):
        phase_dict = {}
        phases = self.get_all_phases()
        for phase in phases:
            phase_dict[phase] = {}
            sub_df = self.df.iloc[self.df['phases'].map(lambda x: phase in x).to_numpy()]
            phase_dict[phase]['Dwell'] = sub_df['Dwell'].to_list()
            phase_dict[phase]['Tpeak'] = sub_df['Tpeak'].to_list()
            phase_dict[phase]['refined_lps'] = sub_df['refined_lps'].to_list()
            if self.cations:
                for _, cation in enumerate(self.cations):
                    phase_dict[phase][cation] = sub_df[cation].to_list()

        return phase_dict
        

    def get_dict_for_lp_plot(self, phase):
        indices = self.get_index_with_phase(phase)
        lp_dict = {}
        lp_dict["phases"] =self.df.loc[indices, "phases"].to_list()
        lp_dict["Tpeak"] = self.df.loc[indices, 'Tpeak'].to_list()
        lp_dict["Dwell"] = self.df.loc[indices, 'Dwell'].to_list()
        lp_dict["x"] = self.df.loc[indices, 'x'].to_list()
        lp_dict["y"] = self.df.loc[indices, 'y'].to_list()
        cations = self.get_cations()
        for _, cation in enumerate(cations):
            lp_dict[cation] = self.df.loc[indices, cation].to_list()
        lp_dict["refined_lps"] = self.df.loc[indices, 'refined_lps'].to_list()
        lp_dict["refined_lps_uncer"] = self.df.loc[indices, 'refined_lps_uncer'].to_list()
        lp_dict["act"] = self.df.loc[indices, 'act'].to_list()
        lp_dict["act_uncer"] = self.df.loc[indices, 'act_uncer'].to_list()
        lp_dict["width"] = self.df.loc[indices, 'width'].to_list()
        lp_dict["width_uncer"] = self.df.loc[indices, 'width_uncer'].to_list()
        return lp_dict


    def save_lp(self, filename, phase):
        lp_dict = self.get_dict_for_lp_plot(phase)

        with open(filename, 'w') as f:
            json.dump(lp_dict, f)


    def get_index_with_phase(self, phase):
        if hasattr(self, 'df'):
            return [i for i, p in enumerate(self.df['phases']) if phase in p]
        return []

    def get_labeled_phases(self, idx):
        return [self.df['phases'][i] for i in idx]

    def get_refined_lp_of_phase(self, phase_name):
        indices = self.get_index_with_phase(phase_name)
        phases_ls = self.get_labeled_phases(indices)
        refined_lp = []

        for i, phases in zip(indices, phases_ls):
            for j, phase in enumerate(phases):
                if phase_name == phase:
                    refined_lp.append(self.refined_lps[i][j])

        return refined_lp


    def update_temp_profile_for_stored_labels(self):
        # Go throgh labeled data, find the ones that belongs to 
        for label in self.labeldata:
            if label.sample_num == self.current_ind:
                label.tpeak = self.temp_func(np.array(self.transform_data_idx_to_x(label.x_idx)))
        

    def get_temp_profile_at(self, ind):
        """
        Returns xaxis, temperature function, and  center position (in data index)
        Temperature profile can be obtained by temp_func(xaxis)
        """
        temp_func = self.temperature_profile_func_dict[self.temperature_profile_year](self.current_dwell,
                                                                                      self.current_tpeak)
        width_func = self.width_func_dict[self.temperature_profile_year]

        center = self.df['center_idx'][ind]
        if center is None:
            # FIXME: Hiddne side effect
            center = get_center_asym(
                    self.df['data'][ind],
                    *width_func(self.df['Tpeak'][ind], self.df['Dwell'][ind]))
            self.df.at[ind, 'center_idx'] = center 
       
        
        xaxis = self.get_xaxis(ind)
        self.xaxis = xaxis           
        self.temp_func = temp_func           
        return xaxis, temp_func, center

    def get_current_temp_profile(self):
        return self.get_temp_profile_at(self.ind)

    def get_current_phases_bw_x_range(self, x_min_ind, x_max_ind):
        phases = set()
        x_range = range(x_min_ind, x_max_ind+1)
        if self.df['center_idx'][self.ind] in x_range:
            for phase in self.df['phases'][self.ind]:
                phases.add(phase)
        for x_idx in x_range:
            phase_names = self.labeldata.get_phase(self.ind, x_idx)
            if phase_names is not None:
                for phase_name in phase_names:
                    phases.add(phase_name)
        return list(phases)

    def get_xaxis(self, ind):
        if 'xx' in self.df:
            xmin = self.df['xx'][ind][0]
            xmax = self.df['xx'][ind][-1]
            xaxis = np.arange(xmax-xmin)# (xmax-xmin)/2
            xaxis -= xaxis[int(len(xaxis)*self.df['center_idx'][self.ind]/len(self.df['xx'][self.ind]))]
        else:
            xmin = 0
            xmax = self.df['data'][self.ind].shape[1]
            xaxis = np.arange(xmax-xmin) * 10 # 10 um per column
            xaxis -= xaxis[self.df['center_idx'][self.ind]]
        self.xaxis = xaxis
        return xaxis

    def get_current_xaxis(self):
        if self.current_xx is not None:
            xmin = self.current_xx[0]
            xmax = self.current_xx[-1]
            xaxis = np.arange(xmax-xmin)# (xmax-xmin)/2
            xaxis -= xaxis[int(len(xaxis)*self.df['center_idx'][self.ind]/len(self.current_xx))]
        else:
            xmin = 0
            xmax = self.df['data'][self.ind].shape[1]
            xaxis = np.arange(xmax-xmin) * 10 # 10 um per column
            xaxis -= xaxis[self.df['center_idx'][self.ind]]
        self.xaxis = xaxis
        return xaxis #self.current_xx 


    def get_stripe_xlabel(self):
        if self.current_xx is None: return "Index"
        return "Location (um)"

    def get_center_idx(self):
        center =  self.df['center_idx'].tolist()
        rc = []
        for c in center:
            if c is not None:
                rc.append(int(c))
            else:
                rc.append(None)
        return rc

    def load_center_idx(self, center_idx):
        self.df['center_idx'] = [None for _ in range(self.size)]
        for i, center in enumerate(center_idx):
            if center is not None:
                self.df.at[i, "center_idx"] = center

    @property
    def current_ind(self):
        return self._ind

    @property
    def labeled(self):
        return [bool(phase) for phase in self.df['phases']]


    @property
    def current(self):
        tmp = np.zeros(len(self.df))
        tmp[self.ind] = 1
        tmp = tmp.astype(bool)
        return tmp


    @property
    def current_data(self):
        return self.__getitem__(self.ind)


    @property
    def current_xx(self):
        if 'xx' in self.df:
            return self.df['xx'][self.ind]
        return None

    @property
    def current_center(self):
        return self.df['center_idx'][self.ind]


    @property
    def ind(self):
        return self._ind


    @ind.setter
    def ind(self, new_ind):
        if new_ind == len(self.conds):
            return
        self._ind = new_ind
        self.current_dwell, self.current_tpeak = get_condition(
            self.conds[new_ind])
        if self.cations:
            self.current_composition = {c: f for c, f in zip(self.df.iloc[new_ind]['cations'],
                                          self.df.iloc[new_ind]['fracs'])}
        else:
            self.current_composition = dict()

    @property
    def current_x(self):
        return self.x[self._ind]


    @property
    def current_y(self):
        return self.y[self._ind]


    @property
    def current_filename(self):
        filename = self.conds[self._ind]
        if 'cations' in self.h5[filename].attrs:
            for idx, cation in enumerate(
                    self.h5[self.conds[self._ind]].attrs['cations']):
                frac = self.h5[self.conds[self._ind]].attrs['fracs'][idx][0]
                filename += f'_{cation}_{frac:.3f}'
        return filename


    def transform_data_idx_to_x(self, x):
        return int(self.xaxis[0] + x / self.df['data'][self.ind].shape[1] * len(self.xaxis))

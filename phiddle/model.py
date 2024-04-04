import re
import json

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import h5py

from util import collect_data_and_q, collect_conditions, collect_positions, get_condition


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
        self.df_data["phases"] = [[] for i in range(self.size)]
        self.df_data["refined_lps"] = [[] for _ in range(self.size)]
        self.df_data["refined_lps_uncer"] = [[] for _ in range(self.size)] 
        self.df_data["act"] = [[] for _ in range(self.size)]
        self.df_data["act_uncer"] = [[] for _ in range(self.size)]
        self.df_data["width"] = [[] for _ in range(self.size)]
        self.df_data["width_uncer"] = [[] for _ in range(self.size)]
        self.df_data["is_refined"] = [False for _ in range(self.size)]

        self.df = pd.DataFrame(self.df_data)

        
    def read_h5(self, file_path, q_path=None):
        self.df_data = {}
        self.h5_path = file_path  # keep for now, may be redundent
        try:
            self.h5 = h5py.File(file_path, 'r')['exp']
        except FileNotFoundError:
            print("h5 path not found")
            return

        self._ind = 0 # State holder.. Should model has state..?
        self.conds = sorted(list(self.h5))
        self.df_data['data'], self.df_data['q'] = collect_data_and_q(self.h5, self.conds)
        self.size = len(self.conds)

        self.df_data['Dwell'], self.df_data['Tpeak'] = collect_conditions(list(self.h5))
        
        self.df_data['x'], self.df_data['y'] = collect_positions(self.h5, self.conds)
        if 'fracs' in list(self.h5[self.conds[0]].attrs):
            self.cations = self.h5[self.conds[0]].attrs['cations']

            for i, cation in enumerate(self.cations):
                _frac = []
                for j, cond in enumerate(self.conds):
                    _frac.append(self.h5[cond].attrs['fracs'][i][0])
                self.df_data[cation] = _frac
            self.df_data['cations'] = [self.cations for _ in range(self.size)]


        if 'xx' in list(self.h5[self.conds[0]].attrs):
            xx = []
            for idx, cond in enumerate(self.conds):
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

        self.current_dwell, self.current_tpeak = get_condition(self.conds[self._ind])
        self.df = pd.DataFrame(self.df_data)


    def get_cations(self):
        if hasattr(self, 'cations'):
            return self.cations
        return []


    def update(self, ind):
        self.ind = ind


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


    def update_phases(self, phases): 
        self.df['phases'] = phases


    def is_refined(self, phase_name):
        indicies = self.get_index_with_phase(phase_name)
        return np.all([self.df['is_refined'][ind] for ind in indicies])


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
            if hasattr(self, 'cations'):
                for j, cation in enumerate(self.cations):
                    phase_dict[phase][cation] = sub_df[cation].to_list()

        return phase_dict
        

    def get_dict_for_lp_plot(self, phase):
        indicies = self.get_index_with_phase(phase)
        lp_dict = {}
        lp_dict["Tpeak"] = self.df.loc[indicies, 'Tpeak'].to_list()
        lp_dict["Dwell"] = self.df.loc[indicies, 'Dwell'].to_list()
        lp_dict["x"] = self.df.loc[indicies, 'x'].to_list()
        lp_dict["y"] = self.df.loc[indicies, 'y'].to_list()
        cations = self.get_cations()
        for j, cation in enumerate(cations):
            lp_dict[cation] = self.df.loc[indicies, cation].to_list()
        lp_dict["refined_lps"] = self.df.loc[indicies, 'refined_lps'].to_list()
        lp_dict["refined_lps_uncer"] = self.df.loc[indicies, 'refined_lps_uncer'].to_list()
        lp_dict["act"] = self.df.loc[indicies, 'act'].to_list()
        lp_dict["act_uncer"] = self.df.loc[indicies, 'act_uncer'].to_list()
        lp_dict["width"] = self.df.loc[indicies, 'width'].to_list()
        lp_dict["width_uncer"] = self.df.loc[indicies, 'width_uncer'].to_list()
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
    def ind(self):
        return self._ind


    @ind.setter
    def ind(self, new_ind):
        if new_ind == len(self.conds):
            return
        self._ind = new_ind
        self.current_dwell, self.current_tpeak = get_condition(
            self.conds[new_ind])


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

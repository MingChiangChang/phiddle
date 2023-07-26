import re

import numpy as np
from scipy.interpolate import CubicSpline
import h5py

from PyQt6.QtCore import Slot
from util import collect_data_and_q, collect_conditions, collect_positions, get_condition


class datamodel():

    def __init__(self):
        pass
        
    def read_h5(self, file_path, q_path=None):
        self.h5_path = file_path  # keep for now, may be redundent
        try:
            self.h5 = h5py.File(file_path, 'r')['exp']
        except FileNotFoundError:
            print("h5 path not found")
            return

        self._ind = 0
        self.conds = sorted(list(self.h5))
        self.data = collect_data_and_q(self.h5, self.conds)
        self.size = len(self.conds)

        self.dwells, self.tpeaks = collect_conditions(list(self.h5))
        self.current_dwell, self.current_tpeak = get_condition(
            self.conds[self._ind])

        self.x, self.y = collect_positions(self.h5, self.conds)
        if 'cations' in list(self.h5[self.conds[0]].attrs):
            self.cations = self.h5[self.conds[0]].attrs['cations']

            self.fractions = np.zeros((len(self.conds), len(self.cations)))
            for i, _ in enumerate(self.cations):
                for j, cond in enumerate(self.conds):
                    self.fractions[j,
                                   i] = self.h5[cond].attrs['fracs'][i][0]

        if 'xx' in list(self.h5[self.conds[0]].attrs):
            self.xx = []
            for idx, cond in enumerate(self.conds):
                self.xx.append(self.h5[cond].attrs['xx'])
        self.phases = ["" for i in range(len(self.conds))]


    def read_udi(self, udi):
        with open(udi, 'r') as f:
            d = f.read()

        d = d.split()

        for line in d:
            if line.startswith('N='):
                self.size = int(line[line.index('=')+1:])
            if line.startswith("Composition="):
                self.cations = line[line.index('=')+1:].split(',')
            if line.startswith('X='):
                self.x = np.array(list(map(float, line[line.index('=')+1:].split(','))))
            if line.startswith('Y='):
                self.y = np.array(list(map(float, line[line.index('=')+1:].split(','))))
            if line.startswith('Q='):
                self.q = np.array(list(map(float, line[line.index('=')+1:].split(','))))

        # self.data = np.zeros((self.size, len(self.q)))
        self.fractions = np.zeros((self.size, len(self.cations)))
        for idx, cation in enumerate(self.cations):
            for line in d:
                if line.startswith(cation + '='):
                    fractions = line[line.index('=')+1:].split(',')
                    self.fractions[:,idx] = np.array(list(map(float, fractions)))

        self.data = {}
        for line in d:
            if re.match(r'^I\d+=', line):
                ind = int(line[1:line.index('=')]) - 1
                self.data[ind] ={} 
                d = np.array(list(map(float, line[line.index('=')+1:].split(',')))) 
                q = self.q
                if len(d) > 2048:
                    q = np.linspace(self.q[0], self.q[-1], 2048)
                    cs = CubicSpline(self.q, d)
                    d = cs(q)
                self.data[ind]['data'] = np.tile(d, (151, 1)).T
                self.data[ind]['q'] = q
                self.data[ind]['cond'] = "tau_0_T_0"
                self.data[ind]['x'] = self.x[ind]
                self.data[ind]['y'] = self.y[ind]
                self.data[ind]['fracs'] = [[self.fractions[ind, j]] for j, _ in enumerate(self.cations)] 
                self.data[ind]['cations'] = self.cations

        self.phases = ["" for _ in range(self.size)]
        self.dwells = np.array([0 for _ in range(self.size)])
        self.tpeaks = np.array([0 for _ in range(self.size)])
        self.current_dwell = 0
        self.current_tpeak = 0
        self.conds = ["tau_0_T_0" for _ in range(self.size)]


                    

    def update(self, ind):
        self.ind = ind

    def add_to_phase_diagram(self, phase_names):
        self.phases[self._ind] = phase_names

    def __getitem__(self, ind):
        if not ind in self.data: 
            data = self.data[self.conds[ind]]
        else:
            data = self.data[ind]
        return data

    def get_dict_for_phase_diagram(self):

        phase_dict = {}
        for idx, phases in enumerate(self.phases):
            if phases:
                for phase in phases:
                    if phase not in phase_dict:
                        phase_dict[phase] = {}
                        phase_dict[phase]['dwell'] = []
                        phase_dict[phase]['tpeak'] = []
                    phase_dict[phase]['dwell'].append(self.dwells[idx])
                    phase_dict[phase]['tpeak'].append(self.tpeaks[idx])
        return phase_dict

    @property
    def labeled(self):
        return [bool(phase) for phase in self.phases]

    @property
    def labeled_dwells(self):
        return self.dwells[self.labeled]

    @property
    def labeled_tpeaks(self):
        return self.tpeaks[self.labeled]

    @property
    def labeled_x(self):
        return self.x[self.labeled]

    @property
    def labeled_y(self):
        return self.y[self.labeled]

    @property
    def current_data(self):
        return self.__getitem__(self.ind)

    @property
    def current_xx(self):
        if hasattr(self, 'xx'):
            return self.xx[self.ind]
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

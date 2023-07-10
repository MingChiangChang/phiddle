import numpy as np
import h5py

from PyQt6.QtCore import Slot
from util import collect_data_and_q, collect_conditions, collect_positions, get_condition


class datamodel():

    def __init__(self):  # , file_path, q_path=None):
        pass
        # if file_path.endswith("h5"):
        #    self.h5_path = file_path # keep for now, may be redundent
        #    self.h5 = h5py.File(file_path, 'r')['exp']
        #
        #    self._ind = 0
        #    self.conds = sorted(list(self.h5))
        #    self.data = collect_data_and_q(self.h5, self.conds)
        #    self.size = len(self.conds)

        #    self.dwells, self.tpeaks = collect_conditions(list(self.h5))
        #    self.current_dwell, self.current_tpeak = get_condition(self.conds[self._ind])

        #    self.x, self.y = collect_positions(self.h5, self.conds)

        # if file_path.endswith("npy"):
        #    self.data = np.load(file_path)
        #    self.q = np.load(q_path)
        #
        # self.phases = ["" for i in range(len(self.conds))]

    def read_h5(self, file_path, q_path=None):
        if file_path.endswith("h5"):
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

        if file_path.endswith("npy"):
            self.data = np.load(file_path)
            self.q = np.load(q_path)

        self.phases = ["" for i in range(len(self.conds))]

    # @Slot(int)

    def update(self, ind):
        self.ind = ind

    def add_to_phase_diagram(self, phase_names):
        self.phases[self._ind] = phase_names

    def __getitem__(self, ind):
        data = self.data[self.conds[ind]]
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

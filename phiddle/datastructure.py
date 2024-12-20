from dataclasses import dataclass
from typing import List, Dict


@dataclass
class PatternLabel():
    tpeak: float 
    dwell: float 
    composition: Dict[str, float] 
    phase: List[str] 
    sample_num: int   # To back track the XRD pattern in the datamodel
    x_idx: int        # This too 

# Just going to maitain a list of PatternLabel
# Try doing seperate data class for h5 files and phase label information

class LabelData():
    """ This data structure keep track of finer phase labeling details """

    def __init__(self):
        self.labels = []
        self.sample_nums = []
        self.x_indices = {} 
        
    def __getitem__(self, idx):
        return self.labels[idx]

    def update(self, tpeak, dwell, composition, phase, sample_num, x_idx):
        label = PatternLabel(tpeak, dwell, composition, phase, sample_num, x_idx)
        self.labels.append(label)
        if not sample_num in self.sample_nums:
            self.sample_nums.append(sample_num)
            self.x_indices[str(sample_num)] = [x_idx] # official python json automatically makes key to be strings
            return 

        if not x_idx in self.x_indices[str(sample_num)]:
            self.x_indices[str(sample_num)].append(x_idx) 
            return

        # If its entry is already in there, do O(n) search 
        # (Could later do O(nlogn) binary search if it is too slow
        # We will have to maintain the list to be sorted

        for label in self.labels:
            if label.sample_num == sample_num and label.x_idx == x_idx:
                label.phase = phase

    def remove(self, sample_num, x_idx):
        for idx, label in enumerate(self.labels):
            if label.sample_num == sample_num and label.x_idx == x_idx:
                del self.labels[idx]

        if str(sample_num) in self.x_indices:
            if x_idx in self.x_indices[str(sample_num)]:
                self.x_indices[str(sample_num)].remove(x_idx)
            if not self.x_indices[str(sample_num)]:
                del self.x_indices[str(sample_num)]
                self.sample_nums.remove(sample_num)


    def load_stored_label_data(self, data):
        if 'labels' not in data:
            return
        self.labels = [PatternLabel(**d) for d in data['labels']]
        self.sample_nums = data['sample_nums']
        self.x_indices = data['x_indices']

    def serialize_data(self):
        d = {}
        d['labels'] = [label.__dict__ for label in self.labels]
        d['sample_nums'] = self.sample_nums 
        d['x_indices'] = self.x_indices
        return d 

    def get_dict_for_phase_diagram(self):

        phase_dict = dict()
        phases = self.get_all_phases()
        cation = self.get_cation()

        for phase in phases:
            phase_dict[phase] = {}
            phase_dict[phase]['Dwell'] = []
            phase_dict[phase]['Tpeak'] = []
            for c in cation:
                phase_dict[phase][c] = []

        for label in self.labels:
            for p in label.phase: 
                phase_dict[p]['Dwell'].append(label.dwell) 
                phase_dict[p]['Tpeak'].append(label.tpeak) 
                for c in cation:
                    phase_dict[p][c].append(label.composition[c]) 
        return phase_dict 
        
    def get_cation(self):
        if self.labels:
            return list(self.labels[0].composition)
        else:
            return []

    def get_labeled_indices(self, sample_num):
        # if str(sample_num) in self.x_indices
        # return self.x_indices[str(sample_num)]
        return self.x_indices.get(str(sample_num), [])


    def get_all_phases(self):
        all_phases = set()
        for label in self.labels:
            for phase in label.phase:
                all_phases.add(phase)
        return all_phases

    def get_phase(self, n_sample, x_idx):
        """ O(nm) if there is n labels and m queries"""
        if not n_sample in self.sample_nums:
            return

        for label in self.labels:
            if label.sample_num == n_sample:
                if label.x_idx == x_idx:
                    return label.phase 
        return 


    def get_phase_in_x_range(self, n_sample, x_range):
        """ O(n) regardless of x range size """
        rp = set()
        if not n_sample in self.sample_nums:
            return rp

        for label in self.labels:
            if label.sample_num == n_sample:
                if label.x_idx in x_range:
                    for p in label.phase:
                        rp.add(p)
        return rp

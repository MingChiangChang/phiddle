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
        if not sample_num in self.sample_nums:
            self.labels.append(PatternLabel(tpeak, dwell, composition, phase, sample_num, x_idx))
            self.sample_nums.append(sample_num)
            self.x_indices[str(sample_num)] = [x_idx] # official python json automatically makes key to be strings
            return 

        if not x_idx in self.x_indices[str(sample_num)]:
            self.labels.append(PatternLabel(tpeak, dwell, composition, phase, sample_num, x_idx))
            self.x_indices[str(sample_num)].append(x_idx) 
            return

        # If its entry is already in there, do O(n) search 
        # (Could later do O(nlogn) binary search if it is too slow
        # We will have to maintain the list to be sorted

        for label in self.labels:
            if label.sample_num == sample_num: 
                if label.x_idx == x_idx:
                    label.phase = phase

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



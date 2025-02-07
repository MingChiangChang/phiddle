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
    """
    A data structure for managing and storing phase labeling information related
    to sample data.

    Attributes:
        labels (list): A list to store instances of PatternLabel that represent
                       phase labels.
        sample_nums (list): A list of sample numbers.
        x_indices (dict): A dictionary mapping sample numbers (as strings) to 
                          their corresponding x indices.
    """

    def __init__(self):
        """
        Initializes the LabelData instance with empty lists and dictionaries.
        """
        self.labels = []
        self.sample_nums = []
        self.x_indices = {} 
        
    def __getitem__(self, idx):
        """ allow indexing to get labels at required index """
        return self.labels[idx]

    def update(self, tpeak, dwell, composition, phase, sample_num, x_idx):
        """
        Update the label data with new information or modify an existing label.

        If a label with the same sample number and x index already exists,
        its phase is updated.
        If it does not exist, a new label is created and added.

        parameters:
            tpeak: float
                The peak temperature value associated with the label.
            dwell: float
                The dwell time associated with the label.
            composition: dict
                A dictionary representing the composition of the sample,
                where keys are elements and values are their concentrations.
            phase: list
                A list of phases for the sample.
            sample_num: int
                The sample number for this label.
            x_idx: int
                The x index for this label.
        """
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
        # (Could later change to O(logn) binary search if it is too slow
        # We will have to maintain the list to be sorted

        for label in self.labels:
            if label.sample_num == sample_num and label.x_idx == x_idx:
                label.phase = phase

    def remove(self, sample_num, x_idx):
        """
        Remove the label associated with a specific sample number and x index.

        parameters:
            sample_num: int
                The sample number whose label to remove.
            x_idx: int
                The x index whose label to remove.
        """
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
        """
        Load previously stored label data into the current LabelData instance.

        parameters:
            data: dict
                A dictionary containing previously serialized label data,
                which includes 'labels', 'sample_nums', and 'x_indices' keys.
        """
        if 'labels' not in data:
            return
        self.labels = [PatternLabel(**d) for d in data['labels']]
        self.sample_nums = data['sample_nums']
        self.x_indices = data['x_indices']

    def serialize_data(self):
        """
        Serialize the label data into a dictionary format for storage or transmission.

        return:
            dict: A dictionary containing the serialized label data,
                  including 'labels', 'sample_nums', and 'x_indices'.
        """
        d = {}
        d['labels'] = [label.__dict__ for label in self.labels]
        d['sample_nums'] = self.sample_nums 
        d['x_indices'] = self.x_indices
        return d 

    def get_dict_for_phase_diagram(self):
        """
        Prepare and return a dictionary formatted for use in plotting the phase diagram.

        return:
            dict: A dictionary containing phase data,
                  including Dwell, Tpeak, and composition for each phase.
        """

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
        """
        Retrieve the list of cation elements present in the composition of the labels.

        return:
            list: A list of cation elements in the composition of the first label.
                  Returns an empty list if no labels are present.
        """
        if self.labels:
            return list(self.labels[0].composition)
        else:
            return []

    def get_labeled_indices(self, sample_num):
        """
        Get the labeled x indices associated with a specific sample number.

        parameters:
            sample_num: int
                The sample number to search for.

        return:
            list: A list of labeled x indices associated with the specified sample number.
                  Returns an empty list if no indices are found.
        """
        return self.x_indices.get(str(sample_num), [])


    def get_all_phases(self):
        """
        Get all unique phases across all labels.

        return:
            set: A set of all unique phases across all labels.
        """
        all_phases = set()
        for label in self.labels:
            for phase in label.phase:
                all_phases.add(phase)
        return all_phases

    def get_phase(self, n_sample, x_idx):
        """
        Retrieve the phase for a specific sample number and x index.
        O(nm) if there are n labels and m queries

        parameters:
            n_sample: int
                The sample number to search for.
            x_idx: int
                The x index to search for.

        return:
            list or None: The list of phases for the given sample number and x index,
                          or None if no label is found.
        """
        if not n_sample in self.sample_nums:
            return

        for label in self.labels:
            if label.sample_num == n_sample:
                if label.x_idx == x_idx:
                    return label.phase 
        return 


    def get_phase_in_x_range(self, n_sample, x_range):
        """
        Retrieve all unique phases for a given sample number within a specified x index range.
        O(n) regardless of x range size

        parameters:
            n_sample: int
                The sample(stripe) number to search for.
            x_range: range or list
                A range or list of x indices to check.

        return:
            set: A set of phases associated with the specified sample number and x range.
        """
        rp = set()
        if not n_sample in self.sample_nums:
            return rp

        for label in self.labels:
            if label.sample_num == n_sample:
                if label.x_idx in x_range:
                    for p in label.phase:
                        rp.add(p)
        return rp

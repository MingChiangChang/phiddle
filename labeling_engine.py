from copy import deepcopy
import numpy as np

from pyPhaseLabel import PhaseModel, CrystalPhase, EQ, BackgroundModel, FixedPseudoVoigt
from pyPhaseLabel import create_phases, evaluate_obj, optimize_phase, Lorentz, PseudoVoigt
from julia.Main import Wildcard, Lazytree, search, get_probabilities, get_fraction 
from julia.BackgroundSubtraction import mcbl



class labeler():
    

    def __init__(self):
        pass

        #with open(cifs, 'r') as f:
        #    t = f.read()
        #self.phases = create_phases(t, .2, FixedPseudoVoigt(0.1))

        self.std_noise = 0.05
        self.mean_θ = [1., .5, .1]
        self.std_θ = [0.001, 0.5, 0.05]

        self.max_phase = 1
        self.expand_k = 2
        self.use_background = False
        self.background_option = "MCBL"
        self.optimize_mode = "Simple"
        self.background_length = 8. 
        self.max_iter = 512

        self.results = []
        self.bg = []
        self.has_labeled = False
        #self.phase_names = [phase.name for phase in self.phases]

    def read_cifs(self, cifs):
        with open(cifs, 'r') as f:
            t = f.read()
        self.phases = create_phases(t, .2, FixedPseudoVoigt(0.1))

        self.phase_names = [phase.name for phase in self.phases]

    def fit(self, q, d):
        data = deepcopy(d)
        tree = Lazytree(self.phases, q)

        if self.background_option == "MCBL":
            self.bg = mcbl(data, q, self.background_length)
            data -= self.bg 
            data[data<0] = 1E-5
        if self.background_option == "Default":
            self.use_background = True
        if self.background_option == "None":
            self.use_background = False
        #data -= np.min(data)
        norm = np.max(data)
        data /= norm 
        result = search(tree, q, data,
                        self.max_phase, self.expand_k, .1, False,
                        self.use_background, self.background_length,
                        self.std_noise, self.mean_θ, self.std_θ,
                        optimize_mode="Simple",
                        em_loop_num=5,
                        maxiter=self.max_iter,
                        regularization=True)

        results = [r for subresults in result[1:] for r in subresults]
        t = get_probabilities(results, q, data, self.std_noise, self.mean_θ, self.std_θ)
        t[t==0] = -1          # Trick to remove identical elements
        m = np.min(np.abs(t)) # Replace with tiny random number 
        t[t==-1] = m*np.random.rand(np.sum(t==-1))
        ind = np.argmax(t)
        self.t, self.results = zip(*sorted(zip(t, results), reverse=True))
        print(results[0].phase_model.CPs)

        if self.background_option in ["None", "Default"]:
            self.bg = np.zeros(q.shape)

        self.has_labeled = True
        self.label_ind = 0
        #return results[:5], bg # bg is only nonzero when constant background is used (MCBL0

    def fit_phases(self, q, d, phase_names):
        data = deepcopy(d)
        phases = self.get_phase_w_phase_names(phase_names)

        if self.background_option == "MCBL":
            self.bg = mcbl(data, q, self.background_length)
            data -= self.bg
            data[data<0] = 1E-5
            bg = None
        elif self.background_option == "Default":
            bg = BackgroundModel(x, EQ(), 8., 10.)
        else:
            bg = None
        #data -= np.min(data)
        norm = np.max(data)
        data /= norm

        pm = PhaseModel(phases, None, bg)
        result = optimize_phase(pm, q, data,
                                self.std_noise, self.mean_θ, self.std_θ, 
                                optimize_mode=self.optimize_mode, 
                                maxiter=self.max_iter)
        self.results = [result]
        self.has_labeled = True
        self.label_ind = 0
        if self.background_option in ["None", "Default"]:
            self.bg = np.zeros(q.shape)

    def next_label_result(self):
        self.label_ind += 1
        if self.label_ind == len(self.results):
            self.label_ind = 0
        return self.label_ind+1, self.results[self.label_ind], self.bg


    def previous_label_result(self):
        self.label_ind -= 1
        if self.label_ind < 0:
            self.label_ind = len(self.results)-1
        return self.label_ind+1, self.results[self.label_ind], self.bg


    def get_phase_names(self, isChecked_ls):
        phase_names = [self.phases[idx].name for idx, check in enumerate(isChecked_ls[:-2])
                                        if check]
        if isChecked_ls[-1]:
            phase_names.append("Amorphous")
        return phase_names

    @property
    def params(self):
        return (self.std_noise, self.mean_θ, self.std_θ, self.max_phase,
               self.expand_k, self.background_length, self.max_iter,
               self.optimize_mode, self.background_option)

    def set_hyperparams(self, std_noise, mean, std, max_phase,
                        expand_k, background_length, max_iter,
                        optimize_mode, background_option):

        self.std_noise = std_noise
        self.mean_θ = mean
        self.std_θ = std
        self.max_phase = max_phase
        self.expand_k = expand_k
        self.background_length = background_length
        self.max_iter = max_iter
        self.optimize_mode = optimize_mode
        self.background_option = background_option

    def get_phase_w_phase_names(self, phase_names):
        return [self.phases[self.phase_names.index(phase_name)] for phase_name in phase_names] 

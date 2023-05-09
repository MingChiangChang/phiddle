import numpy as np

from pyPhaseLabel import PhaseModel, CrystalPhase, EQ, BackgroundModel, FixedPseudoVoigt
from pyPhaseLabel import create_phases, evaluate_obj, optimize_phase, Lorentz, PseudoVoigt
from julia.Main import Wildcard, Lazytree, search, get_probabilities, get_fraction 
from julia.BackgroundSubtraction import mcbl



class labeler():
    

    def __init__(self, cifs):

        with open(cifs, 'r') as f:
            t = f.read()
        self.phases = create_phases(t, .2, FixedPseudoVoigt(0.1))

        self.std_noise = 0.05
        self.mean_θ = [1., .5, .1]
        self.std_θ = [0.001, 0.05, 0.05]


    def fit(self, q, data):
        tree = Lazytree(self.phases, q)
        data -= mcbl(data, q, 8.) 
        data[data<0] = 1E-5
        data -= np.min(data)
        data /= np.max(data)
        result = search(tree, q, data, 2, 2, .1, False, False, 8.,
                self.std_noise, self.mean_θ, self.std_θ, optimize_mode="EM",
                em_loop_num=5,
                maxiter=512, regularization=True)

        results = [r for subresults in result[1:] for r in subresults]
        t = get_probabilities(results, q, data, self.std_noise, self.mean_θ, self.std_θ)
        ind = np.argmax(t)
        print(results[ind].phase_model.CPs)
        return results[ind]

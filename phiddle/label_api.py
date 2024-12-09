import operator
from copy import deepcopy

import requests
import numpy as np
from scipy.interpolate import CubicSpline

from pyPhaseLabel import PhaseModel, CrystalPhase, EQ, BackgroundModel, FixedPseudoVoigt
from pyPhaseLabel import create_phases, evaluate_obj, optimize_phase, Lorentz, PseudoVoigt
from julia.Main import Wildcard, Lazytree, search_b, get_probabilities, optimize_b
from julia.BackgroundSubtraction import mcbl
from julia import CrystalShift as CS

from util import minmax_norm

class labeler():

    def __init__(self, address="127.0.0.1", port="8080"):
        self.std_noise = .05
        self.mean_θ = [1., .5, .1]
        self.std_θ = [.05, .5, .1]
        self.address = address
        self.port = port

        self.max_phase = 1
        self.expand_k = 2
        self.use_background = False
        self.background_option = "MCBL"
        self.optimize_mode = "Simple"
        self.background_length = 8.
        self.max_iter = 512
        self.csv_file = ""

        self.results = []
        self.bg = []
        self.has_labeled = False

    def construct_setting_dict(self):
        _dict = {}
        _dict['std_noise'] = self.std_noise
        _dict['mean_θ'] = self.mean_θ
        _dict['std_θ'] = self.std_θ
        _dict["maxiter"] = self.max_iter
        _dict["regularization"] = True 
        _dict["method"] = "LM" # Method is fixed for now
        _dict["objective"] = "LS" # Objective is fixed for now
        _dict["optimization_mode"] = self.optimize_mode 
        _dict["em_loop_num"] = 5
        _dict["λ"] = 1.  # Means nothing for now
        _dict["verbose"] = False
        _dict["tol"] = 1E-6
        _dict["depth"] = self.max_phase 
        _dict["k"] = self.expand_k 
        _dict["amorphous"] = False
        _dict["background"] = self.use_background
        _dict["background_length"] = self.background_length 
        _dict["background_option"] = self.background_option
        _dict["csv_file"] = self.csv_file
        return _dict


    def read_csv(self, csv):
        self.csv_file = csv
        _dict = {"csv_file": csv}
        # TODO: Write a specific module to interact with the Julia server
        # TODO: For checking return code
        response = requests.put(f"{self.server_ip}/set_csv", json=_dict)
        response = requests.get(f"{self.server_ip}/phase_names")
        self.phase_names = response.json()#[phase.name for phase in self.phases]

    def fit(self, q, d, selected_phase_names = None):
        _dict = self.construct_setting_dict()
        self.q = q
        self.data = d 
        _dict["q"] = q.tolist()
        _dict["data"] = d.tolist()

        if self.optimize_mode == "With Uncertainty":
            print("Uncertainty return is not availble for tree seach labeling .. yet.")
            print("Defaulting to Simple optimization")
            _dict["optimize_mode"] = "Simple"
        else:
            _dict["optimize_mode"] = self.optimize_mode

        if selected_phase_names is None:
            response = requests.post(f"{self.server_ip}/label", json=_dict)
        else:
            _dict["names"] = selected_phase_names
            response = requests.post(f"{self.server_ip}/label_with_phase",
                                     json=_dict)
        results, self.probs, self.bg = response.json()
        self.print_fitted_result(results[0]["phase_model"])
        self.results = [result["phase_model"] for result in results]
        self.has_labeled = True
        self.label_ind = 0

    def fit_phases(self, q, d, phase_names):
        _dict = self.construct_setting_dict()
        self.q = q
        self.data = d 
        _dict["q"] = q.tolist()
        _dict["data"] = d.tolist()
        _dict["names"] = phase_names
        
        response = requests.post(f"{self.server_ip}/fit_with_phase",
                                 json=_dict)
        results, self.probs, self.bg = response.json()

        if self.optimize_mode == "With Uncertainty":
            print("Uncertainty not supported in server mode yet")
            print("default to simple result")
            uncer = np.zeros(8*len(results["CPs"]))

            # TODO: Implement this
            # response = requests.post("http://127.0.0.1:8080/eval_uncertainty", 
            #                      json=_dict)
            # var = response.json()
            # TODO: Move this process to backend
            # all_params = np.ravel(np.array([CS.get_eight_params(cp) for cp in results.CPs])) 
            # if np.all(var >= 0):
            #     log_uncer = np.sqrt(var)
            #     uncer = np.maximum(np.abs(np.exp(np.log(all_params + log_uncer)) - all_params),
            #                        np.abs(np.exp(np.log(all_params + log_uncer)) - all_params))
            # else:
            #     print("CAUTION: There is negative hessian value, meaning the optimization is not successsful.")
            #     print("         The uncertainty is default to 50% of the parameters")
            #     uncer = np.array([ 0.5*param if v != 0 else 0 for v, param in zip(var, all_params) ])
        else:
            uncer = np.zeros(8*len(results["CPs"]))
            

        self.print_fitted_result(results, uncer)

        # print("############## Output ################")
        # print("")
        # print("Refinement result:")
        # for i, cp in enumerate(result.CPs): 
        #     if self.optimize_mode == "With Uncertainty": 
        #         self.print_refined_result_w_uncer(cp, uncer[i*8:(i+1)*8])
        #     else:
        #         self.print_refined_result(cp)

        # print("")
        # print("Fractions:")
        # for i, xi in enumerate(fractions):
        #     print(f"    {result.CPs[i].name}: {xi}")
        # print("")
        # print("#####################################")

        self.results = [results]
        self.has_labeled = True
        self.label_ind = 0
        if self.background_option in ["None", "Default"]:
            self.bg = np.zeros(q.shape)
        return results, uncer 

    def get_peaks_at(self, idx: int):
        return requests.get(f"{self.server_ip}/peak_info/{idx}").json()
        

    @property
    def residual(self):
        spectrum = evaluate_obj(self.results[self.label_ind], self.q)
        d = self.data - spectrum
        d -= self.bg
        return d

    @property
    def fit_result(self):
        return self.evaluate_fitted(self.label_ind, self.q)

    @property
    def server_ip(self):
        return f"http://{self.address}:{self.port}"


    def next_label_result(self):
        self.label_ind += 1
        if self.label_ind == len(self.results):
            self.label_ind = 0
        return (self.label_ind + 1,
                self.probs[self.label_ind],
                self.fit_result,
                self.get_fraction(self.results[self.label_ind]["CPs"]),
                self.bg)

    def previous_label_result(self):
        self.label_ind -= 1
        if self.label_ind < 0:
            self.label_ind = len(self.results) - 1
        return (self.label_ind + 1,
                self.probs[self.label_ind],
                self.fit_result,
                self.get_fraction(self.results[self.label_ind]["CPs"]),
                self.bg)

    @property
    def params(self):
        return (self.std_noise, self.mean_θ, self.std_θ, self.max_phase,
                self.expand_k, self.background_length, self.max_iter,
                self.optimize_mode, self.background_option)

    def get_subset_phases_idx_from_names(self, phase_names):
        idx = []
        for i, phase_name in enumerate(self.phase_names):
            if phase_name in phase_names:
                idx.append(i)

        return idx

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



    def get_dict_for_storing(self):
        datadict = {}
        datadict['q'] = self.q.tolist()
        datadict['XRD'] = self.data.tolist()

        res = self.results[self.label_ind]
        if res.background is not None:
            datadict['background'] = (
                evaluate_obj(
                    res.background,
                    self.q) +
                self.bg).tolist()
        else:
            datadict['background'] = self.bg.tolist()

        fractions = get_fraction(res.CPs)
        phase_dict = {}
        for idx, phase in enumerate(res.CPs):
            phase_dict[phase.name] = {}
            phase_dict[phase.name]["lattice"] = self.get_dict_from_cl(phase.cl)
            phase_dict[phase.name]["ref_lattice"] = self.get_dict_from_cl(
                phase.origin_cl)
            phase_dict[phase.name]["pattern"] = evaluate_obj(
                phase, self.q).tolist()
            phase_dict[phase.name]["fraction"] = fractions[idx]
            phase_dict[phase.name]["act"] = phase.act
            phase_dict[phase.name]["width"] = phase.σ
        datadict['phase'] = phase_dict
        return datadict

    def get_dict_from_cl(self, cl):
        c_dict = {}
        c_dict["a"] = cl.a
        c_dict["b"] = cl.b
        c_dict["c"] = cl.c
        c_dict["α"] = cl.α * 180 / np.pi
        c_dict["β"] = cl.β * 180 / np.pi
        c_dict["γ"] = cl.γ * 180 / np.pi
        return c_dict

    def get_optimize_enum(self, optimize_mode_str):
        if optimize_mode_str == "Simple":
            return CS.Simple
        if optimize_mode_str == "EM":
            return CS.EM
        if optimize_mode_str == "With Uncertainty":
            return CS.WithUncer

    def print_refined_result(self, CP):
        a_strain, b_strain, c_strain, α_strain, β_strain, γ_strain = CS.get_strain(CP)
        print(f"Phase name: {CP.name}, ID: {CP.id}")
        print(f"Optimization parameters:")
        print(f"Activation: {CP.act}, Peak width: {CP.σ}")
        print(f"Normalization: {CP.norm_constant}")
        print(f"Lattice information:")
        print(f"a: {CP.cl.a:.6f}, b: {CP.cl.b:.6f}, c: {CP.cl.c:.6f}")
        print(f"α: {CP.cl.α/np.pi*180:.2f}, β: {CP.cl.β/np.pi*180:2f}, γ: {CP.cl.γ/np.pi*180:.2f}")
        print(f"Strain: a:{a_strain:.4f}, b:{b_strain:.4f}, c:{c_strain:.4f}, α:{α_strain:.4f}, β:{β_strain:.4f}, γ:{γ_strain:.4f}")
        print("")


    def print_refined_result_w_uncer(self, CP, uncer):
        a_strain, b_strain, c_strain, α_strain, β_strain, γ_strain = CS.get_strain(CP)
        print(f"Phase name: {CP.name}, ID: {CP.id}")
        print(f"Optimization parameters:")
        print(f"Activation: {CP.act:.4f}±{uncer[6]:.4f}, Peak width: {CP.σ:.4f}±{uncer[7]:.4f}")
        print(f"Normalization: {CP.norm_constant}")
        print(f"Lattice information:")
        print(f"a: {CP.cl.a:.6f}±{uncer[0]:.6f}, b: {CP.cl.b:.6f}±{uncer[1]:.6f}, c: {CP.cl.c:.6f}±{uncer[2]:.6f}")
        print(f"α: {CP.cl.α/np.pi*180:.2f}±{uncer[3]/np.pi*180:.2f}, β: {CP.cl.β/np.pi*180:.2f}±{uncer[4]/np.pi*180:.2f}, γ: {CP.cl.γ/np.pi*180:.2f}±{uncer[5]/np.pi*180:.2f}")
        print(f"Strain: a:{a_strain:.4f}, b:{b_strain:.4f}, c:{c_strain:.4f}, α:{α_strain:.4f}, β:{β_strain:.4f}, γ:{γ_strain:.4f}")
        print("")

    def get_fraction(self, CPs):
        fraction = []
        for CP in CPs:
            fraction.append(self.get_moles(CP))
        fraction = np.array(fraction)
        fraction /= np.sum(fraction)
        return fraction

    def get_strain(self, CP):
        a = (CP["cl"]["a"] - CP["origin_cl"]["a"]) / CP["origin_cl"]["a"]
        b = (CP["cl"]["b"] - CP["origin_cl"]["b"]) / CP["origin_cl"]["b"]
        c = (CP["cl"]["c"] - CP["origin_cl"]["c"]) / CP["origin_cl"]["c"]
        α = (CP["cl"]["α"] - CP["origin_cl"]["α"]) / CP["origin_cl"]["α"]
        β = (CP["cl"]["β"] - CP["origin_cl"]["β"]) / CP["origin_cl"]["β"]
        γ = (CP["cl"]["γ"] - CP["origin_cl"]["γ"]) / CP["origin_cl"]["γ"]
        return a, b, c, α, β, γ

    def get_moles(self, CP):
        return CP["act"] * self.get_n(CP["profile"]["α"], CP["σ"]) / CP["norm_constant"]

    def get_n(self, α, σ) :
        return α*np.pi*σ + (1 - α)*σ*np.sqrt(2*np.pi)


    ######### API call helpers #######
    def evaluate(self, idx, q):
        _dict = {'q': q.tolist()}
        return np.array(requests.get(f"{self.server_ip}/evaluate/{idx}", json=_dict).json())

    def evaluate_fitted(self, idx, q):
        _dict = {'q': q.tolist()}
        return requests.get(f"{self.server_ip}/evaluate_fitted/{idx}", json=_dict).json()

    def evaluate_fitted_background(self, idx, q):
        _dict = {'q': q.tolist()}
        return np.array(requests.get(f"{self.server_ip}/evaluate_fitted_background/{idx}", json=_dict).json())

    ######### Print utils ########
    def print_current_CPs(self, uncer=None):
        self.print_CPs(self.results[self.label_ind]["CPs"], uncer)

    def print_CPs(self, CPs, uncer=None):
        i = 0
        for _, CP in enumerate(CPs):
            a, b, c, α, β, γ = self.get_strain(CP)
            print(f"Phase name: {CP['name']}, ID: {CP['id']}")
            print("Optimization parameters:")

            if uncer is not None:
                print(f"Activation: {CP['act']:.4f}±{uncer[i+6]:.4f}, Peak width: {CP['σ']:.4f}±{uncer[i+7]:.4f}")
            else:
                print(f"Activation: {CP['act']:.4f}, Peak width: {CP['σ']:.4f}")
            print(f"Normalization: {CP['norm_constant']}")
            print("Lattice information:")
            if uncer is not None:
                print(f"a: {CP['cl']['a']}±{uncer[i]:.6f}, b: {CP['cl']['b']}±{uncer[i+1]:.6f}, c: {CP['cl']['c']}±{uncer[i+2]:.6f}")
                print(f"α: {CP['cl']['α']/np.pi*180:.2f}±{uncer[i+3]/np.pi*180:.2f}, β: {CP['cl']['β']/np.pi*180:.2f}±{uncer[i+4]/np.pi*180:.2f}, γ: {CP['cl']['γ']/np.pi*180:.2f}±{uncer[i+5]/np.pi*180:.2f}")
                print(f"Strain: a:{a:.4f}, b:{b:.4f}, c:{c:.4f}, α:{α:.4f}, β:{β:.4f}, γ:{γ:.4f}")
                print("")
                i += 8
                continue

            print(f"a: {CP['cl']['a']}, b: {CP['cl']['b']}, c: {CP['cl']['c']}")
            print(f"α: {CP['cl']['α']/np.pi*180}, β: {CP['cl']['β']/np.pi*180}, γ: {CP['cl']['γ']/np.pi*180}")
            print(f"Strain: a:{a:.4f}, b:{b:.4f}, c:{c:.4f}, α:{α:.4f}, β:{β:.4f}, γ:{γ:.4f}")
            print("")


    def print_fitted_result(self, result, uncer=None):
        fractions = self.get_fraction(result["CPs"])
        print("############## Output ################")
        print("")
        print("Labeling result:")
        self.print_CPs(result["CPs"], uncer)
        print("")
        print(f"Probability: {self.probs[0]}")
        print("Fractions:")
        for i, xi in enumerate(fractions):
            print(f"    {result['CPs'][i]['name']}: {xi}")
        print("") 
        print("#######################################")
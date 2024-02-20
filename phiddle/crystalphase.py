import time

from pymatgen.io.cif import CifParser
from xrayutilities.materials.cif import CIFFile
from xrayutilities.materials.material import Crystal
from xrayutilities.simpack import PowderDiffraction
import numpy as np
import matplotlib.pyplot as plt

from peak_profile import Gauss
from crystal import get_crystal
from util import _get_strain


class CrystalPhase():

    def __init__(self, cl, origin_cl, name, peaks, peakprofile, norm_constant,
                 act = 1.0, width = 0.1, _id=0):
        """
        Input parameter:

        cl: any object inherited the BaseCrystal class
        peaks: array-like object that as the dimension of nx4
        """

        self.cl = cl # Might need to do deep copy here
        self.origin_cl = origin_cl
        self.name = name
        self.peaks = peaks
        self.peakprofile = peakprofile
        self.norm_constant = norm_constant
        self.act = act 
        self.width = width 
        self.id = _id


    def __repr__(self):
        return (f"Phase name: {self.name}, ID: {self.id}\n"
                "Optimization parameters:\n"
                f"Activation: {self.act}, Peak width: {self.width}\n"
                f"Normalization: {self.norm_constant}\n"
                "Lattice information:\n"
                f"{self.cl}\n")
        

    def __call__(self, x):

        s = np.zeros_like(x)
        for peak in self.peaks: # vectorize this later
            h, k, l, I = peak
            peak_loc = self.cl.get_peak(h, k, l)*10
            s += self.act * I * self.peakprofile((x - peak_loc)/self.width)

        return s

    @property
    def free_param_num(self):
        return self.cl.free_param_num + self.peakprofile.free_param_num + 2

    def get_free_params(self):
        p = self.cl.get_free_params()
        p.extend([self.act, self.width])
        return p


    def set_params(self, params):
        assert len(params) == self.free_param_num
        idx = self.cl.free_param_num
        self.cl.set_params(*params[:idx])
        self.act = params[idx]
        self.width = params[idx+1]
        if self.peakprofile.free_param_num > 0:
            self.peakprofile.set_param(params[idx+2])



    def get_strain(self):

        return [ _get_strain(self.cl.a, self.origin_cl.a),
                 _get_strain(self.cl.b, self.origin_cl.b),
                 _get_strain(self.cl.c, self.origin_cl.c),
                 _get_strain(self.cl.α, self.origin_cl.α),
                 _get_strain(self.cl.β, self.origin_cl.β),
                 _get_strain(self.cl.γ, self.origin_cl.γ) ]


def get_crystalphase_from_cif(cif_path,
                              q_range,
                              peakprofile,
                              act = 1.0,
                              width = 0.1,
                              _id=0,
                              wvlen=1.2782,
                              int_thresh=1e-4):

    lattice = CIFFile(cif_path).SGLattice()
    crystal = Crystal('test', lattice)
    xrd = PowderDiffraction(crystal, wl=wvlen).data # wvlen is used to compute debye-waller
    cif = CifParser(cif_path)
    cif_dict = cif.as_dict()
    phase_name = _get_phase_name(cif_dict[list(cif_dict)[0]])

    # hkls = []
    # Is = []
    peaks = []

    for hkl in xrd:
        q = xrd[hkl]['qpos']*10 # Default to use nm-1 unit
        I = xrd[hkl]['r']
        if q_range[0] < q < q_range[1] and I>int_thresh:
            hkl = list(hkl)
            hkl.append(I)
            peaks.append(hkl)
            # hkls.append(hkl)
            # Is.append(I)

    lps = [crystal.a, crystal.b, crystal.c, crystal.alpha, crystal.beta, crystal.gamma]
    peaks = np.array(peaks)
    norm_constant = np.max(peaks[:,-1])
    peaks[:,-1] /= norm_constant

    crystal = get_crystal(*lps)
    crystal_origin = get_crystal(*lps)

    return CrystalPhase(crystal, crystal_origin, phase_name, peaks, peakprofile, norm_constant, act, width, _id)


def _get_phase_name(info_dict):

    if "_chemical_formula_structural" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_structural"])
    elif "_chemical_formula_moiety" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_moiety"])
    elif "_chemical_formula_sum" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_sum"])
    else:
        print("cannot find phase name in cif")
    try:
        space_group = remove_blank(info_dict["_space_group_name_H-M_alt"])
    except KeyError:
        print(f"No _space_group_name_H-M_alt for {phase_name}")
        space_group = ""

    return phase_name + "_" + space_group


def remove_blank(string):
    while ' ' in string:
        string = string.replace(' ', '')
    return string
        

def collect_free_params(cs):
    params = []
    for i, c in enumerate(cs):
        params.extend(c.get_free_params())

    return np.array(params)

def evaluate(cs: list, x):
    y = np.zeros_like(x)
    return evaluate_in_place(y, cs, x)

def evaluate_in_place(y, cs:list, x):
    for i, c in enumerate(cs):
        y += c(x)
    return y


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    wvlen = 1.5406
    _id = 0
    cif_path = '/Users/ming/Downloads/cif/Fe2O3_mp-19770_computed.cif'
    peakprofile = Gauss()
    cs = get_crystalphase_from_cif(cif_path, (8, 60), peakprofile, 1.0, 0.3, _id, wvlen)
    print(cs.peaks.shape)
    print(cs.cl)

    x = np.linspace(8, 60, 1024)

    start = time.time()
    plt.plot(x, cs(x))
    plt.show()
    print(time.time()-start)

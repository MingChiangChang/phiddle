# Module to turn cif file into input text file that
# Can be turned into CrystalPhase object with the constructor
# CrystalPhase(path_of_cif)
# File start with the crystal information
# Crystal system, a=?, b=?, c=?, alpha=?, beta=?, gamma=?
# and follow by
# the information of the sticks as follows
# h k l intensity
# For each sticks

import glob
import os
from pathlib import Path

from pymatgen.io.cif import CifParser
from xrayutilities.materials.cif import CIFFile
from xrayutilities.materials.material import Crystal
from xrayutilities.simpack import PowderDiffraction
import numpy as np


def cif_to_input(cif_paths, output_path, q_range,
                 wvlen=0.15406, _type=float):
    '''
    cif_to_input(cif_path, output_path, q_range, output_name='sticks')

    Main function that you should be interfacing with this module
    '''
    if not output_path.endswith('.csv'):
        output_path += '.csv'
    with open(output_path, 'w') as f:
        for idx, cif_path in enumerate(cif_paths):
            cif = CifParser(cif_path)
            cif_dict = cif.as_dict()
            key = _get_key(cif_dict)
            info_dict = cif_dict[key]
            info_dict['file_name'] = os.path.basename(cif_path)[:-4]
            lattice = CIFFile(cif_path).SGLattice()
            write_cif(f, idx, info_dict, lattice, _type, q_range, wvlen)


def write_cif(f, idx, info_dict, lattice, _type, q_range, wvlen):
    f.write(f'{idx},')
    write_crystal_info(f, info_dict, _type)
    write_peaks_info(f, lattice, q_range, wvlen)


def write_crystal_info(f, info_dict, _type):
    f.write(_get_phase_name(info_dict, default=info_dict['file_name']))
    f.write(',')
    f.write(_get_crystal_system(info_dict))
    f.write(',')
    f.write(_get_lattice_parameters(info_dict, _type))


def write_peaks_info(f, lattice, q_range, wvlen):
    crystal = Crystal('test', lattice)
    xrd = PowderDiffraction(crystal).data

    qs = []
    Is = []
    for peak in xrd:
        q = xrd[peak]['qpos'] * 10
        I = xrd[peak]['r']
        qs.append(q)
        Is.append(I)
    Is = np.array(Is)
    Is /= np.max(Is)

    for peak, q, I in zip(xrd, qs, Is):
        if q_range[0] < q < q_range[1] and I > 0.0001:
            f.write(f'\n{peak[0]},{peak[1]},{peak[2]},{q},{I}')
    f.write('#\n')


def q_to_two_theta(wvlen, *args):
    two_thetas = []
    for q in args:
        two_thetas.append(360 * np.arcsin(wvlen * q / (4 * np.pi)) / np.pi)
    return tuple(two_thetas)


def _get_key(cif_dict):
    return list(cif_dict.keys())[0]


def _get_phase_name(info_dict, default=None):
    # TODO: Generalize this
    print(info_dict)
    if "_chemical_formula_structural" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_structural"])
    elif "_chemical_formula_moiety" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_moiety"])
    elif "_chemical_formula_sum" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_sum"])
    elif default is not None:
        phase_name = default
        print("cannot find phase name in cif, default to file name")

    if "_symmetry_space_group_name_H-M" in info_dict:
        space_group = remove_blank(info_dict["_symmetry_space_group_name_H-M" ])
    elif "_space_group_name_H-M_alt" in info_dict:
        space_group = remove_blank(info_dict["_space_group_name_H-M_alt"])
    else: 
        print(f"No _space_group_name_H-M_alt for {phase_name}")
        space_group = ""
    return phase_name + "_" + space_group


def _get_lattice_parameters(info_dict, _type):
    a, b, c = _get_cell_length(info_dict)
    alpha, beta, gamma = _get_cell_angle(info_dict)
    print(a, b, c, alpha, beta, gamma)
    return (f"{_type(a)},{_type(b)},{_type(c)},"
            f"{_type(alpha)},{_type(beta)},{_type(gamma)}")


def _get_cell_length(info_dict):
    return (remove_parentheses(info_dict['_cell_length_a']),
            remove_parentheses(info_dict['_cell_length_b']),
            remove_parentheses(info_dict['_cell_length_c']))


def remove_parentheses(string):
    if '(' in string:
        return string[:string.index('(')]
    return string


def remove_blank(string):
    while ' ' in string:
        string = string.replace(' ', '')
    return string


def _get_cell_angle(info_dict):
    return (remove_parentheses(info_dict['_cell_angle_alpha']),
            remove_parentheses(info_dict['_cell_angle_beta']),
            remove_parentheses(info_dict['_cell_angle_gamma']))


def _get_crystal_system(info_dict):
    try:
        sg_num = int(info_dict['_space_group_IT_number'])
    except KeyError:
        print("No space group info in cif. Default to triclinic")
        return "triclinic"
    if sg_num in [1, 2]:
        return "triclinic"
    elif 3 <= sg_num <= 15:
        return "monoclinic"
    elif 16 <= sg_num <= 74:
        return "orthohombic"
    elif 75 <= sg_num <= 142:
        return "tetragonal"
    elif 143 <= sg_num <= 167:
        return "trigonal"
    elif 168 <= sg_num <= 194:
        return "hexagonal"
    elif 195 <= sg_num <= 230:
        return "cubic"


if __name__ == "__main__":
    home = Path.home()
    # path = home / 'Desktop' / 'github' /\
    #        'Crystallography_based_shifting' / 'data'
    path = home / 'Downloads' / 'CIF-3'
    path = home / "Downloads" / "CrFeV_toCornell" / "icdd"
    path = home / 'Desktop' / 'Code' / 'CrystalShift.jl' / 'data' / 'calibration'
    path = home / "Downloads" / "tio2"
    path = home / "Downloads" / "test"
    path = home / "Downloads" / "YourCustomFileName"
    path = home / "Downloads" / "CIFs-2" / "LaOx"
    path = home / "Desktop" / "All_CIFs"
    path = home / "Downloads" / "ino"
    path = home / "Desktop" / "Code" / "SARA.jl" / "BiTiO" / "cifs"
    path = home / "Downloads" / "toCornell_Ming" / "cifs"
    path = home / "Downloads" / "al2o3"
    # path = home / "Desktop" / "TaSnCoO" / "cifs"
    # path = home / "Downloads" / "igzo"
    # path = home / "Downloads" / "AlLiFeO copy"
    cif_paths = list(path.glob('*.cif'))
    # cif_paths = path.glob("Ta-Sn-O/*/*.cif")

    # cif_paths = [str(path / 'Bi2Ti2O7_ICSD.cif') , str(path / 'Delta.cif')]
    out_path = path  # / 'Ta-Sn-O'
    print(cif_paths)
    print(out_path)
    cif_to_input(cif_paths, str(out_path), (10, 80))

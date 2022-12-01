from scipy.constants import c, e, m_e, epsilon_0, pi
from glob import glob
import numpy as np


def get_f1(element, keV):
    file_path = glob(f'./rsc/f1/*{element}*')[0]
    _, element_keV, element_f1 = np.loadtxt(
        file_path, skiprows=1, delimiter=',').T
    return np.interp([keV],  element_keV, element_f1)[0]


def xray_scat_length(f1):
    return e**2 / (4 * pi * epsilon_0 * m_e * c**2) * f1


def get_cell_volume(cif_file):
    """
    -----
    ToDo:
    if _cell_volume not found, calculate from lat params
    """
    tag = '_cell_volume'
    with open(cif_file, 'r') as f:
        lines = f.read().split('\n')
    cell_volume = [i for i in lines if i.startswith(tag)]
    cell_volume = cell_volume[0].split()[-1].split('(')[0]
    return cell_volume


def get_sym_constraints(cif_file):
    with open(cif_file, 'r') as f:
        lines = f.read().split('\n')
    xyz_index = [i for i, l in enumerate(lines) if l.startswith('loop')]
    xyz_slice = slice(xyz_index[0]+2, xyz_index[1])
    xyz_lines = lines[xyz_slice]
    # clean up the xyz lines remove tabs, commas, and quotes
    for ch in [' ', '\t', "'"]:
        xyz_lines = [line.replace(ch, '') for line in xyz_lines]
    transforms = [lambda x, y, z, coef=i: eval(f'[{coef}]') for i in xyz_lines]
    return transforms


def get_scale(recipe):
    """
    ------
    get Scale from recipe
    """
    ...


def main():
    print(get_sym_constraints('/home/cipmin/5_Felix/GitHub_EzFit/cifs/NiO.cif'))


if __name__ == '__main__':
    main()

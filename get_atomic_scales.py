from scipy.constants import c, e, m_e, epsilon_0, pi
from glob import glob
import numpy as np
import pandas as pd


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
    xyz_index1 = [i for i, l in enumerate(lines) if '_xyz' in l]
    xyz_index2 = [i for i, l in enumerate(lines) if i > xyz_index1[0]
                  and 'loop_' in l]
    xyz_lines = lines[xyz_index1[0]+1:xyz_index2[0]]

    if '_space_group_symop_id' in lines:
        xyz_lines = [xyz[len(f'{i+1}'):] for i, xyz in enumerate(xyz_lines)]

    # clean up the xyz lines remove tabs, commas, and quotes
    for ch in [' ', '\t', "'"]:
        xyz_lines = [line.replace(ch, '') for line in xyz_lines] 
    out = [[], [], []]
    for xyz in xyz_lines:
        for comp in xyz.split(','):
            comp_ind = [(ind, comp) for ind, char in enumerate(['x', 'y', 'z'])
                        if char in comp]
            for ind, comp in comp_ind:
                out[ind].append(comp)
    xyz_lines = [','.join(comp) for comp in np.array(out).T]
    transforms = [lambda x, y, z, coef=i: eval(f'[{coef}]') for i in xyz_lines]
    return transforms


def get_fract_cord(cif_file):
    with open(cif_file, 'r') as f:
        lines = f.read().split('\n')
    site_ind = [
        ind for ind, line in enumerate(lines) if line.startswith('_atom_site_')
    ]
    site_slice = slice(site_ind[0], site_ind[-1]+1)
    col_names = [line.split('_atom_site_')[-1] for line in lines[site_slice]]
    l_site = site_ind[-1] + 1
    number_symb_ind = [l_site + i for i, line in enumerate(lines[l_site + 1:])
                       if '#' in line]
    if number_symb_ind:
        data = [line.split() for line in lines[site_ind[-1]+1:
                number_symb_ind[0]+1]]
    else:
        data = [line.split() for line in lines[site_ind[-1]+1:]]
    df = pd.DataFrame(data, columns=col_names)
    return df


def gen_site_pos(df, transforms):
    frac_label = {}
    for lab in df.label.unique():  
        row = df.loc[df.label == lab, :]
        x, y, z = [row[xyz].to_numpy()[0]
                   for xyz in ['fract_x', 'fract_y', 'fract_z']]

        x, y, z = [float(i.split('(')[0]) for i in [x, y, z]]
        frac = np.unique([t(x, y, z) for t in transforms], axis=0)
        frac = frac % 1
        frac = np.unique(frac, axis=0)
        frac_label[lab] = frac
    return frac_label


def count_atom(frac_label, df):
    atom_multiplicity = {}
    for label, pos in frac_label.items():
        atom_multiplicity[label] = np.array(pos).shape[0]
    for label, occ in zip(df.label, df.occupancy):
        atom_multiplicity[label] *= float(occ)
    return atom_multiplicity


def get_scale(phases, atom_multiplicity, cell_volume, f1):
    phase_molar_contrib = {}
    for phase, scale in phases.keys():
        molar_contrib = 0
        for label, mult in atom_multiplicity.items():
            molar_contrib += (1 / (mult / cell_volume * f1[label]**2))
        molar_contrib *= scale
        phase_molar_contrib[phase] = molar_contrib
    phase_molar_contrib = {k: v / sum(phase_molar_contrib.values()) for k, v in phase_molar_contrib.items()}
    return phase_molar_contrib


def get_atomic_scales(cif_file, phases, keV):
    cell_volume = float(get_cell_volume(cif_file))
    transforms = get_sym_constraints(cif_file)
    df = get_fract_cord(cif_file)
    frac_label = gen_site_pos(df, transforms)
    atom_multiplicity = count_atom(frac_label, df)
    f1 = {label: get_f1(label, keV) for label in df.label.unique()}
    get_scale(phases, atom_multiplicity, cell_volume, f1)

def main():
    cif_file = '/home/cipmin/5_Felix/GitHub_EzFit/cifs/NiO.cif'
    df = get_fract_cord(cif_file)
    transforms = get_sym_constraints(cif_file)
    frac_label = gen_site_pos(df, transforms)
    count_atom(frac_label)
    # get scale for eacj phase


if __name__ == '__main__':
    main()

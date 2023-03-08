import pandas as pd
import numpy as np
from scipy.constants import c, e, m_e, epsilon_0, pi
from glob import glob
from molmass import Formula

class GetScales():
    '''
    Class

    '''

    def calc_scale(self):
        self.get_cifs()
        self.get_scales()
        self.get_cell_volume()
        self.get_sym_constraints()
        self.get_fract_cord()
        self.gen_site_pos()
        self.count_atom()
        self.get_number_density()
        self.get_xray_scat_len()
        self.get_avg_scat_len()
        return self.get_real_scales(), self.get_weight_percent()

    def get_scales(self):
        '''
        Function
        '''
        self.scales = {}
        for phase in self.phases:
            scale = getattr(self.recipe, phase + '_scale').value
            self.scales[phase] = scale
        return self.scales

    def get_cifs(self):
        self.cif_contents = {}
        for phase in self.phases:
            with open(self.cif_files[phase], 'r') as f:
                lines = f.readlines()
            lines = [line.strip() for line in lines]
            self.cif_contents[phase] = lines
        return self.cif_contents

    def get_cell_volume(self):
        """
        -----
        ToDo:
        if _cell_volume not found, calculate from lat params
        """
        self.cell_volume = {}
        tag = '_cell_volume'
        for phase in self.phases:
            lines = self.cif_contents[phase]
            cell_line = [i for i in lines if i.startswith(tag)]
            self.cell_volume[phase] = float(cell_line[0].split()[-1].split('(')[0])
            self.cell_volume[phase] = self.cell_volume[phase] * 10**-30
        return self.cell_volume

    def get_sym_constraints(self):
        """
        args:
        cif_file
        returns:
        symmetrie operations
        """
        self.transforms = {}
        for phase in self.phases:
            lines = self.cif_contents[phase]
            xyz_index1 = [i for i, l in enumerate(lines) if '_xyz' in l]
            xyz_index2 = [i for i, l in enumerate(lines) if i > xyz_index1[0] and 'loop_' in l]
            xyz_lines = lines[xyz_index1[0]+1:xyz_index2[0]]
            if '_space_group_symop_id' in lines:
                xyz_lines = [xyz[len(f'{i+1}'):] for i, xyz in enumerate(xyz_lines)]
            for ch in [' ', '\t', "'"]:
                xyz_lines = [line.replace(ch, '') for line in xyz_lines] 
            self.transforms[phase] = [lambda x, y, z, coef=i: eval(f'[{coef}]') for i in xyz_lines] 
        return self.transforms

    def get_fract_cord(self):
        """
        args:
        cif_file
        returns:
        dataframe with fractional coordinates
        """
        self.fract_cord = {}
        for phase in self.phases:
            lines = self.cif_contents[phase]
            site_ind = [
                ind for ind, line in enumerate(lines) if '_atom_site_' in line
            ]
            site_slice = slice(site_ind[0], site_ind[-1]+1)
            col_names = [line.split('_atom_site_')[-1] for line in lines[site_slice]]
            l_site = site_ind[-1] + 1
            number_symb_ind = [l_site + i for i, line in enumerate(lines[l_site + 1:])
                            if '#' in line]
            if number_symb_ind:
                data = [line.split() for line in lines[site_ind[-1]+1:number_symb_ind[0]] if line != '']
            else:
                data = [line.split() for line in lines[site_ind[-1]+1:] if line != '']
            self.fract_cord[phase] = pd.DataFrame(data, columns=col_names)
        return self.fract_cord

    def gen_site_pos(self):
        """
        args:
        df: dataframe with fractional coordinates,
        transforms: symmetrie operations
        returns:
        site_pos: dictionary with phase as keys, labels as keys and fractional coordinates as values
        """
        self.site_pos = {}
        for phase in self.phases:
            self.site_pos[phase] = {}
            for lab in self.fract_cord[phase].label.unique(): 
                row = self.fract_cord[phase].loc[self.fract_cord[phase].label == lab, :]
                x, y, z = [row[xyz].to_numpy()[0]
                        for xyz in ['fract_x', 'fract_y', 'fract_z']]
                x, y, z = [float(i.split('(')[0]) for i in [x, y, z]]
                frac = np.unique([t(x, y, z) for t in self.transforms[phase]], axis=0)
                frac = frac % 1
                for i, (x, y, z) in enumerate(frac):
                    x = round(x, 4)
                    y = round(y, 4)
                    z = round(z, 4)
                    frac[i] = [x, y, z]
                frac = np.unique(frac, axis=0)
                self.site_pos[phase][lab] = frac
        return self.site_pos

    def count_atom(self):
        """
        args: 
        df: dataframe with fractional coordinates,
        site_pos: dictionary with labels as keys and fractional coordinates as values
        returns:
        atom_mult: dictionary with phase as key, labels as keys and number of atoms as values
        """
        self.atom_mult = {}
        for phase in self.phases:
            self.atom_mult[phase] = {}
            for label, pos in self.site_pos[phase].items():
                self.atom_mult[phase][label] = np.array(pos).shape[0]
            for label, occ in zip(self.fract_cord[phase].label, self.fract_cord[phase].occupancy):
                self.atom_mult[phase][label] *= float(occ)
        return self.atom_mult

    def get_number_density(self):
        """
        args:
        cell_volume
        atom_multiplicity
        returns:
        n_x: number density of element x
        """
        self.number_dens = {}
        for phase in self.phases:
            n = 0  #n = total amount of atoms
            n_x = list(self.atom_mult[phase].values())  #number of element x in unit cell
            dens = sum(n_x) / self.cell_volume[phase]
            self.number_dens[phase] = dens
        return self.number_dens

    def get_xray_scat_len(self):
        """
        args: 
        atom_multiplicity: we get the atom label and use it to track the atomic form factor
        xray_energy: still needs to be put in manually
        returns:
        f1: dictionary, keys = phases, values = list of real part of atomic form factor for each atom site
        b1: dictionary, keys = phases, values = list of scattering length for each atom site
        """
        f1 = {}
        self.b1 = {}
        elements = {}
        for phase in self.phases:
            elements[phase] = []
            f1[phase] = []
            self.b1[phase] = []
            for keys in self.atom_mult[phase]:
                element = str()
                for symbol in keys:
                    if symbol.isalpha():
                        element += symbol
                elements[phase].append(element)
            for element in elements[phase]:
                file_path = glob(f'/home/feilix/Git_Fit/ezfit/rsc/f1*{element}*')[0]
                _, element_keV, element_f1 = np.loadtxt(
                    file_path, skiprows=1, delimiter=',').T
                f = np.interp([self.config['Measurement']['keV']],  element_keV, element_f1)[0]
                bi = e**2 / (4 * pi * epsilon_0 * m_e * c**2) * f
                f1[phase].append(f)
                self.b1[phase].append(bi)
        return self.b1

    def get_avg_scat_len(self):
        """
        args:
        b1
        atom_mult
        returns:
        avg_scat_len: dictionary with phase as keys and average scatter length of unit cell as values
        """
        self.avg_scat_len = {}
        n = {}  #n = total amount of atoms
        f_sum = {}  #sum of scattering length
        for phase in self.phases:
            self.avg_scat_len[phase] = 0
            n[phase] = 0
            f_sum[phase] = 0
            n_x = list(self.atom_mult[phase].values()) #number of element x in unit cell
            for i in range(len(self.b1[phase])):
                f_sum[phase] += self.b1[phase][i] * n_x[i]
                n[phase] += n_x[i]
            #print('f_sum', f_sum, 'n', n)
            self.avg_scat_len[phase] = f_sum[phase] / n[phase]
        return self.avg_scat_len
        
    def get_real_scales(self):
        """
        args: 
        scale
        avg_scat_len
        number_dens
        returns:
        real_scale
        """
        self.real_scales = {}
        s = {}
        for phase in self.phases:
            s[phase] = self.scales[phase] / (self.number_dens[phase] * self.avg_scat_len[phase]**2)
        sum_scales = sum(s.values())
        for phase in self.phases:
            self.real_scales[phase] = s[phase] / sum_scales
        return self.real_scales

    def get_weight_percent(self):
        self.weight_percent = {}
        total = 0
        for phase in self.phases:
            formula = self.formulas[phase]
            total += Formula(formula).isotope.mass * self.real_scales[phase]
        for phase in self.phases:
            formula = self.formulas[phase]
            self.weight_percent[phase] = Formula(formula).isotope.mass * self.real_scales[phase] / total
        return self.weight_percent



import diffpy.srfit.pdf.characteristicfunctions as CF
from diffpy.srfit.fitbase import FitResults
from itertools import count
from pathlib import Path
from typing import List
import numpy as np
from . import diffpy_wrap as dw
import toml
from ezfit.get_scales import GetScales



def _read_config(config_location: str = None):
    if config_location is None:
        cwd = Path().resolve()
        config_path = list(Path(cwd).glob('*.toml'))[0]
    else:
        config_path = Path(config_location).expanduser().resolve()
    config: dict = toml.load(config_path)
    return config


def validate_file_locations(config: dict):
    for key, val in config['files'].items():
        if not Path(val).exists():
            raise FileNotFoundError(f"File {val} does not exist.")


def _phase_counter(self, phase):
    if not hasattr(self, f'{phase}_count'):
        setattr(self, f'{phase}_count', count(1))
        return 0
    else:
        return next(getattr(self, f'{phase}_count'))


def _parse_phases(self, phases):
    phases_cp = phases.copy()
    for i, phase in enumerate(phases):
        ocur = phases.count(phase)
        if ocur > 1:
            cnt = _phase_counter(self, phase)
            phases_cp[i] = f'{phase}Γ{cnt}'

    return phases_cp


def _fetch_function(phase, function):
    func_param = {
        'sphericalCF':
        (CF.sphericalCF, ['r', f'{phase}_psize']),
        'spheroidalCF':
        (CF.spheroidalCF, ['r', f'{phase}_erad', f'{phase}_prad']),
        'spheroidalCF2':
        (CF.spheroidalCF2, ['r', f'{phase}_psize', f'{phase}_axrat']),
        'lognormalSphericalCF':
        (CF.lognormalSphericalCF, ['r', f'{phase}_psize', f'{phase}_sig']),
        'sheetCF':
        (CF.sheetCF, ['r', f'{phase}_sthick']),
        'shellCF':
        (CF.shellCF, ['r', f'{phase}_radius', f'{phase}_thickness']),
        'shellCF2':
        (CF.shellCF, ['r', f'{phase}_a', f'{phase}_delta']),
        'bulfkCF':
        (lambda r: r * 1, ['r']),
        }
    return func_param[function]


def create_cif_files_string(phases, config):
    cif_files = {}
    for phase in list(phases):
        cif_files[f'{phase}'] = f'{config["files"]["cifs"]}{phase.split("Γ")[0]}.cif'

    return cif_files


def create_equation_string(phases, nanoparticle_shapes):
    equation_list = []
    for phase, function in zip(phases, nanoparticle_shapes):
        equation_list.append(f'{phase} * {phase}{function}')
    equation = ' + '.join(equation_list)

    return equation


def create_functions(phases, nanoparticle_shapes):
    functions = {}
    for phase, function in zip(phases, nanoparticle_shapes):
        function_definition = _fetch_function(phase, function)
        functions[f'{phase}{function}'] = function_definition

    return functions


class FitPDF(GetScales):
    def __init__(
            self,
            file: str,
            phasesetup: dict,
            config_location: str = None,
    ):
        self.phases = list(phasesetup.keys())
        self.nanoparticle_shapes = np.array(list(phasesetup.values())).T[0]
        self.formulas = dict(zip(self.phases, np.array(list(phasesetup.values())).T[1]))
        self.file = file
        self.phases = _parse_phases(self, self.phases)

        self.config = _read_config(config_location)
        validate_file_locations(self.config)

        self.cif_files = create_cif_files_string(self.phases, self.config)
        self.equation = create_equation_string(self.phases, self.nanoparticle_shapes)
        self.functions = create_functions(self.phases, self.nanoparticle_shapes)

    def update_recipe(self):
        self.recipe, self.pg = dw.create_recipe_from_files(
             data_file=self.file,
             meta_data=self.config['PDF'],
             equation=self.equation,
             cif_files=self.cif_files,
             functions=self.functions
        )

    def apply_restraints(self):
        recipe = self.recipe
        for phase in self.phases:

            delta2 = getattr(self.recipe, f'{phase}_delta2')
            recipe.restrain(delta2, lb=1, ub=5, sig=1e-3)
            delta2.value = 3

            scale = getattr(self.recipe, f'{phase}_scale')
            recipe.restrain(scale, lb=0.01, ub=2, sig=1e-3)
            scale.value = 1

            for abc in ['a', 'b', 'c']:
                try:
                    lat = getattr(self.recipe, f'{phase}_{abc}')
                    lb_lat = lat.value - 0.5
                    ub_lat = lat.value + 0.5
                    recipe.restrain(lat, lb=lb_lat, ub=ub_lat, sig=1e-3)
                except AttributeError:
                    pass

            for func in self.functions.values():
                params = func[1][1:]
                for p in params:
                    param = getattr(self.recipe, p)
                    recipe.restrain(param, lb=10, ub=100, sig=1e-3)
                    param.value = 50

    def create_param_order(self):
        nCF = []
        for phase in self.phases:
            for func in self.functions.values():
                for varn in func[1][1:]:
                    nCF.append(varn)
        nCF = [n for n in nCF if n]
        self.param_order = [
            ['free', 'lat', 'scale'],
            ['free', *nCF],
            ['free', 'delta2', 'adp'],
        ]

    def run_fit(self):
        #self.update_recipe()
        self.apply_restraints()
        self.create_param_order()
        dw.optimize_params_manually(
            self.recipe,
            self.param_order,
            rmin=self.config['R_val']['rmin'],
            rmax=self.config['R_val']['rmax'],
            rstep=self.config['R_val']['rstep'],
            ftol=1e-5,
            print_step=self.config['Verbose']['step'],
        )
        res = FitResults(self.recipe)
        if self.config['Verbose']['results']:
            res.printResults()
        molscale, weighscale = self.calc_scale()
        print('Mol Scales:\n', [f'{k} = {v:1.3}' for k, v in molscale.items()])
        print('Weight Scales:\n', [f'{k} = {v:1.3}' for k, v in weighscale.items()])
        
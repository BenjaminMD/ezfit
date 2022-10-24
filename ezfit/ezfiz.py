from distutils.command.config import config
import diffpy.srfit.pdf.characteristicfunctions as CF
from itertools import count
from pathlib import Path
import ezfit
from typing import List
import toml



def _read_config(config_location):
    if config_location is None:
        print('No config file location provided. Finding in "./"')
        cwd = Path().resolve()
        config_path = Path(cwd / 'FitPDF_config.toml')
    else:
        config_path = Path(config_location).expanduser().resolve()

    config: dict = toml.load(config_path)

    for key, val in config['files'].items():
        if not Path(val).exists():
            raise FileNotFoundError(f"File {val} does not exist.")

    return config


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
        'blukCF':
        (lambda r: r * 1, ['r']),
        }
    return func_param[function]


def create_cif_files(phases):
    cif_files = {}
    for phase in list(phases):
        cif_files[f'{phase}'] = f'./CIFS/{phase.split("Γ")[0]}.cif'

    return cif_files


def create_equation(phases, nanoparticle_shapes):
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


class FitPDF():
    def __init__(
            self,
            phases: List[str],
            nanoparticle_shapes: List[str],
            config_location: str = None,
    ):
        phases = _parse_phases(self, phases)

        self.config = _read_config(config_location)
        self.cif_files = create_cif_files(phases)
        self.equation = create_equation(phases, nanoparticle_shapes)
        self.functions = create_functions(phases, nanoparticle_shapes)

    def update_recipe(self):
        print(ezfit)
        # self.recipe = ezfit.dw.create_recipe_from_files(
        #     data_file=self.file,
        #     meta_data=config['PDF'],
        #     equation=self.equation,
        #     cif_files=self.cif_files,
        #     functions=self.functions
        # )
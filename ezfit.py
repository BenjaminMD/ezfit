import diffpy.srfit.pdf.characteristicfunctions as CF
from diffpy.srfit.fitbase import FitResults
from itertools import count
from pathlib import Path
from typing import List
import numpy as np
from . import diffpy_wrap as dw
from .contribution import Contribution
import toml






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
        'bulkCF':
        (lambda r: 1, ['r']),
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


class FitPDF():
    def __init__(
            self,
            file: str,
            contributions: List[Contribution],
            config_location: str = None,
    ):
        
        self.phases = [contribution.cif_name for contribution in contributions]
        self.nanoparticle_shapes = [contribution.cf_name for contribution in contributions]
        self.formulas = {contribution.cif_name: contribution.formula for contribution in contributions}
        
        
        self.file = file
        self.phases = _parse_phases(self, self.phases)

        self.config = self.load_toml_config(config_location)

        self.cif_files = create_cif_files_string(self.phases, self.config)
        self.equation = create_equation_string(self.phases, self.nanoparticle_shapes)
        self.functions = create_functions(self.phases, self.nanoparticle_shapes)
        self.dw = dw

    def load_toml_config(self, config_location: str = None):
        if config_location is None:
            cwd = Path().resolve()
            config_path = list(Path(cwd).glob('*.toml'))[0]
        else:
            config_path = Path(config_location).expanduser().resolve()
        config: dict = toml.load(config_path)
        return config


    def update_recipe(self):

        self.recipe, self.pgs = dw.create_recipe_from_files(
             data_file=self.file,
             meta_data=self.config['PDF'],
             equation=self.equation,
             cif_files=self.cif_files,
             functions=self.functions
        )
        if not self.config['PDF']:
            self.add_instr_params()

    def add_instr_params(self) -> None:
        print('attempting to fit instrumental parameters')
        if len(self.pgs) != 1:
            raise ValueError('determining instrument param only one pg allowed')
        pg = list(self.pgs.values())[0]
        self.recipe.addVar(
                pg.qdamp,
                name='qdamp',
                value=0.1,
                fixed=True,
                tags='qdamp'
            ).boundRange(0.)

        self.recipe.addVar(
            pg.qbroad,
            name='qbroad',
            value=0.1,
            fixed=True,
            tags='qbroad'
            ).boundRange(0.)
    
        self.config['param_order'][-1]['free'].extend(['qdamp', 'qbroad'])

    def apply_restraints(self):
        recipe = self.recipe
        for phase in self.phases:

            delta2 = getattr(self.recipe, f'{phase}_delta2')
            recipe.restrain(delta2, lb=0, ub=5, sig=1e-3)
            delta2.value = 2

            
            # biso stuff --------
            fc = recipe.PDF
            pg = getattr(fc, phase)
            atoms: typing.List[ParameterSet] = pg.phase.getScatterers()

            for atom in atoms:
                site = atom.name 
                Biso = getattr(self.recipe, f'{phase}_{site}_Biso')
                recipe.restrain(Biso, lb=0.1, ub=2, sig=1e-3)
            
            # biso stuff --------

            scale = getattr(self.recipe, f'{phase}_scale')
            recipe.restrain(scale, lb=0.01, ub=2, sig=1e-3)
            scale.value = 1
            
            
            recipe.fix('all')
            recipe.free('occ')

            for occ_name in recipe.getNames():
                occ = getattr(recipe, occ_name)
                recipe.restrain(occ, lb=0.0, ub=1, sig=1e-3)
            recipe.fix('all')
            
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
        for func in self.functions.values():
            for varn in func[1][1:]:
                nCF.append(varn)
        nCF = [n for n in nCF if n]

        for order in self.config['param_order']:
            if "cfs" in order["free"]:
                id = order["free"].index("cfs")
                order["free"].pop(id)
                order["free"].extend(nCF)


    def run_fit(self):
        self.apply_restraints()
        self.create_param_order()
        dw.optimize_params(
            self.recipe,
            self.config['param_order'],
            rmin=self.config['R_val']['rmin'],
            rmax=self.config['R_val']['rmax'],
            rstep=self.config['R_val']['rstep'],
            ftol=1e-5,
            print_step=self.config['Verbose']['step'],
        )
        self.res = FitResults(self.recipe)
        if self.config['Verbose']['results']:
            self.res.printResults()
        return self.res

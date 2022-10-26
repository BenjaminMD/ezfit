import diffpy.srfit.pdf.characteristicfunctions as CF
from diffpy.srfit.fitbase import FitResults
from dataclasses import dataclass
from itertools import count
import rsc.diffpy_wrap as dw
import re


@dataclass(frozen=True)
class MetaData:
    qbroad: float = 1.59e-02
    qdamp: float = 1.07e-02

    def __call__(self) -> str:
        return(self.__dict__)
#hallo#
    def fetch_function(self, phase, function):
        print(f"Fetching {function} for {phase}")
        func_param = {
            'sphericalCF':
            (CF.sphericalCF, ['r', f'{phase}_psize']),
            'spheroidalCF':
            (CF.spheroidalCF, ['r', f'{phase}_erad', f'{phase}_prad']),
            'spheroidalCF2':
            (CF.spheroidalCF2, ['r', f'{phase}_psize', f'{phase}_axrat']),
            'lognormalSphericalCF':
            (CF.lognormalSphericalCF, ['r', f'{phase}_psize', f'{phase}_psig']),
            'sheetCF':
            (CF.sheetCF, ['r', f'{phase}_sthick']),
            'shellCF':
            (CF.shellCF, ['r', f'{phase}_radius', f'{phase}_thickness']),
            'shellCF2':
                (CF.shellCF, ['r', f'{phase}_a', f'{phase}_delta']),
            'bulkCF':
                (lambda r: r * 1, ['r']),
            }

        return func_param[function]


class RecipeClass(MetaData):
    def __init__(self, file, phases, nano_shape):
        self.file = file
        self.phases = self._parse_phases(phases)
        self.nano_shape = nano_shape

        self.create_cif_files()
        self.create_equation()
        self.create_functions()

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
                cnt = self._phase_counter(phase)
                phases_cp[i] = f'{phase}Γ{cnt}'

        return phases_cp

    def create_cif_files(self):
        self.cif_files = {}
        for phase in list(self.phases):
            self.cif_files[f'{phase}'] = f'./CIFS/{phase.split("Γ")[0]}.cif'

    def create_equation(self):
        equation_list = []
        for phase, function in zip(self.phases, self.nano_shape):
            equation_list.append(f'{phase} * {phase}{function}')
        self.equation = ' + '.join(equation_list)

    def create_functions(self):
        self.functions = {}
        for phase, function in zip(self.phases, self.nano_shape):
            function_definition = MetaData().fetch_function(phase, function)
            self.functions[f'{phase}{function}'] = function_definition

    def update_recipe(self):
        self.recipe = dw.create_recipe_from_files(
            data_file=self.file,
            meta_data=MetaData()(),
            equation=self.equation,
            cif_files=self.cif_files,
            functions=self.functions
        )

    def apply_restraints(self):
        recipe = self.recipe
        for phase in self.phases:

            delta2 = getattr(self.recipe, f'{phase}_delta2')
            recipe.restrain(delta2, lb=1, ub=5, sig=1e-3)
            delta2.value = 3.0

            scale = getattr(self.recipe, f'{phase}_scale')
            recipe.restrain(scale, lb=0.01, ub=2, sig=1e-3)
            scale.value = 0.5

            for abc in ['a', 'b', 'c']:
                try:
                    lat = getattr(self.recipe, f'{phase}_{abc}')
                    lb_lat = lat.value - 0.2
                    ub_lat = lat.value + 0.2
                    recipe.restrain(lat, lb=lb_lat, ub=ub_lat, sig=1e-3)
                except AttributeError:
                    pass

            for func in self.functions.values():
                params = func[1][1:]
                for p in params:
                    param = getattr(self.recipe, p)


    def create_param_order(self):
        ns = []
        for phase in self.phases:
            for func in self.functions.values():
                for varn in func[1][1:]:
                    ns.append(varn)
        ns = [n for n in ns if n]
        self.param_order = [
            ['free', 'lat', 'scale'],
            ['free', *ns],
            ['free', 'adp', 'delta2'],
            ['free', 'all']
        ]
    
    def run_fit(self):
        self.update_recipe()
        self.apply_restraints()
        self.create_param_order()
        dw.optimize_params_manually(
            self.recipe,
            self.param_order,
            rmin=1,
            rmax=8,
            rstep=0.01,
            ftol=1e-5,
            print_step=True
        )
        res = FitResults(self.recipe)
        res.printResults()


def get_gr(recipe):
    """
    Get the gr of a recipe and for each phase contribution
    returns:
    - r: list of floats
    - gobs: list of floats
    - gcalc: list of floats
    - gdiff: list of floats
    - baseline: float 
    - gr_composition: dict of list of floats
    """
    def remove_consecutive_duplicates(string, char):
        indices = [m.start() for m in re.finditer(char * 2, string)]
        if indices:
            for i in indices:
                string = string[:i] + string[i+1:]
            return remove_consecutive_duplicates(string, char)
        else:
            return string

    equation = recipe.PDF.getEquation()        
    for char in ['\)', '\(']:
        equation = (remove_consecutive_duplicates(equation, char))
    

    prof = recipe._contributions['PDF'].profile
    r = prof.x
    gobs = prof.y
    gcalc = recipe._contributions['PDF'].evaluate()
    baseline = 1.35 * gobs.min()
    gdiff = gobs - gcalc


    gr_composition = {}
    for eq in equation.split(' + '):
        gr = recipe.PDF.evaluateEquation(eq[1:])
        gr_composition[eq[1:]] = gr

    return r, gobs, gcalc, gdiff, baseline, gr_composition



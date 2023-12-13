import diffpy.srfit.pdf.characteristicfunctions as CF
from diffpy.srfit.fitbase import FitResults, initializeRecipe
from itertools import count
from pathlib import Path
from typing import List
from . import diffpy_wrap as dw
from .contribution import Contribution
import toml
from .get_scales import GetScales
from .ezconstraints import Ezrestraint


def _phase_counter(self, phase):
    if not hasattr(self, f"{phase}_count"):
        setattr(self, f"{phase}_count", count(1))
        return 0
    else:
        return next(getattr(self, f"{phase}_count"))


def _parse_phases(self, phases):
    phases_cp = phases.copy()
    for i, phase in enumerate(phases):
        ocur = phases.count(phase)
        if ocur > 1:
            cnt = _phase_counter(self, phase)
            phases_cp[i] = f"{phase}Γ{cnt}"

    return phases_cp


def _fetch_function(phase, function):
    func_param = {
        "sphericalCF":
            (CF.sphericalCF, ["r", f"{phase}_psize"]),
        "spheroidalCF":
            (CF.spheroidalCF, ["r", f"{phase}_erad", f"{phase}_prad"]),
        "spheroidalCF2":
            (CF.spheroidalCF2, ["r", f"{phase}_psize", f"{phase}_axrat"]),
        "lognormalsphericalCF":
            (CF.lognormalSphericalCF, ["r", f"{phase}_psize", f"{phase}_sig"]),
        "sheetCF":
            (CF.sheetCF, ["r", f"{phase}_sthick"]),
        "shellCF":
            (CF.shellCF, ["r", f"{phase}_radius", f"{phase}_thickness"]),
        "shellCF2":
            (CF.shellCF, ["r", f"{phase}_a", f"{phase}_delta"]),
        "bulkCF":
            (lambda r: 1, ["r"]),
    }
    return func_param[function]


def create_cif_files_string(phases, config):
    cif_files = {}
    for phase in list(phases):
        cif = config["files"]["cifs"]
        cif_files[f"{phase}"] = f'{cif}{phase.split("Γ")[0]}.cif'

    return cif_files


def create_equation_string(phases, nanoparticle_shapes):
    equation_list = []
    for phase, function in zip(phases, nanoparticle_shapes):
        equation_list.append(f"{phase} * {phase}{function}")
    equation = " + ".join(equation_list)

    return equation


def create_functions(phases, nanoparticle_shapes):
    functions = {}
    for phase, function in zip(phases, nanoparticle_shapes):
        function_definition = _fetch_function(phase, function)
        functions[f"{phase}{function}"] = function_definition

    return functions


class FitPDF(Ezrestraint, GetScales):
    def __init__(
        self,
        file: str,
        contributions: List[Contribution],
        config_location: str = "",
    ):

        self.phases = [contribution.cif_name for contribution in contributions]
        self.nanoparticle_shapes = [
            contribution.cf_name for contribution in contributions
        ]
        self.formulas = {
            contribution.cif_name: contribution.formula
            for contribution in contributions
        }

        self.file = file
        self.phases = _parse_phases(self, self.phases)

        self.config = self.load_toml_config(config_location)

        self.cif_files = create_cif_files_string(self.phases, self.config)
        self.equation = create_equation_string(
            self.phases,
            self.nanoparticle_shapes
        )
        self.functions = create_functions(
            self.phases,
            self.nanoparticle_shapes
        )
        self.dw = dw

    def load_toml_config(self, config_location: str = ""):
        if config_location:
            config_path = Path(config_location).expanduser().resolve()
        else:
            cwd = Path().resolve()
            config_path = list(Path(cwd).glob("*.toml"))[0]
        config: dict = toml.load(config_path)
        return config

    def update_recipe(self):

        self.recipe, self.pgs = dw.create_recipe_from_files(
            data_file=self.file,
            meta_data=self.config["PDF"],
            equation=self.equation,
            cif_files=self.cif_files,
            functions=self.functions,
        )
        if not self.config["PDF"]:
            self.add_instr_params()

        self.fc = self.recipe.PDF

    def add_instr_params(self) -> None:
        print("attempting to fit instrumental parameters")
        if len(self.pgs) != 1:
            raise ValueError("instrument param only one pg allowed")
        pg = list(self.pgs.values())[0]
        self.recipe.addVar(
            pg.qdamp, name="qdamp", value=0.1, fixed=True, tags="qdamp"
        ).boundRange(0.0)  # type: ignore

        self.recipe.addVar(
            pg.qbroad, name="qbroad", value=0.1, fixed=True, tags="qbroad"
        ).boundRange(0.0)  # type: ignore

        self.config["param_order"][-1]["free"].extend(["qdamp", "qbroad"])

    def apply_restraints(self):
        for param in self.config["Restraints"].keys():
            self.restrain_param(param, self.config)
        for phase in self.phases:
            self.shared_occ(phase)


    def create_param_order(self):
        nCF = []
        for func in self.functions.values():
            for varn in func[1][1:]:
                nCF.append(varn)
        nCF = [n for n in nCF if n]

        for order in self.config["param_order"]:
            if "cfs" in order["free"]:
                id = order["free"].index("cfs")
                order["free"].pop(id)
                order["free"].extend(nCF)
                
    def LoadResFromFile(self, path_to_results: str):
        initializeRecipe(self.recipe, path_to_results)
        

    def run_fit(self):
        self.apply_restraints()
        self.create_param_order()
        dw.optimize_params(
            self.recipe,
            self.config["param_order"],
            rmin=self.config["R_val"]["rmin"],
            rmax=self.config["R_val"]["rmax"],
            rstep=self.config["R_val"]["rstep"],
            ftol=1e-5,
            print_step=self.config["Verbose"]["step"],
        )
        self.res = FitResults(self.recipe)
#        self.molscale, self.weighscale = self.calc_scale()
#        self.all_scales = {'mol_scale': self.molscale, 'wt_scale': self.weighscale}
#        print('Mol Scales:\n', [f'{k} = {v:1.3}' for k, v in self.molscale.items()])
#        print('Weight Scales:\n', [f'{k} = {v:1.3}' for k, v in self.weighscale.items()])
        if self.config["Verbose"]["results"]:
            self.res.printResults()
        return self.res

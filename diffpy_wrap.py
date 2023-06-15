import typing
from pathlib import Path
import numpy as np
from diffpy.srfit.fitbase import FitContribution, FitRecipe, Profile
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.fitbase import FitResults
from ezpdf import get_gr
from pyobjcryst import loadCrystal
from pyobjcryst.crystal import Crystal
from scipy.optimize import least_squares


def _create_recipe(
        equation: str,
        crystals: typing.Dict[str, Crystal],
        functions: typing.Dict[
                str, typing.Tuple[typing.Callable, typing.List[str]]
            ],
        profile: Profile,
        fc_name: str = "PDF"
) -> FitRecipe:
    pgs = {}
    fr = FitRecipe()
    fc = FitContribution(fc_name)
    for name, crystal in crystals.items():
        pg = PDFGenerator(name)
        pg.setStructure(crystal, periodic=True)
        pg.parallel(32)
        #pg._calc.evaluatortype = 'OPTIMIZED'
        #pg._calc.evaluatortype = 'BASIC'
        fc.addProfileGenerator(pg)

        pgs[name] = pg

    if functions:
        for name, (f, argnames) in functions.items():
            fc.registerFunction(f, name=name, argnames=argnames)
    fc.setEquation(equation)
    fc.setProfile(profile, xname="r", yname="G", dyname="dG")
    fr.addContribution(fc)
    return fr, pgs


def _get_tags(phase: str, param: str) -> typing.List[str]:
    return [param, phase, "{}_{}".format(phase, param)]


def _get_name(*args: str) -> str:
    return "_".join(args)


def _rename_par(name: str, atoms: list) -> str:
    parts = name.split("_")
    np = len(parts)
    na = len(atoms)
    if np > 1 and parts[1].isdigit() and -1 < int(parts[1]) < na:
        parts[1] = atoms[int(parts[1])].name
        parts = parts[::-1]
    return "_".join(parts)


def _add_params_in_pg(recipe: FitRecipe, pg: PDFGenerator, meta_data) -> None:
    name: str = pg.name
    if not meta_data:
        print('No meta data found for {}'.format(name))

    recipe.addVar(
        pg.scale,
        name=_get_name(name, "scale"),
        value=0.,
        fixed=True,
        tags=_get_tags(name, "scale")
    ).boundRange(0.)
    recipe.addVar(
        pg.delta2,
        name=_get_name(name, "delta2"),
        value=0.,
        fixed=True,
        tags=_get_tags(name, "delta2")
    ).boundRange(0.)
    latpars = pg.phase.sgpars.latpars
    for par in latpars:
        recipe.addVar(
            par,
            name=_get_name(name, par.name),
            fixed=True,
            tags=_get_tags(name, "lat")
        ).boundRange(0.)
    atoms: typing.List[ParameterSet] = pg.phase.getScatterers()
    for atom in atoms:
        par = atom.Biso
        recipe.addVar(
            par,
            name=_get_name(name, atom.name, "Biso"),
            value=0.02,
            fixed=True,
            tags=_get_tags(name, "adp")
        ).boundRange(0.)
    xyzpars = pg.phase.sgpars.xyzpars

    for par in xyzpars:
        par_name = _rename_par(par.name, atoms)
        try:
            recipe.addVar(
                par,
                name=_get_name(name, par_name),
                fixed=True,
                tags=_get_tags(name, "xyz")
            )
        except ValueError:
            print(name, par_name, 'already constrained')
    for atom in atoms:
        atom_type = filter(lambda x: x.isalpha(), atom.name)
        atom_type = ''.join(atom_type)
        tag_list = [
                *_get_tags(name, "occ"),
                *_get_tags(name, f"occ_{atom_type}")
            ]
        par = atom.occ
        recipe.addVar(
            par,
            name=_get_name(name, atom.name, "occ"),
            fixed=True,
            tags=np.unique(tag_list)
        ).boundRange(0.)
    return


def _add_params_in_fc(
        recipe: FitRecipe,
        fc: FitContribution,
        names: typing.List[str],
        tags: typing.List[str]
) -> None:
    for name in names:
        par = getattr(fc, name)
        recipe.addVar(
            par,
            value=10.,
            fixed=True,
            tags=tags
        )
    return


def _initialize_recipe(
        recipe: FitRecipe,
        functions: typing.Dict[
                str, typing.Tuple[typing.Callable, typing.List[str]]
            ],
        crystals: typing.Dict[str, Crystal],
        fc_name: str = "PDF",
        meta_data=None
) -> None:

    fc: FitContribution = getattr(recipe, fc_name)
    if functions:
        for name, (_, argnames) in functions.items():
            _add_params_in_fc(recipe, fc, argnames[1:], tags=[name, "cfs"])
    for name in crystals.keys():
        pg: PDFGenerator = getattr(fc, name)
        _add_params_in_pg(recipe, pg, meta_data)
    recipe.clearFitHooks()
    return


def create_recipe_from_files(
        equation: str,
        cif_files: typing.Dict[str, str],
        data_file: str,
        functions: typing.Dict[
                str, typing.Tuple[typing.Callable, typing.List[str]]
            ] = {},
        meta_data: typing.Dict[str, typing.Union[str, int, float]] = None,
        fc_name: str = "PDF"
) -> typing.Tuple[FitRecipe, typing.Dict[str, PDFGenerator]]:

    if meta_data is None:
        meta_data = {}
    crystals = {n: loadCrystal(f) for n, f in cif_files.items()}
    pp = PDFParser()
    pp.parseFile(data_file)
    profile = Profile()
    profile.loadParsedData(pp)
    profile.meta.update(meta_data)
    recipe, pgs = _create_recipe(
        equation, crystals, functions, profile, fc_name=fc_name
    )
    _initialize_recipe(
        recipe, functions, crystals, fc_name=fc_name, meta_data=meta_data
    )
    return recipe, pgs


def optimize_params(
    recipe: FitRecipe,
    steps: typing.List[typing.List[str]],
    rmin: float = None,
    rmax: float = None,
    rstep: float = None,
    print_step: bool = True,
    fc_name: str = "PDF",
    **kwargs
) -> None:

    n = len(steps)
    free_steps = [order["free"] for order in steps]
    fix_steps = [order["fix"] for order in steps]

    fc: FitContribution = getattr(recipe, fc_name)
    p: Profile = fc.profile
    p.setCalculationRange(xmin=rmin, xmax=rmax, dx=rstep)
    for step in free_steps:
        recipe.fix(*step)
    for i, (free_step, fix_step) in enumerate(zip(free_steps, fix_steps)):
        if fix_step:
            recipe.fix(*fix_step)
        if free_step:
            recipe.free(*free_step)
        if print_step:
            print(
                "Step {} / {}: params {}".format(
                    i + 1, n, ", ".join(recipe.getNames())
                ),
                end="\r"
            )
        least_squares(
            recipe.residual,
            recipe.getValues(),
            bounds=recipe.getBounds2(),
            **kwargs
        )
    return


def save_results(
        recipe: FitRecipe,
        footer: str,
        directory: str,
        file_stem: str,
        pg_names: typing.List[str] = None,
        fc_name: str = "PDF"
) -> None:
    """Save the parameters, fits and structures in the FitRecipe object.

    Parameters
    ----------
    recipe :
        The FitRecipe object.
    directory :
        The directory to output the files.
    file_stem :
        The stem of the filename.
    pg_names :
        The name of the PDFGenerators to save. If None, not to save.
    fc_name
        The name of the FitContribution in the FitRecipe. Default "PDF".
    Returns
    -------
    None.
    """
    r, gobs, gcalc, gdiff, baseline, gr_composition = get_gr(recipe)
    phases = list(gr_composition.keys())
    header = ['r', 'g(r)', 'g(r)_calc', 'g(r)_diff', *phases]
    data = np.vstack([r, gobs, gcalc, gdiff, *gr_composition.values()]).T

    d_path = Path(directory)
    d_path.mkdir(parents=True, exist_ok=True)
    f_path = d_path.joinpath(file_stem)
    fr = FitResults(recipe)
    fr.saveResults(str(f_path.with_suffix(".res")), footer=f'{footer}')
    fc: FitContribution = getattr(recipe, fc_name)
    # profile: Profile = fc.profile

    np.savetxt(
        str(f_path.with_suffix(".fgr")),
        data, header=';'.join(header),
        comments='', delimiter=';'
    )
    # profile.savetxt(str(f_path.with_suffix(".fgr")))
    if pg_names is not None:
        for pg_name in pg_names:
            pg: PDFGenerator = getattr(fc, pg_name)
            stru: Crystal = pg.stru
            cif_path = f_path.with_name(
                "{}_{}".format(f_path.stem, pg_name)
            ).with_suffix(".cif")
            with cif_path.open("w") as f:
                stru.CIFOutput(f)
    return

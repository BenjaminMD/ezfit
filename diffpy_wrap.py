import typing
from pathlib import Path
import numpy as np
from diffpy.srfit.fitbase import FitContribution, FitRecipe, Profile
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.fitbase import FitResults
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
    fr = FitRecipe()
    fc = FitContribution(fc_name)
    for name, crystal in crystals.items():
        pg = PDFGenerator(name)
        pg.setStructure(crystal, periodic=True)
        pg.parallel(32)

        fc.addProfileGenerator(pg)
    if functions:
        for name, (f, argnames) in functions.items():
            fc.registerFunction(f, name=name, argnames=argnames)
    fc.setEquation(equation)
    fc.setProfile(profile, xname="r", yname="G", dyname="dG")
    fr.addContribution(fc)
    return fr, pg


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
            pg.qdamp,
            name='qdamp',
            value=0.,
            fixed=True,
            tags='qdamp'
        ).boundRange(0.)

        recipe.addVar(
            pg.qbroad,
            name='qbroad',
            value=0.,
            fixed=True,
            tags='qbroad'
            ).boundRange(0.)

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
        recipe.addVar(
            par,
            name=_get_name(name, par_name),
            fixed=True,
            tags=_get_tags(name, "xyz")
        )
    for atom in atoms:
        atom_type = filter(lambda x: x.isalpha(), atom.name)
        atom_type = ''.join(atom_type)
        tag_list = [
                *_get_tags(name, "occ"),
                *_get_tags(name, f"occ{atom_type}")
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
            _add_params_in_fc(recipe, fc, argnames[1:], tags=[name])
    for name in crystals.keys():
        pg: PDFGenerator = getattr(fc, name)
        _add_params_in_pg(recipe, pg, meta_data)
    recipe.clearFitHooks()
    return pg


def create_recipe_from_files(
        equation: str,
        cif_files: typing.Dict[str, str],
        data_file: typing.Dict[str, str],
        functions: typing.Dict[
                str, typing.Tuple[typing.Callable, typing.List[str]]
            ] = {},
        meta_data: typing.Dict[str, typing.Union[str, int, float]] = None,
        fc_name: str = "PDF"
) -> FitRecipe:
    if meta_data is None:
        meta_data = {}
    crystals = {n: loadCrystal(f) for n, f in cif_files.items()}
    pp = PDFParser()
    pp.parseFile(data_file)
    profile = Profile()
    profile.loadParsedData(pp)
    profile.meta.update(meta_data)
    recipe, pg = _create_recipe(
        equation, crystals, functions, profile, fc_name=fc_name
    )
    pg = _initialize_recipe(
        recipe, functions, crystals, fc_name=fc_name, meta_data=meta_data
    )
    return recipe  , pg


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
    fc: FitContribution = getattr(recipe, fc_name)
    p: Profile = fc.profile
    p.setCalculationRange(xmin=rmin, xmax=rmax, dx=rstep)
    for step in steps:
        recipe.fix(*step)
    for i, step in enumerate(steps):
        recipe.free(*step)
        if print_step:
            print(
                "Step {} / {}: refine {}".format(
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


def optimize_params_manually(
    recipe: FitRecipe,
    steps: typing.List[typing.List[str]],
    rmin: float = None,
    rmax: float = None,
    rstep: float = None,
    print_step: bool = True,
    fc_name: str = "PDF",
    **kwargs
) -> None:

    from warnings import warn
    warn("Felix Fabio? Do we really need the functionality to fix a step during the fit or do we just want to free one by one", DeprecationWarning)
    n = len(steps)
    fc: FitContribution = getattr(recipe, fc_name)
    p: Profile = fc.profile
    p.setCalculationRange(xmin=rmin, xmax=rmax, dx=rstep)

    for i, step in enumerate(steps):
        eval(f'recipe.{step[0]}(*{step[1:]})')
        if print_step:
            print(
                "Step {} / {}: refine {}".format(
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
    d_path = Path(directory)
    d_path.mkdir(parents=True, exist_ok=True)
    f_path = d_path.joinpath(file_stem)
    fr = FitResults(recipe)
    fr.saveResults(str(f_path.with_suffix(".res")))
    fc: FitContribution = getattr(recipe, fc_name)
    profile: Profile = fc.profile
    profile.savetxt(str(f_path.with_suffix(".fgr")))
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


def visualize_fits(
            ax,
            recipe: FitRecipe,
            xlim: typing.Tuple = None,
            fc_name: str = "PDF"
        ) -> None:
    """Visualize the fits in the FitRecipe object.

    Parameters
    ----------
    recipe :
        The FitRecipe object.
    xlim :
        The boundary of the x to show in the plot.
    fc_name :
        The name of the FitContribution in the FitRecipe. Default "PDF".

    Returns
    -------
    None.
    """
    # get data
    fc = getattr(recipe, fc_name)
    r = fc.profile.x
    g = fc.profile.y
    gcalc = fc.profile.ycalc
    if xlim is not None:
        sel = np.logical_and(r >= xlim[0], r <= xlim[1])
        r = r[sel]
        g = g[sel]
        gcalc = gcalc[sel]
    gdiff = g - gcalc
    diffzero = -0.8 * np.max(g) * np.ones_like(g)
    # plot figure
    ax.plot(r, g, 'bo', label="G(r) Data")
    ax.plot(r, gcalc, 'r-', label="G(r) Fit")
    ax.plot(r, gdiff + diffzero, 'g-', label="G(r) Diff")
    ax.plot(r, diffzero, 'k-')
    ax.set_xlabel(r"$r (\AA)$")
    ax.set_ylabel(r"$G (\AA^{-2})$")
    ax.legend(loc=1)
    return

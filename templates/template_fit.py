from diffpy.srfit.fitbase import FitResults
from matplotlib import pyplot as plt
from ezfit import FitPDF, Contribution
from ezpdf import plot_PDF
from pathlib import Path
from glob import glob


def get_atoms(fit, phase, element):
    fc = fit.recipe.PDF
    phase_atoms = {}
    pg = getattr(fc, phase)
    atoms = pg.phase.getScatterers()
    phase_atoms[f"{phase}"] = [atom for atom in atoms if atom.name.startswith(element)]
    return phase_atoms


def constrain_Bisos(recipe, Bisos_dict):
    for phase, atom_list in Bisos_dict.items():
        if atom_list:
            constraints = [f"{phase}_{atom.name}_Biso" for atom in atom_list]
            first_Biso = atom_list[0].Biso
            for constraint in constraints[1:]:
                recipe.constrain(constraint, first_Biso)


data_file = glob("./gr/SBa200.gr")[0]
name = Path(data_file).stem

Alox = "delta4_Al2O3"
Al2O3 = Contribution(cif_name=Alox, cf_name="sphericalCF", formula="Al2O3")
NiO = Contribution(cif_name="NiO", cf_name="sphericalCF", formula="NiO")
Ni = Contribution(cif_name="Ni", cf_name="sphericalCF", formula="Ni")

fit = FitPDF(data_file, contributions=[Al2O3])

fit.update_recipe()
fit.shared_param(Alox, "Al", "Biso")
fit.shared_param(Alox, "O", "Biso")

fit.run_fit()
fig, ax = plot_PDF(fit, fit.recipe, fit.res)
ax.set_title(f"{Alox} fit of {name}")
fit.dw.save_results(
    recipe=fit.recipe,
    footer=f"Fit of {name}",
    directory="./res_single/",
    file_stem=f"{name}{Alox}",
    pg_names=fit.phases,
)
plt.savefig(f"./res_single/{name}_{Alox}.pdf")

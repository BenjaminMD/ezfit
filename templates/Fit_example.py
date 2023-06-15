from pathlib import Path
from ezfit import FitPDF, Contribution

path_CeO2 = Path('/home/ben/DESY_PDF/0_data/config/mask_MP/CeO2_Standard_abs0_10s_00026.gr')


CeO2 = Contribution(cif_name='CeO2', cf_name='bulkCF', formula='CeO2')

fit = FitPDF(path_CeO2, contributions=[CeO2])

fit.update_recipe()
fit.apply_restraints()
fit.run_fit()
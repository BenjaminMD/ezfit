class Ezrestraint:
    def __init__(self, fit):
        self.fit = fit
        self.recipe = fit.recipe
        self.fc = fit.recipe.PDF

    def get_atoms(self, phase, element):
        phase_atoms = {}
        pg = getattr(self.fc, phase)
        atoms = pg.phase.getScatterers()
        phase_atoms[f"{phase}"] = [
            atom for atom in atoms if atom.name.startswith(element)
        ]
        return phase_atoms

    def shared_param(self, phase, element, param):
        atoms = self.get_atoms(phase, element)
        atom_list = atoms.get(phase, [])
        if not atom_list:
            return

        recipe = self.recipe
        constraints = [f"{phase}_{atom.name}_{param}" for atom in atom_list]
        first_param = getattr(atom_list[0], param)
        for constraint in constraints[1:]:
            recipe.constrain(constraint, first_param)

    def restrain_param(self, param, lb=None, ub=None, lr=None, initial=None):
        recipe = self.recipe
        recipe.fix("all")
        recipe.free(param)
        for param_name in recipe.getNames():
            param = getattr(recipe, param_name)
            if lr:
                lb = param.value - lr
                ub = param.value + lr
            recipe.restrain(param, lb=lb, ub=ub, sig=1e-3)
            if initial:
                param.value = initial
        recipe.fix("all")
        return self

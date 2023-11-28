from collections import defaultdict



class Ezrestraint:
    def __init__(self, fit):
        self.fit = fit
        self.recipe = fit.recipe
        self.fc = fit.recipe.PDF

    def _get_atoms(self, phase, element):
        phase_atoms = {}
        pg = getattr(self.fc, phase)
        atoms = pg.phase.getScatterers()
        phase_atoms[f"{phase}"] = [
            atom for atom in atoms if atom.name.startswith(element)
        ]
        return phase_atoms

    def shared_param(self, phase, element, param):
        atoms = self._get_atoms(phase, element)
        atom_list = atoms.get(phase, [])
        if not atom_list:
            return

        recipe = self.recipe
        constraints = [f"{phase}_{atom.name}_{param}" for atom in atom_list]
        first_param = getattr(atom_list[0], param)
        for constraint in constraints[1:]:
            recipe.constrain(constraint, first_param)

    def shared_occ(self, phase):
        pg = getattr(self.fc, phase)
        atoms = pg.phase.getScatterers()
        atom_xyz = {i: (i.x.value, i.y.value, i.z.value) for i in atoms}
        duplicate_keys = Ezrestraint.find_duplicate_keys(atom_xyz)
        for _, shared_atoms in duplicate_keys:
            occ_eq = []
            for atom in shared_atoms:
                occ_eq.append(f"{phase}_{atom.name}_occ")
            occ_eq = " + ".join(occ_eq)
            self.recipe.newVar(f"{phase}_{atom.name}_Occ_sum")
            self.recipe.constrain(f"{phase}_{atom.name}_Occ_sum", f'{occ_eq}')
            self.recipe.restrain(f"{phase}_{atom.name}_Occ_sum", lb=0.0, ub=1.0, sig=1e-3)

    @staticmethod   
    def find_duplicate_keys(dictionary):
        value_keys = defaultdict(list)
        duplicate_keys = []

        for key, value in dictionary.items():
            value_keys[value].append(key)

        for value, keys in value_keys.items():
            if len(keys) > 1:
                duplicate_keys.append((value, keys))

        return duplicate_keys

    def restrain_param(self, param, config, initial=None):
        lr = None
        lb = None
        ub = None
        recipe = self.recipe
        lb_ub_ini = config["Restraints"][param]
        print(lb_ub_ini)
        if type(lb_ub_ini) != float:
            if len(lb_ub_ini) == 3:
                if type(lb_ub_ini) != float:
                    lb, ub, initial = lb_ub_ini
                else:
                    lr = lb_ub_ini
            if len(lb_ub_ini) == 2:
                if type(lb_ub_ini) != float:
                    lb, ub = lb_ub_ini
                else:
                    lr = lb_ub_ini
        else:
            lr = lb_ub_ini
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

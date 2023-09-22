from dataclasses import dataclass


@dataclass
class Contribution:
    cif_name: str
    cf_name: str
    formula: str
    names: str = None
from dataclasses import dataclass


@dataclass
class Contribution:
    cif_name: str
    cf_name: str
    formula: str
    name: str = None

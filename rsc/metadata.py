import diffpy.srfit.pdf.characteristicfunctions as CF
from dataclasses import dataclass


@dataclass(frozen=True)
class MetaData:
    qbroad: float = 1.12563856e-02
    qdamp: float = 1.52097999e-02

    def __call__(self) -> str:
        return(self.__dict__)

    def fetch_function(self, phase, function):
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
            'blukCF': 
                (lambda r: r  * 1, ['r']),
            }
        return func_param[function]

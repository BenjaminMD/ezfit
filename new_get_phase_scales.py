from ezfit import ezfit
from diffpy.srfit.fitbase import FitResults
from ezfit.get_gr import get_gr
from scipy.constants import c, e, m_e, epsilon_0, pi
from glob import glob
import numpy as np
import pandas as pd


cifs = fit_pdf.cif_files
PDF = fit_pdf.recipe.PDF

class GetScales:
    '''
    Class
    '''
    def __init__(self, cifs, PDF, xray_energy):
        self.cifs = cifs
        self.PDF = PDF
        self.phases = cifs.keys()
        #self.f1 = get_f1('Ni', xray_energy)

    def get_scale_from_PDF(self):
            '''
            method in a class
            '''   
            self.scales = {}
            for phase in self.phases:
                PDF_phase = getattr(self.PDF, phase)
                scale = PDF_phase.scale.getValue()
                self.scales[phase] = scale
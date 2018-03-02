##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from Functions import *
import os
import pickle
import naima
import astropy.units as units
from naima.models import PionDecay, TableModel

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##--------------##
# Path for files #
##--------------##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files')

##-----------##
# Computation #
##-----------##

Emin = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV
number_bin_E = 50

pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/Gamma_emission/'

with open('spectra', 'rb') as spectra_read:

        # Loading data
    spectrum_energy = pickle.load(spectra_read)
    spectrum_unit = pickle.load(spectra_read)
    spectrum = spectrum_energy*spectrum_unit
    #print(spectrum_unit.unit)
    sed_PD = pickle.load(spectra_read)
    sed_PD_unit = pickle.load(spectra_read)
    lum_gamma = sed_PD * sed_PD_unit

    figure_number = 1

    n = len(spectrum_energy)

    for i in range (n):
        #lum = SPD[i,0] * units.erg*1/units.s*1/units.GeV
        fig = plt.figure(figsize=(8,5))
        plt.rc('font', family='sans')
        plt.rc('mathtext', fontset='custom')
        plt.loglog(spectrum[i,0], lum_gamma[i,0], lw=2, c=naima.plot.color_cycle[0])
        plt.title('Production of photons by Pion Decay')
        plt.xlabel('Photon energy [{0}]'.format(spectrum.unit.to_string('latex_inline')))
        #plt.ylabel('$L_E$ [{0}]'.format(spectrum_unit.unit.to_string('latex_inline')))
        plt.ylabel('$L_E$ [{0}]'.format(lum_gamma.unit.to_string('latex_inline')))
        plt.tight_layout()
        fig.savefig(pathfigure+'T%.2f.eps'%i)
        figure_number += 1

    plt.show()

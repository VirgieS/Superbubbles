##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
import scipy.integrate as integrate
import astropy.units as units
import os
import pickle
import naima
from naima.models import PionDecay, TableModel
from Functions import *
from Functions_CR import *
from Functions_SB import *
from Functions_gamma import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_system import *

ngas = [1, 10] * 1/units.cm**3 # density of target protons (cm^-3)

N_E = power_law_distribution(ECR) * 1/(units.GeV)
model = TableModel(E_CR, N_E, amplitude = 1)

    # For the first density
PD1 = PionDecay(model, nh = ngas[0], nuclear_enhancement = True, useLUT = False)
sed_PD1 = PD1.sed(spectrum_energy, distance = 0 * units.pc)

    # For the second density
PD2 = PionDecay(model, nh = ngas[1], nuclear_enhancement = True, useLUT = False)
sed_PD2 = PD2.sed(spectrum_energy, distance = 0 * units.pc)

    # Plot
figure_number = 1
number_of_plot = 2
x = spectrum
y = [numpy.asarray(sed_PD1), numpy.asarray(sed_PD2)]
title = ''
label_name = ['1 $cm^{-3}$', '10 $cm^{-3}$']
xlabel = 'E$_\gamma$ [GeV]'
ylabel = 'E$^2_\gamma$ dN/dE [erg s$^{-1}$]'
symbol = ['', '']
linestyle = ['-.', ':']
color = ['cornflowerblue', 'orangered']
text = ''
log_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig('/Users/stage/Documents/Virginie/Superbubbles/figures/Examples.pdf')
plt.show()

"""
It plots the gamma-ray luminosities in the HE and VHE energy ranges and the spectral index in both energy ranges for one or two iterations.

All the parameters must be given in the Parameters_system.

Make sure that you have already run the code 'Plotting.py' to have the concatenisation of the iterations.
"""


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

##====##
# Path #
##====##

    # you need to change it
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parametric_studies/diffusion/delta0_33/Gamma_emission/'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

    # Initialization
figure_number = 1

    ## ------- ##
    # Load data #
    ## ------- ##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parametric_studies/diffusion/delta0_33/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS = pickle.load(iteration_write)
    Lum_Fermi = pickle.load(iteration_write)
    Lum = pickle.load(iteration_write)
    Gamma_HESS = pickle.load(iteration_write)
    Gamma_GeV = pickle.load(iteration_write)
    Gamma_MeV = pickle.load(iteration_write)
    Lum_pwn = pickle.load(iteration_write)
    Lum_psr = pickle.load(iteration_write)

    # One iteration
j = 20  # one iteration

xlabel = 'Time [Myr]'

xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5

        # Gamma emission comparison
ymin = 1e29
ymax = 1e36

sym_mean = ['', '']
linestyle_mean = ['-', '-']
sym_pwn = ['', '']
linestyle_pwn = [':', ':']

color_mean = ['navy', 'orangered']
color_pwn = ['green', 'darkviolet']

            # H.E.S.S. range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

y_mean = [Lum_HESS[20], Lum_HESS[98]]
y_pwn = [Lum_pwn[20], Lum_pwn[98]]

semilog_plot(figure_number, 2, t6, y_mean, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, xmin, xmax, ymin, ymax)
semilog_plot(figure_number, 2, t6, y_pwn, xlabel, ylabel_HESS, sym_pwn, linestyle_pwn, color_pwn, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_HESS_tot.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_mean, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_HESS.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_pwn, xlabel, ylabel_HESS, sym_pwn, linestyle_pwn, color_pwn, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_pwn.pdf')
figure_number += 1

            # Fermi range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

y_mean = [Lum_Fermi[20], Lum_Fermi[98]]
y_psr = [Lum_psr[20], Lum_psr[98]]

semilog_plot(figure_number, 2, t6, y_mean, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, xmin, xmax, ymin, ymax)
semilog_plot(figure_number, 2, t6, y_psr, xlabel, ylabel_HESS, sym_pwn, linestyle_pwn, color_pwn, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_Fermi_tot.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_mean, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_Fermi.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_psr, xlabel, ylabel_HESS, sym_pwn, linestyle_pwn, color_pwn, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_psr.pdf')
figure_number += 1

label = ['VHE', 'HE']
y_mean = [Lum_HESS[20], Lum_Fermi[20]]

semilog_plot(figure_number, 2, t6, y_mean, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, xmin, xmax, ymin, ymax, label_name = label)
plt.savefig(pathfigure_gamma+'Gamma_emission_comparison.pdf')

plt.show()

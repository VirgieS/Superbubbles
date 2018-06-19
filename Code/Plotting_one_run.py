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
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/100/Total/Gamma_emission/'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

    # Initialization
figure_number = 1

    ## ------- ##
    # Load data #
    ## ------- ##

        # Nob = 30

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_30 = pickle.load(iteration_write)
    Lum_Fermi_30 = pickle.load(iteration_write)
    Lum_30 = pickle.load(iteration_write)
    Gamma_HESS_30 = pickle.load(iteration_write)
    Gamma_GeV_30 = pickle.load(iteration_write)
    Gamma_MeV_30 = pickle.load(iteration_write)
    Lum_pwn_30 = pickle.load(iteration_write)
    Lum_psr_30 = pickle.load(iteration_write)

    # One iteration
j = 20  # one iteration

xlabel = 'Time [Myr]'

text = ''
Title = ''

xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5

        # Gamma emission comparison
label_mean = 'none'
ymin = 1e29
ymax = 1e36

sym_mean = ['', '']
linestyle_mean = ['-', '-']

color_mean = ['navy', 'orangered']
color_pwn = ['green', 'darkviolet']

            # H.E.S.S. range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

y_mean = [Lum_HESS_30[20], Lum_HESS_30[98]]
y_pwn = [Lum_pwn_30[20], Lum_pwn_30[98]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
semilog_plot(figure_number, 2, t6, y_pwn, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_pwn, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_HESS_tot.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_HESS.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_pwn, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_pwn, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_pwn.pdf')
figure_number += 1

            # Fermi range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

y_mean = [Lum_Fermi_30[20], Lum_Fermi_30[98]]
y_psr = [Lum_psr_30[20], Lum_psr_30[98]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
semilog_plot(figure_number, 2, t6, y_psr, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_pwn, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_Fermi_tot.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_Fermi.pdf')
figure_number += 1

semilog_plot(figure_number, 2, t6, y_psr, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_pwn, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_psr.pdf')
figure_number += 1

        # Comparison energy range
label_comparison = ['HE CRs', 'VHE CRs']

            # Low diffusion coefficient
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$]'

y_mean = [Lum_Fermi_30[20], Lum_HESS_30[20]]

semilog_plot(figure_number, 2, t6, y_mean, label_comparison, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Gamma_emission_one_run.pdf')
figure_number += 1

    # spectral index
ymin = 0.0
ymax = 3.5

        # H.E.S.S. range
ylabel_HESS = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'

y_mean = [Gamma_HESS_30[20], Gamma_HESS_30[98]]

plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Photon_index_HESS.pdf')
figure_number += 1

        # Fermi range
ylabel_HESS = '$\Gamma_{ph}$ (1 GeV - 10 GeV)'

y_mean = [Gamma_GeV_30[20], Gamma_GeV_30[98]]

plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Photon_index_Fermi.pdf')
figure_number += 1

            # High diffusion coefficient
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$]'

y_mean = [Lum_Fermi_300[j], Lum_HESS_300[j]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'One_run_300.pdf')

plt.show()

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

## NEED TO WRITE CLEARLY WHAT I DO

    # IRAP
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/Total/Gamma_emission/'
    # Home
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/Parameters/stars/100/Gamma_emission/bis_300'
#pathfigure_remain = '/home/vivi/Documents/Master_2/Superbubbles/figures/Parameters/stars/100/Remain/bis_300'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

"""
Plot the graphics for all iterations
"""

    # Total number of iterations
nit_tot = 100                                                          #you need to change it for your simulations

    # Fix time array (yr)
tmin = 3/yr26yr         # yr
tmax = 10/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1

    ## ------- ##
    # Load data #
    ## ------- ##

        # Nob = 30

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/30/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_30 = pickle.load(iteration_write)
    Lum_Fermi_30 = pickle.load(iteration_write)
    Lum_30 = pickle.load(iteration_write)
    Gamma_HESS_30 = pickle.load(iteration_write)
    Gamma_GeV_30 = pickle.load(iteration_write)
    Gamma_MeV_30 = pickle.load(iteration_write)
    Lum_pwn_30 = pickle.load(iteration_write)
    Lum_psr_30 = pickle.load(iteration_write)
"""
        # Nob = 100

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_100 = pickle.load(iteration_write)
    Lum_Fermi_100 = pickle.load(iteration_write)
    Lum_100 = pickle.load(iteration_write)
    Gamma_HESS_100 = pickle.load(iteration_write)
    Gamma_GeV_100 = pickle.load(iteration_write)
    Gamma_MeV_100 = pickle.load(iteration_write)
    Lum_pwn_100 = pickle.load(iteration_write)
    Lum_psr_100 = pickle.load(iteration_write)

"""
        # Nob = 300

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/300/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_300 = pickle.load(iteration_write)
    Lum_Fermi_300 = pickle.load(iteration_write)
    Lum_300 = pickle.load(iteration_write)
    Gamma_HESS_300 = pickle.load(iteration_write)
    Gamma_GeV_300 = pickle.load(iteration_write)
    Gamma_MeV_300 = pickle.load(iteration_write)
    Lum_pwn_300 = pickle.load(iteration_write)
    Lum_psr_300 = pickle.load(iteration_write)


    # One iteration
j = 10  # one iteration

xlabel = 'Time [Myr]'

text = ''
Title = ''

xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5
ymin = 1e29
ymax = 1e36
        # Gamma emission comparison
label_mean = ['30', '300']

sym_mean = ['', '']
linestyle_mean = ['-', '-']

color_mean = ['purple', 'green']

            # H.E.S.S. range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

y_mean = [Lum_HESS_30[j], Lum_HESS_300[j]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'One_run_HESS.pdf')
figure_number += 1

            # Fermi range
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

y_mean = [Lum_Fermi_30[j], Lum_Fermi_300[j]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'One_run_Fermi.pdf')
figure_number += 1

        # Comparison energy range
label_mean = ['VHE CRs', 'HE CRs']

sym_mean = ['', '']
linestyle_mean = ['-', '-']

color_mean = ['orangered', 'cornflowerblue']

            # Low diffusion coefficient
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$]'

y_mean = [Lum_Fermi_30[j], Lum_HESS_30[j]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'One_run_30.pdf')
figure_number += 1

            # High diffusion coefficient
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$]'

y_mean = [Lum_Fermi_300[j], Lum_HESS_300[j]]

semilog_plot(figure_number, 2, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'One_run_300.pdf')

plt.show()

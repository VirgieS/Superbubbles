"""
It plots a comparison of each gamma-ray luminosities in the HE and VHE energy ranges and the spectral index in both energy ranges.

All the parameters must be given in the Parameters_system.

Make sure that you have already run the code 'Iterations.py' and 'Plotting.py' to have the concatenisation of the data set.
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
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/diffusion/Total/Gamma_emission/'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

    # Initialization
figure_number = 1

    ## ------- ##
    # Load data #
    ## ------- ##

        # First set for comparison

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/diffusion/1e26/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_30 = pickle.load(iteration_write)
    Lum_Fermi_30 = pickle.load(iteration_write)
    Lum_30 = pickle.load(iteration_write)
    Gamma_HESS_30 = pickle.load(iteration_write)
    Gamma_GeV_30 = pickle.load(iteration_write)
    Gamma_MeV_30 = pickle.load(iteration_write)
    Lum_pwn_30 = pickle.load(iteration_write)
    Lum_psr_30 = pickle.load(iteration_write)

##---------------------------##
# Mean and standard deviation #
##---------------------------##

    # Initialization

        # Mean
Lum_HESS_mean_30 = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_Fermi_mean_30 = numpy.zeros(number_bin_t)              # from 100 MeV to 100 GeV
Gamma_HESS_mean_30 = numpy.zeros(number_bin_t)             # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_mean_30 = numpy.zeros(number_bin_t)              # photon spectral index from 1 GeV to 10 GeV
Lum_pwn_mean_30 = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean_30 = numpy.zeros(number_bin_t)                # GeV emission of PSRs

        # Standard deviation
Lum_HESS_std_30 = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std_30 = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Gamma_HESS_std_30 = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_std_30 = numpy.zeros(number_bin_t)               # photon spectral index from 1 GeV to 10 GeV

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_30[:, j]

    Lum_Fermi = Lum_Fermi_30[:, j]

    Gamma_HESS = Gamma_HESS_30[:, j]

    Gamma_GeV = Gamma_GeV_30[:, j]

    Lum_pwn = Lum_pwn_30[:, j]
    Lum_pwn = Lum_pwn[numpy.where(Lum_pwn > 0.0)[0]]

    Lum_psr = Lum_psr_30[:, j]
    Lum_psr = Lum_psr[numpy.where(Lum_psr > 0.0)[0]]

    Lum_HESS_mean_30[j] = numpy.mean(Lum_HESS)
    Lum_Fermi_mean_30[j] = numpy.mean(Lum_Fermi)
    Gamma_HESS_mean_30[j] = numpy.mean(Gamma_HESS)
    Gamma_GeV_mean_30[j] = numpy.mean(Gamma_GeV)
    Lum_pwn_mean_30[j] = 10**(numpy.mean(numpy.log10(Lum_pwn)))
    Lum_psr_mean_30[j] = 10**(numpy.mean(numpy.log10(Lum_psr)))

    Lum_HESS_std_30[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std_30[j] = numpy.std(Lum_Fermi)
    Gamma_HESS_std_30[j] = numpy.std(Gamma_HESS)
    Gamma_GeV_std_30[j] = numpy.std(Gamma_GeV)

Lum_HESS_mean_30 = numpy.nan_to_num(Lum_HESS_mean_30)
Lum_HESS_std_30 = numpy.nan_to_num(Lum_HESS_std_30)
Lum_HESS_pst_30 = Lum_HESS_mean_30 + Lum_HESS_std_30
Lum_HESS_pst_30 = Lum_HESS_mean_30 + Lum_HESS_std_30
Lum_HESS_mst_30 = Lum_HESS_mean_30 - Lum_HESS_std_30
Lum_HESS_mst_30 = numpy.nan_to_num(Lum_HESS_mst_30)
ind0 = numpy.where(Lum_HESS_mst_30 < 0)[0]
Lum_HESS_mst_30[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_mean_30 = numpy.nan_to_num(Lum_Fermi_mean_30)
Lum_Fermi_std_30 = numpy.nan_to_num(Lum_Fermi_std_30)
Lum_Fermi_pst_30 = Lum_Fermi_mean_30 + Lum_Fermi_std_30
Lum_Fermi_pst_30 = Lum_Fermi_mean_30 + Lum_Fermi_std_30
Lum_Fermi_mst_30 = Lum_Fermi_mean_30 - Lum_Fermi_std_30
Lum_Fermi_mst_30 = numpy.nan_to_num(Lum_Fermi_mst_30)
ind0 = numpy.where(Lum_Fermi_mst_30 < 0)[0]
Lum_Fermi_mst_30[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_pst_30 = Gamma_HESS_mean_30 + Gamma_HESS_std_30
Gamma_HESS_mst_30 = Gamma_HESS_mean_30 - Gamma_HESS_std_30
ind0 = numpy.where(Gamma_HESS_mst_30 < 0)[0]
Gamma_HESS_mst_30[ind0] = numpy.zeros(len(ind0))

Gamma_GeV_pst_30 = Gamma_GeV_mean_30 + Gamma_GeV_std_30
Gamma_GeV_mst_30 = Gamma_GeV_mean_30 - Gamma_GeV_std_30
ind0 = numpy.where(Gamma_GeV_mst_30 < 0)[0]
Gamma_GeV_mst_30[ind0] = numpy.zeros(len(ind0))

        # Second data set for comparison

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/diffusion/1e27/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_100 = pickle.load(iteration_write)
    Lum_Fermi_100 = pickle.load(iteration_write)
    Lum_100 = pickle.load(iteration_write)
    Gamma_HESS_100 = pickle.load(iteration_write)
    Gamma_GeV_100 = pickle.load(iteration_write)
    Gamma_MeV_100 = pickle.load(iteration_write)
    Lum_pwn_100 = pickle.load(iteration_write)
    Lum_psr_100 = pickle.load(iteration_write)

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

    # Initialization

        # Mean
Lum_HESS_mean_100 = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_Fermi_mean_100 = numpy.zeros(number_bin_t)              # from 100 MeV to 100 GeV
Gamma_HESS_mean_100 = numpy.zeros(number_bin_t)             # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_mean_100 = numpy.zeros(number_bin_t)              # photon spectral index from 1 GeV to 10 GeV
Lum_pwn_mean_100 = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean_100 = numpy.zeros(number_bin_t)                # GeV emission of PSRs

        # Standard deviation
Lum_HESS_std_100 = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std_100 = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Gamma_HESS_std_100 = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_std_100 = numpy.zeros(number_bin_t)               # photon spectral index from 1 GeV to 10 GeV

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_100[:, j]

    Lum_Fermi = Lum_Fermi_100[:, j]

    Gamma_HESS = Gamma_HESS_100[:, j]

    Gamma_GeV = Gamma_GeV_100[:, j]

    Lum_pwn = Lum_pwn_100[:, j]
    Lum_pwn = Lum_pwn[numpy.where(Lum_pwn > 0.0)[0]]

    Lum_psr = Lum_psr_100[:, j]
    Lum_psr = Lum_psr[numpy.where(Lum_psr > 0.0)[0]]

    Lum_HESS_mean_100[j] = numpy.mean(Lum_HESS)
    Lum_Fermi_mean_100[j] = numpy.mean(Lum_Fermi)
    Gamma_HESS_mean_100[j] = numpy.mean(Gamma_HESS)
    Gamma_GeV_mean_100[j] = numpy.mean(Gamma_GeV)
    Lum_pwn_mean_100[j] = 10**(numpy.mean(numpy.log10(Lum_pwn)))
    Lum_psr_mean_100[j] = 10**(numpy.mean(numpy.log10(Lum_psr)))

    Lum_HESS_std_100[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std_100[j] = numpy.std(Lum_Fermi)
    Gamma_HESS_std_100[j] = numpy.std(Gamma_HESS)
    Gamma_GeV_std_100[j] = numpy.std(Gamma_GeV)

Lum_HESS_mean_100 = numpy.nan_to_num(Lum_HESS_mean_100)
Lum_HESS_std_100 = numpy.nan_to_num(Lum_HESS_std_100)
Lum_HESS_pst_100 = Lum_HESS_mean_100 + Lum_HESS_std_100
Lum_HESS_mst_100 = Lum_HESS_mean_100 - Lum_HESS_std_100
Lum_HESS_mst_100 = numpy.nan_to_num(Lum_HESS_mst_100)
ind0 = numpy.where(Lum_HESS_mst_100 < 0)[0]
Lum_HESS_mst_100[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_mean_100 = numpy.nan_to_num(Lum_Fermi_mean_100)
Lum_Fermi_std_100 = numpy.nan_to_num(Lum_Fermi_std_100)
Lum_Fermi_pst_100 = Lum_Fermi_mean_100 + Lum_Fermi_std_100
Lum_Fermi_mst_100 = Lum_Fermi_mean_100 - Lum_Fermi_std_100
Lum_Fermi_mst_100 = numpy.nan_to_num(Lum_Fermi_mst_100)
ind0 = numpy.where(Lum_Fermi_mst_100 < 0)[0]
Lum_Fermi_mst_100[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_pst_100 = Gamma_HESS_mean_100 + Gamma_HESS_std_100
Gamma_HESS_mst_100 = Gamma_HESS_mean_100 - Gamma_HESS_std_100
ind0 = numpy.where(Gamma_HESS_mst_100 < 0)[0]
Gamma_HESS_mst_100[ind0] = numpy.zeros(len(ind0))

Gamma_GeV_pst_100 = Gamma_GeV_mean_100 + Gamma_GeV_std_100
Gamma_GeV_mst_100 = Gamma_GeV_mean_100 - Gamma_GeV_std_100
ind0 = numpy.where(Gamma_GeV_mst_100 < 0)[0]
Gamma_GeV_mst_100[ind0] = numpy.zeros(len(ind0))

        # Third data set for comparison

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/')

with open('Total', 'rb') as iteration_write:

    Lum_HESS_300 = pickle.load(iteration_write)
    Lum_Fermi_300 = pickle.load(iteration_write)
    Lum_300 = pickle.load(iteration_write)
    Gamma_HESS_300 = pickle.load(iteration_write)
    Gamma_GeV_300 = pickle.load(iteration_write)
    Gamma_MeV_300 = pickle.load(iteration_write)
    Lum_pwn_300 = pickle.load(iteration_write)
    Lum_psr_300 = pickle.load(iteration_write)

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

    # Initialization

        # Mean
Lum_HESS_mean_300 = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_Fermi_mean_300 = numpy.zeros(number_bin_t)              # from 100 MeV to 100 GeV
Gamma_HESS_mean_300 = numpy.zeros(number_bin_t)             # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_mean_300 = numpy.zeros(number_bin_t)              # photon spectral index from 1 GeV to 10 GeV
Lum_pwn_mean_300 = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean_300 = numpy.zeros(number_bin_t)                # GeV emission of PSRs

        # Standard deviation
Lum_HESS_std_300 = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std_300 = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Gamma_HESS_std_300 = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_std_300 = numpy.zeros(number_bin_t)               # photon spectral index from 1 GeV to 10 GeV

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_300[:, j]

    Lum_Fermi = Lum_Fermi_300[:, j]

    Gamma_HESS = Gamma_HESS_300[:, j]

    Gamma_GeV = Gamma_GeV_300[:, j]

    Lum_pwn = Lum_pwn_300[:, j]
    Lum_pwn = Lum_pwn[numpy.where(Lum_pwn > 0.0)[0]]

    Lum_psr = Lum_psr_300[:, j]
    Lum_psr = Lum_psr[numpy.where(Lum_psr > 0.0)[0]]

    Lum_HESS_mean_300[j] = numpy.mean(Lum_HESS)
    Lum_Fermi_mean_300[j] = numpy.mean(Lum_Fermi)
    Gamma_HESS_mean_300[j] = numpy.mean(Gamma_HESS)
    Gamma_GeV_mean_300[j] = numpy.mean(Gamma_GeV)
    Lum_pwn_mean_300[j] = 10**(numpy.mean(numpy.log10(Lum_pwn)))
    Lum_psr_mean_300[j] = 10**(numpy.mean(numpy.log10(Lum_psr)))

    Lum_HESS_std_300[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std_300[j] = numpy.std(Lum_Fermi)
    Gamma_HESS_std_300[j] = numpy.std(Gamma_HESS)
    Gamma_GeV_std_300[j] = numpy.std(Gamma_GeV)

Lum_HESS_mean_300 = numpy.nan_to_num(Lum_HESS_mean_300)
Lum_HESS_std_300 = numpy.nan_to_num(Lum_HESS_std_300)
Lum_HESS_pst_300 = Lum_HESS_mean_300 + Lum_HESS_std_300
Lum_HESS_mst_300 = Lum_HESS_mean_300 - Lum_HESS_std_300
Lum_HESS_mst_300 = numpy.nan_to_num(Lum_HESS_mst_300)
ind0 = numpy.where(Lum_HESS_mst_300 < 0)[0]
Lum_HESS_mst_300[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_mean_300 = numpy.nan_to_num(Lum_Fermi_mean_300)
Lum_Fermi_std_300 = numpy.nan_to_num(Lum_Fermi_std_300)
Lum_Fermi_pst_300 = Lum_Fermi_mean_300 + Lum_Fermi_std_300
Lum_Fermi_mst_300 = Lum_Fermi_mean_300 - Lum_Fermi_std_300
Lum_Fermi_mst_300 = numpy.nan_to_num(Lum_Fermi_mst_300)
ind0 = numpy.where(Lum_Fermi_mst_300 < 0)[0]
Lum_Fermi_mst_300[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_pst_300 = Gamma_HESS_mean_300 + Gamma_HESS_std_300
Gamma_HESS_mst_300 = Gamma_HESS_mean_300 - Gamma_HESS_std_300
ind0 = numpy.where(Gamma_HESS_mst_300 < 0)[0]
Gamma_HESS_mst_300[ind0] = numpy.zeros(len(ind0))

Gamma_GeV_pst_300 = Gamma_GeV_mean_300 + Gamma_GeV_std_300
Gamma_GeV_mst_300 = Gamma_GeV_mean_300 - Gamma_GeV_std_300
ind0 = numpy.where(Gamma_GeV_mst_300 < 0)[0]
Gamma_GeV_mst_300[ind0] = numpy.zeros(len(ind0))


    # Plot
label_mean = ['$D_0 = 10^{26}$ [cm$^2$ s$^{-1}$]', '$D_0 = 10^{27}$ [cm$^2$ s$^{-1}$]', '$D_0 = 10^{28}$ [cm$^2$ s$^{-1}$]']
label_std = 'none'

sym_mean = ['', '', '']
linestyle_mean = ['-', '-', '-']

color_mean = ['purple', 'blue', 'green']

xlabel = 'Time [Myr]'
text = ''
Title = 'none'

xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5



        # Gamma luminosity of the superbubble
ymin = 1e29
ymax = 1e36
            # HESS energy range

y_mean = [Lum_HESS_mean_30, Lum_HESS_mean_100, Lum_HESS_mean_300]
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

semilog_plot(figure_number, 3, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_HESS.pdf')

figure_number += 1

            # Fermi energy range
y_mean = [Lum_Fermi_mean_30, Lum_Fermi_mean_100, Lum_Fermi_mean_300]
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

semilog_plot(figure_number, 3, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_Fermi.pdf')

figure_number += 1

        # Spectral index
ymin = 0.0
ymax = 3.5

            # HESS energy range
y_mean = [Gamma_HESS_mean_30, Gamma_HESS_mean_100, Gamma_HESS_mean_300]
ylabel_HESS = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'

plot(figure_number, 3, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Photon_index_HESS.pdf')

figure_number += 1

            # 1 GeV to 10 GeV
y_mean = [Gamma_GeV_mean_30, Gamma_GeV_mean_100, Gamma_GeV_mean_300]
ylabel_HESS = '$\Gamma_{ph}$ (1 GeV - 10 GeV)'

plot(figure_number, 3, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym_mean, linestyle_mean, color_mean, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Photon_index_GeV.pdf')

figure_number += 1

plt.show()

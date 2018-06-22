"""
It plots the gamma-ray luminosities in the HE and VHE energy ranges and the spectral index in both energy ranges for all iterations and makes the statistic analyses of the results.

All the parameters must be given in the Parameters_system.

Make sure that you have already run the code 'Iterations.py' to have the data set.
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

    # You need to change it
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/100/'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/100/'
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/100/'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

    # Number of iterations per files
nit = 100                                                                        #you need to change it for your simulations (depends on the number of iterations per files)

    # Number of Files
nfiles = 1                                                                     #you need to change it for your simulations (depends on the number of files (paralelization you have done))

    # Total number of iterations
nit_tot = nit * nfiles                                                          #you need to change it for your simulations

    # Initialization
figure_number = 1
k = 0       # for the concatenisation of the all iterations

Lum_HESS_it = numpy.zeros((nit_tot, number_bin_t))          # gamma-ray luminosity in the H.E.S.S. energy range
Lum_Fermi_it = numpy.zeros((nit_tot, number_bin_t))         # gamma-ray luminosity in the Fermi energy range
Lum_it = numpy.zeros((nit_tot, number_bin_t))               # gamma-ray luminosity in the whole energy range
Gamma_HESS_it = numpy.zeros((nit_tot, number_bin_t))        # spectral index in the H.E.S.S. energy range
Gamma_GeV_it = numpy.zeros((nit_tot, number_bin_t))         # spectral index in the HE energy range (1 GeV - 10 GeV)
Gamma_MeV_it = numpy.zeros((nit_tot, number_bin_t))         # spectral index in the HE energy range (100 MeV - 100 GeV)
Lum_pwn_it = numpy.zeros((nit_tot, number_bin_t))           # TeV emission of PWNe
Lum_psr_it = numpy.zeros((nit_tot, number_bin_t))           # GeV emission of PSRs
tsn_it = []                                                 # SN explosion times (yr)
nsn_it = numpy.zeros((nit_tot, number_bin_t))               # number of supernova per iterations

    ## ------- ##
    # Load data #
    ## ------- ##
for i in range (nfiles):

    nfile = i + 1
    file = '%d'%nfile

        # You need to change it
    if nfiles > 1:
        os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/' + file)

    else:
        os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/')

    with open('SB', 'rb') as SB_load:

            # SN explosions time
        tsn = pickle.load(SB_load)  # yrs
        tsn = tsn * yr26yr          # Myrs

            # number of SN
        nsn = pickle.load(SB_load)

    with open('General', 'rb') as data_load:

            # Gamma-ray luminosities
        Lum_HESS = pickle.load(data_load)   # erg/s
        Lum_Fermi = pickle.load(data_load)  # erg/s
        Lum = pickle.load(data_load)        # erg/s

            # Spectral index
        Gamma_HESS = pickle.load(data_load)
        Gamma_GeV = pickle.load(data_load)
        Gamma_MeV = pickle.load(data_load)

            # Gamma-ray luminosity of PSRs and PWNe
        Lum_pwn = pickle.load(data_load)    # erg/s
        Lum_psr = pickle.load(data_load)    # erg/s

            # Concatenisation of all iterations
        for j in range (nit):

            Lum_HESS_it[j + k] = Lum_HESS[j]
            Lum_Fermi_it[j + k] = Lum_Fermi[j]
            Lum_it[j + k] = Lum[j]
            Gamma_HESS_it[j + k] = Gamma_HESS[j]
            Gamma_GeV_it[j + k] = Gamma_GeV[j]
            Gamma_MeV_it[j + k] = Gamma_MeV[j]
            Lum_pwn_it[j + k] = Lum_pwn[j]
            Lum_psr_it[j + k] = Lum_psr[j]
            tsn_it.append(tsn[j])
            nsn_it[j + k] = nsn[j]

        k += nit

    # Recording of the concatenisation (you need to change it)
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/100/')

        # For the others
with open('Total', 'wb') as iteration_write:
    pickle.dump(Lum_HESS_it, iteration_write)
    pickle.dump(Lum_Fermi_it, iteration_write)
    pickle.dump(Lum_it, iteration_write)
    pickle.dump(Gamma_HESS_it, iteration_write)
    pickle.dump(Gamma_GeV_it, iteration_write)
    pickle.dump(Gamma_MeV_it, iteration_write)
    pickle.dump(Lum_pwn_it, iteration_write)
    pickle.dump(Lum_psr_it, iteration_write)

        # For 30 Dor C
"""
with open('General', 'wb') as iteration_write:

    pickle.dump(Lum_HESS_it, iteration_write)
    pickle.dump(Lum_Fermi_it, iteration_write)
    pickle.dump(Lum_it, iteration_write)
    pickle.dump(Gamma_HESS_it, iteration_write)
    pickle.dump(Gamma_GeV_it, iteration_write)
    pickle.dump(Gamma_MeV_it, iteration_write)
    pickle.dump(Lum_pwn_it, iteration_write)
    pickle.dump(Lum_psr_it, iteration_write)

with open('SB', 'wb') as iteration_write:

    pickle.dump(tsn_it, iteration_write)
    pickle.dump(nsn_it, iteration_write)
"""
    ##-------------------------------------------##
    # Histogramme of the sn in our time interval  #
    ##-------------------------------------------##
if nit_tot > 1:

        # Computation of the probability to get nsn
    title = ''
    xlabel = '$n_{sn}$'
    ylabel = 'counts'
    figure = figure_number
    label = 'none'
    len_bins = 1
    histmin = numpy.min(nsn_it)
    histmax = numpy.max(nsn_it)

    histogramme(figure, nsn_it, label, title, xlabel, ylabel, len_bins)
    plt.savefig(pathfigure+'Histogramme_nsn.pdf')

    figure_number = figure + 1

    # Computation of the probability to get nsn
title = ''
xlabel = '$t_{sn}$'
ylabel = 'counts'
figure = figure_number
label = 'none'
len_bins = 1    # Myr

for i in range (nfiles):
    histogramme(figure, tsn_it[i], label, title, xlabel, ylabel, len_bins)

plt.savefig(pathfigure+'Histogramme_tsn.pdf')

figure_number = figure + 1

    ## ---------------------------------------------------------------- ##
    # Histogramme of the probability to have one luminosity/photon index #
    ## ---------------------------------------------------------------- ##
indt = [1000, 2500]
if nit_tot > 1:

            # choosen time
    ind_hist = [indt[0], indt[-1]]

            # Computation of the probability to get a luminosity L
    title = ''
    xlabel_LH = '$L_\gamma$ (1 TeV - 10 TeV)'
    xlabel_GH = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'
    xlabel_LF = '$L_\gamma$ (100 MeV - 100 GeV)'
    xlabel_GG = '$\Gamma_{ph}$ (1 GeV - 10 GeV)'
    xlabel_GM = '$\Gamma_{ph}$ (100 MeV - 1 GeV)'
    xlabel_PSR = '$Lum_{\gamma, psr}$ (100 MeV - 10 GeV)'
    xlabel_PWN = '$Lum_{\gamma, psr}$ (1 TeV - 10 TeV)'
    ylabel = 'counts'
    figure_LH = figure_number
    figure_GH = figure_LH + 1
    figure_LF = figure_GH + 1
    figure_L = figure_LF + 1
    figure_GG = figure_L + 1
    figure_GM = figure_GG + 1
    figure_PSR = figure_GM + 1
    figure_PWN = figure_PSR + 1

    for j in (ind_hist):

        label = 't = %.2e yr'%t_fix[j]

        Lum_HESS = Lum_HESS_it[:, j]
        len_bins = 1e32

        histogramme(figure_LH, Lum_HESS, label, title, xlabel_LH, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_HESS.pdf')

        Gamma_HESS = Gamma_HESS_it[:, j]
        len_bins = 0.001

        histogramme(figure_GH, Gamma_HESS, label, title, xlabel_GH, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Gamma_HESS.pdf')

        Lum_Fermi = Lum_Fermi_it[:, j]
        len_bins = 1e33

        histogramme(figure_L, Lum_Fermi, label, title, xlabel_LF, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_Fermi.pdf')

        Lum = Lum_it[:, j]
        len_bins = 1e33

        histogramme(figure_LF, Lum, label, title, xlabel_LF, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_.pdf')

        Gamma_GeV = Gamma_GeV_it[:, j]
        len_bins = 0.001

        histogramme(figure_GG, Gamma_GeV, label, title, xlabel_GG, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Gamma_GeV.pdf')

        Gamma_MeV = Gamma_MeV_it[:, j]
        len_bins = 0.001

        histogramme(figure_GM, Gamma_MeV, label, title, xlabel_GM, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Gamma_MeV.pdf')

        Lum_PSR = Lum_psr_it[:, j]
        len_bins = 1e33

        histogramme(figure_PSR, Lum_PSR, label, title, xlabel_PSR, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_PSR.pdf')

        Lum_PWN = Lum_pwn_it[:, j]
        len_bins = 1e33

        histogramme(figure_PWN, Lum_PWN, label, title, xlabel_PWN, ylabel, len_bins)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_PWN.pdf')

    figure_number = figure_PWN + 1

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

        # Initialization

            # Mean
Lum_HESS_mean = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_Fermi_mean = numpy.zeros(number_bin_t)              # from 100 MeV to 100 GeV
Lum_mean = numpy.zeros(number_bin_t)                    # from 100 MeV to 100 TeV
Gamma_HESS_mean = numpy.zeros(number_bin_t)             # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_mean = numpy.zeros(number_bin_t)              # photon spectral index from 1 GeV to 10 GeV
Gamma_MeV_mean = numpy.zeros(number_bin_t)              # photon spectral index from 100 MeV to 1 GeV
Lum_pwn_mean = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean = numpy.zeros(number_bin_t)                # GeV emission of PSRs

            # Standard deviation
Lum_HESS_std = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Lum_std = numpy.zeros(number_bin_t)                     # from 100 MeV to 100 TeV
Gamma_HESS_std = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_std = numpy.zeros(number_bin_t)               # photon spectral index from 1 GeV to 10 GeV
Gamma_MeV_std = numpy.zeros(number_bin_t)               # photon spectral index from 100 MeV to 1 GeV
Lum_pwn_std = numpy.zeros(number_bin_t)                 # TeV emission of PWNe
Lum_psr_std = numpy.zeros(number_bin_t)                 # GeV emission of PSRs

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_it[:, j]

    Lum_Fermi = Lum_Fermi_it[:, j]

    Lum = Lum_it[:, j]

    Gamma_HESS = Gamma_HESS_it[:, j]

    Gamma_GeV = Gamma_GeV_it[:, j]

    Gamma_MeV = Gamma_MeV_it[:, j]

    Lum_pwn = Lum_pwn_it[:, j]
    Lum_pwn = Lum_pwn[numpy.where(Lum_pwn > 0.0)[0]]

    Lum_psr = Lum_psr_it[:, j]
    Lum_psr = Lum_psr[numpy.where(Lum_psr > 0.0)[0]]

    Lum_HESS_mean[j] = numpy.mean(Lum_HESS)
    Lum_Fermi_mean[j] = numpy.mean(Lum_Fermi)
    Lum_mean[j] = numpy.mean(Lum)
    Gamma_HESS_mean[j] = numpy.mean(Gamma_HESS)
    Gamma_GeV_mean[j] = numpy.mean(Gamma_GeV)
    Gamma_MeV_mean[j] = numpy.mean(Gamma_MeV)
    Lum_pwn_mean[j] = 10**(numpy.mean(numpy.log10(Lum_pwn)))
    Lum_psr_mean[j] = 10**(numpy.mean(numpy.log10(Lum_psr)))

    Lum_HESS_std[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std[j] = numpy.std(Lum_Fermi)
    Lum_std[j] = numpy.std(Lum)
    Gamma_HESS_std[j] = numpy.std(Gamma_HESS)
    Gamma_GeV_std[j] = numpy.std(Gamma_GeV)
    Gamma_MeV_std[j] = numpy.std(Gamma_MeV)

Lum_HESS_mean = numpy.nan_to_num(Lum_HESS_mean)
Lum_HESS_std = numpy.nan_to_num(Lum_HESS_std)
Lum_HESS_pst = Lum_HESS_mean + Lum_HESS_std
Lum_HESS_mst = Lum_HESS_mean - Lum_HESS_std
Lum_HESS_mst = numpy.nan_to_num(Lum_HESS_mst)
ind0 = numpy.where(Lum_HESS_mst < 0)[0]
Lum_HESS_mst[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_mean = numpy.nan_to_num(Lum_Fermi_mean)
Lum_Fermi_std = numpy.nan_to_num(Lum_Fermi_std)
Lum_Fermi_pst = Lum_Fermi_mean + Lum_Fermi_std
Lum_Fermi_mst = Lum_Fermi_mean - Lum_Fermi_std
Lum_Fermi_mst = numpy.nan_to_num(Lum_Fermi_mst)
ind0 = numpy.where(Lum_Fermi_mst < 0)[0]
Lum_Fermi_mst[ind0] = numpy.zeros(len(ind0))

Lum_mean = numpy.nan_to_num(Lum_mean)
Lum_std = numpy.nan_to_num(Lum_std)
Lum_pst = Lum_mean + Lum_std
Lum_mst = Lum_mean - Lum_std
Lum__mst = numpy.nan_to_num(Lum_mst)
ind0 = numpy.where(Lum_mst < 0)[0]
Lum_mst[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_pst = Gamma_HESS_mean + Gamma_HESS_std
Gamma_HESS_mst = Gamma_HESS_mean - Gamma_HESS_std
ind0 = numpy.where(Gamma_HESS_mst < 0)[0]
Gamma_HESS_mst[ind0] = numpy.zeros(len(ind0))

Gamma_GeV_pst = Gamma_GeV_mean + Gamma_GeV_std
Gamma_GeV_mst = Gamma_GeV_mean - Gamma_GeV_std
ind0 = numpy.where(Gamma_GeV_mst < 0)[0]
Gamma_GeV_mst[ind0] = numpy.zeros(len(ind0))

Gamma_MeV_pst = Gamma_MeV_mean + Gamma_MeV_std
Gamma_MeV_mst = Gamma_MeV_mean - Gamma_MeV_std
ind0 = numpy.where(Gamma_MeV_mst < 0)[0]
Gamma_MeV_mst[ind0] = numpy.zeros(len(ind0))

    # Plot
label = 'none'
sym = ''
linestyle = '-'
xlabel = 'Time [Myr]'
text = ''
Title = ''
xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5

        # Gamma luminosity of the superbubble
ymin = 1e29
ymax = 1e36

            # HESS energy range
y_mean = Lum_HESS_mean
label_mean = 'CRs only'
y_pwn = Lum_pwn_mean
label_pwn = 'PWNe'
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'
color_mean = 'cornflowerblue'
color_pwn = 'green'

semilog_plot(figure_number, 1, t6, y_mean, label_mean, Title, xlabel, ylabel_HESS, sym, linestyle, color_mean, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Lum_HESS_pst, Lum_HESS_mst, color = color_mean, alpha = '0.15')
semilog_plot(figure_number, 1, t6, y_pwn, label_pwn, Title, xlabel, ylabel_HESS, sym, linestyle, color_pwn, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_HESS.pdf')

figure_number += 1

            # Fermi energy range
y_mean = Lum_Fermi_mean
label_mean = 'CRs only'
ylabel_Fermi = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'
color_mean = 'orangered'
y_psr = Lum_psr_mean
label_psr = 'PSRs'
color_psr = 'orange'

semilog_plot(figure_number, 1, t6, y_mean, label_mean, Title, xlabel, ylabel_Fermi, sym, linestyle, color_mean, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Lum_Fermi_pst, Lum_Fermi_mst, color = color_mean, alpha = '0.15')
semilog_plot(figure_number, 1, t6, y_psr, label_psr, Title, xlabel, ylabel_Fermi, sym, linestyle, color_psr, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_Fermi.pdf')

figure_number += 1

            # whole energy range
y = Lum_mean
ylabel = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 TeV)'
color = 'cyan'

semilog_plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Lum_pst, Lum_mst, color = color, alpha = '0.15')
plt.savefig(pathfigure_gamma+'Mean_gamma_emission.pdf')

figure_number += 1

        # TeV emission of PWN
y = Lum_pwn_mean
ylabel = '$\log(L_{\gamma, pwn})$ [erg s$^{-1}$] (1 TeV - 10 TeV)'
color = 'green'

semilog_plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_remain+'Mean_luminosity_pwn.pdf')

figure_number += 1

        # GeV emission of PWN
y = Lum_psr_mean
ylabel = '$\log(L_{\gamma, psr})$ [erg s$^{-1}$] (100 MeV - 100 GeV)'
color = 'orange'

semilog_plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_remain+'Mean_luminosity_psr.pdf')

figure_number += 1

        # Spectral index
ymin = 0.0
ymax = 3.5

            # HESS energy range
y = Gamma_HESS_mean
ylabel_HESS = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'
color = 'cornflowerblue'

plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel_HESS, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Gamma_HESS_pst, Gamma_HESS_mst, color = color, alpha = '0.15')
plt.savefig(pathfigure_gamma+'Photon_index_HESS.pdf')

figure_number += 1

            # 1 GeV to 10 GeV
y = Gamma_GeV_mean
ylabel_GeV = '$\Gamma_{ph}$ (1 GeV - 10 GeV)'
color = 'orangered'

plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel_GeV, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Gamma_GeV_pst, Gamma_GeV_mst, color = color, alpha = '0.15')
plt.savefig(pathfigure_gamma+'Photon_index_GeV.pdf')

figure_number += 1

            # 100 MeV to 1 GeV
y = Gamma_MeV_mean
ylabel_MeV = '$\Gamma_{ph}$ (100 MeV - 1 GeV)'
color = 'green'

plot(figure_number, 1, t6, y, label, Title, xlabel, ylabel_MeV, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.fill_between(t6, Gamma_MeV_pst, Gamma_MeV_mst, color = color, alpha = '0.15')
plt.savefig(pathfigure_gamma+'Photon_index_MeV.pdf')

figure_number += 1

plt.show()

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
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/percentage/20/Gamma_emission/'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/percentage/20/Remain/'
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/percentage/20/SB/'
    # Home
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/Parameters/stars/100/Gamma_emission/bis_300'
#pathfigure_remain = '/home/vivi/Documents/Master_2/Superbubbles/figures/Parameters/stars/100/Remain/bis_300'

## ======================================= ##
# Statistic for a high number of iterations #
## ======================================= ##

"""
Plot the graphics for all iterations
"""

    # Number of iterations per files
nit = 10                                                                       #you need to change it for your simulations

    # Number of Files
nfiles = 10                                                                     #you need to change it for your simulations

    # Total number of iterations
nit_tot = nit * nfiles                                                          #you need to change it for your simulations

    # Fix time array (yr)
tmin = 3/yr26yr         # yr
tmax = 10/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1

Lum_HESS_it = numpy.zeros((nit_tot, number_bin_t))              # total gamma luminosity for the energy range of HESS
Lum_Fermi_it = numpy.zeros((nit_tot, number_bin_t))             # total gamma luminosity for the energy range of Fermi
Lum_it = numpy.zeros((nit_tot, number_bin_t))                   # total gamma luminosity for the whole energy range
Gamma_HESS_it = numpy.zeros((nit_tot, number_bin_t))            # spectral index in the energy range of HESS
Gamma_GeV_it = numpy.zeros((nit_tot, number_bin_t))             # spectral index in the energy range 1 GeV to 10 Gev
Gamma_MeV_it = numpy.zeros((nit_tot, number_bin_t))             # spectral index in the energy range 100 MeV to 1 GeV
Lum_pwn_it = numpy.zeros((nit_tot, number_bin_t))               # TeV emission of PWNe
Lum_psr_it = numpy.zeros((nit_tot, number_bin_t))               # GeV emission of PSRs
nob_it = numpy.zeros((nit_tot, number_bin_t))                   # total of remained OB stars
nsn_it = numpy.zeros((nit_tot))                                 # number of supernovae
tsn_it = []                                                     # supernova explosion times
Rsb_t = numpy.zeros(number_bin_t)
Vsb_t = numpy.zeros(number_bin_t)
ns_t = numpy.zeros(number_bin_t)
Ms_t = numpy.zeros(number_bin_t)
k = 0

    ## ------- ##
    # Load data #
    ## ------- ##
for i in range (nfiles):

    nfile = i + 1
    file = '%d'%nfile

        # IRAP
    os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/percentage/20/' + file)

        # HOME
    #os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Parameters/stars/30/'+ file)

    with open('SB', 'rb') as SB_load:

        t0 = pickle.load(SB_load)
        t0 = t0 * yr26yr
        nsn = pickle.load(SB_load)

    with open('General', 'rb') as data_load:

        Lum_HESS = pickle.load(data_load)
        Lum_Fermi = pickle.load(data_load)
        Lum = pickle.load(data_load)
        Gamma_HESS = pickle.load(data_load)
        Gamma_GeV = pickle.load(data_load)
        Gamma_MeV = pickle.load(data_load)
        Lum_pwn = pickle.load(data_load)
        Lum_psr = pickle.load(data_load)
        nob = pickle.load(data_load)
        Rsb = pickle.load(data_load)
        Vsb = pickle.load(data_load)
        Ms = pickle.load(data_load)
        ns = pickle.load(data_load)
        ind = numpy.where(Rsb > 0.0)[0]
        Rsb_t[ind] = Rsb[ind]
        Vsb_t[ind] = Vsb[ind]
        ns_t[ind] = ns[ind]
        Ms_t[ind] = Ms[ind]
        """
            # Fermi energy range
        Emin = 100 * MeV2GeV                # 100 MeV (GeV)
        Emax = 100                          # 100 GeV
        indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]
                # Gamma luminosity
        spectrum_Fermi = spectrum[indE]
        spectrum_erg = spectrum_Fermi * 1.0/erg2GeV     # only in the energy range (erg)
        spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
        lum_Fermi = Flux[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
        Lum_Fermi = luminosity(lum_Fermi, spectrum_ev) # erg s^-1
                # Spectral photon index
        Fluxmin = Flux[:, :, indE[0]]
        Fluxmax = Flux[:, :, indE[-1]]
        Gamma_Fermi = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
        Gamma_Fermi = numpy.nan_to_num(Gamma_Fermi)
        """
        for j in range (nit):

            Lum_HESS_it[j + k] = Lum_HESS[j]
            Lum_Fermi_it[j + k] = Lum_Fermi[j]
            Lum_it[j + k] = Lum[j]
            Gamma_HESS_it[j + k] = Gamma_HESS[j]
            Gamma_GeV_it[j + k] = Gamma_GeV[j]
            Gamma_MeV_it[j + k] = Gamma_MeV[j]
            Lum_pwn_it[j + k] = Lum_pwn[j]
            Lum_psr_it[j + k] = Lum_psr[j]
            nob_it[j + k] = nob[j]
            tsn_it.append(t0[j])
            nsn_it[j + k] = nsn[j]

        k += nit

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/percentage/20/')

with open('Total', 'wb') as iteration_write:

    pickle.dump(Lum_HESS_it, iteration_write)
    pickle.dump(Lum_Fermi_it, iteration_write)
    pickle.dump(Lum_it, iteration_write)
    pickle.dump(Gamma_HESS_it, iteration_write)
    pickle.dump(Gamma_GeV_it, iteration_write)
    pickle.dump(Gamma_MeV_it, iteration_write)
    pickle.dump(Lum_pwn_it, iteration_write)
    pickle.dump(Lum_psr_it, iteration_write)
    pickle.dump(nob_it, iteration_write)
    pickle.dump(tsn_it, iteration_write)
    pickle.dump(nsn_it, iteration_write)

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
#indt = numpy.where((t6 >= 5) & (t6 <= 6))[0]
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
        len_bins = 1e33

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
nob_mean = numpy.zeros(number_bin_t)                    # number of remained OB stars
Lum_pwn_mean = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean = numpy.zeros(number_bin_t)                # GeV emission of PSRs

            # Standard deviation
Lum_HESS_std = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Lum_std = numpy.zeros(number_bin_t)                     # from 100 MeV to 100 TeV
Gamma_HESS_std = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_GeV_std = numpy.zeros(number_bin_t)               # photon spectral index from 1 GeV to 10 GeV
Gamma_MeV_std = numpy.zeros(number_bin_t)               # photon spectral index from 100 MeV to 1 GeV
nob_std = numpy.zeros(number_bin_t)                     # number of remained OB stars
Lum_pwn_std = numpy.zeros(number_bin_t)                 # TeV emission of PWNe
Lum_psr_std = numpy.zeros(number_bin_t)                 # GeV emission of PSRs

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_it[:, j]
    Lum_HESS = Lum_HESS[numpy.where(Lum_HESS > 0.0)[0]]

    Lum_Fermi = Lum_Fermi_it[:, j]
    Lum_Fermi = Lum_Fermi[numpy.where(Lum_Fermi > 0.0)[0]]

    Lum = Lum_it[:, j]
    Lum = Lum[numpy.where(Lum > 0.0)[0]]

    Gamma_HESS = Gamma_HESS_it[:, j]
    #Gamma_HESS = Gamma_HESS[numpy.where(Gamma_HESS > 0.0)[0]]

    Gamma_GeV = Gamma_GeV_it[:, j]
    #Gamma_GeV = Gamma_GeV[numpy.where(Gamma_GeV > 0.0)[0]]

    Gamma_MeV = Gamma_MeV_it[:, j]
    #Gamma_MeV = Gamma_MeV[numpy.where(Gamma_MeV > 0.0)[0]]

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
    nob_mean[j] = numpy.mean(nob_it[:, j])
    Lum_pwn_mean[j] = numpy.mean(numpy.log10(Lum_pwn))
    Lum_psr_mean[j] = numpy.mean(numpy.log10(Lum_psr))

    Lum_HESS_std[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std[j] = numpy.std(Lum_Fermi)
    Lum_std[j] = numpy.std(Lum)
    Gamma_HESS_std[j] = numpy.std(Gamma_HESS)
    Gamma_GeV_std[j] = numpy.std(Gamma_GeV)
    Gamma_MeV_std[j] = numpy.std(Gamma_MeV)
    nob_std[j] = numpy.std(nob_it[:, j])
    Lum_pwn_std[j] = numpy.std(numpy.log10(Lum_pwn))
    Lum_psr_std[j] = numpy.std(numpy.log10(Lum_psr))

Lum_HESS_mean = numpy.nan_to_num(Lum_HESS_mean)
Lum_HESS_std = numpy.nan_to_num(Lum_HESS_std)
Lum_HESS_pst = Lum_HESS_mean + Lum_HESS_std

if Nob == 30:

    Lum_HESS_mst = numpy.zeros(number_bin_t)

else:

    Lum_HESS_mst = Lum_HESS_mean - Lum_HESS_std
    Lum_HESS_mst = numpy.nan_to_num(Lum_HESS_mst)
    ind0 = numpy.where(Lum_HESS_mst < 0)[0]
    Lum_HESS_mst[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_mean = numpy.nan_to_num(Lum_Fermi_mean)
Lum_Fermi_std = numpy.nan_to_num(Lum_Fermi_std)
Lum_Fermi_pst = Lum_Fermi_mean + Lum_Fermi_std

if Nob == 30:

    Lum_Fermi_mst = numpy.zeros(number_bin_t)

else:

    Lum_Fermi_mst = Lum_Fermi_mean - Lum_Fermi_std
    Lum_Fermi_mst = numpy.nan_to_num(Lum_Fermi_mst)
    ind0 = numpy.where(Lum_Fermi_mst < 0)[0]
    Lum_Fermi_mst[ind0] = numpy.zeros(len(ind0))

Lum_mean = numpy.nan_to_num(Lum_mean)
Lum_std = numpy.nan_to_num(Lum_std)
Lum_pst = Lum_mean + Lum_std

if Nob == 30:

    Lum_mst = numpy.zeros(number_bin_t)

else:

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

nob_pst = nob_mean + nob_std
nob_mst = nob_mean - nob_std
indNob = numpy.where(nob_pst >= Nob)[0]
nob_pst[indNob] = Nob * numpy.ones(len(indNob))

Lum_pwn_pst = Lum_pwn_mean + Lum_pwn_std
Lum_pwn_mst = Lum_pwn_mean - Lum_pwn_std
Lum_pwn_mst = numpy.nan_to_num(Lum_pwn_mst)
ind0 = numpy.where(Lum_pwn_mst < 0)[0]
Lum_pwn_mst[ind0] = numpy.zeros(len(ind0))

Lum_psr_pst = Lum_psr_mean + Lum_psr_std
Lum_psr_mst = Lum_psr_mean - Lum_psr_std
Lum_psr_mst = numpy.nan_to_num(Lum_psr_mst)
ind0 = numpy.where(Lum_psr_mst < 0)[0]
Lum_psr_mst[ind0] = numpy.zeros(len(ind0))

    # Plot
label = 'none'
sym = ['', '', '']
linestyle = ['-.', ':', ':']
xlabel = 'Time [Myr]'
text = ''
Title = ''
color = ['cornflowerblue', 'green', 'green']
xmin = tmin * yr26yr - 0.5
xmax = tmax * yr26yr + 0.5

        # Gamma luminosity of the superbubble

            # HESS energy range
y = [Lum_HESS_mean, Lum_HESS_pst, Lum_HESS_mst]
print('Luminosity in the H.E.S.S. range')
print(Lum_HESS_mean[2999])
print(Lum_HESS_std[2999])
ind0 = 1000
ymin = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] - 3)
ymax = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] + 2)
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'
#ind0 = numpy.where(Lum_HESS_mst > 0)[0]

semilog_plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_HESS, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_Mean_gamma_emission_HESS.pdf')

figure_number += 1

            # Fermi energy range
y = [Lum_Fermi_mean, Lum_Fermi_pst, Lum_Fermi_mst]
print('Luminosity in the Fermi range')
print(Lum_Fermi_mean[2999])
print(Lum_Fermi_std[2999])
ind0 = 1000
ymin = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] - 3)
ymax = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] + 2)
ylabel_Fermi = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'
#ind0 = numpy.where(Lum_Fermi_mst > 0)[0]

semilog_plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_Fermi, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_Mean_gamma_emission_Fermi.pdf')

figure_number += 1

            # whole energy range
y = [Lum_mean, Lum_pst, Lum_mst]
ind0 = 1000
ymin = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] - 3)
ymax = 10**(numpy.asarray(numpy.log10(y[0]))[ind0] + 2)
ylabel = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 TeV)'
#ind0 = numpy.where(Lum_mst > 0)[0]

semilog_plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_Mean_gamma_emission.pdf')

figure_number += 1

        # Spectral index

            # HESS energy range
y = [Gamma_HESS_mean, Gamma_HESS_pst, Gamma_HESS_mst]
print('Photon spectral index in the H.E.S.S. range')
print(Gamma_HESS_mean[2999])
print(Gamma_HESS_std[2999])
y = numpy.nan_to_num(y)
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 0.02
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 0.02
ylabel_HESS = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'
#ind0 = numpy.where(Gamma_HESS_mst > 0)[0]

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_HESS, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_stat_Photon_index_HESS.pdf')

figure_number += 1

            # 1 GeV to 10 GeV
y = [Gamma_GeV_mean, Gamma_GeV_pst, Gamma_GeV_mst]
print('Photon spectral index in the GeV range')
print(Gamma_GeV_mean[2999])
print(Gamma_GeV_std[2999])
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 0.02
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 0.02
ylabel_GeV = '$\Gamma_{ph}$ (1 GeV - 10 GeV)'
#ind0 = numpy.where(Gamma_GeV_mst > 0)[0]

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_GeV, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_stat_Photon_index_GeV.pdf')

figure_number += 1

            # 100 MeV to 1 GeV
y = [Gamma_MeV_mean, Gamma_MeV_pst, Gamma_MeV_mst]
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 0.02
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 0.02
ylabel_MeV = '$\Gamma_{ph}$ (100 MeV - 1 GeV)'
#ind0 = numpy.where(Gamma_MeV_mst > 0)[0]

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_MeV, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'new_stat_Photon_index_MeV.pdf')

figure_number += 1

        # TeV emission of PWN
y = [Lum_pwn_mean, Lum_pwn_pst, Lum_pwn_mst]
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 1
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 1
ylabel = '$\log(L_{\gamma, pwn})$ [erg s$^{-1}$] (1 TeV - 10 TeV)'
#ind0 = numpy.where(Lum_pwn_mst > 0)[0]

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_remain+'Mean_luminosity_pwn.pdf')

figure_number += 1

        # GeV emission of PWN
y = [Lum_psr_mean, Lum_psr_pst, Lum_psr_mst]
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 1
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 1
ylabel = '$\log(L_{\gamma, psr})$ [erg s$^{-1}$] (100 MeV - 100 GeV)'
#ind0 = numpy.where(Lum_psr_mst > 0)[0]

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_remain+'Mean_luminosity_psr.pdf')

figure_number += 1

        # Number of remained OB stars in the association
y = [nob_mean, nob_pst, nob_mst]
ind0 = numpy.where(y[0] > 0.0)[0]
ymin = numpy.min(numpy.asarray(y[2])[ind0]) - 1
ymax = numpy.max(numpy.asarray(y[1])[ind0]) + 1
ylabel_ob = '$n_{ob}$'

plot(figure_number, 3, t6, y, label, Title, xlabel, ylabel_ob, sym, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_remain+'Mean_ob.pdf')

figure_number += 1

    # Parameters of the superbubble
ind = numpy.where(Rsb_t > 0.0)[0]
x = t6
color = 'cornflowerblue'
symbol = 'x'
linestyle = ''

    # Radius
y = Rsb_t/100.
ymin = numpy.min(y) - 0.1
ymax = numpy.max(y) + 0.1
ylabel = '$R_{sb}$/(100 [pc])'
color_2 = ['cornflowerblue', 'orange']
symbol_2 = ['x', '+']
linestyle_2 = ['', '']

plot(figure_number, 1, x, y, label, Title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure+'Radius.pdf')
print(Rsb_t[2999])
figure_number += 1

    # Velocity
y = Vsb_t/10.
ymin = numpy.min(y) - 0.1
ymax = numpy.max(y) + 0.1
ylabel = '$V_{sb}$/(10 [km/s])'

plot(figure_number, 1, x, y, label, Title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure+'Velocity.pdf')
print(Vsb_t[2999])
figure_number += 1

    # Mass
y = Ms_t/1000000.
ymin = numpy.min(y) - 0.1
ymax = numpy.max(y) + 0.1
ylabel = '$M_{shell}$/(100000 [M$_\odot$])'

plot(figure_number, 1, x, y, label, Title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure+'Mass.pdf')
print(Ms_t[2999])
figure_number += 1

    # Density
y = ns_t/10.
ymin = numpy.min(y) - 0.01
ymax = numpy.max(y) + 0.01
ylabel = '$n_{shell}$/(1 [cm$^{-3}$])'

plot(figure_number, 1, x, y, label, Title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure+'Density.pdf')
print(ns_t[2999])
figure_number += 1

plt.show()

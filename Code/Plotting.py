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
from Parameters_system_30_Dor_C import *

##====##
# Path #
##====##

## NEED TO WRITE CLEARLY WHAT I DO

    # IRAP
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/30/Gamma_emission/1500'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/Parameters/stars/30/Remain/1500'

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
nit = 150                                                                       #you need to change it for your simulations

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

Lum_HESS_it = numpy.zeros((nit_tot, number_bin_t))            # total gamma luminosity for the energy range of HESS
Lum_Fermi_it = numpy.zeros((nit_tot, number_bin_t))           # total gamma luminosity for the energy range of Fermi
Lum_it = numpy.zeros((nit_tot, number_bin_t))                 # total gamma luminosity for the whole energy range
Gamma_HESS_it = numpy.zeros((nit_tot, number_bin_t))          # spectral index in the energy range of HESS
Gamma_Fermi_it = numpy.zeros((nit_tot, number_bin_t))         # spectral index in the energy range of Fermi
n_pwn_it = numpy.zeros((nit_tot, number_bin_t))               # total number of pulsar wind nebula
Lum_pwn_it = numpy.zeros((nit_tot, number_bin_t))             # TeV emission of PWNe
Lum_psr_it = numpy.zeros((nit_tot, number_bin_t))             # TeV emission of PSRs
nob_it = numpy.zeros((nit_tot, number_bin_t))                 # total of remained OB stars
t0_it = numpy.zeros((nit_tot, number_bin_t))                  # SN explosion times (yr)
k = 0

    ## ------- ##
    # Load data #
    ## ------- ##
for i in range (nfiles):

    nfile = i + 1
    file = '%d'%nfile

        # IRAP
    os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/30/'+ file)

        # HOME
    #os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Parameters/stars/30/'+ file)

    with open('SB', 'rb') as SB_load:

        t0 = pickle.load(SB_load)

    with open('General', 'rb') as data_load:

        Lum_HESS = pickle.load(data_load)
        Lum_Fermi = pickle.load(data_load)
        Lum = pickle.load(data_load)
        #Flux = pickle.load(data_load)
        Gamma_HESS = pickle.load(data_load)
        Gamma_Fermi = pickle.load(data_load)
        n_pwn = pickle.load(data_load)
        Lum_pwn = pickle.load(data_load)
        Lum_psr = pickle.load(data_load)
        nob = pickle.load(data_load)
        spectrum = pickle.load(data_load)
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
            Gamma_Fermi_it[j + k] = Gamma_Fermi[j]
            n_pwn_it[j + k] = n_pwn[j]
            Lum_pwn_it[j + k] = Lum_pwn[j]
            Lum_psr_it[j + k] = Lum_psr[j]
            nob_it[j + k] = nob[j]

        k += nit

    ##-----------------------------------------------------##
    # Histogramme of the probability to have one luminosity #
    ##-----------------------------------------------------##

if nit_tot > 1:
            # choosen time
    ind_hist = [1, 299, 599, 899, 1199, 1499, 1799, 2099, 2399, 2699]
            # Computation of the probability to get a luminosity L
    title = ''
    xlabel_LH = '$L_\gamma$ (1 TeV - 10 TeV)'
    xlabel_GH = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'
    xlabel_LF = '$L_\gamma$ (100 MeV - 100 GeV)'
    xlabel_GF = '$\Gamma_{ph}$ (100 MeV - 100 GeV)'
    ylabel = 'counts'
    figure_LH = figure_number
    figure_GH = figure_LH + 1
    figure_LF = figure_GH + 1
    figure_GF = figure_LF + 1

    for j in (ind_hist):

        label = 't = %.2e yr'%t_fix[j]

        Lum_HESS = Lum_HESS_it[:, j]
        #Lum_HESS = Lum_HESS[numpy.where(Lum_HESS > 0.0)[0]]

        histogramme(figure_LH, Lum_HESS, label, title, xlabel_LH, ylabel)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_HESS.pdf')

        Gamma_HESS = Gamma_HESS_it[:, j]
        #Gamma_HESS = Gamma_HESS[numpy.where(Gamma_HESS > 0.0)[0]]

        histogramme(figure_GH, Gamma_HESS, label, title, xlabel_GH, ylabel)
        plt.savefig(pathfigure_gamma+'Histogramme_Gamma_HESS.pdf')

        Lum_Fermi = Lum_Fermi_it[:, j]
        #Lum_Fermi = Lum_Fermi[numpy.where(Lum_Fermi > 0.0)[0]]

        histogramme(figure_LF, Lum_Fermi, label, title, xlabel_LF, ylabel)
        plt.savefig(pathfigure_gamma+'Histogramme_Lum_Fermi.pdf')

        Gamma_Fermi = Gamma_Fermi_it[:, j]
        #Gamma_Fermi = Gamma_Fermi[numpy.where(Gamma_Fermi > 0.0)[0]]

        histogramme(figure_GF, Gamma_Fermi, label, title, xlabel_GF, ylabel)
        plt.savefig(pathfigure_gamma+'Histogramme_Gamma_Fermi.pdf')

    figure_number = figure_GF + 1

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

        # Initialization
indt = numpy.where(t6 >= 3 + 0.01)[0]
t6 = t6[indt]
number_bin_t = len(t6)

            # Mean
Lum_HESS_mean = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_Fermi_mean = numpy.zeros(number_bin_t)              # from 100 MeV to 100 GeV
Lum_mean = numpy.zeros(number_bin_t)                    # from 100 MeV to 100 TeV
Gamma_HESS_mean = numpy.zeros(number_bin_t)             # photon spectral index from 1 TeV to 10 TeV
Gamma_Fermi_mean = numpy.zeros(number_bin_t)            # photon spectral index from 100 MeV to 100 GeV
n_pwn_mean = numpy.zeros(number_bin_t)                  # number of pulsar wind nebula
nob_mean = numpy.zeros(number_bin_t)                    # number of remained OB stars
Lum_pwn_mean = numpy.zeros(number_bin_t)                # TeV emission of PWNe
Lum_psr_mean = numpy.zeros(number_bin_t)                # GeV emission of PSRs

            # Standard deviation
Lum_HESS_std = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_Fermi_std = numpy.zeros(number_bin_t)               # from 100 MeV to 100 GeV
Lum_std = numpy.zeros(number_bin_t)                     # from 100 MeV to 100 TeV
Gamma_HESS_std = numpy.zeros(number_bin_t)              # photon spectral index from 1 TeV to 10 TeV
Gamma_Fermi_std = numpy.zeros(number_bin_t)             # photon spectral index from 100 MeV to 100 GeV
n_pwn_std = numpy.zeros(number_bin_t)                   # number of pulsar wind nebula
nob_std = numpy.zeros(number_bin_t)                     # number of remained OB stars
Lum_pwn_std = numpy.zeros(number_bin_t)                 # TeV emission of PWNe
Lum_psr_std = numpy.zeros(number_bin_t)                 # GeV emission of PSRs

for j in range (number_bin_t):

    Lum_HESS = Lum_HESS_it[:, j + indt[0]]
    Lum_HESS = Lum_HESS[numpy.where(Lum_HESS > 0.0)[0]]

    Lum_Fermi = Lum_Fermi_it[:, j + indt[0]]
    Lum_Fermi = Lum_Fermi[numpy.where(Lum_Fermi > 0.0)[0]]

    Lum = Lum_it[:, j + indt[0]]
    Lum = Lum[numpy.where(Lum > 0.0)[0]]

    Gamma_HESS = Gamma_HESS_it[:, j + indt[0]]
    Gamma_HESS = Gamma_HESS[numpy.where(Gamma_HESS > 0.0)[0]]

    Gamma_Fermi = Gamma_Fermi_it[:, j + indt[0]]
    Gamma_Fermi = Gamma_Fermi[numpy.where(Gamma_Fermi > 0.0)[0]]

    Lum_pwn = Lum_pwn_it[:, j + indt[0]]
    Lum_pwn = Lum_pwn[numpy.where(Lum_pwn > 0.0)[0]]

    Lum_psr = Lum_psr_it[:, j + indt[0]]
    Lum_psr = Lum_psr[numpy.where(Lum_psr > 0.0)[0]]

    Lum_HESS_mean[j] = numpy.mean(Lum_HESS)
    Lum_Fermi_mean[j] = numpy.mean(Lum_Fermi)
    Lum_mean[j] = numpy.mean(Lum)
    Gamma_HESS_mean[j] = numpy.mean(Gamma_HESS)
    Gamma_Fermi_mean[j] = numpy.mean(Gamma_Fermi)
    n_pwn_mean[j] = numpy.mean(n_pwn_it[:, j])
    nob_mean[j] = numpy.mean(nob_it[:, j])
    Lum_pwn_mean[j] = numpy.mean(Lum_pwn)
    Lum_psr_mean[j] = numpy.mean(Lum_psr)

    Lum_HESS_std[j] = numpy.std(Lum_HESS)
    Lum_Fermi_std[j] = numpy.std(Lum_Fermi)
    Lum_std[j] = numpy.std(Lum)
    Gamma_HESS_std[j] = numpy.std(Gamma_HESS)
    Gamma_Fermi_std[j] = numpy.std(Gamma_Fermi)
    n_pwn_std[j] = numpy.std(n_pwn_it[:, j])
    nob_std[j] = numpy.std(nob_it[:, j])
    Lum_pwn_std[j] = numpy.std(Lum_pwn)
    Lum_psr_std[j] = numpy.std(Lum_psr)

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
Lum_mst = numpy.nan_to_num(Lum_mst)
ind0 = numpy.where(Lum_mst < 0)[0]
Lum_mst[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_mean = numpy.nan_to_num(Gamma_HESS_mean)
Gamma_HESS_std = numpy.nan_to_num(Gamma_HESS_std)
Gamma_HESS_pst = Gamma_HESS_mean + Gamma_HESS_std
Gamma_HESS_mst = Gamma_HESS_mean - Gamma_HESS_std
Gamma_HESS_mst = numpy.nan_to_num(Gamma_HESS_mst)
ind0 = numpy.where(Gamma_HESS_mst < 0)[0]
Gamma_HESS_mst[ind0] = numpy.zeros(len(ind0))

Gamma_Fermi_mean = numpy.nan_to_num(Gamma_Fermi_mean)
Gamma_Fermi_std = numpy.nan_to_num(Gamma_Fermi_std)
Gamma_Fermi_pst = Gamma_Fermi_mean + Gamma_Fermi_std
Gamma_Fermi_mst = Gamma_Fermi_mean - Gamma_Fermi_std
Gamma_Fermi_mst = numpy.nan_to_num(Gamma_Fermi_mst)
ind0 = numpy.where(Gamma_Fermi_mst < 0)[0]
Gamma_Fermi_mst[ind0] = numpy.zeros(len(ind0))

n_pwn_pst = n_pwn_mean + n_pwn_std
n_pwn_mst = n_pwn_mean - n_pwn_std
ind0 = numpy.where(n_pwn_mst < 0)[0]
n_pwn_mst[ind0] = numpy.zeros(len(ind0))

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
#text = r'$D_0$ = %.2e $cm^2 s^{-1}$, $\delta$ = %.2f'u'\n'r'$p_0$ =%.2e $GeV c^{-1}$, $\alpha$ = %.2f'u'\n'r' $n_0$ = %.2e $cm^{-3}$' u'\n' r'$n_{SN}$ = %d, $n_{it}$ = %d'%(D0, delta, p0, alpha, n0, nsn, nit)
color = ['cornflowerblue', 'green', 'green']

        # Gamma luminosity of the superbubble

            # HESS energy range
figure_HESS = figure_number

#Title_HESS = 'Mean gamma emission of the SB in the energy range [1 TeV, 10 TeV]\n'
Title_HESS = ''
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

plot(figure_HESS, 3, t6, [Lum_HESS_mean, Lum_HESS_pst, Lum_HESS_mst], label, Title_HESS, xlabel, ylabel_HESS, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_HESS.pdf')

            # Fermi energy range
figure_Fermi = figure_HESS + 1

#Title_Fermi = 'Mean gamma emission of the SB in the energy range [100 MeV, 100 GeV]\n'
Title_Fermi = ''
ylabel_Fermi = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

plot(figure_Fermi, 3, t6, [Lum_Fermi_mean, Lum_Fermi_pst, Lum_Fermi_mst], label, Title_Fermi, xlabel, ylabel_Fermi, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_Fermi.pdf')

            # whole energy range
figure = figure_Fermi + 1

#Title = 'Mean gamma emission in the energy range [100 MeV, 100 TeV]\n'
Title = ''
ylabel = '$L_\gamma$ [erg s$^{-1}$] (100 MeV - 100 TeV)'

plot(figure, 3, t6, [Lum_mean, Lum_pst, Lum_mst], label, Title, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission.pdf')

figure_number = figure + 1

        # Spectral index

            # HESS energy range
figure_Gamma_HESS = figure_number
#Title_Gamma_HESS = 'Photon index in the energy range [1 TeV, 10 TeV]'
Title_Gamma_HESS = ''
ylabel_HESS = '$\Gamma_{ph}$ (1 TeV - 10 TeV)'

plot(figure_Gamma_HESS, 3, t6, [Gamma_HESS_mean, Gamma_HESS_pst, Gamma_HESS_mst], label, Title_Gamma_HESS, xlabel, ylabel_HESS, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Photon_index_HESS.pdf')

            # Fermi energy range
figure_Gamma_Fermi = figure_Gamma_HESS + 1
#Title_Gamma_Fermi = 'Photon index in the energy range [100 MeV, 100 GeV]'
Title_Gamma_Fermi = ''
ylabel_Fermi = '$\Gamma_{ph}$ (100 MeV - 100 GeV)'

plot(figure_Gamma_Fermi, 3, t6, [Gamma_Fermi_mean, Gamma_Fermi_pst, Gamma_Fermi_mst], label, Title_Gamma_Fermi, xlabel, ylabel_Fermi, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Photon_index_Fermi.pdf')

figure_number = figure_Gamma_Fermi + 1

        # Number of pulsar wind nebula
figure_pwn = figure_number
label = 'none'
#Title_pwn = 'Number of pulsar wind nebula in the superbubble\n'
Title_pwn = ''
ylabel = '$n_{pwn}$'

plot(figure_pwn, 3, t6, [n_pwn_mean, n_pwn_pst, n_pwn_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_pwn.pdf')

figure_number = figure_pwn + 1

        # TeV emission of PWN
figure_pwn = figure_number
#Title_pwn = 'Mean gamma luminosity of the pulsar wind nebulae inside the SB\n'
Title_pwn = ''
ylabel = '$L_{\gamma, pwn}$ [erg s$^{-1}$] (1 TeV - 10 TeV)'

plot(figure_pwn, 3, t6, [Lum_pwn_mean, Lum_pwn_pst, Lum_pwn_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_luminosity_pwn.pdf')

figure_number = figure_pwn + 1

        # GeV emission of PWN
figure_psr = figure_number
#Title_pwn = 'Mean gamma luminosity of the pulsar wind nebulae inside the SB\n'
Title_psr = ''
ylabel = '$L_{\gamma, psr}$ [erg s$^{-1}$] (100 MeV - 100 GeV)'

plot(figure_psr, 3, t6, [Lum_psr_mean, Lum_psr_pst, Lum_psr_mst], label, Title_psr, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_luminosity_psr.pdf')

figure_number = figure_psr + 1

        # Number of remained OB stars in the association
figure_ob = figure_number
label = 'none'
#Title_ob = 'Number of remained Ob stars\n'
Title_ob = ''
ylabel = '$n_{ob}$'
plot(figure_ob, 3, t6, [nob_mean, nob_pst, nob_mst], label, Title_ob, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_ob.pdf')

figure_number = figure_ob + 1

plt.show()

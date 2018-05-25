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
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/10/5')
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/10/Gamma_emission/'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/Parameters/stars/10/Remain/'

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Test/')
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/Test/Gamma_emission'
#pathfigure_remain = '/home/vivi/Documents/Master_2/Superbubbles/figures/Test/Remain/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 100                                                                      #you need to change it for your simulations

    # Which zone for the Computation
zones = [2]                                                                     #you need to change it for your simulations

    # Correction factor

need_correction = False

if need_correction:

    t_end_6 = 4.5           # Myr
    t_end = t_end_6/yr26yr  # yr
    Rsb = 47.0                          # observed radius (pc)                  #you need to change it for your simulations
    Rw = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
    correction_factor = Rsb/Rw

else:
    correction_factor = 1
    t_end_6 = 7

    # Fix time array (yr)
t0min = 3           # Myr
t0max = 7t_end_6    # Myr
tmin = t0min/yr26yr # yr
tmax = (t0max+1)/yr26yr         # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1
time_sn = (t0max - t0min) * 1.0/Nob * 1.0/yr26yr # yr
mean_sn = int(tmax/time_sn)                 # mean value of SN explosions already happen        #you need to change it for your simulations

Lum_it = []                 # total gamma luminosity for the whole energy range
Flux_it = []                # photon flux for the whole energy range
n_pwn_it = []               # total number of pulsar wind nebula
Lum_pwn_it = []             # TeV emission of PWNe
Lum_psr_it = []             # GeV emission of PSRs
nob_it = []                 # total of remained OB stars
t0_it = []                  # SN explosion times (yr)

    ##----------##
    # Iterations #
    ##----------##
print('For %d SN'%Nob)
print('For %d iterations' %nit)

        # SN explosions time (yr)

with open('SB', 'wb') as SB_write:

    for i in range (nit):

        #nsn = numpy.random.poisson(lam = mean_sn)       # random number of SN which mean value equals to mean_sn
        nsn = mean_sn
        t0 = random_uniform(t0min, t0max, nsn)/yr26yr   # random SN explosion times with an uniform distribution from t0min to t0max
        t0 = sorted(t0)

        t0_it.append(t0)

    t0_it = numpy.asarray(t0_it)

    pickle.dump(t0_it, SB_write)

with open('SB', 'rb') as SB_load:

    t0_it = pickle.load(SB_load)

with open('General', 'wb') as data_write:

    for i in range (nit):

        t0 = t0_it[i]

        Lum, Flux, spectrum, n_pwn, Lum_pwn, Lum_psr, nob = data(correction_factor, t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, zones)

        Lum_it.append(Lum)
        Flux_it.append(Flux)
        n_pwn_it.append(n_pwn)
        Lum_pwn_it.append(Lum_pwn)
        Lum_psr_it.append(Lum_psr)
        nob_it.append(nob)

        print('end of the iteration %d' %i)

    Lum_it = numpy.asarray(Lum_it)
    Flux_it = numpy.asarray(Flux_it)

        # in the H.E.S.S. energy range
    Emin = 1 * TeV2GeV                  # 1 TeV (GeV)
    Emax = 10 * TeV2GeV                 # 10 TeV (GeV)
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

            # Gamma luminosity
    spectrum_HESS = spectrum[indE]
    spectrum_erg = spectrum_HESS * 1.0/erg2GeV     # only in the energy range (erg)
    spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
    lum_HESS = Flux_it[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
    Lum_HESS_it = luminosity(lum_HESS, spectrum_ev) # erg s^-1

            # Spectral photon index
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Gamma_HESS_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_HESS_it = numpy.nan_to_num(Gamma_HESS_it)
    #Gamma_HESS_it = numpy.asarray(Gamma_HESS_it)

        # in the Fermi energy range
    Emin = 100 * MeV2GeV                # 100 MeV (GeV)
    Emax = 100                          # 100 GeV
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

            # Gamma luminosity
    spectrum_Fermi = spectrum[indE]
    spectrum_erg = spectrum_Fermi * 1.0/erg2GeV     # only in the energy range (erg)
    spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
    lum_Fermi = Flux_it[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
    Lum_Fermi_it = luminosity(lum_Fermi, spectrum_ev) # erg s^-1

            # Spectral photon index
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Gamma_Fermi_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_Fermi_it = numpy.nan_to_num(Gamma_Fermi_it)
    #Gamma_Fermi_it = numpy.asarray(Gamma_Fermi_it)

    n_pwn_it = numpy.asarray(n_pwn_it)
    Lum_pwn_it = numpy.asarray(Lum_pwn_it)
    Lum_psr_it = numpy.asarray(Lum_psr_it)
    nob_it = numpy.asarray(nob_it)

    pickle.dump(Lum_HESS_it, data_write)
    pickle.dump(Lum_Fermi_it, data_write)
    pickle.dump(Lum_it, data_write)
    pickle.dump(Flux_it, data_write)
    pickle.dump(Gamma_HESS_it, data_write)
    pickle.dump(Gamma_Fermi_it, data_write)
    pickle.dump(n_pwn_it, data_write)
    pickle.dump(Lum_pwn_it, data_write)
    pickle.dump(Lum_psr_it, data_write)
    pickle.dump(nob_it, data_write)
    pickle.dump(spectrum, data_write)

    ##---------##
    # Load data #
    ##---------##

with open('General', 'rb') as data_load:
    Lum_HESS_it = pickle.load(data_load)
    Lum_Fermi_it = pickle.load(data_load)
    Lum_it = pickle.load(data_load)
    Flux_it = pickle.load(data_load)
    Gamma_HESS_it = pickle.load(data_load)
    Gamma_Fermi_it = pickle.load(data_load)
    n_pwn_it = pickle.load(data_load)
    Lum_pwn_it = pickle.load(data_load)
    Lum_psr_it = pickle.load(data_load)
    nob_it = pickle.load(data_load)
    spectrum = pickle.load(data_load)

    ##-----------------------------------------------------##
    # Histogramme of the probability to have one luminosity #
    ##-----------------------------------------------------##

if nit > 1:
            # choosen time
    ind_hist = [0, 299, 599, 899, 1199, 1499, 1799, 2099, 2399, 2699]
            # Computation of the probability to get a luminosity L
    title = 'Probability of the luminosity in the range of energy [1 TeV, 10 TeV]'
    xlabel = '$L_\gamma$'
    ylabel = 'counts'
    for j in (ind_hist):
        L_hist = Lum_HESS_it[:,j]
        label = 't = %.2e yr'%t_fix[j]
        histogramme(figure_number, L_hist, label, title, xlabel, ylabel)
        plt.legend(loc = 'best')
    plt.savefig(pathfigure_gamma+'Histogramme_L_gamma.pdf')
    figure_number += 1

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

        # Initialization
#number_bin_E = len(spectrum_HESS)
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

    Lum_HESS_mean[j] = numpy.mean(Lum_HESS_it[:, j])
    Lum_Fermi_mean[j] = numpy.mean(Lum_Fermi_it[:, j])
    Lum_mean[j] = numpy.mean(Lum_it[:, j])
    Gamma_HESS_mean[j] = numpy.mean(Gamma_HESS_it[:, j])
    Gamma_Fermi_mean[j] = numpy.mean(Gamma_Fermi_it[:, j])
    n_pwn_mean[j] = numpy.mean(n_pwn_it[:, j])
    nob_mean[j] = numpy.mean(nob_it[:, j])
    Lum_pwn_mean[j] = numpy.mean(Lum_pwn_it[:, j])
    Lum_psr_mean[j] = numpy.mean(Lum_psr_it[:, j])

    Lum_HESS_std[j] = numpy.std(Lum_HESS_it[:, j])
    Lum_Fermi_std[j] = numpy.std(Lum_Fermi_it[:, j])
    Lum_std[j] = numpy.std(Lum_it[:, j])
    Gamma_HESS_std[j] = numpy.std(Gamma_HESS_it[:, j])
    Gamma_Fermi_std[j] = numpy.std(Gamma_Fermi_it[:, j])
    n_pwn_std[j] = numpy.std(n_pwn_it[:, j])
    nob_std[j] = numpy.std(nob_it[:, j])
    Lum_pwn_std[j] = numpy.std(Lum_pwn_it[:, j])
    Lum_psr_std[j] = numpy.std(Lum_pwn_it[:, j])

Lum_HESS_pst = Lum_HESS_mean + Lum_HESS_std
Lum_HESS_mst = Lum_HESS_mean - Lum_HESS_std
ind0 = numpy.where(Lum_HESS_mst <= 0)[0]
Lum_HESS_mst[ind0] = numpy.zeros(len(ind0))

Lum_Fermi_pst = Lum_Fermi_mean + Lum_Fermi_std
Lum_Fermi_mst = Lum_Fermi_mean - Lum_Fermi_std
ind0 = numpy.where(Lum_Fermi_mst <= 0)[0]
Lum_Fermi_mst[ind0] = numpy.zeros(len(ind0))

Lum_pst = Lum_mean + Lum_std
Lum_mst = Lum_mean - Lum_std
ind0 = numpy.where(Lum_mst <= 0)[0]
Lum_mst[ind0] = numpy.zeros(len(ind0))

Gamma_HESS_pst = Gamma_HESS_mean + Gamma_HESS_std
Gamma_HESS_mst = Gamma_HESS_mean - Gamma_HESS_std
ind0 = numpy.where(Gamma_HESS_mst <= 0)[0]
Gamma_HESS_mst[ind0] = numpy.zeros(len(ind0))

Gamma_Fermi_pst = Gamma_Fermi_mean + Gamma_Fermi_std
Gamma_Fermi_mst = Gamma_Fermi_mean - Gamma_Fermi_std
ind0 = numpy.where(Gamma_Fermi_mst <= 0)[0]
Gamma_Fermi_mst[ind0] = numpy.zeros(len(ind0))

n_pwn_pst = n_pwn_mean + n_pwn_std
n_pwn_mst = n_pwn_mean - n_pwn_std
ind0 = numpy.where(n_pwn_mst <= 0)[0]
n_pwn_mst[ind0] = numpy.zeros(len(ind0))

nob_pst = nob_mean + nob_std
nob_mst = nob_mean - nob_std
ind0 = numpy.where(nob_mst <= 0)[0]
nob_mst[ind0] = numpy.zeros(len(ind0))
indNob = numpy.where(nob_pst >= Nob)[0]
nob_pst[indNob] = Nob * numpy.ones(len(indNob))

Lum_pwn_pst = Lum_pwn_mean + Lum_pwn_std
Lum_pwn_mst = Lum_pwn_mean - Lum_pwn_std
ind0 = numpy.where(Lum_pwn_mst <= 0)[0]
Lum_pwn_mst[ind0] = numpy.zeros(len(ind0))

Lum_psr_pst = Lum_psr_mean + Lum_psr_std
Lum_psr_mst = Lum_psr_mean - Lum_psr_std
ind0 = numpy.where(Lum_pwn_mst <= 0)[0]
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

    ##===============================================================##
    # Correlation between the gamma luminosity and the spectral index #
    ##===============================================================##

        # choosen time (yr)
indt = numpy.where(t6 >= 4.5)[0]
indt = indt[0]

        # choosen bin of spectral index
Gamma_ph = numpy.zeros(nit)

for i in range (nit):
    Gamma_ph[i] = Gamma_it[i, indt]

        # Initialization
Gamma_min = min(Gamma_ph)
Gamma_max = max(Gamma_ph)
number_bin_Gamma = 5
delta_gamma = (Gamma_max - Gamma_min)/number_bin_Gamma
Gamma_inf = Gamma_min
Lum_HESS_Gamma = numpy.zeros((number_bin_Gamma, number_bin_t))
label = []
sym = []
linestyle = []

for k in range (number_bin_Gamma):

    Gamma_sup = Gamma_inf + delta_gamma
    indGamma = numpy.where((Gamma_ph >= Gamma_inf) & (Gamma_sup <= Gamma_sup))[0]
    label.append('$\Gamma_{ph} \in$ [%.1f, %.1f]'%(Gamma_inf, Gamma_sup))
    linestyle.append('')

    for j in range (number_bin_t):

        Lum_HESS_Gamma[k, j] = numpy.mean(Lum_HESS_it[indGamma, j])

    Gamma_inf = Gamma_sup

figure_correlation = figure_number
Title_correlation = ''
ylabel = '$L_{\gamma}$ [erg s$^{-1}$]'
color = ['brown', 'orange', 'forestgreen', 'cornflowerblue', 'purple']
sym = ['+', 'x', '*', '^', '.']
text = ''

plot(figure_correlation, number_bin_Gamma, t6, Lum_HESS_Gamma, label, Title_correlation, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Correlation.pdf')

plt.show()

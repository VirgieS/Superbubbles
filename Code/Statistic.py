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
os.chdir('/Users/stage/Dropbox/Superbubbles/Files/30_Dor_C/1000_iterations/')#10_SN/')#
pathfigure_gamma = '/Users/stage/Dropbox/Superbubbles/figures/30_Dor_C/Bons/1000_iterations/Poisson/Gamma_emission/'#10_SN/Gamma_emission/'#
pathfigure_remain = '/Users/stage/Dropbox/Superbubbles/figures//30_Dor_C/Bons/1000_iterations/Poisson/Remain/'#10_SN/Remain/'#

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Test/')
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/Test/Gamma_emission'
#pathfigure_remain = '/home/vivi/Documents/Master_2/Superbubbles/figures/Test/Remain/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 1000                                                                      #you need to change it for your simulations

    # Which zone for the Computation
zones = [2]                                                                     #you need to change it for your simulations

    # Correction factor
t_end_6 = 4.5           # Myr
t_end = t_end_6/yr26yr  # yr
Rsb = 47.0                          # observed radius (pc)              #you need to change it for your simulations
Rw = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
correction_factor = Rsb/Rw

    # Fix time array (yr)
t0min = 3              # Myr
t0max = t_end_6 + 1    # Myr
tmin = t0min/yr26yr         # yr
tmax = (t0max + 1)/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1
mean_sn = 5                 # mean value of SN explosions already happen        #you need to change it for your simulations

Lum_it = []                 # total gamma luminosity for the whole energy range
Flux_it = []                # photon flux for the whole energy range
n_pwn_it = []               # total number of pulsar wind nebula
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

        nsn = numpy.random.poisson(lam = mean_sn)       # random number of SN which mean value equals to mean_sn
        t0 = random_uniform(t0min, t0max, nsn)/yr26yr   # random SN explosion times with an uniform distribution from t0min to t0max
        t0 = sorted(t0)

        t0_it.append(t0)

    t0_it = numpy.asarray(t0_it)

    pickle.dump(t0_it, SB_write)
"""
with open('SB', 'rb') as SB_load:

    t0_it = pickle.load(SB_load)
"""
with open('General', 'wb') as data_write:

    for i in range (nit):

        t0 = t0_it[i]

        Lum, Flux, spectrum, n_pwn, nob, figure_number = data(correction_factor, t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, Esep, zones, pathfigure_gamma, i, figure_number)

        Lum_it.append(Lum)
        Flux_it.append(Flux)
        n_pwn_it.append(n_pwn)
        nob_it.append(nob)

        print('end of the iteration %d' %i)

    Lum_it = numpy.asarray(Lum_it)
    Flux_it = numpy.asarray(Flux_it)

        # Gamma luminosity in a precise energy range
    indE = numpy.where((spectrum >= Esep[0]) & (spectrum <= Esep[-1]))[0]
    E1min = Esep[0]
    E1max = Esep[-1]
    spectrum_HESS = spectrum[indE]
    spectrum_erg = spectrum_HESS * 1.0/erg2GeV     # only in the energy range (erg)
    print(spectrum_HESS)
    spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
    lum_HESS = Flux_it[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
    Lum_HESS_it = luminosity(lum_HESS, spectrum_ev) # erg s^-1

        # Spectral photon index
    Emin = Esep[0]
    Emax = Esep[-1]
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Gamma_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_it = numpy.nan_to_num(Gamma_it)
    Gamma_it = numpy.asarray(Gamma_it)

    n_pwn_it = numpy.asarray(n_pwn_it)
    nob_it = numpy.asarray(nob_it)

    pickle.dump(Lum_HESS_it, data_write)
    pickle.dump(Lum_it, data_write)
    pickle.dump(Flux_it, data_write)
    pickle.dump(Gamma_it, data_write)
    pickle.dump(n_pwn_it, data_write)
    pickle.dump(nob_it, data_write)
    pickle.dump(spectrum, data_write)
"""
    ##---------##
    # Load data #
    ##---------##

with open('General', 'rb') as data_load:
    E1min = Esep[0]
    E1max = Esep[-1]
    Lum_HESS_it = pickle.load(data_load)
    Lum_it = pickle.load(data_load)
    Flux_it = pickle.load(data_load)
    Gamma_it = pickle.load(data_load)
    n_pwn_it = pickle.load(data_load)
    nob_it = pickle.load(data_load)
    spectrum = pickle.load(data_load)
    indE = numpy.where((spectrum >= Esep[0]) & (spectrum <= Esep[-1]))
    spectrum_HESS = spectrum[indE]
"""
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
number_bin_E = len(spectrum_HESS)
            # Mean
Lum_HESS_mean = numpy.zeros(number_bin_t)               # from 1 TeV to 10 TeV
Lum_mean = numpy.zeros(number_bin_t)                    # from 100 MeV to 100 TeV
Flux_mean = numpy.zeros((number_bin_t, number_bin_E))   # photon flux from 1 TeV to 10 TeV
Gamma_mean = numpy.zeros(number_bin_t)                  # photon spectral index
n_pwn_mean = numpy.zeros(number_bin_t)                  # number of pulsar wind nebula
nob_mean = numpy.zeros(number_bin_t)                    # number of remained OB stars

            # Standard deviation
Lum_HESS_std = numpy.zeros(number_bin_t)                # from 1 TeV to 10 TeV
Lum_std = numpy.zeros(number_bin_t)                     # from 100 MeV to 100 TeV
Flux_std = numpy.zeros((number_bin_t, number_bin_E))    # photon flux from 1 TeV to 10 TeV
Gamma_std = numpy.zeros(number_bin_t)                   # photon spectral index
n_pwn_std = numpy.zeros(number_bin_t)                   # number of pulsar wind nebula
nob_std = numpy.zeros(number_bin_t)                     # number of remained OB stars

for j in range (number_bin_t):

    Lum_HESS_mean[j] = numpy.mean(Lum_HESS_it[:, j])
    Lum_mean[j] = numpy.mean(Lum_it[:, j])
    Gamma_mean[j] = numpy.mean(Gamma_it[:, j])
    n_pwn_mean[j] = numpy.mean(n_pwn_it[:, j])
    nob_mean[j] = numpy.mean(nob_it[:, j])

    Lum_HESS_std[j] = numpy.std(Lum_HESS_it[:, j])
    Lum_std[j] = numpy.std(Lum_it[:, j])
    Gamma_std[j] = numpy.std(Gamma_it[:, j])
    n_pwn_std[j] = numpy.std(n_pwn_it[:, j])
    nob_std[j] = numpy.std(nob_it[:, j])

    for k in range (number_bin_E):

        Flux_mean[j, k] = numpy.mean(Flux_it[:, j, k])
        Flux_std[j, k] = numpy.std(Flux_it[:, j, k])

Lum_HESS_pst = Lum_HESS_mean + Lum_HESS_std
Lum_HESS_mst = Lum_HESS_mean - Lum_HESS_std
ind0 = numpy.where(Lum_HESS_mst <= 0)[0]
Lum_HESS_mst[ind0] = numpy.zeros(len(ind0))

Lum_pst = Lum_mean + Lum_std
Lum_mst = Lum_mean - Lum_std
ind0 = numpy.where(Lum_mst <= 0)[0]
Lum_mst[ind0] = numpy.zeros(len(ind0))

#Flux_pst = Flux_mean + Flux_std
#Flux_mst = Flux_mean - Flux_std

Gamma_pst = Gamma_mean + Gamma_std
Gamma_mst = Gamma_mean - Gamma_std
ind0 = numpy.where(Gamma_mst <= 0)[0]
Gamma_mst[ind0] = numpy.zeros(len(ind0))

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

    # Plot
label = 'none'
sym = ['', '', '']
linestyle = ['-.', ':', ':']
xlabel = 'Time [Myr]'
text = ''
#text = r'$D_0$ = %.2e $cm^2 s^{-1}$, $\delta$ = %.2f'u'\n'r'$p_0$ =%.2e $GeV c^{-1}$, $\alpha$ = %.2f'u'\n'r' $n_0$ = %.2e $cm^{-3}$' u'\n' r'$n_{SN}$ = %d, $n_{it}$ = %d'%(D0, delta, p0, alpha, n0, nsn, nit)
color = ['cornflowerblue', 'green', 'green']

        # Gamma luminosity of the superbubble
figure_HESS = figure_number
figure = figure_HESS + 1

#Title_HESS = 'Mean gamma emission of the SB in the energy range [1 TeV, 10 TeV]\n'
#Title = 'Mean gamma emission in the energy range [100 MeV, 100 TeV]\n'
Title_HESS = ''
Title = ''
ylabel_HESS = '$L_\gamma$ [erg s$^{-1}$] (%.2f GeV - %.2f GeV)'%(E1min, E1max)
ylabel = '$L_\gamma$ [erg s$^{-1}$] (%.2f GeV - %.2f GeV)'%(spectrum[0], spectrum[-1])

plot(figure_HESS, 3, t6, [Lum_HESS_mean, Lum_HESS_pst, Lum_HESS_mst], label, Title_HESS, xlabel, ylabel_HESS, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_range1.pdf')

plot(figure, 3, t6, [Lum_mean, Lum_pst, Lum_mst], label, Title, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_all.pdf')

figure_number = figure + 1

        # spectral index
figure_Gamma = figure_number
#Title_Gamma = 'Photon index in the energy range [1 TeV, 10 TeV]'
Title_Gamma = ''
ylabel = '$\Gamma_{ph}$'

plot(figure_Gamma, 3, t6, [Gamma_mean, Gamma_pst, Gamma_mst], label, Title_Gamma, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_gamma+'Photon_index.pdf')
figure_number = figure_Gamma + 1

        # Number of pulsar wind nebula
figure_pwn = figure_number
label = 'none'
#Title_pwn = 'Number of pulsar wind nebula in the superbubble\n'
Title_pwn = ''
ylabel = '$n_{pwn}$'

plot(figure_pwn, 3, t6, [n_pwn_mean, n_pwn_pst, n_pwn_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_pwn.pdf')

figure_number = figure_pwn + 1

        # gamma luminosity of the pulsar wind nebula
L_pwn = 1e35    # erg s^-1
L_pwn_mean = n_pwn_mean * L_pwn
L_pwn_pst = n_pwn_pst * L_pwn
L_pwn_mst = n_pwn_mst * L_pwn

figure_pwn = figure_number
#Title_pwn = 'Mean gamma luminosity of the pulsar wind nebulae inside the SB\n'
Title_pwn = ''
ylabel = '$L_{\gamma, pwn}$ [erg s$^{-1}$]'

plot(figure_pwn, 3, t6, [L_pwn_mean, L_pwn_pst, L_pwn_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, color, text)
plt.savefig(pathfigure_remain+'Mean_luminosity_pwn.pdf')

figure_number = figure_pwn + 1

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

"""
        # Flux
figure_flux = figure_number
Title_flux = 'Intrinsic luminosity for all the energy range of HESS\n'
xlabel = 'E [TeV]'
ylabel = '$\Phi_{ph}$ [ph s$^{-1}$ eV$^{-1}$]'
label = ['simulations', 'fit']
sym = ['+', '']
linestyle = ['', '-.']
color = ['cornflowerblue', 'orange']
text = ''

indt = numpy.where(Flux_it[0, :, -1] != 0)[0]

spectrum_HESS_ev = spectrum_HESS * GeV2eV   # in eV
spectrum_HESS_TeV = spectrum_HESS_ev/TeV2eV # in TeV

b_it = Flux_it[0, indt[0], -1] * (spectrum_HESS_ev[-1])**(Gamma_it[0, indt[0]])
fit_it = b_it * spectrum_HESS_ev**(-Gamma_it[0, indt[0]])

y = [Flux_it[0, indt[0]], fit_it]

log_plot(figure_flux, 2, spectrum_HESS_TeV, y, label, Title_flux, xlabel, ylabel, sym, linestyle, color, text)

plt.savefig(pathfigure_gamma+'Photon_flux_all.pdf')
figure_number = figure_flux + 1

#print('for %d iteration(s)'%nit)
#print('for a random number of sn with mean value = %d' %mean)
"""
plt.show()

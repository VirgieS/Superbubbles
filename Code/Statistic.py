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

    # IRAP
#os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/30_Dor_C/')
#pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Gamma_emission/Test/'
#pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Remain/Test/'
#pathfigure_CR = '/Users/stage/Documents/Virginie/Superbubbles/figures/stat_SN/CR/1/'

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/30_Dor_C/Test/')
pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/30_Dor_C/Verif/Gamma_emission/'
pathfigure_remain = '/home/vivi/Documents/Master_2/Superbubbles/figures/30_Dor_C/Verif/Remain/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 1                                                                         #you need to change for your simulations

    # Correction factor
t_end = 4.5e6               # time at which R = Robs                            #you need to change for your simulations
t_end_6 = t_end * yr26yr    # in Myr
Robs = 47.0                 # observed radius (pc)                              #you need to change for your simulations
Rsb = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
correction_factor = Robs/Rsb

    # Which zone for the Computation
zones = [2]                                                                     #you need to change for your simulations

    # Fix time array (yr)
t0min = 3                   # Myr                                               #you need to change for your simulations
t0max = t_end_6             # Myr
tmin = t0min/yr26yr         # yr
tmax = (t0max + 1)/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1
Lum_it_HESS = []    # total gamma luminosity from 1 TeV to 10 TeV
Lum_it = []         # total gamma luminosity from 100 MeV to 100 TeV
Gamma_it = []       # spectral index from 1 TeV to 10 TeV
n_pwn_it = []       # total number of pulsar wind nebula
nob_it = []         # total of remained OB stars
t0_it = []          # SN explosion times (yr)

    ##----------##
    # Iterations #
    ##----------##

with open('SN', 'wb') as SN_write:

    for i in range (nit):

            # SN explosion time
        n = 5             # number of massive stars in the OB association       #you need to change for your simulations
        t0 = random_SN(t0min, t0max, n)/yr26yr
        t0 = sorted(t0)

        t0_it.append(t0)

    t0_it = numpy.asarray(t0_it)

    pickle.dump(t0_it, SN_write)

with open('SN', 'rb') as SN_load:

    t0_it = pickle.load(SN_load)
    n = len(t0_it[0])
    t0 = t0_it[0,0]

for i in range (nit):

    t0 = t0_it[i]

    Lum_HESS, Lum, Gamma, Lum_units, n_pwn, nob, figure_number = data(correction_factor, t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, Esep, p0, alpha, D0, delta, zones, pathfigure_gamma, i, figure_number)
    indHESS = numpy.where(Lum_HESS != 0)
    indGamma = numpy.where(Gamma != 0)
    Lum_it_HESS.append(Lum_HESS)
    Lum_it.append(Lum)
    Gamma_it.append(Gamma)
    n_pwn_it.append(n_pwn)
    nob_it.append(nob)

    print('end of the iteration %d' %i)

Lum_it_HESS = numpy.asarray(Lum_it_HESS)
ind = numpy.where(Lum_HESS > 0)[0]
Lum_it = numpy.asarray(Lum_it)
Gamma_it = numpy.asarray(Gamma_it)
n_pwn_it = numpy.asarray(n_pwn_it)
nob_it = numpy.asarray(nob_it)
"""
    ##-----------------------------------------------------##
    # Histogramme of the probability to have one luminosity #
    ##-----------------------------------------------------##
        # choosen time
ind_hist = [0, 299, 599, 899, 1199, 1499, 1799, 2099, 2399, 2699]

        # Computation of the probability to get a luminosity L
title = 'Probability of the luminosity in the range of energy [1 TeV, 10 TeV]'
xlabel = '$L_\gamma$'
ylabel = 'counts'

for j in (ind_hist):
    L_hist = Lum_it_HESS[:,j]
    label = 't = %.2e yr'%t_fix[j]
    histogramme(figure_number, L_hist, label, title, xlabel, ylabel)
    plt.legend(loc = 'best')

plt.savefig(pathfigure_gamma+'Histogramme_L_gamma.pdf')
figure_number += 1
"""
    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

        # Initialization
            # Mean
Lum_HESS_mean = numpy.zeros(number_bin_t)   # from 1 TeV to 10 TeV
Lum_mean = numpy.zeros(number_bin_t)        # from 100 MeV to 100 TeV
Gamma_mean = numpy.zeros(number_bin_t)      # spectral index
n_pwn_mean = numpy.zeros(number_bin_t)      # number of pulsar wind nebula
nob_mean = numpy.zeros(number_bin_t)        # number of remained OB stars

            # Standard deviation
Lum_HESS_std = numpy.zeros(number_bin_t)    # from 1 TeV to 10 TeV
Lum_std = numpy.zeros(number_bin_t)         # from 100 MeV to 100 TeV
Gamma_std = numpy.zeros(number_bin_t)       # spectral index
n_pwn_std = numpy.zeros(number_bin_t)       # number of pulsar wind nebula
nob_std = numpy.zeros(number_bin_t)         # number of remained OB stars

for j in range (number_bin_t):

    Lum_HESS_mean[j] = numpy.mean(Lum_it_HESS[:,j])
    Lum_mean[j] = numpy.mean(Lum_it[:,j])
    Gamma_mean[j] = numpy.mean(Gamma_it[:,j])
    n_pwn_mean[j] = numpy.mean(n_pwn_it[:,j])
    nob_mean[j] = numpy.mean(nob_it[:, j])

    Lum_HESS_std[j] = numpy.std(Lum_it_HESS[:,j])
    Lum_std[j] = numpy.std(Lum_it[:,j])
    Gamma_std[j] = numpy.std(Gamma_it[:,j])
    n_pwn_std[j] = numpy.std(n_pwn_it[:,j])
    nob_std[j] = numpy.std(nob_it[:, j])

Lum_HESS_pst = Lum_HESS_mean + Lum_HESS_std
Lum_HESS_mst = Lum_HESS_mean - Lum_HESS_std
ind = numpy.where(Lum_HESS_mean > 0)[0]

Lum_pst = Lum_mean + Lum_std
Lum_mst = Lum_mean - Lum_std

Gamma_pst = Gamma_mean + Gamma_std
Gamma_mst = Gamma_mean - Gamma_std

n_pwn_pst = n_pwn_mean + n_pwn_std
n_pwn_mst = n_pwn_mean - n_pwn_std

nob_pst = nob_mean + nob_std
nob_mst = nob_mean - nob_std

    # Plot
label = 'none'
sym = ['', '', '']
linestyle = ['-.', ':', ':']
xlabel = 'Time [Myr]'
text = r'$D_0$ = %.2e $cm^2 s^{-1}$, $\delta$ = %.2f'u'\n'r'$p_0$ =%.2e $GeV c^{-1}$, $\alpha$ = %.2f'u'\n'r' $n_0$ = %.2e $cm^{-3}$' u'\n' r'$n_{SN}$ = %d, $n_{it}$ = %d'%(D0, delta, p0, alpha, n0, n, nit)

        # Gamma luminosity
figure_HESS = figure_number
figure = figure_HESS + 1

Title_HESS = 'Mean gamma emission in the energy range [1 TeV, 10 TeV]'
Title = 'Mean gamma emission in the energy range [100 MeV, 100 TeV]'
ylabel = '$L_\gamma$ [erg s$^{-1}$]'

plot(figure_HESS, 3, t6, [Lum_HESS_mean, Lum_HESS_pst, Lum_HESS_mst], label, Title_HESS, xlabel, ylabel, sym, linestyle, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_range1.pdf')

plot(figure, 3, t6, [Lum_mean, Lum_pst, Lum_mst], label, Title, xlabel, ylabel, sym, linestyle, text)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_all.pdf')

figure_number = figure + 1

        # spectral index
figure_Gamma = figure_number
Title_Gamma = 'Photon index in the energy range [1 TeV, 10 TeV]'
ylabel = '$\Gamma_{ph}$'

plot(figure_Gamma, 3, t6, [Gamma_mean, Gamma_pst, Gamma_mst], label, Title_Gamma, xlabel, ylabel, sym, linestyle, text)
plt.savefig(pathfigure_gamma+'Photon_index.pdf')
figure_number = figure_Gamma + 1

        # Number of pulsar wind nebula
figure_pwn = figure_number
label = 'none'
Title_pwn = 'Number of pulsar wind nebula in the superbubble'
ylabel = '$n_{pwn}$'

plot(figure_pwn, 3, t6, [n_pwn_mean, n_pwn_pst, n_pwn_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, text)
plt.savefig(pathfigure_remain+'Mean_pwn.pdf')

figure_number = figure_pwn + 1

        # Number of remained OB stars in the association
figure_ob = figure_number
label = 'none'
Title_pwn = 'Number of remained Ob stars'
ylabel = '$n_{ob}$'
plot(figure_ob, 3, t6, [nob_mean, nob_pst, nob_mst], label, Title_pwn, xlabel, ylabel, sym, linestyle, text)
plt.savefig(pathfigure_remain+'Mean_ob.pdf')

figure_number = figure_ob + 1

plt.show()

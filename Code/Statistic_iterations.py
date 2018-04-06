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
from statistics import mean, stdev
from naima.models import PionDecay, TableModel
from Functions import *
from Functions_CR import *
from Functions_SB import *
from Functions_gamma import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_SB import *

    # IRAP
#os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/stat_SN/Test/')
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/Test/'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/Test/'
#pathfigure_CR = '/Users/stage/Documents/Virginie/Superbubbles/figures/stat_SN/CR/1/'

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/stat_SN/Test/')
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/Gamma_emission/Test/'
#pathfigure_CR = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/CR/2/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 1000

    # Parameters for the cosmic rays

        # Energy (GeV)
Emin_CR = 1            # minimum kinetic energy: Emin = 1GeV
Emax_CR = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV

        # Power-law distribution of the cosmic rays (GeV^-1 cm^-3)
p0 = 10             # normalization constant (GeV/c)
alpha = 2.0         # exponent of the power-law distribution

        # Diffusion coefficient of the cosmic rays (cm^2 s^-1)
delta = 1.0/2       # exponent of the power-law of the diffusion coefficient
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1

    # Which zone for the Computation
zones = [2]

    # Fix time array (yr)
t0min = 3                   # Myr
t0max = 4                  # Myr
tmin = t0min/yr26yr         # yr
tmax = (t0max + 1)/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Energy of the gamma photons (GeV)
Emin_gamma = 100 * MeV2GeV      # 100 MeV (GeV)
Emax_gamma = 100 * TeV2GeV      # 100 TeV (GeV)
Esep = 1*TeV2GeV                # 1 TeV (GeV)

    # Initialization
figure_number = 1
Lum_it_1 = []       # total gamma luminosity from 100 MeV to 100 GeV
Lum_it_2 = []       # total gamma luminosity from 100 GeV to 100 TeV
Lum_it = []         # total gamma luminosity from 100 MeV to 100 TeV
n_pwn_it = []       # total number of pulsar wind nebula
nob_it = []         # total of remained OB stars

    ##----------##
    # Iterations #
    ##----------##

for i in range (nit):

        # SN explosion time
    n = 5             # number of massive stars in the OB association
    t0 = random_SN(t0min, t0max, n)/yr26yr
    t0 = sorted(t0)

    Lum_1, Lum_2, Lum, Lum_units, figure_number, n_pwn, nob = data(t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, Esep, p0, alpha, D0, delta, zones, pathfigure_gamma, i, figure_number)

    Lum_it_1.append(Lum_1)
    Lum_it_2.append(Lum_2)
    Lum_it.append(Lum)
    n_pwn_it.append(n_pwn)
    nob_it.append(nob)
    """
        # spectral index
    indE = 10       # which energy that the spectral index would be computed
    spectral_index(indE, zones, pathfigure_gamma, i, figure_number)
    """
    print('end of the iteration %d' %i)

Lum_it_1 = numpy.asarray(Lum_it_1)
Lum_it_2 = numpy.asarray(Lum_it_2)
Lum_it = numpy.asarray(Lum_it)
n_pwn_it = numpy.asarray(n_pwn_it)
nob_it = numpy.asarray(nob_it)

    ##---------------------------##
    # Mean and standard deviation #
    ##---------------------------##

        # Initialization
            # Mean
Lum_1_mean = numpy.zeros(number_bin_t)      # from 100 MeV to 100 GeV
Lum_2_mean = numpy.zeros(number_bin_t)      # from 100 GeV to 100 TeV
Lum_mean = numpy.zeros(number_bin_t)        # from 100 MeV to 100 TeV
n_pwn_mean = numpy.zeros(number_bin_t)      # number of pulsar wind nebula
nob_mean = numpy.zeros(number_bin_t)        # number of remained OB stars

            # Standard deviation
Lum_1_stdev = numpy.zeros(number_bin_t)     # from 100 MeV to 100 GeV
Lum_2_stdev = numpy.zeros(number_bin_t)     # from 100 GeV to 100 TeV
Lum_stdev = numpy.zeros(number_bin_t)       # from 100 MeV to 100 TeV
n_pwn_stdev = numpy.zeros(number_bin_t)     # number of pulsar wind nebula
nob_stdev = numpy.zeros(number_bin_t)        # number of remained OB stars

for j in range (number_bin_t):

    Lum_1_mean[j] = mean(Lum_it_1[:,j])
    Lum_2_mean[j] = mean(Lum_it_2[:,j])
    Lum_mean[j] = mean(Lum_it[:,j])
    n_pwn_mean[j] = mean(n_pwn_it[:,j])
    nob_mean[j] = mean(nob_it[:, j])

    Lum_1_stdev[j] = stdev(Lum_it_1[:,j])
    Lum_2_stdev[j] = stdev(Lum_it_2[:,j])
    Lum_stdev[j] = stdev(Lum_it[:,j])
    n_pwn_stdev[j] = stdev(n_pwn_it[:,j])
    nob_stdev[j] = stdev(nob_it[:, j])

Lum_1_pst = Lum_1_mean + Lum_1_stdev
Lum_1_mst = Lum_1_mean - Lum_1_stdev

Lum_2_pst = Lum_2_mean + Lum_2_stdev
Lum_2_mst = Lum_2_mean - Lum_2_stdev

Lum_pst = Lum_mean + Lum_stdev
Lum_mst = Lum_mean - Lum_stdev

n_pwn_pst = n_pwn_mean + n_pwn_stdev
n_pwn_mst = n_pwn_mean - n_pwn_stdev

nob_pst = nob_mean + nob_stdev
nob_mst = nob_mean - nob_stdev

    # Plot
        # Gamma luminosity
figure_1 = figure_number
figure_2 = figure_1 + 1
figure = figure_2 + 1

label = 'none'
sym = ['-.', ':', ':']
Title_1 = 'Mean gamma emission for %d SN explosions (%d iterations) from 100 MeV to 1 TeV'%(n, nit)
Title_2 = 'Mean gamma emission for %d SN explosions (%d iterations) from 1 TeV to 100 TeV'%(n, nit)
Title = 'Mean gamma emission for %d SN explosions (%d iterations) from 100 MeV to 100 TeV'%(n, nit)
xlabel = 'Time [Myr]'
ylabel = '$L_\gamma$ [erg s$^{-1}$]'

log_plot(figure_1, 3, t6, [Lum_1_mean, Lum_1_pst, Lum_1_mst], label, Title_1, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_range1.eps')

log_plot(figure_2, 3, t6, [Lum_2_mean, Lum_2_pst, Lum_2_mst], label, Title_2, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_range2.eps')

log_plot(figure, 3, t6, [Lum_mean, Lum_pst, Lum_mst], label, Title, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_all.eps')

        # Number of pulsar wind nebula
figure_pwn = figure + 1
label = 'none'
sym_stdev = ['-.', ':', ':']
Title_pwn = 'Number of pulsar wind nebula in the superbubble'
xlabel = 'Time [Myr]'
ylabel = '$n_{pwn}$'

log_plot(figure_pwn, 3, t6, [n_pwn_mean, n_pwn_pst, n_pwn_mst], label, Title_pwn, xlabel, ylabel, sym)
plt.savefig(pathfigure_remain+'Mean_pwn.eps')

        # Number of remained OB stars in the association
figure_ob = figure_pwn + 1
label = 'none'
sym_stdev = ['-.', ':', ':']
Title_pwn = 'Number of remained Ob stars'
xlabel = 'Time [Myr]'
ylabel = '$n_{ob}$'
log_plot(figure_ob, 3, t6, [nob_mean, nob_pst, nob_mst], label, Title_pwn, xlabel, ylabel, sym)
plt.savefig(pathfigure_remain+'Mean_ob.eps')

plt.show()

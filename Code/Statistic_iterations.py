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
from Parameters_SB import *

    # IRAP
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/30_Dor_C/')
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Gamma_emission/Test/'
pathfigure_remain = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Remain/Test/'
#pathfigure_CR = '/Users/stage/Documents/Virginie/Superbubbles/figures/stat_SN/CR/1/'

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/stat_SN/Test/')
#pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/Gamma_emission/Test/'
#pathfigure_CR = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/CR/2/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 1

    # Correction factor
t_end = 4.5e6               # time at which R = Robs
t_end_6 = t_end * yr26yr    # in Myr
Robs = 47.0                 # observed radius (pc)
Rsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t_end_6)[0] # from Weaver's model (pc and km/s)
correction_factor = Robs/Rsb

    # Parameters for the cosmic rays

        # Energy (GeV)
Emin_CR = 1            # minimum kinetic energy: Emin = 1GeV
Emax_CR = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV

        # Power-law distribution of the cosmic rays (GeV^-1 cm^-3)
p0 = 10             # normalization constant (GeV/c)
alpha = 2.0         # exponent of the power-law distribution

        # Diffusion coefficient of the cosmic rays (cm^2 s^-1)
delta = 1.0/2       # exponent of the power-law of the diffusion coefficient
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1  ==> prendre *10 et /10

    # Which zone for the Computation
zones = [2]

    # Fix time array (yr)
t0min = 3                   # Myr
t0max = t_end_6             # Myr
tmin = t0min/yr26yr         # yr
tmax = (t0max + 1)/yr26yr   # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Energy of the gamma photons (GeV)
Emin_gamma = 100 * MeV2GeV      # 100 MeV (GeV)
Emax_gamma = 100 * TeV2GeV      # 100 TeV (GeV)
Esep = numpy.array([100, 1*TeV2GeV, 10*TeV2GeV]) # ranges of energy (GeV)

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
        n = 5             # number of massive stars in the OB association
        t0 = random_SN(t0min, t0max, n)/yr26yr
        t0 = sorted(t0)

        t0_it.append(t0)

    t0_it = numpy.asarray(t0_it)

    pickle.dump(t0_it, SN_write)

with open('SN', 'rb') as SN_load:

    t0_it = pickle.load(SN_load)
    n = len(t0_it[0])

for i in range (nit):

    t0 = t0_it[i]

    Lum_HESS, Lum, Gamma, Lum_units, figure_number, n_pwn, nob = data(correction_factor, t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, Esep, p0, alpha, D0, delta, zones, pathfigure_gamma, i, figure_number)

    Lum_it_HESS.append(Lum_HESS)
    print(Lum_HESS)
    Lum_it.append(Lum)
    Gamma_it.append(Gamma)
    n_pwn_it.append(n_pwn)
    nob_it.append(nob)

    print('end of the iteration %d' %i)

Lum_it_HESS = numpy.asarray(Lum_it_HESS)
ind = numpy.where(Lum_HESS > 0)[0]
print(Lum_HESS[ind])
Lum_it = numpy.asarray(Lum_it)
Gamma_it = numpy.asarray(Gamma_it)
n_pwn_it = numpy.asarray(n_pwn_it)
nob_it = numpy.asarray(nob_it)
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
print(Lum_HESS_pst[ind])

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
sym = ['-.', ':', ':']
xlabel = 'Time [Myr]'

        # Gamma luminosity
figure_HESS = figure_number
figure = figure_HESS + 1

Title_HESS = 'Mean gamma emission for %d SN explosions (%d iterations) from 1 TeV to 10 TeV'%(n, nit)
Title = 'Mean gamma emission for %d SN explosions (%d iterations) from 100 MeV to 100 TeV'%(n, nit)
ylabel = '$L_\gamma$ [erg s$^{-1}$]'

log_plot(figure_HESS, 3, t6, [Lum_HESS_mean, Lum_HESS_pst, Lum_HESS_mst], label, Title_HESS, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_range1.eps')

log_plot(figure, 3, t6, [Lum_mean, Lum_pst, Lum_mst], label, Title, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_all.eps')

figure_number = figure + 1

        # spectral index
figure_Gamma = figure_number
Title_Gamma = 'Photon index for %d SN explosions (%d iterations) from 1 TeV to 10 TeV'%(n, nit)
ylabel = '$\Gamma_{ph}$'

log_plot(figure, 3, t6, [Gamma_mean, Gamma_pst, Gamma_mst], label, Title_Gamma, xlabel, ylabel, sym)
plt.savefig(pathfigure_gamma+'Mean_gamma_emission_all.eps')
figure_number = figure_Gamma + 1

        # Number of pulsar wind nebula
figure_pwn = figure_number
label = 'none'
Title_pwn = 'Number of pulsar wind nebula in the superbubble'
ylabel = '$n_{pwn}$'

log_plot(figure_pwn, 3, t6, [n_pwn_mean, n_pwn_pst, n_pwn_mst], label, Title_pwn, xlabel, ylabel, sym)
plt.savefig(pathfigure_remain+'Mean_pwn.eps')

figure_number = figure_pwn + 1

        # Number of remained OB stars in the association
figure_ob = figure_number
label = 'none'
Title_pwn = 'Number of remained Ob stars'
ylabel = '$n_{ob}$'
log_plot(figure_ob, 3, t6, [nob_mean, nob_pst, nob_mst], label, Title_pwn, xlabel, ylabel, sym)
plt.savefig(pathfigure_remain+'Mean_ob.eps')

figure_number = figure_ob + 1

plt.show()
"""

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
#os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/stat_SN/iterations/1/')
#pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/stat_SN/Gamma_emission/iterations/1/'
#pathfigure_CR = '/Users/stage/Documents/Virginie/Superbubbles/figures/stat_SN/CR/1/'

    # Home
os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/stat_SN/iterations/3/')
pathfigure_gamma = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/Gamma_emission/iterations/3/'
#pathfigure_CR = '/home/vivi/Documents/Master_2/Superbubbles/figures/stat_SN/CR/2/'

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 1

    # Parameters for the cosmic rays

        # Energy
Emin_CR = 1            # minimum kinetic energy: Emin = 1GeV
Emax_CR = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV

        # Power-law distribution of the cosmic rays (GeV^-1 cm^-3)
p0 = 10             # normalization constant (GeV/c)
alpha = 2.0         # exponent of the power-law distribution

        # Diffusion coefficient of the cosmic rays (cm^2 s^-1)
delta = 1.0/3       # exponent of the power-law of the diffusion coefficient
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1

for i in range (nit):

            # SN explosion time
    with open('SN', 'wb') as t0_write:
        n = 3               # number of massive stars in the OB association
        t0min = 3           # 3 Myr
        t0max = 37          # 37 Myr
        t0 = random_SN(3, 37, n)/yr26yr     # Sn explosions time (yr)
        t0 = sorted(t0)
        pickle.dump(t0, t0_write)

        # Data computation
    data(Emin_CR, Emax_CR, p0, alpha, D0, delta)

        # Spectral energy distribution computation
    Emin_gamma = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
    Emax_gamma = 100 * TeV2GeV            # 100 TeV = 100 000 GeV

    spectrum(Emin_gamma, Emax_gamma)#, t0, ECR, ECR_unit, t, t_unit, Ntotsn, Ntot_unit, ngastot, ngas_unit, radius)

        # gamma luminosity computation
    Esep = 100      # in GeV
    gamma_luminosity_it(Esep, pathfigure_gamma, i)#, t0, t, t_unit, radius, spectrum, spectrum_unit, sed_PD_sn, sed_PD_unit)
    print('end of the iteration')

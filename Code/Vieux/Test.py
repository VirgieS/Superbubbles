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

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/')
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/'
##----##
# Test #
##----##

Emin_CR = 1            # minimum kinetic energy: Emin = 1GeV
Emax_CR = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV

p0 = 10             # normalization constant (GeV/c)
alpha = 2.0         # exponent of the power-law distribution

delta = 1.0/3       # exponent of the power-law of the diffusion coefficient
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1

t0 = [3e6, 4e6]     # Sn explosions time

data(Emin_CR, Emax_CR, p0, alpha, D0, delta, t0)

Emin_gamma = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
Emax_gamma = 100 * TeV2GeV            # 100 TeV = 100 000 GeV

spectrum(Emin_gamma, Emax_gamma)

Esep = 100      # in GeV
gamma_luminosity(Esep, pathfigure)

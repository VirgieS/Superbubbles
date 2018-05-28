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
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/Parameters/stars/30/10')

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Parameters/stars/10/5')

## ======= ##
# Statistic #
## ======= ##

"""
Compute the important parameters for each iterations
"""

    # Number of iterations
nit = 150                                                                      #you need to change it for your simulations

    # Which zone for the Computation
zones = [2]                                                                     #you need to change it for your simulations

    # Correction factor

need_correction = False

if need_correction:

    t_end_6 = 4.0           # Myr
    t_sn = t_end_6 + 1      # Myr
    t_end = t_end_6/yr26yr  # yr
    Rsb = 47.0                          # observed radius (pc)                  #you need to change it for your simulations
    Rw = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
    correction_factor = Rsb/Rw

else:
    correction_factor = 1
    t_sn = 37

    # Fix time array (yr)
t0min = 3           # Myr
t0max = t_sn
tmin = t0min/yr26yr     # yr
tmax = 10/yr26yr        # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr

    # Initialization
figure_number = 1

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

        if need_correction:

            mean_sn = 5                                                         # you need to change it for your simulations
            nsn = numpy.random.poisson(lam = mean_sn)       # random number of SN which mean value equals to mean_sn

        else:
            nsn = 1

        t0 = random_uniform(t0min, t0max, nsn)/yr26yr   # random SN explosion times with an uniform distribution from t0min to t0max
        t0 = sorted(t0)
        t0 = numpy.asarray(t0)
        indsn = numpy.where(t0 <= tmax)[0]
        t0 = t0[indsn]

        t0_it.append(t0)

    t0_it = numpy.asarray(t0_it)

    pickle.dump(t0_it, SB_write)

        # Computation

with open('General', 'wb') as data_write:

    for i in range (nit):

        t0 = t0_it[i]

        Lum, Flux, spectrum, n_pwn, Lum_pwn, Lum_psr, nob = data(correction_factor, t0, t_fix, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, zones)
        print(Lum)

        Lum_it.append(Lum)
        Flux_it.append(Flux)
        n_pwn_it.append(n_pwn)
        Lum_pwn_it.append(Lum_pwn)
        Lum_psr_it.append(Lum_psr)
        nob_it.append(nob)

        print('end of the iteration %d' %i)

    Lum_it = numpy.asarray(Lum_it)
    Flux_it = numpy.asarray(Flux_it)

        # H.E.S.S. energy range

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

        # Fermi energy range

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

    n_pwn_it = numpy.asarray(n_pwn_it)
    Lum_pwn_it = numpy.asarray(Lum_pwn_it)
    Lum_psr_it = numpy.asarray(Lum_psr_it)
    nob_it = numpy.asarray(nob_it)

    pickle.dump(Lum_HESS_it, data_write)
    pickle.dump(Lum_Fermi_it, data_write)
    pickle.dump(Lum_it, data_write)
    #pickle.dump(Flux_it, data_write)
    pickle.dump(Gamma_HESS_it, data_write)
    pickle.dump(Gamma_Fermi_it, data_write)
    pickle.dump(n_pwn_it, data_write)
    pickle.dump(Lum_pwn_it, data_write)
    pickle.dump(Lum_psr_it, data_write)
    pickle.dump(nob_it, data_write)
    pickle.dump(spectrum, data_write)

"""
It computes the gamma-ray luminosities in the HE and VHE energy ranges and the spectral index in both energy ranges for all iterations.

All the parameters must be given in the Parameters_system
"""

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

    # You need to change it
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/30_Dor_C/new/1e25_22_050/2')

## ======= ##
# Statistic #
## ======= ##

    # Number of iterations
nit = 100                                                                      #you need to change it for your simulations

    # Which zone for the Computation
zones = [2]                                                                     #you need to change it for your simulations

    # Correction factor

need_correction = True

if need_correction:             # if any correction factor must be used

    t_end_6 = 4.0                       # Myrs
    Rsb = 47.0                          # observed radius (pc)                  #you need to change it for your simulations
    Rw = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
    correction_factor = Rsb/Rw

else:
    correction_factor = 1

    # Initialization
figure_number = 1

        # For the gamma-ray emission of the superbubble
Lum_it = []                 # total gamma luminosity for the whole energy range
Flux_it = []                # photon flux for the whole energy range
Lum_pwn_it = []             # TeV emission of PWNe
Lum_psr_it = []             # GeV emission of PSRs
nob_it = []                 # total of remained OB stars
tsn_it = []                 # SN explosion times (yr)
nsn_it = []                 # number of supernova per iterations

        # For the parameters of the SB
Rsb = numpy.zeros(number_bin_t)     # size of the superbubble (pc)
Vsb = numpy.zeros(number_bin_t)     # velocity of the shock (km/s)
Ms = numpy.zeros(number_bin_t)      # mass of the supershell (solar mass)
ns = numpy.zeros(number_bin_t)      # density of the supershell (cm^-3)

    ##----------##
    # Iterations #
    ##----------##
print('For %d SN'%Nob)
print('For %d iterations' %nit)

        # SN explosions time (yr)

with open('SB', 'wb') as SB_write:

    for i in range (nit):

        nsn = Nob

        tsn = numpy.random.uniform(tsnmin, tsnmax, nsn)/yr26yr  # random SN explosion times with an uniform distribution from t0min to t0max
        tsn = sorted(tsn)                                       # sorted random tsn
        tsn = numpy.asarray(tsn)
        indsn = numpy.where(tsn <= tmax)[0]
        tsn = tsn[indsn]
        nsn = len(tsn)

        tsn_it.append(tsn)
        nsn_it.append(nsn)

    tsn_it = numpy.asarray(tsn_it)
    nsn_it = numpy.asarray(nsn_it)

    pickle.dump(tsn_it, SB_write)
    pickle.dump(nsn_it, SB_write)

        # Computation

with open('General', 'wb') as data_write:

    for i in range (nit):

        tsn = tsn_it[i]

        Lum, Flux, Lum_pwn, Lum_psr, nob, R_sb, V_sb, M_s, n_s = data(correction_factor, tsn, t_fix, zones)
        ind = numpy.where(R_sb > 0.0)[0]

        for j in (ind):

            Rsb[j] = R_sb[j]
            Vsb[j] = V_sb[j]
            Ms[j] = M_s[j]
            ns[j] = n_s[j]

        Lum_it.append(Lum)
        Flux_it.append(Flux)
        Lum_pwn_it.append(Lum_pwn)
        Lum_psr_it.append(Lum_psr)
        nob_it.append(nob)

        print('end of the iteration %d' %i)

    Lum_it = numpy.asarray(Lum_it)
    Flux_it = numpy.asarray(Flux_it)

        # In each energy range

            # VHE range
    Emin = 1 * TeV2GeV                  # 1 TeV (GeV)
    Emax = 10 * TeV2GeV                 # 10 TeV (GeV)
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

                # Gamma-ray luminosity
    spectrum_HESS = spectrum[indE]
    spectrum_erg = spectrum_HESS * 1.0/erg2GeV     # only in the energy range (erg)
    spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
    lum_HESS = Flux_it[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
    Lum_HESS_it = luminosity(lum_HESS, spectrum_ev) # erg s^-1

                # Spectral photon index
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Emin = spectrum[indE[0]]
    Emax = spectrum[indE[-1]]
    Gamma_HESS_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_HESS_it = numpy.nan_to_num(Gamma_HESS_it)

            # HE range
    Emin = 100 * MeV2GeV                # 100 MeV (GeV)
    Emax = 100                          # 100 GeV
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

                # Gamma-ray luminosity
    spectrum_Fermi = spectrum[indE]
    spectrum_erg = spectrum_Fermi * 1.0/erg2GeV     # only in the energy range (erg)
    spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
    lum_Fermi = Flux_it[:, :, indE] * spectrum_erg   # erg s^-1 eV^-1
    Lum_Fermi_it = luminosity(lum_Fermi, spectrum_ev) # erg s^-1

                # Spectral photon index (1 GeV to 10 GeV)
    Emin = 1        # 1 GeV
    Emax = 10       # 10 GeV
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Emin = spectrum[indE[0]]
    Emax = spectrum[indE[-1]]
    Gamma_GeV_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_GeV_it = numpy.nan_to_num(Gamma_GeV_it)

                # Spectral photon index (100 MeV to 1 GeV)
    Emin = 100 * MeV2GeV    # 100 MeV
    Emax = 1                # 1 GeV
    indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]
    Fluxmin = Flux_it[:, :, indE[0]]
    Fluxmax = Flux_it[:, :, indE[-1]]
    Emin = spectrum[indE[0]]
    Emax = spectrum[indE[-1]]
    Gamma_MeV_it = spectral_index(Emin, Emax, Fluxmin, Fluxmax)
    Gamma_MeV_it = numpy.nan_to_num(Gamma_MeV_it)

    Lum_pwn_it = numpy.asarray(Lum_pwn_it)
    Lum_psr_it = numpy.asarray(Lum_psr_it)
    nob_it = numpy.asarray(nob_it)

    pickle.dump(Lum_HESS_it, data_write)
    pickle.dump(Lum_Fermi_it, data_write)
    pickle.dump(Lum_it, data_write)
    pickle.dump(Gamma_HESS_it, data_write)
    pickle.dump(Gamma_GeV_it, data_write)
    pickle.dump(Gamma_MeV_it, data_write)
    pickle.dump(Lum_pwn_it, data_write)
    pickle.dump(Lum_psr_it, data_write)
    pickle.dump(nob_it, data_write)
    pickle.dump(Rsb, data_write)
    pickle.dump(Vsb, data_write)
    pickle.dump(Ms, data_write)
    pickle.dump(ns, data_write)


    # CHECKING
print('number of SN: %d' %Nob)
print('diffusion coefficient: %.2e cm^2 s^-1' %D0)
print('injection index: %.2f' %alpha)
print('energy dependence: %.2f'%delta)
print('density: %.2f'%n0)

##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
import os
import pickle
import naima
import astropy.units as units
from naima.models import PionDecay, TableModel
from Functions import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##--------------##
# Path for files #
##--------------##

    # At IRAP
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/2_SN/')

    # At home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/n_SN')

##-----------##
# Computation #
##-----------##

with open('spectra', 'wb') as spectra_write:
    with open('gas', 'rb') as gas_load:
        with open('energy', 'rb') as energy_load:
            with open('data', 'rb') as data_load:
                with open('distance', 'rb') as distance_load:
                    with open('time', 'rb') as time_load:

                            # Loading of data: ngas, energy, nCR, radius and time

                        ngassn = pickle.load(gas_load)
                        ngas_unit = pickle.load(gas_load)

                        ECR = pickle.load(energy_load)
                        ECR_unit = pickle.load(energy_load)
                        ECR = ECR * ECR_unit

                        Ntotsn = pickle.load(data_load)
                        Ntot_unit = pickle.load(data_load)

                        r = pickle.load(distance_load)

                        t0 = pickle.load(time_load)
                        t = pickle.load(time_load)
                        t_unit = pickle.load(time_load)

                            # Energy spectrum (GeV)
                        Emin = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
                        Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV
                        number_bin_E = 20
                        spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV

                        nt = len(t)
                        nt0 = len(t0)

                        sedsn = []

                        for i in range (nt0):       # for each SN explosion

                            sed = []

                            for j in range (nt):    # for each time step
                                """
                                if t[j] < t0[i]:
                                    #sed.append(0)
                                    continue
                                else:
                                """
                                    # Initialization
                                distance = r[j]
                                nr = len(distance) + 2

                                    # Important quantities
                                ngas = ngassn[j] * ngas_unit
                                ntot = numpy.asarray(Ntotsn[j, i]) * Ntot_unit

                                    # Recording
                                spectrum_r = []
                                sed_r = []

                                for k in range (nr):     # for each radius step
                                    model = TableModel(ECR, ntot[k], amplitude = 1)
                                    PD = PionDecay(model, nh = ngas[k], nuclear_enhancement = True)

                                    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                    sed_r.append(sed_PD)

                                sed.append(sed_r)

                            sedsn.append(sed)

                        spectrum = numpy.asarray(spectrum_energy)
                        spectrum_unit = spectrum_energy.unit
                        sedsn = numpy.asarray(sedsn)
                        sed_unit = sed_PD.unit
                        pickle.dump(spectrum, spectra_write)
                        pickle.dump(spectrum_unit, spectra_write)
                        pickle.dump(sedsn, spectra_write)
                        pickle.dump(sed_unit, spectra_write)

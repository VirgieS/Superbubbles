##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from Functions import *
import os
import pickle
import naima
import astropy.units as units
from naima.models import PionDecay, TableModel

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##--------------##
# Path for files #
##--------------##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files')

##-----------##
# Computation #
##-----------##

with open('spectra', 'wb') as spectra_write:
    with open('gas', 'rb') as gas_load:
        with open('energy', 'rb') as energy_load:
            with open('data', 'rb') as data_load:
                with open('distance', 'rb') as distance_load:
                    with open('time', 'rb') as time_load:

                        #my_spectra_write = pickle.Pickler(spectra_read)

                            # Loading of data: ngas, energy, nCR, radius and time
                        my_gas_load = pickle.Unpickler(gas_load)
                        ngas = my_gas_load.load()

                        my_energy_load = pickle.Unpickler(energy_load)
                        E = my_energy_load.load()

                        my_data_load = pickle.Unpickler(data_load)
                        Ntot = my_data_load.load()

                        my_distance_load = pickle.Unpickler(distance_load)
                        r = my_distance_load.load()

                        my_time_load = pickle.Unpickler(time_load)
                        t = my_time_load.load()

                            # Energy spectrum (GeV)
                        Emin = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
                        Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV
                        number_bin_E = 20
                        spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV

                        n = len(t)

                        spectrum = []
                        sed = []

                        for i in range (n):
                            m = len(r[i])+2
                            spectrum_r = []
                            sed_r = []

                            for j in range (m):

                                #if numpy.all(Ntot[i, j, :] == 0): # Condition if there are relativistic particles
                                #    continue

                                #else:
                                #    model = TableModel(E * units.GeV, Ntot[i, j, :] * 1/units.GeV, amplitude = 1)
                                #    PD = PionDecay(model, nh = ngas[i,j+1] * 1/units.cm**3, nuclear_enhancement = True)

                                #    spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV
                                #    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.kpc)

                                #    spectrum_r.append(spectrum_energy)
                                    #print(numpy.asarray(spectrum_r).shape)
                                #    sed_r.append(sed_PD)

                                model = TableModel(E * units.GeV, Ntot[i, j, :] * 1/units.GeV, amplitude = 1)
                                PD = PionDecay(model, nh = ngas[i,j] * 1/units.cm**3, nuclear_enhancement = True)

                                #spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV

                                sed_PD = PD.sed(spectrum_energy, distance = 0 * units.kpc)

                                #spectrum_r.append(spectrum_energy)
                                sed_r.append(sed_PD)

                                """
                                if ((i == 2) and (j == 20)):

                                        # Computation of the Pion Decay
                                    model = TableModel(E * units.GeV, Ntot[i, j, :] * 1/units.GeV, amplitude = 1)
                                    PD = PionDecay(model, nh = ngas[i,j] * 1/units.cm**3, nuclear_enhancement = True)

                                    spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV
                                    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.kpc)
                                    Lum = numpy.asarray(sed_PD)
                                    unit_luminosity = units.erg*1/units.s*1/units.GeV

                                        # Plot
                                    plt.figure(figsize=(8,5))
                                    plt.rc('font', family='sans')
                                    plt.rc('mathtext', fontset='custom')
                                    plt.loglog(spectrum_energy,sed_PD,lw=2,c=naima.plot.color_cycle[0])
                                    plt.title('Production of photons by Pion Decay')
                                    plt.xlabel('Photon energy [{0}]'.format(spectrum_energy.unit.to_string('latex_inline')))
                                    plt.ylabel('$L_E$ [{0}]'.format(unit_luminosity.unit.to_string('latex_inline')))
                                    #plt.ylabel('$E^2 dN/dE$ [{0}]'.format(sed_PD.unit.to_string('latex_inline')))
                                    plt.tight_layout()

                                    # Flux
                                        # First zone
                                    E1 = 100        # 100 GeV (GeV)
                                    ind1 = numpy.where(E <= E1)[0]
                                    Flux1 = integration_log(E[ind1], Lum[ind1])
                                    print('Photon flux for the first zone (100 MeV to 100 GeV): %.5e'%Flux1)

                                        # Second zone
                                    ind2 = numpy.where(E > E1)[0]
                                    Flux2 = integration_log(E[ind2], Lum[ind2])
                                    print('Photon flux for the second zone (100 GeV to 100 TeV): %.5e'%Flux2)

                                    plt.show()
                            """
                            #spectrum.append(spectrum_r)
                            sed.append(sed_r)

                        spectrum = numpy.asarray(spectrum_energy)
                        spectrum_unit = spectrum_energy.unit
                        sed = numpy.asarray(sed)
                        sed_unit = sed_PD.unit
                        #lum_unit = units.erg*1/units.s*1/units.GeV
                        pickle.dump(spectrum, spectra_write)
                        pickle.dump(spectrum_unit, spectra_write)
                        pickle.dump(sed, spectra_write)
                        pickle.dump(sed_unit, spectra_write)

                        #print(numpy.asarray(spectrum).shape)
                        #print(numpy.asarray(sed).shape)

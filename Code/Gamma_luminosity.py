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

Emin = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV

with open('gas', 'rb') as gas_load:
    with open('energy', 'rb') as energy_load:
        with open('data', 'rb') as data_load:
            with open('distance', 'rb') as distance_load:
                with open('time', 'rb') as time_load:

                        # Loading of data: ngas, energy, nCR, radius and time
                    my_gas_load = pickle.Unpickler(gas_load)
                    ngas = my_gas_load.load()

                    my_energy_load = pickle.Unpickler(energy_load)
                    E = my_energy_load.load()

                    my_data_load = pickle.Unpickler(data_load)
                    Ntot = my_data_load.load()

                    my_distance_load = pickle.Unpickler(distance_load)
                    r = my_distance_load.load()
                    print(r.shape)

                    my_time_load = pickle.Unpickler(time_load)
                    t = my_time_load.load()

                    n = len(t)
                    l = len(E)

                    for i in range (n):
                        m = len(r[i])

                        for j in range (m):

                            if numpy.all(Ntot[i, j, :] == 0): # Condition if there are relativistic particles
                                continue

                    """
                    for i in range (len(Ne[0, :, 0])):
                        for j in range (len(Ne[0, i, :])):

                            if numpy.all(Ne[:, i, j] == 0): # Condition if there are relativistic particles
                                continue
                            else:
                                #model = TableModel(E * units.GeV, Ne[:, i, j] * 1/units.GeV, amplitude = 1)
                                #IC = InverseCompton(model, seed_photon_fields=['CMB'])

                                if ((i == 2000) and (j == 20)):

                                    model = TableModel(E * units.GeV, Ne[:, i, j] * 1/units.GeV, amplitude = 1)
                                    PD = PionDecay(model, nh = 1 * 1/units.cm**3, nuclear_enhancement = True)

                                    spectrum_energy = numpy.logspace(log10(Emin), log10(Emax), 1000) * units.GeV
                                    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.kpc)

                            #model = TableModel(E * units.GeV, Ne[:, 2000, 20] * 1/units.GeV, amplitude = 1)
                            #IC = InverseCompton(model, seed_photon_fields=['CMB'])
                            #spectrum_energy = numpy.logspace(log10(Emin), log10(Emax), 1000) * units.GeV
                            #sed_IC = IC.sed(spectrum_energy, distance = 0 * units.kpc)
                                    plt.loglog(spectrum_energy,sed_PD,lw=2, label='PD',c=naima.plot.color_cycle[0])
                                    plt.legend(loc='lower left')
                                    plt.show()
                                    """

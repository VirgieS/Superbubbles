##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from math import *
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from pylab import *
from Functions import *
import os
import pickle
import naima
import astropy.units as units
from naima.models import InverseCompton, TableModel

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

Emin = 100 * MeV2GeV          # 100 MeV = 0.1 GeV
Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV

with open('energy.dat', 'rb') as energy_load:

    with open('CR.dat', 'rb') as CR_load:

        my_energy_load = pickle.Unpickler(energy_load)
        E = my_energy_load.load()

        my_CR_load = pickle.Unpickler(CR_load)
        Ne = my_CR_load.load()

        for i in range (len(Ne[0, :, 0])):
            for j in range (len(Ne[0, i, :])):

                if numpy.all(Ne[:, i, j] == 0): # Condition if there are relativistic particles
                    continue
                else:
                    #model = TableModel(E * units.GeV, Ne[:, i, j] * 1/units.GeV, amplitude = 1)
                    #IC = InverseCompton(model, seed_photon_fields=['CMB'])

                    if ((i == 2000) and (j == 20)):

                        model = TableModel(E * units.GeV, Ne[:, i, j] * 1/units.GeV, amplitude = 1)
                        IC = InverseCompton(model, seed_photon_fields=['CMB'])

                        spectrum_energy = numpy.logspace(log10(Emin), log10(Emax), 1000) * units.GeV
                        sed_IC = IC.sed(spectrum_energy, distance = 0 * units.kpc)

                #model = TableModel(E * units.GeV, Ne[:, 2000, 20] * 1/units.GeV, amplitude = 1)
                #IC = InverseCompton(model, seed_photon_fields=['CMB'])
                #spectrum_energy = numpy.logspace(log10(Emin), log10(Emax), 1000) * units.GeV
                #sed_IC = IC.sed(spectrum_energy, distance = 0 * units.kpc)
                        plt.loglog(spectrum_energy,sed_IC,lw=2, label='IC (total)',c=naima.plot.color_cycle[0])
                        plt.show()

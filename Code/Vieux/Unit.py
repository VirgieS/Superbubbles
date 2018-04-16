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

with open('spectra', 'rb') as spectra_read:

        # Loading data
    spectrum_energy = pickle.load(spectra_read)
    spectrum_unit = pickle.load(spectra_read)
    sed_PD = pickle.load(spectra_read)
    sed_PD_unit = pickle.load(spectra_read)

    print(sed_PD_unit)

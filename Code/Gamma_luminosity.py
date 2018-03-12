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
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/n_SN')

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

                        rsn = pickle.load(distance_load)

                        t0 = pickle.load(time_load)
                        t = pickle.load(time_load)
                        t_unit = pickle.load(time_load)

                            # Energy spectrum (GeV)
                        Emin = 100 * MeV2GeV            # 100 MeV = 0.1 GeV
                        Emax = 100 * TeV2GeV            # 100 TeV = 100 000 GeV
                        number_bin_E = 20
                        spectrum_energy = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E) * units.GeV

                        #nt = len(t)
                        nt0 = len(t0)

                        sedsn = []

                        for k in range (nt0):       # for each SN explosion

                                # Initialization
                            print(k)
                            r = rsn[k]
                            ngas = ngassn[k] * ngas_unit
                            Ntot = Ntotsn[k] * Ntot_unit
                            indsn = numpy.where(t >= t0[k])[0]      # only time where the SN already explodes
                            nt = len(t[indsn])

                            sed = []

                            for i in range (nt):    # for each time step

                                nr = len(r[i])+2
                                spectrum_r = []
                                sed_r = []

                                for j in range (nr):     # for each radius step
                                    model = TableModel(ECR, Ntot[i, j], amplitude = 1)
                                    PD = PionDecay(model, nh = ngas[i,j], nuclear_enhancement = True)

                                    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.kpc)
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
                        print(numpy.asarray(sedsn).shape)

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
from scipy.integrate import trapz

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
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/Gamma_emission/'

with open('time', 'rb') as time_load:
    with open('distance', 'rb') as distance_load:
        with open('spectra', 'rb') as spectra_read:

                # Loading data
                    # time
            t = pickle.load(time_load)

                    # distance
            rsb = pickle.load(distance_load)

                    # Energ (GeV) and spectrum (erg s^-1 GeV^-1)
            spectrum_energy = pickle.load(spectra_read)
            spectrum_unit = pickle.load(spectra_read)
            spectrum = spectrum_energy*spectrum_unit
            sed_PD = pickle.load(spectra_read)
            sed_PD_unit = pickle.load(spectra_read)
            #sed_PD = sed_PD * sed_PD_unit
            lum_gamma = sed_PD/spectrum_energy
            #print(lum_gamma.unit)
            #lum_gamma = numpy.asarray(lum_gamma)
            #lum_gamma = lum_gamma * sed_PD_unit * 1/spectrum_unit
            sed_PD = sed_PD * sed_PD_unit

            figure_number = 1

            n = len(t)
            Esep = 100          # The separation of the two spectral zones (GeV)

                # Computation of the gamma luminosity for two ranges of energy (from 0.1 GeV to 100 GeV and from 100 GeV to 100 000 GeV)
            ind1 = numpy.where(spectrum_energy < Esep)[0]         # for the first range
            ind2 = numpy.where(spectrum_energy >= Esep)[0]        # for the second range

                    # For recording
            Fluxsbt1 = []
            Fluxsbt2 = []
            Fluxshellt1 = []
            Fluxshellt2 = []
            Fluxoutt1 = []
            Fluxoutt2 = []
            Fluxtott1 = []
            Fluxtott2 = []

            for i in range (n):

                    # First zone: in the cavity (erg s^-1)
                Fluxsb1 = []
                Fluxsb2 = []
                m = len(rsb[i])
                """
                fig = plt.figure(figsize=(8,5))
                plt.rc('font', family='sans')
                plt.rc('mathtext', fontset='custom')
                plt.loglog(spectrum, sed_PD[i,2], lw=2, c=naima.plot.color_cycle[0])
                plt.title('Production of photons by Pion Decay')
                plt.xlabel('Photon energy [{0}]'.format(spectrum.unit.to_string('latex_inline')))
                #plt.ylabel('$E L_E$ [{0}]'.format(spectrum_unit.unit.to_string('latex_inline')))
                plt.ylabel('$E dN/dE$ [{0}]'.format(sed_PD.unit.to_string('latex_inline')))
                plt.tight_layout()
                fig.savefig(pathfigure+'T%d.eps'%i)
                figure_number += 1
                """

                for j in range (m):

                    #Flux1 = integration_log(spectrum_energy[ind1], lum_gamma[i, j, ind1])
                    Flux1= trapz(lum_gamma[i, j, ind1], spectrum_energy[ind1])
                    Fluxsb1.append(Flux1)

                    #Flux2 = integration_log(spectrum_energy[ind2], lum_gamma[i, j, ind2])
                    Flux2 = trapz(lum_gamma[i, j, ind2], spectrum_energy[ind2])
                    Fluxsb2.append(Flux2)

                    # Sum of each contribution (erg s^-1)
                Lumsb1 = numpy.sum(Fluxsb1)
                Lumsb2 = numpy.sum(Fluxsb2)
                #print(Lumsb1)

                Fluxsbt1.append(Lumsb1)
                Fluxsbt2.append(Lumsb2)

                    # Second zone: in the supershell (erg s^-1)
                #Fluxshell1 = integration_log(spectrum_energy[ind1], lum_gamma[i, m, ind1])
                Fluxshell1 = trapz(lum_gamma[i, m, ind1], spectrum_energy[ind1])
                #Fluxshell2 = integration_log(spectrum_energy[ind2], lum_gamma[i, m, ind2])
                Fluxshell2 = trapz(lum_gamma[i, m, ind2], spectrum_energy[ind2])

                Fluxshellt1.append(Fluxshell1)
                Fluxshellt2.append(Fluxshell2)

                    # Third zone: outside the SB (erg s^-1)
                #Fluxout1 = integration_log(spectrum_energy[ind1], lum_gamma[i, m+1, ind1])
                Fluxout1 = trapz(lum_gamma[i, m+1, ind1], spectrum_energy[ind1])
                #Fluxout2 = integration_log(spectrum_energy[ind2], lum_gamma[i, m+1, ind2])
                Fluxout2 = trapz(lum_gamma[i, m+1, ind2], spectrum_energy[ind2])

                Fluxoutt1.append(Fluxout1)
                Fluxoutt2.append(Fluxout2)

                    # Total gamma luminosity (erg s^-1)
                Fluxtot1 = numpy.sum([Lumsb1, Fluxshell1, Fluxout1])
                Fluxtot2 = numpy.sum([Lumsb2, Fluxshell2, Fluxout2])

                Fluxtott1.append(Fluxtot1)
                Fluxtott2.append(Fluxtot2)

                """
                Fluxsbt1 = numpy.asarray(Fluxsbt1)
                Fluxsbt2 = numpy.asarray(Fluxsbt2)

                Fluxshellt1 = numpy.asarray(Fluxshellt1)
                Fluxshellt2 = numpy.asarray(Fluxshellt2)

                Fluxoutt1 = numpy.asarray(Fluxoutt1)
                Fluxoutt2 = numpy.asarray(Fluxoutt2)
                """

            log_plot(figure_number, 4, t, [Fluxsbt1, Fluxshellt1, Fluxoutt1, Fluxtott1], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time (yr)', r'$L_{\gamma}$ (0.1-100 GeV)', ['+', '+', '+', '-'])
            figure_number += 1
            log_plot(figure_number, 4, t, [Fluxsbt2, Fluxshellt2, Fluxoutt2, Fluxtott2], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time (yr)', r'$L_{\gamma}$ (100-100 000 GeV)', ['+', '+', '+', '-'])
            plt.show()
    """
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
    """

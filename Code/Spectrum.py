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

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_SB import *

##--------------##
# Path for files #
##--------------##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/n_SN/')

##-----------##
# Computation #
##-----------##
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/n_SN/Gamma_emission/'

with open('time', 'rb') as time_load:
    with open('distance', 'rb') as distance_load:
        with open('spectra', 'rb') as spectra_read:

                # Loading data
                    # time
            t = pickle.load(time_load)
            t_unit = pickle.load(time_load)
            time = t * t_unit

                    # distance
            rsb = pickle.load(distance_load)

                    # Energy (GeV)
            spectrum_energy = pickle.load(spectra_read)
            spectrum_unit = pickle.load(spectra_read)
            spectrum = spectrum_energy * spectrum_unit

                     # spectal energy distribution (erg s^-1): E^2 dN/dE = E L_E
            sed_PD = pickle.load(spectra_read)
            sed_PD_unit = pickle.load(spectra_read)
            lum_gamma = sed_PD/spectrum_energy      # intrinsic luminosity (erg S^-1)
            sed_PD = sed_PD * sed_PD_unit           # with unit
            lum_gamma_unit = sed_PD/spectrum        # unit

                # Initialization
            figure_number = 1
            n = len(t)
            Esep = 100          # The separation of the two spectral zones (GeV)

                ##================##
                # Gamma luminosity #
                ##================##

            ind1 = numpy.where(spectrum_energy < Esep)[0]         # for the first range
            ind2 = numpy.where(spectrum_energy >= Esep)[0]        # for the second range

                    # For recording
                        # Within the SB
            Fluxsbt1 = []       # 100 MeV to 100 GeV
            Fluxsbt2 = []       # 100 GeV to 100 TeV
            Fluxsbt = []        # 100 MeV to 100 TeV

                        # Within the supershell
            Fluxshellt1 = []    # 100 MeV to 100 GeV
            Fluxshellt2 = []    # 100 GeV to 100 TeV
            Fluxshellt = []     # 100 MeV to 100 TeV

                        # Outside the SB
            Fluxoutt1 = []      # 100 MeV to 100 GeV
            Fluxoutt2 = []      # 100 GeV to 100 TeV
            Fluxoutt = []       # 100 MeV to 100 TeV

                        # Total
            Fluxtott1 = []      # 100 MeV to 100 GeV
            Fluxtott2 = []      # 100 GeV to 100 TeV
            Fluxtott = []       # 100 MeV to 100 TeV

            for i in range (n):         # for each time step

                    ##------------------------------------##
                    # First zone: in the cavity (erg s^-1) #
                    ##------------------------------------##
                Fluxsb1 = []
                Fluxsb2 = []
                Fluxsb = []
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

                for j in range (m):     # for each radius step

                        # From 100 MeV to 100 GeV
                    Flux1= integrate.trapz(lum_gamma[i, j, ind1], spectrum_energy[ind1])
                    Flux1 = numpy.nan_to_num(Flux1)
                    Fluxsb1.append(Flux1)

                        # From 100 GeV to 100 TeV
                    Flux2 = integrate.trapz(lum_gamma[i, j, ind2], spectrum_energy[ind2])
                    Flux2 = numpy.nan_to_num(Flux2)
                    Fluxsb2.append(Flux2)

                        # From 100 MeV to 100 TeV
                    Flux = integrate.trapz(lum_gamma[i,j], spectrum_energy)
                    Flux = numpy.nan_to_num(Flux)
                    Fluxsb.append(Flux)

                        # Total contribution in the SB (erg s^-1)
                Lumsb1 = numpy.sum(Fluxsb1)
                Lumsb2 = numpy.sum(Fluxsb2)
                Lumsb = numpy.sum(Fluxsb)

                Fluxsbt1.append(Lumsb1)
                Fluxsbt2.append(Lumsb2)
                Fluxsbt.append(Lumsb)

                    ##-----------------------------------------##
                    # Second zone: in the supershell (erg s^-1) #
                    ##-----------------------------------------##

                        # From 100 MeV to 100 GeV
                Fluxshell1 = integrate.trapz(lum_gamma[i, m, ind1], spectrum_energy[ind1])
                Fluxshell1 = numpy.nan_to_num(Fluxshell1)
                Fluxshellt1.append(Fluxshell1)

                        # From 100 GeV to 100 TeV
                Fluxshell2 = integrate.trapz(lum_gamma[i, m, ind2], spectrum_energy[ind2])
                Fluxshell2 = numpy.nan_to_num(Fluxshell2)
                Fluxshellt2.append(Fluxshell2)

                        # From 100 MeV to 100 TeV
                Fluxshell = integrate.trapz(lum_gamma[i, m], spectrum_energy)
                Fluxshell = numpy.nan_to_num(Fluxshell)
                Fluxshellt.append(Fluxshell)

                    ##-------------------------------------##
                    # Third zone: outside the SB (erg s^-1) #
                    ##-------------------------------------##

                        # From 100 MeV to 100 GeV
                Fluxout1 = integrate.trapz(lum_gamma[i, m+1, ind1], spectrum_energy[ind1])
                Fluxout1 = numpy.nan_to_num(Fluxout1)
                Fluxoutt1.append(Fluxout1)

                        # From 100 GeV to 100 TeV
                Fluxout2 = integrate.trapz(lum_gamma[i, m+1, ind2], spectrum_energy[ind2])
                Fluxout2 = numpy.nan_to_num(Fluxout2)
                Fluxoutt2.append(Fluxout2)

                        # From 100 MeV to 100 TeV
                Fluxout = integrate.trapz(lum_gamma[i, m+1], spectrum_energy)
                Fluxout = numpy.nan_to_num(Fluxout)
                Fluxoutt.append(Fluxout)

                    # Total gamma luminosity (erg s^-1)
                Fluxtot1 = numpy.sum([Lumsb1, Fluxshell1, Fluxout1])
                Fluxtot2 = numpy.sum([Lumsb2, Fluxshell2, Fluxout2])
                Fluxtot = numpy.sum([Lumsb, Fluxshell, Fluxout])

                Fluxtott1.append(Fluxtot1)
                Fluxtott2.append(Fluxtot2)
                Fluxtott.append(Fluxtot)

                    # Unit
            Flux_unit = lum_gamma_unit * spectrum.unit

                    # Plot
            log_plot(figure_number, 4, t, [Fluxsbt1, Fluxshellt1, Fluxoutt1, Fluxtott1], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Flux_unit.unit.to_string('latex_inline')), ['+', '+', '+', '-'])
            plt.savefig(pathfigure+'Gamma_luminosity_range1.eps')
            figure_number += 1
            log_plot(figure_number, 4, t, [Fluxsbt2, Fluxshellt2, Fluxoutt2, Fluxtott2], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Flux_unit.unit.to_string('latex_inline')), ['+', '+', '+', '-'])
            plt.savefig(pathfigure+'Gamma_luminosity_range2.eps')
            figure_number += 1
            log_plot(figure_number, 4, t, [Fluxsbt, Fluxshellt, Fluxoutt, Fluxtott], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Flux_unit.unit.to_string('latex_inline')), ['+', '+', '+', '-'])
            plt.savefig(pathfigure+'Gamma_luminosity.eps')
            figure_number += 1

                ##===========================================##
                # Energy released by the gamma emission (erg) #
                ##===========================================##

                    # Initialization
            ts = t * yr2s
            Energysb = numpy.zeros(n-1)
            Energyshell = numpy.zeros(n-1)
            Energyout = numpy.zeros(n-1)
            Energytot = numpy.zeros(n-1)

                    # Computation
            for i in range (1,n):
                Energysb[i-1] = integrate.trapz(Fluxsbt[0:i], ts[0:i])
                Energyshell[i-1] = integrate.trapz(Fluxshellt[0:i], ts[0:i])
                Energyout[i-1] = integrate.trapz(Fluxoutt[0:i], ts[0:i])
                Energytot[i-1] = integrate.trapz(Fluxtott[0:i], ts[0:i])

                    # Plot
            log_plot(figure_number, 4, t[1:], [Energysb, Energyshell, Energyout, Energytot], ['SB', 'Shell', 'Out', 'Total'], 'Energy released by the gamma emission', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), 'E (erg)', ['+', '+', '+', '-'])
            plt.savefig(pathfigure+'Energy_released.eps')
            figure_number += 1

                    ##-------------------------------------------------------------------------------------------------##
                    # VERIFICATION OF THE RATIO OF EMITTED GAMMA RAY OVER THE TOTAL ENERGY OF THE ACCELERATED PARTICLES #
                    ##-------------------------------------------------------------------------------------------------##

                        # Computation
            Esnr = eta * Esn        # in erg
            Eratio = Energytot/Esnr * 100

                        # Plot
            log_plot(figure_number, 1, t[1:], Eratio, 'none', 'Percentage of emitted energy of the particles', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$E_\gamma/E_{SN,released}$ (%)', '+')
            plt.savefig(pathfigure+'Ratio.eps')
            plt.show()

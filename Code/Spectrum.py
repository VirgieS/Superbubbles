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

    # At IRAP
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/n_SN/')
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/n_SN/Gamma_emission/'

    # At home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/n_SN')
#pathfigure = '/home/vivi/Documents/Master_2/Superbubbles/figures/n_SN/Gamma_emission/'


##-----------##
# Computation #
##-----------##
#pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/n_SN/Gamma_emission/'

with open('time', 'rb') as time_load:
    with open('distance', 'rb') as distance_load:
        with open('spectra', 'rb') as spectra_read:

                # Loading data
                    # time
            t0 = pickle.load(time_load)
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
            sed_PD_sn = pickle.load(spectra_read)
            sed_PD_unit = pickle.load(spectra_read)

                    # For recording
                        # Within the SB
            Fluxsbsn1 = []       # 100 MeV to 100 GeV
            Fluxsbsn2 = []       # 100 GeV to 100 TeV
            Fluxsbsn = []        # 100 MeV to 100 TeV

                        # Within the supershell
            Fluxshellsn1 = []    # 100 MeV to 100 GeV
            Fluxshellsn2 = []    # 100 GeV to 100 TeV
            Fluxshellsn = []     # 100 MeV to 100 TeV

                        # Outside the SB
            Fluxoutsn1 = []      # 100 MeV to 100 GeV
            Fluxoutsn2 = []      # 100 GeV to 100 TeV
            Fluxoutsn = []       # 100 MeV to 100 TeV

                        # Total
            Fluxtotsn1 = []      # 100 MeV to 100 GeV
            Fluxtotsn2 = []      # 100 GeV to 100 TeV
            Fluxtotsn = []       # 100 MeV to 100 TeV

                # Initialization
            figure_number = 1
            nt0 = len(t0)
            Esep = 100          # The separation of the two spectral zones (GeV)

                ##================##
                # Gamma luminosity #
                ##================##

            ind1 = numpy.where(spectrum_energy < Esep)[0]         # for the first range
            ind2 = numpy.where(spectrum_energy >= Esep)[0]        # for the second range

            for k in range (nt0):           # for each SN explosion

                    # Initialization
                indsn = numpy.where(t >= t0[k])     # only time where the SN already explodes
                nt = len(t[indsn])

                r = rsb[k]
                lum_gamma = sed_PD_sn[k]/spectrum_energy
                lum_gamma_unit = sed_PD_unit/spectrum_unit
                print(lum_gamma_unit)

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

                for i in range (nt):         # for each time step

                        ##------------------------------------##
                        # First zone: in the cavity (erg s^-1) #
                        ##------------------------------------##
                    Fluxsb1 = []
                    Fluxsb2 = []
                    Fluxsb = []
                    m = len(r[i])
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

                    # Total for each SN explosions

                        # In the SB
                Fluxsbsn1.append(Fluxsbt1)
                Fluxsbsn2.append(Fluxsbt2)
                Fluxsbsn.append(Fluxsbt)

                        # In the shell
                Fluxshellsn1.append(Fluxshellt1)
                Fluxshellsn2.append(Fluxshellt2)
                Fluxshellsn.append(Fluxshellt)

                        # Outside the SB
                Fluxoutsn1.append(Fluxoutt1)
                Fluxoutsn2.append(Fluxoutt2)
                Fluxoutsn.append(Fluxoutt)

                        # Sum of each contribution
                Fluxtotsn1.append(Fluxtott1)
                Fluxtotsn2.append(Fluxtott2)
                Fluxtotsn.append(Fluxtott)

                    # Unit
                Flux_unit = lum_gamma_unit * spectrum.unit
                print(Flux_unit)

                    # Plot
                log_plot(figure_number, 4, t[indsn], [Fluxsbt1, Fluxshellt1, Fluxoutt1, Fluxtott1], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[k], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                plt.savefig(pathfigure+'Gamma_luminosity_range1_t0-%d.eps' %k)
                figure_number += 1
                log_plot(figure_number, 4, t[indsn], [Fluxsbt2, Fluxshellt2, Fluxoutt2, Fluxtott2], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[k], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                plt.savefig(pathfigure+'Gamma_luminosity_range2_t0-%d.eps' %k)
                figure_number += 1
                log_plot(figure_number, 4, t[indsn], [Fluxsbt, Fluxshellt, Fluxoutt, Fluxtott], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[k], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                plt.savefig(pathfigure+'Gamma_luminosity_t0-%d.eps' %k)
                figure_number += 1

                    ##===========================================##
                    # Energy released by the gamma emission (erg) #
                    ##===========================================##

                        # Initialization
                ts = t[indsn] * yr2s
                tsn = t[indsn]
                Energysb = numpy.zeros(nt-1)
                Energyshell = numpy.zeros(nt-1)
                Energyout = numpy.zeros(nt-1)
                Energytot = numpy.zeros(nt-1)

                        # Computation
                for i in range (1,nt):
                    Energysb[i-1] = integrate.trapz(Fluxsbt[0:i], ts[0:i])
                    Energyshell[i-1] = integrate.trapz(Fluxshellt[0:i], ts[0:i])
                    Energyout[i-1] = integrate.trapz(Fluxoutt[0:i], ts[0:i])
                    Energytot[i-1] = integrate.trapz(Fluxtott[0:i], ts[0:i])

                        # Plot
                log_plot(figure_number, 4, tsn[1:], [Energysb, Energyshell, Energyout, Energytot], ['SB', 'Shell', 'Out', 'Total'], 'Energy released by the gamma emission (t_{0} = %.2e yr)'%t0[k], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), 'E (erg)', ['+', '+', '+', '-'])
                plt.savefig(pathfigure+'Energy_released_t0-%d.eps'%k)
                figure_number += 1

                        ##=================================================================================================##
                        # VERIFICATION OF THE RATIO OF EMITTED GAMMA RAY OVER THE TOTAL ENERGY OF THE ACCELERATED PARTICLES #
                        ##=================================================================================================##

                            # Computation
                Esnr = eta * Esn        # in erg
                Eratio = Energytot/Esnr * 100

                            # Plot
                log_plot(figure_number, 1, tsn[1:], Eratio, 'none', 'Percentage of emitted energy of the particles ($t_{0}$ = %.2e yr)'%t0[k], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$E_\gamma/E_{SN,released}$ (%)', '+')
                plt.savefig(pathfigure+'Ratio_t0-%d.eps'%k)
                figure_number += 1

            print(numpy.asarray(Fluxsbsn).shape)


# Plot for all SN explosions

nt = len(t)

Fluxtot1_sn_tot = numpy.zeros(nt)
Fluxtot2_sn_tot = numpy.zeros(nt)
Fluxtot_sn_tot = numpy.zeros(nt)
print(nt)

for i in range (nt0):

    indsn = numpy.where(t >= t0[i])

        # From 100 MeV to 100 GeV
    Fluxsb1_sn = Fluxsbsn1[i]
    Fluxshell1_sn = Fluxshellsn1[i]
    Fluxout1_sn = Fluxoutsn1[i]
    Fluxtot1_sn = Fluxtotsn1[i]
    k = 0

    for j in range (nt):

        if t[j] < t0[i]:
            continue
        else:
            Fluxtot1_sn_tot[j] += Fluxtot1_sn[k]
            k += 1

        # From 100 GeV to 100 TeV
    Fluxsb2_sn = Fluxsbsn2[i]
    Fluxshell2_sn = Fluxshellsn2[i]
    Fluxout2_sn = Fluxoutsn2[i]
    Fluxtot2_sn = Fluxtotsn2[i]
    k = 0

    for j in range (nt):

        if t[j] < t0[i]:
            continue
        else:
            Fluxtot2_sn_tot[j] += Fluxtot2_sn[k]
            k += 1

        # From 100 MeV to 100 TeV
    Fluxsb_sn = Fluxsbsn[i]
    Fluxshell_sn = Fluxshellsn[i]
    Fluxout_sn = Fluxoutsn[i]
    Fluxtot_sn = Fluxtotsn[i]
    k = 0

    for j in range (nt):

        if t[j] < t0[i]:
            continue
        else:
            Fluxtot_sn_tot[j] += Fluxtot_sn[k]
            k += 1

    print(numpy.asarray(Fluxtot_sn_tot).shape)
    sym = ['-', '--', ':', '-']

        # Plot
    log_plot(figure_number, 4, t[indsn], [Fluxsb1_sn, Fluxshell1_sn, Fluxout1_sn, Fluxtot1_sn], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), sym)
    log_plot(figure_number + 1, 4, t[indsn], [Fluxsb2_sn, Fluxshell2_sn, Fluxout2_sn, Fluxtot2_sn], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), sym)
    log_plot(figure_number + 2, 4, t[indsn], [Fluxsb_sn, Fluxshell_sn, Fluxout_sn, Fluxtot_sn], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), sym)

log_plot(figure_number, 1, t, Fluxtot1_sn_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), '-')
plt.savefig(pathfigure+'Gamma_luminosity_all_range1.eps')
log_plot(figure_number + 1, 1, t, Fluxtot2_sn_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), '-')
plt.savefig(pathfigure+'Gamma_luminosity_all_range2.eps')
log_plot(figure_number + 2, 1, t, Fluxtot_sn_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Flux_unit.to_string('latex_inline')), '-')
plt.savefig(pathfigure+'Gamma_luminosity_all.eps')

plt.show()

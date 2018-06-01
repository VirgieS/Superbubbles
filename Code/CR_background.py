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

## NEED TO WRITE CLEARLY WHAT I DO

    # IRAP
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/30_Dor_C/')
pathfigure_gamma = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Bons/background/'
print(Nob)

    # Home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/Parameters/stars/10/5')

## =========== ##
# CR background #
## =========== ##

    # Fix time array (yr)
tmin = 3/yr26yr     # yr
tmax = 10/yr26yr    # yr
number_bin_t = 3000
t_fix = numpy.linspace(tmin, tmax, number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr
t7 = t6 * s6yr27yr

    # Correction factor

need_correction = False

if need_correction:

    t_end_6 = 4.0           # Myr
    Rsb = 47.0                          # observed radius (pc)                  #you need to change it for your simulations
    Rw = radius_velocity_SB(t_end_6)[0] # from Weaver's model (pc and km/s)
    correction_factor = Rsb/Rw

else:
    correction_factor = 1

    # Parameters of the SB

Rsb, Vsb = radius_velocity_SB(t6)   # radius of the SB (pc)
correction_factor = 1
Rsb = correction_factor * Rsb       # correction of the radius
Msb, Mswept = masses(t7, Rsb)       # swept-up and inner masses (solar masses)
Ms = Mswept - Msb                   # mass in the shell (solar masses)
ns, hs = density_thickness_shell_percentage(percentage, Rsb, Mswept, Msb)  # thickness (pc) and density (cm^-3) of the shell
#ns, hs = density_thickness_shell(Vsb, Mswept, Msb, Rsb)                     # thickness (pc) and density (cm^-3) of the shell
Vs = (Ms*Msun2g)/(ns * mu * mpg)

with open('CRbackground', 'wb') as CR_write:

        # CR particles density (cm^-3 GeV^-1)
    n_CRb = 0.3 * cosmicray_lis(ECR)

        # Gamma luminosity (erg s^-1)
    Lum_CRb = numpy.zeros(number_bin_t)
    Lum_HESS_CRb = numpy.zeros(number_bin_t)
    Lum_Fermi_CRb = numpy.zeros(number_bin_t)

    for i in range (number_bin_t):

                # Density of gas (cm^-3)
        ngas = ns[i] * 1/units.cm**3

                # Particles distribution (GeV^-1)
        N_CRb = n_CRb * Vs[i] * 1/units.GeV


        if numpy.asarray(N_CRb).all() == 0:
            continue

            # For all the range of energy (100 MeV to 100 TeV)
                # intrisic differential luminosity (eV^-1 s^-1)
        model = TableModel(E_CR, N_CRb, amplitude = 1)
        PD = PionDecay(model, nh = ngas, nuclear_enhancement = True, useLUT = False)
        flux_PD = PD.flux(spectrum_energy, distance = 0 * units.pc)
        flux_PD = numpy.nan_to_num(numpy.asarray(flux_PD))

                # Gamma luminosity (whole energy range) (erg s^-1)
        spectrum_erg = spectrum * 1.0/erg2GeV     # erg
        spectrum_ev = spectrum_erg * 1.0/eV2erg   # eV
        lum_energy = flux_PD * spectrum_erg          # erg s^-1 eV^-1
        Lum_CRb[i] = luminosity(lum_energy, spectrum_ev)

            # H.E.S.S energy range
        Emin = 1 * TeV2GeV                  # 1 TeV (GeV)
        Emax = 10 * TeV2GeV                 # 10 TeV (GeV)
        indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

                # Gamma luminosity
        spectrum_HESS = spectrum[indE]
        spectrum_erg = spectrum_HESS * 1.0/erg2GeV     # only in the energy range (erg)
        spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
        lum_HESS = flux_PD[indE] * spectrum_erg   # erg s^-1 eV^-1
        Lum_HESS_CRb[i] = luminosity(lum_HESS, spectrum_ev) # erg s^-1

            # H.E.S.S energy range
        Emin = 100 * MeV2GeV                # 100 MeV (GeV)
        Emax = 100                          # 100 GeV (GeV)
        indE = numpy.where((spectrum >= Emin) & (spectrum <= Emax))[0]

                # Gamma luminosity
        spectrum_Fermi = spectrum[indE]
        spectrum_erg = spectrum_Fermi * 1.0/erg2GeV     # only in the energy range (erg)
        spectrum_ev = spectrum_erg * 1.0/eV2erg         # eV
        lum_Fermi = flux_PD[indE] * spectrum_erg   # erg s^-1 eV^-1
        Lum_Fermi_CRb[i] = luminosity(lum_Fermi, spectrum_ev) # erg s^-1

        print('t = %.2f Myrs done'%t6[i])

    pickle.dump(Lum_CRb, CR_write)
    pickle.dump(Lum_HESS_CRb, CR_write)
    pickle.dump(Lum_Fermi_CRb, CR_write)
print(Lum_CRb)
print(Lum_HESS_CRb)
print(Lum_Fermi_CRb)
print(Nob)

    # plot

        # whole energy range
figure_number = 1
number_of_plot = 1
x = t6
y = Lum_CRb
label_name = 'none'
title = ''
xlabel = 'Time [Myrs]'
ylabel = '$L_\gamma$ (100 MeV - 100 Tev) [erg s$^{-1}$]'
symbol = ''
linestyle = '-'
color = 'red'
text = ''
xmin = numpy.min(x) - 0.5
xmax = numpy.max(x) + 0.5
ymin = 10**numpy.min(numpy.log10(y) - 3)
ymax = 10**numpy.max(numpy.log10(y) + 2)

semilog_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Background_all.pdf')
figure_number += 1

        # Fermi energy range
y = Lum_Fermi_CRb
ylabel = '$L_\gamma$ (100 MeV - 100 Gev) [erg s$^{-1}$]'
ymin = 10**numpy.min(numpy.log10(y) - 3)
ymax = 10**numpy.max(numpy.log10(y) + 2)

semilog_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Background_Fermi.pdf')
figure_number += 1

        # H.E.S.S. energy range
y = Lum_HESS_CRb
ylabel = '$L_\gamma$ (1 MeV - 10 Tev) [erg s$^{-1}$]'
ymin = 10**numpy.min(numpy.log10(y) - 3)
ymax = 10**numpy.max(numpy.log10(y) + 2)

semilog_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, color, text, xmin, xmax, ymin, ymax)
plt.savefig(pathfigure_gamma+'Background_HESS.pdf')
figure_number += 1

plt.show()

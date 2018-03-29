# In this program, we will define all the important functions to compute the gamma emission of a superbubble

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

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_SB import *

##---------##
# Functions #
##---------##

def data(Emin_CR, Emax_CR, p0, alpha, D0, delta, zones):

    """
    Return 5 files with differents quantities.
    Inputs:
        Emin_CR     :       minimum kinetic energy of the relativistic particles (GeV)
        Emax_CR     :       maximum kinetic energy of the relativistic particles (GeV)
        p0          :       normalization constant of the power-law distribution of the cosmic rays (GeV/c)
        alpha       :       exponent of the power-law distribution of the cosmic rays
        D0          :       diffusion coefficient at p0 (cm^2 s^-1)
        delta       :       exponent of the power-law distribution of the diffusion coefficient
        zones       :       which zone do you want to compute (1: cavity of the SB, 2: supershell and 3: outside)

    Outputs:
        gas         :       file with the density of gas for all time step and radius step (cm^-3)
        energy      :       file with the energy of the relativistic particles (GeV)
        time        :       file with the time array from the first SN explosion to 1 millions years after the last SN explosion (yr)
        radius      :       file with the radius array for each time step (pc)
        data        :       file with the distribution of relativistic particles for each time step, each SN explosion and each zones (GeV^-1)
    """

    # Parameters for the system

        # luminosity
    L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
    L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s

        # in the ISM
    pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
    C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

    with open('SN', 'rb') as t0_load:
        with open('time', 'rb') as time_load:
            with open('data', 'wb') as data_write:                          # for the total gamma luminosity

                                # Energy
                            number_bin_E = 10
                            ECR = numpy.logspace(numpy.log10(Emin_CR), numpy.log10(Emax_CR), number_bin_E)  # GeV
                            ECR_unit = units.GeV      # recording the unit of the energy array

                            pickle.dump(ECR, energy_write)
                            pickle.dump(ECR_unit, energy_write)

                                ## ====================================================== ##
                                # power-law distribution of the cosmic rays (GeV^-1 cm^-3) #
                                ## ====================================================== ##
                                    # N(p) = N0 * (p/p0)^(-alpha)
                                    # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)
                            Esng = Esn*erg2GeV  # in GeV
                            N_E = power_law_distribution(Emin_CR, Emax_CR, ECR, alpha, eta, Esng, p0)

                                ## =============================== ##
                                # diffusion coefficient (cm^2 s^-1) #
                                ## =============================== ##
                                    # D(p) = D0 * (p/p0)^(-delta)
                                    # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
                            D = diffusion_coefficient(p0, D0, ECR, delta)

                                ## ====================================================== ##
                                # Computation of N(E,t,r), particles distribution (GeV^-1) #
                                ## ====================================================== ##

                                    # SN explosion: time (yr)
                            t0 = pickle.load(t0_load)
                            nt0 = len(t0)

                                    # recording
                            sed_PD_sn = []              # spectral energy distribution for each SN, each time, each zone and each energy: len(t0)xlen(t)x(len(r)+2)xlen(E)

                            for i in range (nt0):

                                tmin = t0[i]    # only when the SN occurs
                                t6 = tmin * yr26yr      # 10^6 yr
                                Rsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)[0]           # radius and velocity of the SB

                                if i == 0:

                                    dt = Rsb**2/(6*D[len(D)-1])

                                tau = Rsb**2/(6*D[len(D)-1])
                                tmax = tmin + 3*tau
                                number_bin_t = (tmax - tmin)/dt
                                t = numpy.linspace(tmin, tmax, number_bin_t)

                                    # Initialization
                                figure_number = 1

                                    # Recording
                                sed = []

                                for j in range (number_bin_t):        # for each time step


                                        # Initialization
                                    t6 = t[j] * yr26yr      # 10^6 yr
                                    t7 = t6 * s6yr27yr      # 10^7 yr
                                    delta_t = t[j] - t0[i]
                                    sed_r = []

                                        # Parameters of the SB
                                    Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)           # radius and velocity of the SB
                                    Msb, Mswept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu)   # swept-up and inner masses (solar masses)
                                    hs, ns = density_thickness_shell(mu, n0, Vsb, C02, Mswept, Msb, Rsb)            # thickness and density of the shell (pc)

                                    for zone in (zones):

                                        if zone == 1:                           # in the SB

                                                # distance array (pc)
                                            rmin = 0.01                 # minimum radius (pc)
                                            rmax = Rsb-hs               # maximum radius (pc)
                                            number_bin_r = 15           # number of bin for r from 0 to Rsb-hs
                                            r = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_r)    # position (pc)

                                                # Density of gas (cm^-3)
                                            nsb = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, r, Rsb)[1]
                                            ngas = nsb * 1/units.cm**3

                                                # Particles distribution (GeV^-1)
                                            r_in = 0
                                            r_out = r[0]
                                            N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                                                # Spectral distribution
                                            model = TableModel(ECR, N_part, amplitude = 1)
                                            PD = PionDecay(model, nh = ngas[0], nuclear_enhancement = True)
                                            sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                            sed_r.append(sed_PD)

                                            for k in range (1, number_bin_r):   # for each radius inside the SB

                                                    # Particles distribution (GeV^-1)
                                                r_in = r_out
                                                r_out = r[k]
                                                N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                                                    # Spectral distribution (erg s^-1)
                                                model = TableModel(ECR, N_part, amplitude = 1)
                                                PD = PionDecay(model, nh = ngas[k], nuclear_enhancement = True)
                                                sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                                sed_r.append(sed_PD)

                                        elif zone == 2:                         # in the supershell

                                                # Density of gas
                                            ngas = ns * 1/units.cm**3

                                                # Particles distribution (GeV^-1)
                                            r_in = Rsb - hs     # in pc
                                            r_out = Rsb         # in pc
                                            N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                                                # Spectral disbution (erg s^-1)
                                            model = TableModel(ECR, N_part, amplitude = 1)
                                            PD = PionDecay(model, nh = ngas, nuclear_enhancement = True)
                                            sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                            sed_r.append(sed_PD)

                                        else:                                   # outside the SB

                                                # Density of gas
                                            ngas = n0 * 1/units.cm**3

                                                # Distribution of particles (GeV^-1)
                                            N_part = inf_particles(Rsb, N_E, D, delta_t) * 1/units.GeV

                                                 # Spectral disbution (erg s^-1)
                                            model = TableModel(ECR, N_part, amplitude = 1)
                                            PD = PionDecay(model, nh = ngas, nuclear_enhancement = True)
                                            sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                            sed_r.append(sed_PD)

                                    sed.append(sed_r)

                                sed_PD_sn.append(sed)



                            Ntotsn = numpy.asarray(Ntotsn)
                            Ntot_unit = 1/(units.GeV)               # units of Ntot (GeV^-1)
                            pickle.dump(Ntotsn, data_write)         # Recording the number of particles for each SN explosion time, each time, each zone and each energy
                            pickle.dump(Ntot_unit, data_write)

                            ngastot = numpy.asarray(ngastot)
                            ngas_unit = 1/units.cm**3               # units of ngas (cm^-3)
                            pickle.dump(ngastot, gas_write)         # Recording the density of gas for each SN explosion, each time and each zone
                            pickle.dump(ngas_unit, gas_write)

    return t

def spectrum(Emin_gamma, Emax_gamma, zones):

    """
    Return one file

    Inputs:
        Emin_gamma  :       minimum energy of the gamma photon (GeV)
        Emax_gamma  :       maximum energy of the gamma photon (GeV)
        zones       :       which zone do you want to compute (1: cavity of the SB, 2: supershell and 3: outside)

    Output:
        spectra     :       file with the energy spectrum of the gamma photons and the spectral energy distribution for each SN explosions, time step and zones
    """

    with open('spectra', 'wb') as spectra_write:
        with open('gas', 'rb') as gas_load:
            with open('energy', 'rb') as energy_load:
                with open('data', 'rb') as data_load:
                    with open('radius', 'rb') as radius_load:
                        with open('time', 'rb') as time_load:
                            with open('SN', 'rb') as t0_load:

                                    # Loading of data: ngas, energy, nCR, radius and time

                                ngassn = pickle.load(gas_load)
                                ngas_unit = pickle.load(gas_load)

                                ECR = pickle.load(energy_load)
                                ECR_unit = pickle.load(energy_load)
                                ECR = ECR * ECR_unit

                                Ntotsn = pickle.load(data_load)
                                Ntot_unit = pickle.load(data_load)

                                radius = pickle.load(radius_load)

                                t0 = pickle.load(t0_load)
                                t = pickle.load(time_load)
                                t_unit = pickle.load(time_load)

                                    # Energy spectrum (GeV)
                                number_bin_E = 20
                                spectrum_energy = numpy.logspace(numpy.log10(Emin_gamma), numpy.log10(Emax_gamma), number_bin_E) * units.GeV

                                    # Initialization
                                nt0 = len(t0)

                                    # recording
                                sed_PD_sn = []              # spectral energy distribution for each SN, each time, each zone and each energy: len(t0)xlen(t)x(len(r)+2)xlen(E)

                                for i in range (nt0):       # for each SN explosion

                                    sed = []
                                    indSN = numpy.where(t > t0[i])
                                    nt = len(numpy.asarray(t)[indSN])

                                    for j in range (nt):    # for each time step when the SN explosion occurs or after the SN explosion

                                            # Important quantities
                                        ngas = ngassn[j] * ngas_unit
                                        ntot = numpy.asarray(Ntotsn[i])[j] * Ntot_unit

                                            # Recording
                                        spectrum_r = []
                                        sed_r = []
                                        n = 0

                                        for zone in (zones):

                                            if zone == 1:                       # in the SB

                                                    # Initialization
                                                r = radius[j]
                                                nr = len(r)

                                                for k in range (nr):     # for each radius step
                                                    model = TableModel(ECR, ntot[k], amplitude = 1)
                                                    PD = PionDecay(model, nh = ngas[k], nuclear_enhancement = True)
                                                    sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                                    n += 1

                                                sed_r.append(sed_PD)

                                            elif zone == 2:                     # in the supershell

                                                model = TableModel(ECR, ntot[n], amplitude = 1)
                                                PD = PionDecay(model, nh = ngas[n], nuclear_enhancement = True)
                                                sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                                sed_r.append(sed_PD)
                                                n += 1

                                            else:                               # outside the SB
                                                model = TableModel(ECR, ntot[n], amplitude = 1)
                                                PD = PionDecay(model, nh = ngas[n], nuclear_enhancement = True)
                                                sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)

                                        sed.append(sed_r)

                                    sed_PD_sn.append(sed)

                                spectrum = numpy.asarray(spectrum_energy)
                                spectrum_unit = spectrum_energy.unit
                                sed_PD_sn = numpy.asarray(sed_PD_sn)
                                sed_unit = sed_PD.unit
                                pickle.dump(spectrum, spectra_write)
                                pickle.dump(spectrum_unit, spectra_write)
                                pickle.dump(sed_PD_sn, spectra_write)
                                pickle.dump(sed_unit, spectra_write)
    return

def luminosity(lum_energy, energy):

    """
    Return the luminosity in a range of energy

    Inputs:
        lum_energy  :       intrinsic luminosity array for a range of energy (erg/s/GeV)
        energy      :       range of energy in which we will compute the luminosity (GeV)
    """
    lum = integrate.trapz(lum_energy, energy)
    lum = numpy.nan_to_num(lum)

    return lum

def gamma_luminosity(Esep, pathfigure):

    """
    Return the gamma luminosity in the ranges of energy

    Inputs:
        Esep        :       limit of the two ranges of energy (GeV)
        pathfigure  :       path to save figures
    """
    with open('SN', 'rb') as t0_load:
        with open('time', 'rb') as time_load:
            with open('radius', 'rb') as radius_load:
                with open('spectra', 'rb') as spectra_read:

                        # Loading data
                            # time
                    t0 = pickle.load(t0_load)
                    t = pickle.load(time_load)
                    t_unit = pickle.load(time_load)
                    time = t * t_unit

                            # radius
                    r = pickle.load(radius_load)

                            # Energy (GeV)
                    spectrum = pickle.load(spectra_read)
                    spectrum_unit = pickle.load(spectra_read)
                    spectrum_energy = spectrum * spectrum_unit

                             # spectal energy distribution (erg s^-1): E^2 dN/dE = E L_gamma
                    sed_PD_sn = pickle.load(spectra_read)
                    sed_PD_unit = pickle.load(spectra_read)

                    #spectrum_energy = spectrum * spectrum_unit
                            # For recording
                                # Within the SB
                    Lumsb_sn_1 = []         # 100 MeV to 100 GeV
                    Lumsb_sn_2 = []         # 100 GeV to 100 TeV
                    Lumsb_sn = []           # 100 MeV to 100 TeV

                                # Within the supershell
                    Lumshell_sn_1 = []      # 100 MeV to 100 GeV
                    Lumshell_sn_2 = []      # 100 GeV to 100 TeV
                    Lumshell_sn = []        # 100 MeV to 100 TeV

                                # Outside the SB
                    Lumout_sn_1 = []       # 100 MeV to 100 GeV
                    Lumout_sn_2 = []       # 100 GeV to 100 TeV
                    Lumout_sn = []         # 100 MeV to 100 TeV

                                # Total
                    Lumtot_sn_1 = []       # 100 MeV to 100 GeV
                    Lumtot_sn_2 = []       # 100 GeV to 100 TeV
                    Lumtot_sn = []         # 100 MeV to 100 TeV

                        # Initialization
                    figure_number = 1
                    nt0 = len(t0)

                        ##================##
                        # Gamma luminosity #
                        ##================##

                    ind1 = numpy.where(spectrum < Esep)[0]         # for the first range
                    ind2 = numpy.where(spectrum >= Esep)[0]        # for the second range

                    for i in range (nt0):           # for each SN explosion

                        lum_energy = sed_PD_sn[i]/spectrum
                        lum_energy_unit = sed_PD_unit/spectrum_unit

                            # For recording
                                # Within the SB
                        Lumsb_t_1 = []          # 100 MeV to 100 GeV
                        Lumsb_t_2 = []          # 100 GeV to 100 TeV
                        Lumsb_t = []            # 100 MeV to 100 TeV

                                # Within the supershell
                        Lumshell_t_1 = []       # 100 MeV to 100 GeV
                        Lumshell_t_2 = []       # 100 GeV to 100 TeV
                        Lumshell_t = []         # 100 MeV to 100 TeV

                                # Outside the SB
                        Lumout_t_1 = []         # 100 MeV to 100 GeV
                        Lumout_t_2 = []         # 100 GeV to 100 TeV
                        Lumout_t = []           # 100 MeV to 100 TeV

                                # Total
                        Lumtot_t_1 = []         # 100 MeV to 100 GeV
                        Lumtot_t_2 = []         # 100 GeV to 100 TeV
                        Lumtot_t = []           # 100 MeV to 100 TeV

                            # Initialization
                        indSN = numpy.where(t > t0[i])[0]     # only time where the SN already explodes
                        nt = len(numpy.asarray(t)[indSN])

                        for j in range (nt):         # for each time step

                                ##------------------------------------##
                                # First zone: in the cavity (erg s^-1) #
                                ##------------------------------------##
                                    # Initialization
                            radius = r[j]
                            nr = len(radius)

                                    # Recording
                            Lumsb_1_r = []
                            Lumsb_2_r = []
                            Lumsb_r = []

                            for k in range (nr):     # for each radius step

                                    # From 100 MeV to 100 GeV
                                Lum_1 = luminosity(lum_energy[j, k, ind1], spectrum[ind1])
                                Lumsb_1_r.append(Lum_1)

                                    # From 100 GeV to 100 TeV
                                Lum_2 = luminosity(lum_energy[j, k, ind2], spectrum[ind2])
                                Lumsb_2_r.append(Lum_2)

                                    # From 100 MeV to 100 TeV
                                Lum = luminosity(lum_energy[j,k], spectrum)
                                Lumsb_r.append(Lum)

                                    # Total contribution in the SB (erg s^-1)
                            Lumsb_1 = numpy.sum(Lumsb_1_r)
                            Lumsb_2 = numpy.sum(Lumsb_2_r)
                            Lumsb = numpy.sum(Lumsb_r)

                            Lumsb_t_1.append(Lumsb_1)
                            Lumsb_t_2.append(Lumsb_2)
                            Lumsb_t.append(Lumsb)

                                ##-----------------------------------------##
                                # Second zone: in the supershell (erg s^-1) #
                                ##-----------------------------------------##

                                    # From 100 MeV to 100 GeV
                            Lumshell_1 = luminosity(lum_energy[j, nr, ind1], spectrum[ind1])
                            Lumshell_t_1.append(Lumshell_1)

                                    # From 100 GeV to 100 TeV
                            Lumshell_2 = luminosity(lum_energy[j, nr, ind2], spectrum[ind2])
                            Lumshell_t_2.append(Lumshell_2)

                                    # From 100 MeV to 100 TeV
                            Lumshell = luminosity(lum_energy[j, nr], spectrum)
                            Lumshell_t.append(Lumshell)

                                ##-------------------------------------##
                                # Third zone: outside the SB (erg s^-1) #
                                ##-------------------------------------##

                                    # From 100 MeV to 100 GeV
                            Lumout_1 = luminosity(lum_energy[j, nr+1, ind1], spectrum[ind1])
                            Lumout_t_1.append(Lumout_1)

                                    # From 100 GeV to 100 TeV
                            Lumout_2 = luminosity(lum_energy[j, nr+1, ind2], spectrum[ind2])
                            Lumout_t_2.append(Lumout_2)

                                    # From 100 MeV to 100 TeV
                            Lumout = luminosity(lum_energy[j, nr+1], spectrum)
                            Lumout_t.append(Lumout)

                                # Total gamma luminosity (erg s^-1)
                            #Lumtot1 = numpy.sum([Lumsb_1, Lumshell_1, Lumout_1])
                            #Lumtot2 = numpy.sum([Lumsb_2, Lumshell_2, Lumout_2])
                            #Lumtot = numpy.sum([Lumsb, Lumshell, Lumout])
                            Lumtot1 = numpy.sum([Lumsb_1, Lumshell_1])
                            Lumtot2 = numpy.sum([Lumsb_2, Lumshell_2])
                            Lumtot = numpy.sum([Lumsb, Lumshell])


                            Lumtot_t_1.append(Lumtot1)
                            Lumtot_t_2.append(Lumtot2)
                            Lumtot_t.append(Lumtot)

                            # Total for each SN explosions

                                # In the SB
                        Lumsb_sn_1.append(Lumsb_t_1)
                        Lumsb_sn_2.append(Lumsb_t_2)
                        Lumsb_sn.append(Lumsb_t)

                                # In the shell
                        Lumshell_sn_1.append(Lumshell_t_1)
                        Lumshell_sn_2.append(Lumshell_t_2)
                        Lumshell_sn.append(Lumshell_t)

                                # Outside the SB
                        Lumout_sn_1.append(Lumout_t_1)
                        Lumout_sn_2.append(Lumout_t_2)
                        Lumout_sn.append(Lumout_t)

                                # Sum of each contribution
                        Lumtot_sn_1.append(Lumtot_t_1)
                        Lumtot_sn_2.append(Lumtot_t_2)
                        Lumtot_sn.append(Lumtot_t)

                            # Unit
                        Lum_unit = lum_energy_unit * spectrum_energy.unit

                            # Plot
                        sym = ['-.', '--', ':', '-']
                        log_plot(figure_number, 4, t[indSN], [Lumsb_t_1, Lumshell_t_1, Lumout_t_1, Lumtot_t_1], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
                        plt.savefig(pathfigure+'Gamma_luminosity_range1_t0-%d.eps' %i)
                        figure_number += 1
                        log_plot(figure_number, 4, t[indSN], [Lumsb_t_2, Lumshell_t_2, Lumout_t_2, Lumtot_t_2], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
                        plt.savefig(pathfigure+'Gamma_luminosity_range2_t0-%d.eps' %i)
                        figure_number += 1
                        log_plot(figure_number, 4, t[indSN], [Lumsb_t, Lumshell_t, Lumout_t, Lumtot_t], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
                        plt.savefig(pathfigure+'Gamma_luminosity_t0-%d.eps' %i)
                        figure_number += 1

        Lumtot_sn_1 = numpy.asarray(Lumtot_sn_1)
        Lumtot_sn_2 = numpy.asarray(Lumtot_sn_2)
        Lumtot_sn = numpy.asarray(Lumtot_sn)

        # Plot for all SN explosions

        Lum_sn_1_tot = numpy.zeros(len(t))
        Lum_sn_2_tot = numpy.zeros(len(t))
        Lum_sn_tot = numpy.zeros(len(t))

        for i in range (nt0):

            indSN = numpy.where(t > t0[i])[0]  # consider only time when there is already the SN explosion
            t_sn = numpy.asarray(t)[indSN]* yr26yr   # in Myr
            t_sn_unit = units.Myr
            nt = len(t_sn)
            Lum_sn_1 = numpy.asarray(Lumtot_sn_1[i])
            Lum_sn_2 = numpy.asarray(Lumtot_sn_2[i])
            Lum_sn = numpy.asarray(Lumtot_sn[i])

            for j in range (nt):

                Lum_sn_1_tot[j] += Lum_sn_1[j]
                Lum_sn_2_tot[j] += Lum_sn_2[j]
                Lum_sn_tot[j] += Lum_sn[j]

            sym = ['-.', '--', ':', '-']

                # Plot
                    # First range
            log_plot(figure_number, 1, t_sn, Lumsb_sn_1[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-.')
            log_plot(figure_number, 1, t_sn, Lumshell_sn_1[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '--')
            log_plot(figure_number, 1, t_sn, Lumout_sn_1[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ':')
            log_plot(figure_number, 1, t_sn, Lumtot_sn_1[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')

                    # Second range
            log_plot(figure_number + 1, 1, t_sn, Lumsb_sn_2[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-.')
            log_plot(figure_number + 1, 1, t_sn, Lumshell_sn_2[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '--')
            log_plot(figure_number + 1, 1, t_sn, Lumout_sn_2[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ':')
            log_plot(figure_number + 1, 1, t_sn, Lumtot_sn_2[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')

                    # All energy
            log_plot(figure_number + 2, 1, t_sn, Lumsb_sn[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-.')
            log_plot(figure_number + 2, 1, t_sn, Lumshell_sn[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '--')
            log_plot(figure_number + 2, 1, t_sn, Lumout_sn[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ':')
            log_plot(figure_number + 2, 1, t_sn, Lumtot_sn[i], 'none', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')

        ind0_1 = numpy.where(Lum_sn_1_tot != 0)[0]
        ind0_2 = numpy.where(Lum_sn_2_tot != 0)[0]
        ind0 = numpy.where(Lum_sn_tot != 0)[0]

        log_plot(figure_number, 1, t[ind0_1]* yr26yr, Lum_sn_1_tot[ind0_1], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all_range1.eps')
        log_plot(figure_number + 1, 1, t[ind0_2] * yr26yr, Lum_sn_2_tot[ind0_2], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all_range2.eps')
        log_plot(figure_number + 2, 1, t[ind0] * yr26yr, Lum_sn_tot[ind0], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all.eps')

        #plt.show()
    return

def gamma_luminosity_it(Esep, pathfigure, zones, iteration, figure_number):

    """
    Return the gamma luminosity in the ranges of energy

    Inputs:
        Esep        :       limit of the two ranges of energy (GeV)
        pathfigure  :       path to save figures
        zones       :       which zone do you want to compute (1: cavity of the SB, 2: supershell and 3: outside)
        iteration   :       number of the iteration
        figure_number:      number of the figure
    """
    with open('SN', 'rb') as t0_load:
        with open('time', 'rb') as time_load:
            with open('radius', 'rb') as radius_load:
                with open('spectra', 'rb') as spectra_read:

                        # Loading data
                            # Sn explosion time (yr)
                    t0 = pickle.load(t0_load)

                            # time array (yr)
                    t = pickle.load(time_load)
                    t_unit = pickle.load(time_load)
                    time = t * t_unit

                            # radius
                    r = pickle.load(radius_load)

                            # Energy of the gamma photon (GeV)
                    spectrum = pickle.load(spectra_read)
                    spectrum_unit = pickle.load(spectra_read)
                    spectrum_energy = spectrum * spectrum_unit

                             # spectal energy distribution (erg s^-1): E^2 dN/dE = E L_gamma
                    sed_PD_sn = pickle.load(spectra_read)
                    sed_PD_unit = pickle.load(spectra_read)


                            # Unit
                    lum_energy_unit = sed_PD_unit/spectrum_unit
                    Lum_unit = lum_energy_unit * spectrum_energy.unit

                        # For recording
                    for zone in (zones):

                        if zone == 1:                                           # in the SB

                            Lumsb_sn_1 = []         # 100 MeV to 100 GeV
                            Lumsb_sn_2 = []         # 100 GeV to 100 TeV
                            Lumsb_sn = []           # 100 MeV to 100 TeV

                        elif zone == 2:                                         # in the supershell

                            Lumshell_sn_1 = []      # 100 MeV to 100 GeV
                            Lumshell_sn_2 = []      # 100 GeV to 100 TeV
                            Lumshell_sn = []        # 100 MeV to 100 TeV

                        else:                                                   # outside the SB

                            Lumout_sn_1 = []       # 100 MeV to 100 GeV
                            Lumout_sn_2 = []       # 100 GeV to 100 TeV
                            Lumout_sn = []         # 100 MeV to 100 TeV

                                # Total
                    Lumtot_sn_1 = []       # 100 MeV to 100 GeV
                    Lumtot_sn_2 = []       # 100 GeV to 100 TeV
                    Lumtot_sn = []         # 100 MeV to 100 TeV

                        # Initialization
                    nt0 = len(t0)

                        ##================##
                        # Gamma luminosity #
                        ##================##

                    ind1 = numpy.where(spectrum < Esep)[0]         # for the first range
                    ind2 = numpy.where(spectrum >= Esep)[0]        # for the second range

                    for i in range (nt0):           # for each SN explosion

                        lum_energy = sed_PD_sn[i]/spectrum

                            # For each time step
                        for zone in (zones):

                            if zone == 1:                                       # in the SB
                                Lumsb_t_1 = []          # 100 MeV to 100 GeV
                                Lumsb_t_2 = []          # 100 GeV to 100 TeV
                                Lumsb_t = []            # 100 MeV to 100 TeV

                            elif zone == 2:                                     # in the supershell

                                Lumshell_t_1 = []       # 100 MeV to 100 GeV
                                Lumshell_t_2 = []       # 100 GeV to 100 TeV
                                Lumshell_t = []         # 100 MeV to 100 TeV

                            else:                                               # outside the SB

                                Lumout_t_1 = []         # 100 MeV to 100 GeV
                                Lumout_t_2 = []         # 100 GeV to 100 TeV
                                Lumout_t = []           # 100 MeV to 100 TeV

                                # Total
                        Lumtot_t_1 = []         # 100 MeV to 100 GeV
                        Lumtot_t_2 = []         # 100 GeV to 100 TeV
                        Lumtot_t = []           # 100 MeV to 100 TeV

                            # Initialization
                        indSN = numpy.where(t > t0[i])[0]     # only time where the SN already explodes
                        t_sn = numpy.asarray(t)[indSN] * yr26yr
                        t_sn_unit = units.Myr
                        nt = len(t_sn)

                        for j in range (nt):         # for each time step

                            n = 0
                            SB = False

                            for zone in (zones):

                                if zone == 1:                                   # in the SB

                                            # Initialization
                                    radius = r[j]
                                    nr = len(radius)

                                            # Recording
                                    Lumsb_1_r = []
                                    Lumsb_2_r = []
                                    Lumsb_r = []

                                    for k in range (nr):     # for each radius step

                                            # From 100 MeV to 100 GeV
                                        Lum_1 = luminosity(lum_energy[j, k, ind1], spectrum[ind1])
                                        Lumsb_1_r.append(Lum_1)

                                            # From 100 GeV to 100 TeV
                                        Lum_2 = luminosity(lum_energy[j, k, ind2], spectrum[ind2])
                                        Lumsb_2_r.append(Lum_2)

                                            # From 100 MeV to 100 TeV
                                        Lum = luminosity(lum_energy[j,k], spectrum)
                                        Lumsb_r.append(Lum)

                                            # Total contribution in the SB (erg s^-1)
                                    Lumsb_1 = numpy.sum(Lumsb_1_r)
                                    Lumsb_2 = numpy.sum(Lumsb_2_r)
                                    Lumsb = numpy.sum(Lumsb_r)

                                    Lumsb_t_1.append(Lumsb_1)
                                    Lumsb_t_2.append(Lumsb_2)
                                    Lumsb_t.append(Lumsb)

                                    n = nr
                                    SB = True

                                elif zone == 2:                                 # in the supershell

                                            # From 100 MeV to 100 GeV
                                    Lumshell_1 = luminosity(lum_energy[j, n, ind1], spectrum[ind1])
                                    Lumshell_t_1.append(Lumshell_1)

                                            # From 100 GeV to 100 TeV
                                    Lumshell_2 = luminosity(lum_energy[j, n, ind2], spectrum[ind2])
                                    Lumshell_t_2.append(Lumshell_2)

                                            # From 100 MeV to 100 TeV
                                    Lumshell = luminosity(lum_energy[j, n], spectrum)
                                    Lumshell_t.append(Lumshell)

                                    n += 1

                                else:                                           # outside the SB

                                            # From 100 MeV to 100 GeV
                                    Lumout_1 = luminosity(lum_energy[j, n, ind1], spectrum[ind1])
                                    Lumout_t_1.append(Lumout_1)

                                            # From 100 GeV to 100 TeV
                                    Lumout_2 = luminosity(lum_energy[j, n, ind2], spectrum[ind2])
                                    Lumout_t_2.append(Lumout_2)

                                            # From 100 MeV to 100 TeV
                                    Lumout = luminosity(lum_energy[j, nr+1], spectrum)
                                    Lumout_t.append(Lumout)

                            if SB:
                                    # Total gamma luminosity (erg s^-1)
                                Lumtot1 = numpy.sum([Lumsb_1, Lumshell_1])
                                Lumtot2 = numpy.sum([Lumsb_2, Lumshell_2])
                                Lumtot = numpy.sum([Lumsb, Lumshell])
                            else:
                                Lumtot1 = Lumshell_1
                                Lumtot2 = Lumshell_2
                                Lumtot = Lumshell

                            Lumtot_t_1.append(Lumtot1)
                            Lumtot_t_2.append(Lumtot2)
                            Lumtot_t.append(Lumtot)

                            # Total for each SN explosions + plot
                        figure_1 = figure_number
                        figure_2 = figure_1 + 1
                        figure_tot = figure_2 + 1

                        for zone in (zones):

                            if zone == 1:                                       # in the SB
                                    # recording
                                Lumsb_sn_1.append(Lumsb_t_1)
                                Lumsb_sn_2.append(Lumsb_t_2)
                                Lumsb_sn.append(Lumsb_t)

                                    # plot
                                sym = '--'
                                label = 'SB'

                                        # 100 MeV to 100 GeV
                                log_plot(figure_1, 1, t_sn, Lumsb_t_1, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 GeV to 100 TeV
                                log_plot(figure_2, 1, t_sn, LumSB_t_2, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 MeV to 100 TeV
                                log_plot(figure_tot, 1, t_sn, LumSB_t, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            elif zone == 2:                                     # in the supershell

                                    # recording
                                Lumshell_sn_1.append(Lumshell_t_1)
                                Lumshell_sn_2.append(Lumshell_t_2)
                                Lumshell_sn.append(Lumshell_t)

                                    # plot
                                sym = '--'
                                label = 'Shell'

                                        # 100 MeV to 100 GeV
                                log_plot(figure_1, 1, t_sn, Lumshell_t_1, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 GeV to 100 TeV
                                log_plot(figure_2, 1, t_sn, Lumshell_t_2, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 MeV to 100 TeV
                                log_plot(figure_tot, 1, t_sn, Lumshell_t, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            else:                                               # outside the SB

                                    # recording
                                Lumout_sn_1.append(Lumout_t_1)
                                Lumout_sn_2.append(Lumout_t_2)
                                Lumout_sn.append(Lumout_t)

                                    # plot
                                sym = ':'
                                label = 'Out'

                                        # 100 MeV to 100 GeV
                                log_plot(figure_1, 1, t_sn, Lumout_t_1, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 GeV to 100 TeV
                                log_plot(figure_2, 1, t_sn, Lumout_t_2, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 MeV to 100 TeV
                                log_plot(figure_tot, 1, t_sn, Lumout_t, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)


                                # Sum of each contribution
                        Lumtot_sn_1.append(Lumtot_t_1)
                        Lumtot_sn_2.append(Lumtot_t_2)
                        Lumtot_sn.append(Lumtot_t)

                                    # Plot
                        sym = '-'
                        label = 'Total'

                                        # 100 MeV to 100 GeV
                        log_plot(figure_1, 1, t_sn, Lumtot_t_1, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 GeV to 100 TeV
                        log_plot(figure_2, 1, t_sn, Lumtot_t_2, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                                        # 100 MeV to 100 TeV
                        log_plot(figure_tot, 1, t_sn, Lumtot_t, label, 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
                        figure_number = figure_tot + 1


            # for all SN explosions

        Lumtot_sn_1 = numpy.asarray(Lumtot_sn_1)        # from 100 MeV to 100 GeV
        Lumtot_sn_2 = numpy.asarray(Lumtot_sn_2)        # from 100 GeV to 100 TeV
        Lumtot_sn = numpy.asarray(Lumtot_sn)            # from 100 MeV to 100 TeV

            # Initialization
        Lum_sn_1_tot = numpy.zeros(len(t))
        Lum_sn_2_tot = numpy.zeros(len(t))
        Lum_sn_tot = numpy.zeros(len(t))

        figure_1 = figure_number
        figure_2 = figure_1 + 1
        figure_tot = figure_2 + 1

        for i in range (nt0):

            indSN = numpy.where(t > t0[i])[0]                       # only time when a SN has already occured
            l = indSN[0]                                            # first index of t > t0
            t_sn = numpy.asarray(t)[indSN] * yr26yr                 # in Myr
            t_sn_unit = units.Myr

                # Initialization
            nt = len(t_sn)
            Lum_sn_1 = numpy.asarray(Lumtot_sn_1[i])
            Lum_sn_2 = numpy.asarray(Lumtot_sn_2[i])
            Lum_sn = numpy.asarray(Lumtot_sn[i])

            for j in (indSN):

                Lum_sn_1_tot[j] += Lum_sn_1[j-l]
                Lum_sn_2_tot[j] += Lum_sn_2[j-l]
                Lum_sn_tot[j] += Lum_sn[j-l]

                # Plot

            for zone in (zones):

                label = 'none'

                if zone == 1:               # in the SB

                    sym = '-.'

                            # from 100 MeV to 100 GeV
                    log_plot(figure_1, 1, t_sn, Lumsb_sn_1[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 GeV to 100 TeV
                    log_plot(figure_2, 1, t_sn, Lumsb_sn_2[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 MeV to 100 TeV
                    log_plot(figure_tot, 1, t_sn, Lumsb_sn[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                elif zone == 2:             # in the supershell

                    sym = '--'
                            # from 100 MeV to 100 GeV
                    log_plot(figure_1, 1, t_sn, Lumshell_sn_1[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 GeV to 100 TeV
                    log_plot(figure_2, 1, t_sn, Lumshell_sn_2[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 MeV to 100 TeV
                    log_plot(figure_tot, 1, t_sn, Lumshell_sn[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                else:                       # outside the SB

                    sym = ':'

                            # from 100 MeV to 100 GeV
                    log_plot(figure_1, 1, t_sn, Lumout_sn_1[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 GeV to 100 TeV
                    log_plot(figure_2 + 1, 1, t_sn, Lumout_sn_2[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

                            # from 100 MeV to 100 TeV
                    log_plot(figure_tot, 1, t_sn, Lumout_sn[i], label, 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)


        ind0_1 = numpy.where(Lum_sn_1_tot != 0)[0]
        ind0_2 = numpy.where(Lum_sn_2_tot != 0)[0]
        ind0 = numpy.where(Lum_sn_tot != 0)[0]

        log_plot(figure_1, 1, t[ind0_1] * yr26yr, Lum_sn_1_tot[ind0_1], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all_range1_it%d.eps'%iteration)
        log_plot(figure_2, 1, t[ind0_2] * yr26yr, Lum_sn_2_tot[ind0_2], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all_range2_it%d.eps' %iteration)
        log_plot(figure_tot, 1, t[ind0] * yr26yr, Lum_sn_tot[ind0], 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
        plt.savefig(pathfigure+'Gamma_luminosity_all_it%d.eps'%iteration)
        #plt.show()

        figure_number = figure_tot + 1

    return Lum_sn_1_tot, Lum_sn_2_tot, Lum_sn_tot, figure_number

def spectral_index(indE, zones, pathfigure, iteration, figure_number):

    """
    Return the spectral index at a precise energy_load

    Inputs:
        indE        :       which energy the spectral index is compute (GeV)
        zones       :       which zone do you want to compute (1: cavity of the SB, 2: supershell and 3: outside)
        pathfigure  :       path to save figures
        iteration   :       number of the iteration
        figure_number:      number of the figure
    """

    with open('SN', 'rb') as t0_load:
        with open('time', 'rb') as time_load:
            with open('spectra', 'rb') as spectra_load:
                with open('radius', 'rb') as radius_load:

                        # Sn explosions time (yr)
                    t0 = pickle.load(t0_load)
                    nt0 = len(t0)

                        # time array (yr)
                    t = pickle.load(time_load)
                    t_sn_unit = units.Myr

                        # radius (pc)
                    radius = pickle.load(radius_load)

                        # Energy of the gamma photon (GeV)
                    E_gamma = pickle.load(spectra_load)
                    E_gamma_unit = pickle.load(spectra_load)
                    nE = len(E_gamma)

                        # Spectral distribution of the gamma photons (erg s^-1)
                    sed = pickle.load(spectra_load)          # E^2 dN/dE (erg s^-1)
                    print(sed.shape)
                    sed_unit = pickle.load(spectra_load)


                    alpha_t0 = []

                    for i in range (nt0):

                        indSN = numpy.where(t > t0[i])[0]
                        nt = len(t[indSN])
                        lum_gamma = sed[i]/E_gamma                             # E dN/dE (erg s^- GeV^-1)

                            # spectral index for each time and energy
                        alpha_t = []

                        for j in range (nt):

                            Lum_t = numpy.asarray(lum_gamma[j])
                            n = 0
                            alpha_zone = []

                            for zone in (zones):

                                if zone == 1:                                   # in the SB

                                    r = radius[j]
                                    nr = len(radius)

                                    for k in range (nr):
                                        Lum_r = Lum_t[k]

                                        if indE < nE-1:

                                            alpha = -(numpy.log(Lum_r[indE + 1]) - numpy.log(Lum_r[indE]))/(numpy.log(E_gamma[indE + 1]) - numpy.log(E_gamma[indE]))
                                        else:
                                            alpha = -(numpy.log(Lum_r[indE]) - numpy.log(Lum_r[indE - 1]))/(numpy.log(E_gamma[indE]) - numpy.log(E_gamma[indE - 1]))

                                    alpha_zone.append(alpha)
                                    n = nr

                                elif zone == 2:                                 # in the supershell

                                    Lum_r = Lum_t[n]

                                    if indE < nE-1:

                                        alpha = -(numpy.log(Lum_r[indE + 1]) - numpy.log(Lum_r[indE]))/(numpy.log(E_gamma[indE + 1]) - numpy.log(E_gamma[indE]))
                                    else:
                                        alpha = -(numpy.log(Lum_r[indE]) - numpy.log(Lum_r[indE - 1]))/(numpy.log(E_gamma[indE]) - numpy.log(E_gamma[indE - 1]))

                                    alpha_zone.append(alpha)
                                    n += 1

                                else:                                           # outside the SB

                                    Lum_r = Lum_t[n]

                                    if indE < nE-1:

                                        alpha = -(numpy.log(Lum_r[indE + 1]) - numpy.log(Lum_r[indE]))/(numpy.log(E_gamma[indE + 1]) - numpy.log(E_gamma[indE]))
                                    else:
                                        alpha = -(numpy.log(Lum_r[indE]) - numpy.log(Lum_r[indE - 1]))/(numpy.log(E_gamma[indE]) - numpy.log(E_gamma[indE - 1]))

                                    alpha_zone.append(alpha)

                            alpha_t.append(alpha_zone)

                        alpha_t0.append(alpha_t)

                            # Plot

                        for zone in (zones):

                            n = 0

                            if zone == 1:                                           # in the SB

                                label = 'SB'
                                sym = '-.'

                                for k in range (nr):
                                    alpha = numpy.asarray(alpha_t)[:,k]
                                    log_plot(figure_number, 1, t[indSN] * yr26yr, alpha, label, u'Spectral index of the gamma emission at 'r'$E_\gamma$'u' = %.2e {0}'.format(E_gamma_unit.to_string('latex_inline')) %E_gamma[indE], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), r'$\alpha$', sym)
                                    plt.savefig(pathfigure+'spectral_index_SB_%dit.eps'%iteration)

                                n = nr
                                figure_number += 1

                            elif zone == 2:                                         # in the supershell

                                label = 'Shell'
                                sym = '--'

                                alpha = numpy.asarray(alpha_t)[:, n]

                                log_plot(figure_number, 1, t[indSN] * yr26yr, alpha, label, u'Spectral index of the gamma emission at 'r'$E_\gamma$'u' = %.2e {0}'.format(E_gamma_unit.to_string('latex_inline')) %E_gamma[indE], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), r'$\alpha$', sym)
                                plt.savefig(pathfigure+'spectral_index_shell_%dit.eps'%iteration)

                                n += 1
                                figure_number += 1

                            else:                                                   # outside the SB

                                label = 'Out'
                                sym = ':'

                                alpha = numpy.asarray(alpha_t)[:, n]

                                log_plot(figure_number, 1, t[indSN] * yr26yr, alpha, label, u'Spectral index of the gamma emission at 'r'$E_\gamma$'u' = %.2e {0}'.format(E_gamma_unit.to_string('latex_inline')) %E_gamma[indE], 'Time [{0}]'.format(t_sn_unit.to_string('latex_inline')), r'$\alpha$', sym)
                                plt.savefig(pathfigure+'spectral_index_out_%dit_E%d.eps'%(iteration, E_gamma[indE]))

                    alpha_t0 = numpy.asarray(alpha_t0)
                    #plt.show()
    return

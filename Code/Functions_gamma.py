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

def data(Emin_CR, Emax_CR, p0, alpha, D0, delta, t0):

    """
    Return 5 files with differents quantities.
    Inputs:
        Emin_CR     :       minimum kinetic energy of the relativistic particles (GeV)
        Emax_CR     :       maximum kinetic energy of the relativistic particles (GeV)
        p0          :       normalization constant of the power-law distribution of the cosmic rays (GeV/c)
        alpha       :       exponent of the power-law distribution of the cosmic rays
        D0          :       diffusion coefficient at p0 (cm^2 s^-1)
        delta       :       exponent of the power-law distribution of the diffusion coefficient
        t0          :       array of the Sn explosion time (yr)

    Outputs:
        gas         :       file with the density of gas for all time step and radius step (cm^-3)
        energy      :       file with the energy of the relativistic particles (GeV)
        time        :       file with the time array from the first SN explosion to 1 millions years after the last SN explosion (yr)
        distance    :       file with the radius array for each time step (pc)
        data        :       file with the distribution of relativistic particles for each time step, each SN explosion and each zones (GeV^-1)
    """

    # Parameters for the system

        # luminosity
    L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
    L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s

        # in the ISM
    pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
    C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

    with open('gas', 'wb') as gas_write:                                            # for the density of gas (cm^-3)
        with open('energy', 'wb') as energy_write:                                  # for the energy of the particles (GeV)
            with open('time', 'wb') as time_write:                                  # for the time array (yr)
                with open('distance', 'wb') as distance_write:                      # for the distance array (pc)
                    with open('data', 'wb') as data_write:                          # for the particles distribution (GeV^-1)

                            # Energy
                        number_bin_E = 10
                        E = numpy.logspace(numpy.log10(Emin_CR), numpy.log10(Emax_CR), number_bin_E)  # GeV
                        E_unit = units.GeV      # recording the unit of the energy array

                        pickle.dump(E, energy_write)
                        pickle.dump(E_unit, energy_write)

                            ## ====================================================== ##
                            # power-law distribution of the cosmic rays (GeV^-1 cm^-3) #
                            ## ====================================================== ##
                                # N(p) = N0 * (p/p0)^(-alpha)
                                # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)
                        Esng = Esn*erg2GeV  # in GeV
                        N_E = power_law_distribution(Emin_CR, Emax_CR, E, alpha, eta, Esng, p0)

                            ## =============================== ##
                            # diffusion coefficient (cm^2 s^-1) #
                            ## =============================== ##
                                # D(p) = D0 * (p/p0)^(-delta)
                                # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
                        D = diffusion_coefficient(p0, D0, E, delta)

                            ## ====================================================== ##
                            # Computation of N(E,t,r), particles distribution (GeV^-1) #
                            ## ====================================================== ##

                                # SN explosion: time (yr)
                        nt0 = len(t0)

                        pickle.dump(t0, time_write)     # Recording the SN explosion time (yr)

                                # time vector (yr)
                        tmin = t0[0]       # nothing happens before the first SN explosion (yr)
                        tmax = t0[nt0-1] + 1e6
                        number_bin_t = 200
                        t = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)
                        nt = len(t)
                        t_unit = units.yr

                        pickle.dump(t, time_write)      # Recording the time array (yr)
                        pickle.dump(t_unit, time_write)

                            # Initialization
                        figure_number = 1

                            # Recording
                        Ntotsn = []               # particles distribution for each time, distance and energy: the matrix of dimension len(t0)xlen(t)x(len(r)+2)xlen(E)
                        ngastot = []               # density of gas for each time and distance : matrix of dimension len(t0)xlen(t)x(len(r)+2)
                        distance = []           # distance array for each time: array of dimension len(t0)xlen(t)x(len(r))

                            # For verification
                        #Nverift = []

                        for i in range (nt):        # for each time step
                            Ntot = []               # particles distribution for each time, distance and energy: the matrix of dimension len(t)x(len(r)+1)xlen(E)
                            ngas = []               # density of gas for each time and distance : matrix of dimension len(t)x(len(r)+1)

                                # Initialization
                            t6 = t[i] * yr26yr      # 10^6 yr
                            t7 = t6 * s6yr27yr      # 10^7 yr

                                # r array
                            Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)           # radius and velocity of the SB
                            Msb, Mswept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu)   # swept-up and inner masses (solar masses)
                            hs, ns = density_thickness_shell(mu, n0, Vsb, C02, Mswept, Msb, Rsb)            # thickness and density of the shell (pc)

                            rmin = 0.01                 # minimum radius (pc)
                            rmax = Rsb-hs               # maximum radius (pc)
                            number_bin_r = 15           # number of bin for r from 0 to Rsb-hs
                            r = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_r)    # position (pc)

                            distance.append(r)

                                # Density of gas (cm^-3)
                                    # First zone: in the cavity

                            nsb = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, r, Rsb)[1]
                            n_gas = nsb
                            n_gas = n_gas.tolist()

                                    # Second zone: in the supershell
                            n_gas.append(ns)

                                    # Third zone: outside the superbubble
                            n_gas.append(n0)

                                    # Recording for each time step
                            ngastot.append(n_gas)

                            for j in range (nt0):       # for each SN explosion

                                if t[i] < t0[j]:
                                    n = numpy.zeros((number_bin_r + 2, number_bin_E))
                                    Ntot.append(n)
                                    continue
                                else:

                                    deltaT = t[i]-t0[j]

                                        ## --------------------------------- ##
                                        # First zone: in the cavity of the SB #
                                        ## --------------------------------- ##

                                            # For recording the distribution of particles for each time
                                    Nr = []         # (GeV^-1) : matrix of dimension (len(r)-1) x len(E))

                                            # r array (pc)
                                    r_in = 0
                                    r_out = r[0]
                                    N_part = shell_particles(r_in, r_out, N_E, D, deltaT)
                                    Nr.append(N_part)

                                    for k in range (1, number_bin_r):   # for each radius inside the SB

                                        r_in = r_out
                                        r_out = r[k]
                                        N_part = shell_particles(r_in, r_out, N_E, D, deltaT)
                                        Nr.append(N_part)

                                        ## -------------------------------- ##
                                        # Second zone: inside the supershell #
                                        ## -------------------------------- ##

                                            # For the particle distribution (GeV^-1): N_part = 4 * pi * int_{Rsb-hs}^Rsb Ne
                                    r_in = rmax         # in pc
                                    r_out = Rsb         # in pc
                                    N_part = shell_particles(r_in, r_out, N_E, D, deltaT)
                                    Nr.append(N_part)

                                        ## ------------------------ ##
                                        # Third zone: outside the SB #
                                        ## ------------------------ ##
                                            # For the particle distribution (GeV^-1): N_part = 4 * pi * int_Rsb^inf Ne r^2 dr
                                    N_part = inf_particles(Rsb, N_E, D, deltaT)
                                    Nr.append(N_part)

                                        # --------- #
                                        # Recording
                                    Ntot.append(Nr)

                            Ntotsn.append(Ntot)

                        Ntotsn = numpy.asarray(Ntotsn)
                        Ntot_unit = 1/(units.GeV)               # units of Ntot (GeV^-1)
                        pickle.dump(Ntotsn, data_write)         # Recording the number of particles for each SN explosion time, each time, each zone and each energy
                        pickle.dump(Ntot_unit, data_write)

                        distance = numpy.asarray(distance)      # Recording the r array for each SN explosion time, each time, each zone
                        pickle.dump(distance, distance_write)

                        ngastot = numpy.asarray(ngastot)
                        ngas_unit = 1/units.cm**3           # units of ngas (cm^-3)
                        pickle.dump(ngastot, gas_write)     # Recording the density of gas for each SN explosion, each time and each zone
                        pickle.dump(ngas_unit, gas_write)

    return

def spectrum(Emin_gamma, Emax_gamma):

    """
    Return one file

    Inputs:
        Emin_gamma  :       minimum energy of the gamma photon (GeV)
        Emax_gamma  :       maximum energy of the gamma photon (GeV)

    Output:
        spectra     :       file with the energy spectrum of the gamma photons and the spectral energy distribution for each SN explosions, time step and zones
    """

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

                            r = pickle.load(distance_load)

                            t0 = pickle.load(time_load)
                            t = pickle.load(time_load)
                            t_unit = pickle.load(time_load)

                                # Energy spectrum (GeV)
                            number_bin_E = 20
                            spectrum_energy = numpy.logspace(numpy.log10(Emin_gamma), numpy.log10(Emax_gamma), number_bin_E) * units.GeV

                            nt = len(t)
                            nt0 = len(t0)

                            sedsn = []

                            for i in range (nt0):       # for each SN explosion

                                sed = []

                                for j in range (nt):    # for each time step
                                    """
                                    if t[j] < t0[i]:
                                        #sed.append(0)
                                        continue
                                    else:
                                    """
                                        # Initialization
                                    distance = r[j]
                                    nr = len(distance) + 2

                                        # Important quantities
                                    ngas = ngassn[j] * ngas_unit
                                    ntot = numpy.asarray(Ntotsn[j, i]) * Ntot_unit

                                        # Recording
                                    spectrum_r = []
                                    sed_r = []

                                    for k in range (nr):     # for each radius step
                                        model = TableModel(ECR, ntot[k], amplitude = 1)
                                        PD = PionDecay(model, nh = ngas[k], nuclear_enhancement = True)

                                        sed_PD = PD.sed(spectrum_energy, distance = 0 * units.pc)
                                        sed_r.append(sed_PD)

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
                r = pickle.load(distance_load)

                        # Energy (GeV)
                spectrum_energy = pickle.load(spectra_read)
                spectrum_unit = pickle.load(spectra_read)
                spectrum = spectrum_energy * spectrum_unit

                         # spectal energy distribution (erg s^-1): E^2 dN/dE = E L_gamma
                sed_PD_sn = pickle.load(spectra_read)
                sed_PD_unit = pickle.load(spectra_read)

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
                nt = len(t)

                    ##================##
                    # Gamma luminosity #
                    ##================##

                ind1 = numpy.where(spectrum_energy < Esep)[0]         # for the first range
                ind2 = numpy.where(spectrum_energy >= Esep)[0]        # for the second range

                for i in range (nt0):           # for each SN explosion

                    lum_energy = sed_PD_sn[i]/spectrum_energy
                    lum_energy_unit = sed_PD_unit/spectrum_unit

                        # For recording
                            # Within the SB
                    Lumsb_t_1 = []       # 100 MeV to 100 GeV
                    Lumsb_t_2 = []       # 100 GeV to 100 TeV
                    Lumsb_t = []        # 100 MeV to 100 TeV

                            # Within the supershell
                    Lumshell_t_1 = []    # 100 MeV to 100 GeV
                    Lumshell_t_2 = []    # 100 GeV to 100 TeV
                    Lumshell_t = []     # 100 MeV to 100 TeV

                            # Outside the SB
                    Lumout_t_1 = []      # 100 MeV to 100 GeV
                    Lumout_t_2 = []      # 100 GeV to 100 TeV
                    Lumout_t = []       # 100 MeV to 100 TeV

                            # Total
                    Lumtot_t_1 = []      # 100 MeV to 100 GeV
                    Lumtot_t_2 = []      # 100 GeV to 100 TeV
                    Lumtot_t = []       # 100 MeV to 100 TeV

                        # Initialization
                    indSN = numpy.where(t >= t0[i])[0]     # only time where the SN already explodes

                    for j in range (nt):         # for each time step


                            ##------------------------------------##
                            # First zone: in the cavity (erg s^-1) #
                            ##------------------------------------##
                                # Initialization
                        distance = r[j]
                        nr = len(distance)

                                # Recording
                        Lumsb_1_r = []
                        Lumsb_2_r = []
                        Lumsb_r = []

                        for k in range (nr):     # for each radius step

                                # From 100 MeV to 100 GeV
                            Lum_1 = luminosity(lum_energy[j, k, ind1], spectrum_energy[ind1])
                            Lumsb_1_r.append(Lum_1)

                                # From 100 GeV to 100 TeV
                            Lum_2 = luminosity(lum_energy[j, k, ind2], spectrum_energy[ind2])
                            Lumsb_2_r.append(Lum_2)

                                # From 100 MeV to 100 TeV
                            Lum = luminosity(lum_energy[j,k], spectrum_energy)
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
                        Lumshell_1 = luminosity(lum_energy[j, nr, ind1], spectrum_energy[ind1])
                        Lumshell_t_1.append(Lumshell_1)

                                # From 100 GeV to 100 TeV
                        Lumshell_2 = luminosity(lum_energy[j, nr, ind2], spectrum_energy[ind2])
                        Lumshell_t_2.append(Lumshell_2)

                                # From 100 MeV to 100 TeV
                        Lumshell = luminosity(lum_energy[j, nr], spectrum_energy)
                        Lumshell_t.append(Lumshell)

                            ##-------------------------------------##
                            # Third zone: outside the SB (erg s^-1) #
                            ##-------------------------------------##

                                # From 100 MeV to 100 GeV
                        Lumout_1 = luminosity(lum_energy[j, nr+1, ind1], spectrum_energy[ind1])
                        Lumout_t_1.append(Lumout_1)

                                # From 100 GeV to 100 TeV
                        Lumout_2 = luminosity(lum_energy[j, nr+1, ind2], spectrum_energy[ind2])
                        Lumout_t_2.append(Lumout_2)

                                # From 100 MeV to 100 TeV
                        Lumout = luminosity(lum_energy[j, nr+1], spectrum_energy)
                        Lumout_t.append(Lumout)

                            # Total gamma luminosity (erg s^-1)
                        Lumtot1 = numpy.sum([Lumsb_1, Lumshell_1, Lumout_1])
                        Lumtot2 = numpy.sum([Lumsb_2, Lumshell_2, Lumout_2])
                        Lumtot = numpy.sum([Lumsb, Lumshell, Lumout])

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
                    Lum_unit = lum_energy_unit * spectrum.unit

                        # Plot
                    log_plot(figure_number, 4, t, [Lumsb_t_1, Lumshell_t_1, Lumout_t_1, Lumtot_t_1], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                    plt.savefig(pathfigure+'Gamma_luminosity_range1_t0-%d.eps' %i)
                    figure_number += 1
                    log_plot(figure_number, 4, t, [Lumsb_t_2, Lumshell_t_2, Lumout_t_2, Lumtot_t_2], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_\gamma$ (100 GeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                    plt.savefig(pathfigure+'Gamma_luminosity_range2_t0-%d.eps' %i)
                    figure_number += 1
                    log_plot(figure_number, 4, t, [Lumsb_t, Lumshell_t, Lumout_t, Lumtot_t], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble ($t_{0}$ = %.2e yr)'%t0[i], 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_\gamma$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), ['+', '+', '+', '-'])
                    plt.savefig(pathfigure+'Gamma_luminosity_t0-%d.eps' %i)
                    figure_number += 1

    Lumtot_sn_1 = numpy.asarray(Lumtot_sn_1)
    Lumtot_sn_2 = numpy.asarray(Lumtot_sn_2)
    Lumtot_sn = numpy.asarray(Lumtot_sn)

    # Plot for all SN explosions

    Lum_sn_1_tot = numpy.zeros_like(t)
    Lum_sn_2_tot = numpy.zeros_like(t)
    Lum_sn_tot = numpy.zeros_like(t)

    for i in range (nt0):

        indSN = numpy.where(t >= t0[i])[0]  # consider only time when there is already the SN explosion

        for j in (indSN):

            Lum_sn_1_tot[j] += Lumtot_sn_1[i, j]
            Lum_sn_2_tot[j] += Lumtot_sn_2[i, j]
            Lum_sn_tot[j] += Lumtot_sn[i, j]

        sym = ['-', '--', ':', '-']

            # Plot
        log_plot(figure_number, 4, t, [Lumsb_sn_1[i], Lumshell_sn_1[i], Lumout_sn_1[i], Lumtot_sn_1[i]], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
        log_plot(figure_number + 1, 4, t, [Lumsb_sn_2[i], Lumshell_sn_2[i], Lumout_sn_2[i], Lumtot_sn_2[i]], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)
        log_plot(figure_number + 2, 4, t, [Lumsb_sn[i], Lumshell_sn[i], Lumout_sn[i], Lumtot_sn[i]], ['SB', 'Shell', 'Out', 'Total'], 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), sym)

    log_plot(figure_number, 1, t, Lum_sn_1_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (0.1-100 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
    plt.savefig(pathfigure+'Gamma_luminosity_all_range1.eps')
    log_plot(figure_number + 1, 1, t, Lum_sn_2_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100-100 000 GeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
    plt.savefig(pathfigure+'Gamma_luminosity_all_range2.eps')
    log_plot(figure_number + 2, 1, t, Lum_sn_tot, 'Total SN', 'Gamma emission of a superbubble', 'Time [{0}]'.format(time.unit.to_string('latex_inline')), '$L_E$ (100 MeV - 100 TeV) [{0}]'.format(Lum_unit.to_string('latex_inline')), '-')
    plt.savefig(pathfigure+'Gamma_luminosity_all.eps')

    plt.show()
    return

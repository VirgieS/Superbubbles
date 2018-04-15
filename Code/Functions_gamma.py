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

def luminosity(lum_energy, energy):

    """
    Return the luminosity in a range of energy

    Inputs:
        lum_energy  :       intrinsic luminosity array per energy (erg/s/eV)
        energy      :       range of energy in which we will compute the luminosity (eV)
    """
    lum = integrate.trapz(lum_energy, energy)
    lum = numpy.nan_to_num(lum)

    return lum

def spectral_index(Emin, Emax, lum_ph_min, lum_ph_max):

    """
    Return the spectral index at an energy

    Inputs:
        Emin        :       minimum energy of the range (GeV)
        Emax        :       maximum energy of the range (GeV)
        lum_ph_min  :       intrinsic luminosity in photons at Emin (ph/s)
        lum_ph_max  :       intrinsic luminosity in photons at Emax (ph/s)
    """

    return -(numpy.log(lum_ph_max) - numpy.log(lum_ph_min))/(numpy.log(Emax) - numpy.log(Emin))


def data(correction_factor, t0, t, Emin_CR, Emax_CR, Emin_gamma, Emax_gamma, Esep, p0, alpha, D0, delta, zones, pathfigure, iteration, figure_number):

    """
    Return 5 files with differents quantities.
    Inputs:
        correction_factor   : correction factor for the radius of the SB
        t0          :       (sorted) array of the SN explosion times (yr)
        t           :       time array (yr)
        Emin_CR     :       minimum kinetic energy of the relativistic particles (GeV)
        Emax_CR     :       maximum kinetic energy of the relativistic particles (GeV)
        Emin_gamma  :       minimum energy of the gamma photon (GeV)
        Emax_gamma  :       maximum energy of the gamma photon (GeV)
        Esep        :       at which energy is the separation between the two ranges of energy (GeV)
        p0          :       normalization constant of the power-law distribution of the cosmic rays (GeV/c)
        alpha       :       exponent of the power-law distribution of the cosmic rays
        D0          :       diffusion coefficient at p0 (cm^2 s^-1)
        delta       :       exponent of the power-law distribution of the diffusion coefficient
        zones       :       which zone do you want to compute (1: cavity of the SB, 2: supershell and 3: outside)
        pathfigure  :       path of the file where to save the figure of Gamma emission
        iteration   :       number of the iteration
        figure_number:      number of the figure
    """

        ##===============================##
        # Energy of the cosmic rays (GeV) #
        ##===============================##
    number_bin_E = 10
    ECR = numpy.logspace(numpy.log10(Emin_CR), numpy.log10(Emax_CR), number_bin_E)
    E_CR = ECR * units.GeV

        ## ====================================================== ##
        # power-law distribution of the cosmic rays (GeV^-1 cm^-3) #
        ## ====================================================== ##
            # N(p) = N0 * (p/p0)^(-alpha)
            # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)
    Esng = Esn * erg2GeV    # GeV
    N_E = power_law_distribution(Emin_CR, Emax_CR, ECR, alpha, eta, Esng, p0)

        ## =============================== ##
        # diffusion coefficient (cm^2 s^-1) #
        ## =============================== ##
            # D(p) = D0 * (p/p0)^(-delta)
            # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
    D = diffusion_coefficient(p0, D0, ECR, delta)

        ##===========================##
        # Energy of the gamma photons #
        ##===========================##
    number_bin_E = 20
    spectrum = numpy.logspace(numpy.log10(Emin_gamma), numpy.log10(Emax_gamma), number_bin_E)   # GeV
    spectrum_erg = spectrum * 1.0/erg2GeV   # erg
    spectrum_ev = spectrum * GeV2eV         # eV
    spectrum_energy = spectrum * units.GeV

            # range of energy for the observation
    ind_HESS = numpy.where((spectrum > Esep[0]) & (spectrum <= Esep[1]))[0]      # from 1 TeV to 10 TeV

        ## =============================================================== ##
        # Computation of gamma luminosity in each range of energy(erg s^-1) #
        ## =============================================================== ##

            # Initialization
                # length of the arrays
    nt0 = len(t0)
    nt = len(t)

                # time evolution of the number of OB-stars
    nob = Nob * numpy.ones(nt)

                # total gamma luminosity (erg s^-1)
    for zone in (zones):

        if zone == 1:       # in the SB

            Lumsb_sn_HESS = numpy.zeros(nt)         # 1 TeV to 10 TeV
            Lumsb_sn = numpy.zeros(nt)              # 100 MeV to 100 TeV

        elif zone == 2:     # in the supershell

            Lumshell_sn_HESS = numpy.zeros(nt)      # 1 TeV to 10 TeV
            Lumshell_sn = numpy.zeros(nt)           # 100 MeV to 100 TeV

        else:               # outside the SB

            Lumout_sn_HESS = numpy.zeros(nt)        # 1 TeV to 10 TeV
            Lumout_sn = numpy.zeros(nt)             # 100 MeV to 100 TeV

            # Total
    Lumtot_sn_HESS = numpy.zeros(nt)        # 1 TeV to 10 TeV
    Lumtot_sn = numpy.zeros(nt)             # 100 MeV to 100 TeV

                # total spectral index
    Gamma_sn = numpy.zeros(nt)

        # Pulsar wind nebula
    n_pwn_tot = numpy.zeros(nt)

    for i in range (nt0):                                                       # for each SN explosions

            # Time array (yr)
        tmin = t0[i]        # only when the SN occurs (yr)
        tmax = tmin + 1e6
        number_bin_t = 200
        time = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)
        time6 = time * yr26yr   # in Myr

            # Initialization
                # only time corresponding to the time array
        indt = numpy.where((t >= t0[i]) & (t <= tmax))[0]
        indtob = numpy.where(t >= t0[i])[0]

                # gamma luminosity (erg s^-1)
        for zone in (zones):

            if zone == 1:                       # in the SB

                LumHESS_t_sb = numpy.zeros(number_bin_t)
                Lum_t_sb = numpy.zeros(number_bin_t)

            elif zone == 2:                     # in the supershell

                LumHESS_t_shell = numpy.zeros(number_bin_t)
                Lum_t_shell = numpy.zeros(number_bin_t)

            else:                               # outside the SB

                LumHESS_t_out = numpy.zeros(number_bin_t)
                Lum_t_out = numpy.zeros(number_bin_t)

        LumHESS_t_tot = numpy.zeros(number_bin_t)
        Lum_t_tot = numpy.zeros(number_bin_t)

        Gamma = numpy.zeros(number_bin_t)

                # number of pulsar wind nebula
        t_pwn_min = t0[i]               # in yr
        t_pwn_max = t_pwn_min + 1e5     # in yr
        indt_pwn = numpy.where((t >= t_pwn_min) & (t <= t_pwn_max))[0]

        for j in range (number_bin_t):                                          # for each time step

                # Initialization
            t6 = time6[j]           # 10^6 yr
            t7 = t6 * s6yr27yr      # 10^7 yr
            delta_t = time[j] - time[0]
            SB = False  # we do not compute inside the SB
            fluxmin = 0
            fluxmax = 0

                # Parameters of the SB
            Rsb = radius_velocity_SB(t6) [0]                                                # radius and velocity of the SB
            Rsb = correction_factor * Rsb
            Msb, Mswept = masses(t7, Rsb)                                                   # swept-up and inner masses (solar masses)
            hs, ns = density_thickness_shell_percentage(percentage, Rsb, Mswept, Msb)       # thickness and density of the shell (pc)

                # For each zones
            for zone in (zones):

                if zone == 1:                                                   # inside the SB

                    SB = True   # we compute the interior of the SB

                        # distance array (pc)
                    rmin = 0.01                 # minimum radius (pc)
                    rmax = Rsb-hs               # maximum radius (pc)
                    number_bin_r = 15           # number of bin for r from 0 to Rsb-hs
                    r = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_r)    # position (pc)

                        # Density of gas (cm^-3)
                    nsb = profile_density_temperature(t7, r, Rsb)[1]
                    ngas = nsb * 1/units.cm**3

                        # Particles distribution (GeV^-1)
                    r_in = 0
                    r_out = r[0]
                    N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                        # Differential flux (eV^-1 s^-1)
                    model = TableModel(E_CR, N_part, amplitude = 1)
                    PD = PionDecay(model, nh = ngas[0], nuclear_enhancement = True, useLUT = False)
                    flux_PD = PD.flux(spectrum_energy, distance = 0 * units.pc)

                    fluxmin += numpy.asarray(flux_PD[ind_HESS[0]])
                    fluxmax += numpy.asarray(flux_PD[ind_HESS[-1]])

                        # Gamma luminosity (erg s^-1)
                    lum_energy = numpy.asarray(flux_PD) * spectrum_erg          # erg s^- eV^-1
                    lum_units = flux_PD.units * units.erg * units.eV

                    LumHESS = luminosity(lum_energy[ind_HESS], spectrum_ev[ind_HESS])
                    LumHESS_t[j] += LumHESS
                    Lum = luminosity(lum_energy, spectrum_ev)
                    Lum_t[j] += Lum

                    for k in range (1, number_bin_r):   # for each radius inside the SB

                            # Particles distribution (GeV^-1)
                        r_in = r_out
                        r_out = r[k]
                        N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                            # Differential flux (GeV^-1 s^-1)
                        model = TableModel(E_CR, N_part, amplitude = 1)
                        PD = PionDecay(model, nh = ngas[k], nuclear_enhancement = True, useLUT = False)
                        flux_PD = PD.flux(spectrum_energy, distance = 0 * units.pc)

                        fluxmin += numpy.asarray(flux_PD[ind_HESS[0]])
                        fluxmax += numpy.asarray(flux_PD[ind_HESS[-1]])

                            # Gamma luminosity (erg s^-1)
                        lum_energy = numpy.asarray(flux_PD) * spectrum_erg      # erg s^-1 eV^-1
                        lum_units = flux_Pd.unit * units.erg * units.eV

                        LumHESS = luminosity(lum_energy[ind_HESS], spectrum_ev[ind_HESS])
                        LumHESS_t_sb[j] += LumHESS
                        Lum = luminosity(lum_energy, spectrum_ev)
                        Lum_t_sb[j] += Lum

                elif zone == 2:                                                 # in the supershell

                        # Density of gas (cm^-3)
                    ngas = ns * 1/units.cm**3

                        # Particles distribution (GeV^-1)
                    r_in = Rsb - hs     # in pc
                    r_out = Rsb         # in pc
                    N_part = shell_particles(r_in, r_out, N_E, D, delta_t) * 1/units.GeV

                        # Spectral disbution (erg s^-1)
                    model = TableModel(E_CR, N_part, amplitude = 1)
                    PD = PionDecay(model, nh = ngas, nuclear_enhancement = True, useLUT = False)
                    flux_PD = PD.flux(spectrum_energy, distance = 0 * units.pc)

                    fluxmin += numpy.asarray(flux_PD[ind_HESS[0]])
                    fluxmax += numpy.asarray(flux_PD[ind_HESS[-1]])

                        # Gamma luminosity (erg s^-1)
                    lum_energy = numpy.asarray(flux_PD) * spectrum_erg          # erg s^-1 eV^-1
                    lum_units = flux_PD.unit * units.erg * units.eV

                    LumHESS = luminosity(lum_energy[ind_HESS], spectrum_ev[ind_HESS])
                    LumHESS_t_shell[j] += LumHESS
                    Lum = luminosity(lum_energy, spectrum_ev)
                    Lum_t_shell[j] += Lum

                else:                                                           # outside the SB

                        # Density of gas
                    ngas = n0 * 1/units.cm**3

                        # Distribution of particles (GeV^-1)
                    N_part = inf_particles(Rsb, N_E, D, delta_t) * 1/units.GeV

                         # Differential flux (GeV^-1 s^-1)
                    model = TableModel(E_CR, N_part, amplitude = 1)
                    PD = PionDecay(model, nh = ngas, nuclear_enhancement = True, useLUT = False)
                    flux_PD = PD.flux(spectrum_energy, distance = 0 * units.pc)

                        # Gamma luminosity (erg s^-1)
                    lum_energy = numpy.asarray(flux_PD) * spectrum_erg          # erg s^-1 eV^-1
                    lum_units = flux_PD.units * units.erg * units.eV

                    LumHESS = luminosity(lum_energy[ind_HESS], spectrum_ev[ind_HESS])
                    LumHESS_t_out[j] += LumHESS
                    Lum = luminosity(lum_energy, spectrum_ev)
                    Lum_t_out[j] += Lum

            if SB:  # if we compute what happens inside the SB

                LumHESS_t_tot[j] = LumHESS_t_sb[j] + LumHESS_t_shell[j]
                Lum_t_tot[j] = Lum_t_sb[j] + Lum_t_shell[j]

            else:   # the only relevant gamma luminosity is the one of the supershell

                LumHESS_t_tot[j] = LumHESS_t_shell[j]
                Lum_t_tot[j] = Lum_t_shell[j]

                    # spectral index
            Emin = spectrum[ind_HESS[0]] * GeV2eV       # eV
            Emax = spectrum[ind_HESS[-1]] * GeV2eV      # eV
            Gamma[j] = spectral_index(Emin, Emax, fluxmin, fluxmax)

            # Interpolation + plot

                # number of OB stars
        for l in (indtob):
            nob[l] -= 1

                # Pulsar wind nebula

        for l in (indt_pwn):
            n_pwn_tot[l] += 1

                # Gamma luminosity + spectral index

        Title_HESS = 'Gamma emission from 1 TeV to 10 TeV (t0 = %.2e yr)'%t0[i]
        Title = 'Gamma emission from 100 MeV to 100 TeV (t0 = %.2e yr)'%t0[i]
        xlabel = 'Time [Myr]'
        ylabel = '$L_\gamma$ [erg s$^-1$]'
        figure_HESS = figure_number
        figure = figure_HESS + 1

        for zone in (zones):

            if zone == 1:       # in the SB

                LumHESS_sb = interpolation(time, LumHESS_t_sb)
                Lum_sb = interpolation(time, Lum_t_sb)

                    # plot
                #label = 'SB'
                #sym = '-.'

                #log_plot(figure_HESS, 1, time6, LumHESS_t_sb, label, Title_HESS, xlabel, ylabel, sym)
                #log_plot(figure, 1, time6, Lum_t_sb, label, Title, xlabel, ylabel, sym)

            elif zone == 2:     # in the supershell

                LumHESS_shell = interpolation(time, LumHESS_t_shell)
                Lum_shell = interpolation(time, Lum_t_shell)

                    # plot
                #label = 'shell'
                #sym = '--'

                #log_plot(figure_HESS, 1, time6, LumHESS_t_shell, label, Title_HESS, xlabel, ylabel, sym)
                #log_plot(figure, 1, time6, Lum_t_shell, label, Title, xlabel, ylabel, sym)

            else:               # outside the SB

                LumHESS_out = interpolation(time, LumHESS_t_out)
                Lum_out = interpolation(time, Lum_t_out)

                    # plot
                #label = 'out'
                #sym = ':'

                #log_plot(figure_HESS, 1, time6, LumHESS_t_out, label, Title_HESS, xlabel, ylabel, sym)
                #log_plot(figure, 1, time6, Lum_t_out, label, Title, xlabel, ylabel, sym)

            # plot
        #label = 'total'
        #sym = ':'

        #log_plot(figure_HESS, 1, time6, LumHESS_t_tot, label, Title_HESS, xlabel, ylabel, sym)
        #plt.savefig(pathfigure+'Gamma_luminosity_range1_it%d_t0%d.eps'%(iteration,i))

        #log_plot(figure, 1, time6, Lum_t_tot, label, Title, xlabel, ylabel, sym)
        #plt.savefig(pathfigure+'Gamma_luminosity_all_it%d_t0%d.eps'%(iteration,i))
        #figure_number = figure + 1
        #plt.show()

        LumHESS_tot = interpolation(time, LumHESS_t_tot)
        Lum_tot = interpolation(time, Lum_t_tot)
        indL = numpy.where(Lum_t_tot != 0)

        Gamma_tot = interpolation(time, Gamma)

        for l in (indt):

            for zone in (zones):

                if zone == 1:       # in the SB

                    Lumsb_sn_HESS[l] += LumHESS_sb(t[l])
                    Lumsb_sn[l] += Lum_sb(t[l])

                elif zone == 2:     # in the supershell

                    Lumshell_sn_HESS[l] += LumHESS_shell(t[l])
                    Lumshell_sn[l] += Lum_shell(t[l])

                else:               # outside the SB

                    Lumout_sn_HESS[l] += LumHESS_out(t[l])
                    Lumout_sn[l] += Lum_out(t[l])

            Lumtot_sn_HESS[l] += LumHESS_tot(t[l])
            Lumtot_sn[l] += Lum_tot(t[l])
            Gamma_sn[l] += Gamma_tot(t[l])
            Gamma_sn[l]= numpy.nan_to_num(Gamma_sn[l])

    return Lumtot_sn_HESS, Lumtot_sn, Gamma_sn, lum_units, figure_number, n_pwn_tot, nob

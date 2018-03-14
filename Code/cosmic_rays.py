##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.special import erfc, erf
import astropy.units as units
import os
import pickle
from Functions import *

# Physical constants, conversion factors and parameters fot the SB
from Physical_constants import *
from Conversion_factors import *
from Parameters_SB import *

##--------------##
# Path for files #
##--------------##

    # At IRAP
os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files/5_SN/')
pathfigure = '/Users/stage/Documents/Virginie/Superbubbles/figures/5_SN/CR/'

    # At home
#os.chdir('/home/vivi/Documents/Master_2/Superbubbles/Files/n_SN')
#pathfigure = '/home/vivi/Documents/Master_2/Superbubbles/figures/n_SN/CR/'

##-----------##
# Computation #
##-----------##

# Parameters for the system

    # luminosity
L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s

    # in the ISM
pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

    # Recording

with open('gas', 'wb') as gas_write:                                            # for the density of gas (cm^-3)
    with open('energy', 'wb') as energy_write:                                  # for the energy of the particles (GeV)
        with open('time', 'wb') as time_write:                                  # for the time array (yr)
            with open('distance', 'wb') as distance_write:                      # for the distance array (pc)
                with open('data', 'wb') as data_write:                          # for the particles distribution (GeV^-1)

                        # Energy
                    Emin = 1            # minimum kinetic energy: Emin = 1GeV
                    Emax = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV
                    number_bin_E = 10
                    E = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E)  # GeV
                    E_unit = units.GeV      # recording the unit of the energy array

                    pickle.dump(E, energy_write)
                    pickle.dump(E_unit, energy_write)

                        ## ====================================================== ##
                        # power-law distribution of the cosmic rays (GeV^-1 cm^-3) #
                        ## ====================================================== ##
                            # N(p) = N0 * (p/p0)^(-alpha)
                            # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)
                    Esng = Esn*erg2GeV  # in GeV
                    p0 = 10             # normalization constant (GeV/c)
                    alpha = 2.0         # exponent of the power-law distribution

                    N_E = power_law_distribution(Emin, Emax, E, alpha, eta, Esng, p0)

                        ## =============================== ##
                        # diffusion coefficient (cm^2 s^-1) #
                        ## =============================== ##
                            # D(p) = D0 * (p/p0)^(-delta)
                            # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
                    delta = 1.0/3       # exponent of the power-law of the diffusion coefficient
                    D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1

                    D = diffusion_coefficient(p0, D0, E, delta)

                        ## ====================================================== ##
                        # Computation of N(E,t,r), particles distribution (GeV^-1) #
                        ## ====================================================== ##

                            # SN explosion: time (yr)
                    t0 = [3e6, 4e6, 5e6, 4.5e6]      # in yr
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

                                if deltaT < 1e-9:
                                    print('SN explosion')

                                    print('Radius of the Superbubble: %.2f pc'%Rsb)
                                    print('Diffusion coefficient at %.2e GeV: %.2e cm^2 s^{-1}' %(E[0], D[0]))
                                    tdiff = (Rsb*pc2cm)**2/(6*D[0]) * 1/yr2s
                                    print('Diffusion time: %.2f yr' %tdiff)
                                    print('Diffusion coefficient at %.2e GeV: %.2e cm^2 s^{-1}' %(E[number_bin_E - 1], D[number_bin_E - 1]))
                                    tdiff = (Rsb*pc2cm)**2/(6*D[number_bin_E - 1]) * 1/yr2s
                                    print('Diffusion time: %.2f yr' %tdiff)
                                    deltaT = 0

                                    ## --------------------------------- ##
                                    # First zone: in the cavity of the SB #
                                    ## --------------------------------- ##

                                """
                                if j == 5:
                                    plot(figure_number, 1, r, n_gas, 'none', 'Density of gas in the SB at t=%.2e yr' %t[j], 'radius (pc)', u'n(r) 'r'($cm^{-3}$)', '+', 'r')
                                    figure_number +=1
                                """
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

                                """
                                Nverife = []
                                for k in range (len(E)):

                                    nCR = diffusion_spherical(t[j], r, t0, N_E[k], D[k])
                                    Nverife.append(nCR)

                                Nverift.append(Nverife)
                                """
                                        # VERIFICATION of the density of relavistic particles (cm^-3 GeV^-1): evolution in radius for one time and one energy
                                """
                                if (j == 1):    # chosen of time interval after the SN explosion (yr)
                                    k = 0        # choosen energy
                                    print(t[j])
                                    Nverifr = diffusion_spherical(t[j], r, t0, N_E[k], D[k])
                                    popt, pcov = curve_fit(gauss, r, Nverifr)
                                    fit = gauss(r, *popt)
                                    plot(figure_number, 1, r, Nverifr, 'Simulation at %.2e' %E[k], 'Density of CR in the SB at t=%.2e yr' %t[j], 'radius (pc)', u'N(r) 'r'($GeV^{-1}$)', '+', 'b')
                                    plot(figure_number, 1, r, fit, 'Fit at E = %.2e GeV' %E[k], 'Density of CR in the SB at t=%.2e yr ' %t[j], 'radius (pc)', u'N(r) 'r'($GeV^{-1}$)', '-', 'g')
                                    figure_number += 1

                                        # Verification of the fitting
                                    delta_t = t[j]-t0
                                    print('dt = %.2e yr' %delta_t)
                                    print('The diffusion distance of the CR at dt = %.2e yr and E = %.2e GeV:' %(delta_t, E[k]))
                                    print(numpy.sqrt(6 * D[k] * delta_t * yr2s)/pc2cm)

                                    Dopt = popt[1]
                                    print('The simulated standard deviation of the CR at dt = %.2e and E = %.2e GeV:' %(delta_t, E[k]))
                                    print(numpy.sqrt(6 * Dopt))
                                """
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

# Important quantities
Ninit = numpy.sum(N_E)      # total number of particles injected in the SB (particles)

    ## ============================================================= ##
    # Number of particles at one time for each energy and each radius #
    ## ============================================================= ##

        # Choosen one time and one SN explosion
for j in range (nt0):

        # Initialization
    indSN = numpy.where(t >= t0[j])[0]

        # Parameters
    indt = len(t)-1        # chosen time (yr)

    ntot = Ntotsn[indt, j]
    ntot = numpy.asarray(ntot)

            # Initialization
    nr = len(Nr)     # length of Nr
    nE = len(E)      # length of E
    delta_t = t[indt] - t0[j]

            # Computation
    sum = numpy.zeros(nE)

    for k in range (nE):                # sum on the energy
        sum[k] = numpy.sum(ntot[:,k])

            # Plot
    log_plot(figure_number, nr, E, ntot, 'none', 'Number of CR in the SB at %.2e yr after the SN explosion' %delta_t , 'E [{0}]'.format(E_unit.to_string('latex_inline')), 'N(E) [{0}]'.format(Ntot_unit.unit.to_string('latex_inline')), '+')
    log_plot(figure_number, 1, E, sum, 'Sum', 'Number of CR in the SB at %.2e yr after the SN explosion'%delta_t , 'E [{0}]'.format(E_unit.to_string('latex_inline')), 'N(E) [{0}]'.format(Ntot_unit.unit.to_string('latex_inline')), 'x')
    log_plot(figure_number, 1, E, N_E, 'Injected', 'Number of CR in the SB at %.2e yr after the SN explosion'%delta_t , 'E [{0}]'.format(E_unit.to_string('latex_inline')), 'N(E) [{0}]'.format(Ntot_unit.unit.to_string('latex_inline')), '-')
    plt.savefig(pathfigure+'Particles_T%d_t0%d.eps'%(t[indt],j))
    figure_number +=1

        ## ========================================================= ##
        # VERIFICATION OF THE CONSERVATION OF THE NUMBER OF PARTICLES #
        ## ========================================================= ##

            # Initialization
    Nratio = []                 # ratio of remained particles in time

        # Computation
    for i in (indSN):            # index for time

        ntotSN = Ntotsn[i, j]
        ntotSN = numpy.asarray(ntotSN)
        Nremain = numpy.sum(ntotSN)
        #print(Nremain/Ninit)
        Nratio.append(Nremain/Ninit)

            # Plot
    plot(figure_number, 1, t[indSN], Nratio, 'none', 'Ratio of remained particles in the considered volume', 'Time [{0}]'.format(t_unit.to_string('latex_inline')), '$N_{tot, t}/N_{tot, 0}$', '-')
    plt.savefig(pathfigure+'Conservation_particles_t0%d.eps'%j)
    figure_number += 1

    ## ==================================================== ##
    # VERIFICATION OF THE ESCAPE TIME SCALE OF THE PARTICLES #
    ## ==================================================== ##

    # Initialization
    nr = len(r)
    label_name = []
    symbol = []

        # Superbubble
    NSB = numpy.zeros((nE, nt))

        # Supershell
    indshell = 15       # supershell zone
    Nshell = numpy.zeros((nE, nt))


        # Outside
    indout = 15       # supershell zone
    Nout = numpy.zeros((nE, nt))


        # Computation
    for i in (indSN):          # for each time

        ntot_es = Ntotsn[i, j]
        ntot_es = numpy.asarray(ntot_es)

        for l in range (nE):         # for each energy
            Nsum = 0

            ## ----------- ##
            # Inside the SB #
            ## ----------- ##

            for k in range (nr):    # sum over each radius inside the SB
                Nsum += ntot_es[k, l]

            NSB[l, i] = Nsum/N_E[l]
            NSB[l, i] = numpy.nan_to_num(NSB[l, i])

            ## ------------------- ##
            # Inside the supershell #
            ## ------------------- ##

            Nshell[l, i] = ntot_es[indshell, l]/N_E[l]
            Nshell[l, i] = numpy.nan_to_num(Nshell[l, i])

            ## ------------ ##
            # Outside the SB #
            ## ------------ ##
            Nout[l, i] = ntot_es[indout, l]/N_E[l]
            Nout[l, i] = numpy.nan_to_num(Nout[l, i])

            label_name.append('%.2e GeV'%E[l])
            symbol.append('+')

                # Plot
    log_plot(figure_number, nE, t[indSN], NSB[:,indSN], label_name, 'Time evolution of the number of particles inside the SB', 'Time [{0}]'.format(t_unit.to_string('latex_inline')), r'$N_{tot, r}/N_E$', symbol)
    plt.savefig(pathfigure+'Escape_SB_t0%d.eps'%j)
    figure_number +=1
    log_plot(figure_number, nE, t[indSN], Nshell[:, indSN], label_name, 'Time evolution of the number of particles inside the supershell', 'Time [{0}]'.format(t_unit.to_string('latex_inline')), r'$N_{tot, r}/N_E$', symbol)
    plt.savefig(pathfigure+'Escape_shell_t0%d.eps'%j)
    figure_number += 1
    log_plot(figure_number, nE, t[indSN], Nout[:, indSN], label_name, 'Time evolution of the number of particles outside the SB', 'Time [{0}]'.format(t_unit.to_string('latex_inline')), r'$N_{tot, r}/N_E$', symbol)
    plt.savefig(pathfigure+'Escape_out_t0%d.eps'%j)
    figure_number += 1

"""
    ## ======================================== ##
    # VERIFICATION OF THE PARTICLES DISTRIBUTION #
    ## ======================================== ##

    # Initialization
indr = 7   # choosen one radius
nt = len(t)
Nin = numpy.zeros(nt)

    # Computation
for i in range (nt):
    Nin[i] = numpy.sum(ntot[i, indr])/Ninit

    # plot
log_plot(figure_number, 1, t, Nin, 'none', 'Number of particles at %.2f pc' %r[indr], 'Time [{0}]'.format(t_unit.to_string('latex_inline')), r'$N_{tot, r}$', '-', pathfigure,)
"""

plt.show()

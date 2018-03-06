##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.special import erfc, erf
from Functions import *
import os
import pickle

# Physical constants, conversion factors and parameters fot the SB
from Physical_constants import *
from Conversion_factors import *
from Parameters_SB import *

##--------------##
# Path for files #
##--------------##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files')

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

    # Computation of the diffusion coefficient and the initial density profile
with open('gas', 'wb') as gas_write:
    with open('energy', 'wb') as energy_write:
        with open('time', 'wb') as time_write:
            with open('distance', 'wb') as distance_write:
                with open('data', 'wb') as data_write:

                            # Energies
                    Emin = 1            # minimum kinetic energy: Emin = 1GeV
                    Emax = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV
                    number_bin_E = 10.0
                    E = numpy.logspace(numpy.log10(Emin), numpy.log10(Emax), number_bin_E)  # GeV

                    my_energy_write = pickle.Pickler(energy_write)
                    my_energy_write.dump(E)

                        # power-law distribution of the cosmic rays (GeV^-1 cm^-3)
                            # N(p) = N0 * (p/p0)^(-alpha)
                            # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)

                    eta = 0.1           # efficiency of the cosmic rays acceleration
                    Esn = 1e51          # total kinetic energy released from the SN explosion (erg)
                    Esng = Esn*erg2GeV  # in GeV
                    p0 = 10             # normalization constant (GeV/c)
                    alpha = 2.0         # exponent of the power-law distribution

                    N_E = power_law_distribution(Emin, Emax, E, alpha, eta, Esn, p0)

                        # diffusion coefficient (cm^2 s^-1)
                            # D(p) = D0 * (p/p0)^(-delta)
                            # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
                    delta = 1.0/3       # exponent of the power-law of the diffusion coefficient
                    D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1

                    D = diffusion_coefficient(p0, D0, E, delta)

                        # Computation of N(E,t,r), particles distribution (GeV^-1)

                            # SN explosion: time (yr)
                    t0 = 3e6      # in yr

                            # time vector (yr)
                    tmin = t0       # nothing happens before the first SN explosion (yr)
                    tmax = lifetime
                    number_bin_t = 30
                    t = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)

                    my_time_write = pickle.Pickler(time_write)
                    my_time_write.dump(t)

                        # Initialization
                    figure_number = 1

                        # Recording
                    Ntot = []               # particles distribution for each timey, distance and energy: the matrix of dimension len(t)x(len(r)-1)xlen(E)
                    ngas = []               # density of gas for each time and distance : matrix of dimension len(t)xlen(r)
                    distance = []           # distance array for each time: array of dimension len(t)xlen(r)

                    for j in range (len(t)):

                        t6 = t[j] * yr26yr      # 10^6 yr
                        t7 = t6 * s6yr27yr      # 10^7 yr

                            # First zone: in the cavity of the SB
                                # Computation of the distance array (pc)
                        Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)           # radius and velocity of the SB
                        Msb, Mswept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu)   # swept-up and inner masses (solar masses)
                        hs, ns = density_thickness_shell(mu, n0, Vsb, C02, Mswept, Msb, Rsb)            # thickness and density of the shell (pc)
                        rmin = 0.01                 # minimum radius (pc)
                        rmax = Rsb-hs               # maximum radius (pc)
                        number_bin_r = 15           # number of bin for r from 0 to Rsb-hs
                        r = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_r)    # position (pc)

                        distance.append(r)

                            # Computation of the density of gas in the SB (cm^-3)
                        nsb = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, r, Rsb)[1]
                        n_gas = nsb
                        n_gas = n_gas.tolist()
                        """
                        if j == 5:
                            plot(figure_number, 1, r, n_gas, 'none', 'Density of gas in the SB at t=%.2e yr' %t[j], 'radius (pc)', u'n(r) 'r'($cm^{-3}$)', '+')
                            figure_number +=1
                        """
                            # For recording the distribution of particles for each time
                        Nr = []         # (GeV^-1) : matrix of dimension (len(r)-1) x len(E))

                            # For the verification of the density of particles (cm^-3 GeV^-1)
                        #Nverif = numpy.zeros_like(r)       # density of particles at a chosen energy (cm^-3 GeV^-1)
                        #energy = len(E)-1                  # index of the chosen energy

                        r_in = 0
                        r_out = r[0]
                        N_part = shell_particles(r_in, r_out, N_E, D, t[j])
                        Nr.append(N_part)

                        for i in range (1, number_bin_r):

                            r_in = r_out
                            r_out = r[i]
                            N_part = shell_particles(r_in, r_out, N_E, D, t[j])
                            Nr.append(N_part)

                        #Ne_out = diffusion_spherical(t[j], r, t0, NE[energy], D[energy])

                        #Nverif[len(Nverif)-1] = Ne_out

                        #if (t[j] - t0 == 10*dt):    # chosen of time interval after the SN explosion (yr)
                        #    popt, pcov = curve_fit(gauss, r, Nverif)
                        #    fit = gauss(r, *popt)
                        #    plot(figure_number, 1, r, Nverif, 'Simulation', 'Density of CR in the SB at t=%.2e yr and at %.2e GeV' %(t[j], E[energy]), 'radius (pc)', u'N(r) 'r'($GeV^{-1}$)', '+')
                        #    plot(figure_number, 1, r, fit, 'Fit', 'Density of CR in the SB at t=%.2e yr and at %.2e GeV' %(t[j], E[energy]), 'radius (pc)', u'N(r) 'r'($GeV^{-1}$)', '-')
                        #    figure_number += 1

                            # Verification
                        #    delta_t = t[j]-t0
                        #    print('The diffusion distance of the CR at dt = %.2e yr and E = %.2e GeV:' %(delta_t, E[len(E)-1]))
                        #    print(sqrt(6 * D[len(E)-1] * t[j] * yr2s)/pc2cm)

                        #    Dopt = popt[1]* (pc2cm)**2
                        #    print('The simulated standard deviation of the CR at dt = %.2e and E = %.2e GeV:' %(delta_t, E[len(E)-1]))
                        #    print(sqrt(6 * Dopt * yr2s)/pc2cm)

                            # Second zone: inside the supershell
                                # For the density of gas: n(r) = ns
                        n_gas.append(ns)

                                # For the particle distribution (GeV^-1): N_part = 4 * pi * int_{Rsb-hs}^Rsb Ne
                        r_in = Rsb-hs       # in pc
                        r_out = Rsb          # in pc
                        #a = (rmin*pc2cm)/(numpy.sqrt(4 * D * t[j]*yr2s))
                        #b = (rmax*pc2cm)/(numpy.sqrt(4 * D * t[j]*yr2s))
                        #N_part = N_E/(numpy.sqrt(numpy.pi)) * (numpy.sqrt(numpy.pi) * (erf(b) - erf(a)) + 2 * (numpy.exp(-b**2) * b - numpy.exp(-a**2) * a))
                        N_part = shell_particles(r_in, r_out, N_E, D, t[j])
                        #print(N_part)
                        Nr.append(N_part)

                            # Third zone: outside the SB
                                # For the density of gas (cm^-3): n(r) = n0
                        n_gas.append(n0)

                                # For the particle distribution (GeV^-1): N_part = 4 * pi * int_Rsb^inf Ne r^2 dr
                        #x = (Rsb*pc2cm)/(numpy.sqrt(4 * D * t[j]*yr2s))
                        #N_part = N_E/(numpy.sqrt(numpy.pi)) * (numpy.sqrt(numpy.pi) * erfc(b) + 2 * numpy.exp(-b**2) * b)
                        N_part = inf_particles(Rsb, N_E, D, t[j])
                        Nr.append(N_part)

                            # recording
                        Ntot.append(Nr)
                        ngas.append(n_gas)

                    Ntot = numpy.asarray(Ntot)
                    my_data_write = pickle.Pickler(data_write)
                    my_data_write.dump(Ntot)

                    distance = numpy.asarray(distance)
                    my_distance_write = pickle.Pickler(distance_write)
                    my_distance_write.dump(distance)

                ngas = numpy.asarray(ngas)
                my_gas_write = pickle.Pickler(gas_write)
                my_gas_write.dump(ngas)

    # number of particles at one time for each energy and each radius

        # Choosen one time
ind = 1        # chosen time (yr)
        # Initialization
n = len(Nr)      # length of r
m = len(E)      # length of E
y = numpy.zeros((n, m))

    # Computation
for i in range (n):
    for j in range (m):
        y[i, j] = Ntot[ind, i, j]
sum = numpy.zeros(m)
for k in range (m):
    sum[k] = numpy.sum(y[:,k])

    # Plot
label_name = 'none'
log_plot(figure_number, n, E, y, label_name, 'Number of CR in the SB at %.2e yr after the SN explosion' %t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '+')
log_plot(figure_number, 1, E, sum, 'Sum', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', 'x')
log_plot(figure_number, 1, E, N_E, 'Injected', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '-')
figure_number +=1

    # Verification of the conservation of shell_particles

        # Initialization
deltat = []                 # time interval after the SN explosion (yr)
Ninit = numpy.sum(N_E)      # total number of particles injected in the SB (particles)
Nratio = []                 # ratio of remained particles in time

for i in range (len(t)):            # index for time

    if t[i] <= t0:
        continue

    else:
        deltat.append(t[i]-t0)
        z = numpy.zeros((m, n))

        for j in range (n):         # index for distance
            for k in range (m):     # index for energy
                #if j == 16:
                    #print(Ntot[i, j, k])
                z[k, j] = Ntot[i, j, k]

        Nremain = numpy.sum(z)
        Nratio.append(Nremain/Ninit)
#print(z)
#print(Nremain)
plot(figure_number, 1, deltat, Nratio, 'none', 'Ratio of remained particles in the considered volume', 'Time after the explosion (yr)', r'$N_{tot, t}/N_{tot, 0}$', '-')

plt.show()

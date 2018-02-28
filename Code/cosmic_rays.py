##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from math import *
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from pylab import *
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
        with open('data', 'wb') as data_write:

                    # Energies
            Emin = 1            # minimum kinetic energy: Emin = 1GeV
            Emax = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV
            number_bin_E = 10.0
            E = numpy.logspace(log10(Emin), log10(Emax), number_bin_E)  # GeV

            my_energy_write = pickle.Pickler(energy_write)
            my_energy_write.dump(E)

                # power-law distribution of the cosmic rays (GeV^-1 cm^-3)
                    # N(p) = N0 * (p/p0)^(-alpha)
                    # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2)) d^3(r)
            eta = 0.1           # efficiency of the cosmic rays acceleration
            Esn = 1e51          # total kinetic energy released from the SN explosion (erg)
            Esng = Esn*erg2GeV  # in GeV
            mpgev = mp*MeV2GeV  # mass of the proton in GeV
            p0 = 10             # normalization constant (GeV/c)
            #E0 = sqrt(p0**2 + mpg**2) - mpgev
            alpha = 2.0
            integral_E = integrate.quad(lambda E: (E**2 + 2*mpgev*E)**(-(1 + alpha)/2.0) * (E + mpgev) * E, Emin, Emax)[0]
            N0 = eta * Esng * cl**(1-alpha) * p0**(-alpha) * 1.0/integral_E                             # normalization constant to have 0.1*Esn (GeV^-1 c)
            N_E = N0/cl**(1-alpha) * (E**2 + 2*mpgev*E)**(-(1+alpha)/2.0) * (E + mpgev)/p0**(-alpha)    # GeV^-1

                # diffusion coefficient (cm^2 s^-1)
                    # D(p) = D0 * (p/p0)^(-delta)
                    # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
            delta = 1.0/3
            D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1
            D = D0 * (numpy.sqrt(E**2 + 2*mpgev*E)/p0)**delta

                # Computation of N(E,t,r), particles distribution (GeV^-1)

                    # SN explosion: time (yr)
            t0 = 3e6      # in yr

                    # time vector (yr)
            tmin = t0       # nothing happens before the first SN explosion (yr)
            tmax = lifetime
            number_bin_t = 40
            t = numpy.logspace(log10(tmin), log10(tmax), number_bin_t)

                # Initialization
            figure_number = 1

                # Recording
            Ntot = []               # particles distribution for each timey, distance and energy: the matrix of dimension len(t)x(len(r)-1)xlen(E)
            ngas = []               # desnity of gas for each time and distance : matrix of dimension len(t)xlen(r)-1

            for j in range (len(t)):

                t6 = t[j] * yr26yr      # 10^6 yr
                t7 = t6 * s6yr27yr      # 10^7 yr

                    # Computation of the distance array (pc)
                Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)   # radius and velocity of the SB
                Msb, Mswept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu) # swept-up and inner masses (solar masses)
                hs, ns = density_thickness_shell(mu, n0, Vsb, C02, Mswept, Msb, Rsb)                           # thickness and density of the shell (pc)
                rmin = 0                    # minimum radius (pc)
                rmax = Rsb                  # maximum radius (pc)
                number_bin_r = 50           # number of bin for r from 0 to Rsb
                r = numpy.linspace(rmin, rmax, number_bin_r)    # position in pc

                    # Computation of the density of gas in the SB (cm^-3)
                n_gas = ns * numpy.ones_like(r)
                ind = numpy.where(r < Rsb-hs)
                rsb = r[ind]
                nsb = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, rsb, Rsb)[1]
                n_gas[ind] = nsb
                ngas.append(n_gas)

                if j == 5:
                    plot(figure_number, 1, r, n_gas, 'none', 'Density of gas in the SB at t=%.2e yr' %t[j], 'radius (pc)', u'n(r) 'r'($cm^{-3}$)', '+')
                    figure_number +=1

                    # For recording the distribution of particles for each time
                Nr = []         # (GeV^-1) : matrix of dimension (len(r)-1) x len(E))

                    # For the verification of the density of particles (cm^-3 GeV^-1)
                #Nverif = numpy.zeros_like(r)  # density of particles at a chosen energy (cm^-3 GeV^-1)
                #energy = len(E)-1                # index of the chosen energy

                for i in range (len(r)-1):

                        # Computation of the density of particles (cm^-3 GeV^-1)
                    Ne_in = diffusion_spherical(t[j], r[i], t0, N_E, D)            # density of particles at the inner radius (cm^-3 GeV^-1)
                    Ne_out = diffusion_spherical(t[j], r[i+1], t0, N_E, D)         # density of particles at the outer radius (cm^-3 GeV^-1)

                    #Nverif[i] = Ne_in[energy]

                        # Computation of the particles distribution (GeV^-1): array of dimension len(E)
                    N_part = shell_particles(Ne_in, r[i], Ne_out, r[i+1])

                        # recording
                    Nr.append(N_part)

                #Nverif[len(Nverif)-1] = Ne_out[energy]

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

                    # recording
                Ntot.append(Nr)

            Ntot = numpy.asarray(Ntot)
            my_data_write = pickle.Pickler(data_write)
            my_energy_write.dump(Ntot)

        ngas = numpy.asarray(ngas)
        my_gas_write = pickle.Pickler(gas_write)
        my_gas_write.dump(ngas)

    # number of particles at one time for each energy and each radius

        # Choosen one time
ind = 1        # chosen time (yr)
        # Initialization
n = len(r)-1    # length of r
m = len(E)      # length of E
y = numpy.zeros((n, m))
#label_name = []

    # Computation
for i in range (n):
    #label_name.append(r'$r_{in}$ = %.2f pc'%r[i])
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
                z[k, j] = Ntot[i, j, k]

        Nremain = numpy.sum(z)
        Nratio.append(Nremain/Ninit)
plot(figure_number, 1, deltat, Nratio, 'none', 'Ratio of remained particles in the considered volume', 'Time after the explosion (yr)', r'$N_{tot, t}/N_{tot, 0}$', '-')

plt.show()

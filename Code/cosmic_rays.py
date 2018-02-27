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

    # Initialization
lifetime = 30e6         # lifetime of the OB association (yr)
dt = 1e3                # time interval (yr)

number_bin_r = 20  # number of bin for r from 0 to Rsb

    # Computation of the diffusion coefficient and the initial density profile

with open('energy.dat', 'wb') as energy_write:
    with open('CR.dat', 'wb') as CR_write:
        with open('Gas.dat', 'wb') as gas_write:

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

            L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
            L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s

                    # diffusion coefficient (cm^2 s^-1)
                        # D(p) = D0 * (p/p0)^(-delta)
                        # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
            delta = 1.0/3
            D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1
            D = D0 * (numpy.sqrt(E**2 + 2*mpgev*E)/p0)**delta

                # Computation of N(E,t,r), particles distribution (GeV^-1)

                    # SN explosion: time (yr)
            t0 = 0      # in yr

                    # Initialization
            i = 0
            figure_number = 1

                    # To record
            Ne = []             # particles distribution for each energy, time and distance: the matrix (len(E)xlen(t)xlen(r)) to register the particles distribution

            Rmax = 2000
            dr = 0.2

            while i < len(E):

                    # Initialization
                t = []                  # in yr
                t.append(dt)            # time t = 0 is not interesting
                j = 0
                t6 = t[j] * yr26yr      # 10^6 yr
                t7 = t6 * s6yr27yr      # 10^7 yr

                    # Computation of the density of particles (cm^-3 GeV^-1)
                #Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)   # radius and velocity of the SB
                #dr = Rsb/(number_bin_r - 1)                                             # space interval (pc)
                Nr, r = diffusion_spherical(t[j], Rmax, t0, N_E[i], D[i], dr)          # density of particles (cm^-3 GeV^-1) and r array (pc)

                    # Computation of the particles distribution (GeV^-1)
                Nt = []         # particles distribution for each time (GeV^-1) : matrix (len(t)xlen(r))
                N_part = shell_particles(Nr, r)
                Nt.append(N_part)               # recording

                #if i == 0:
                #    Ngas = []           # to record the desnity of gas for each time and radius: matrix (len(t)xlen(r))

                while t[j] < lifetime:
                        # time (yr)
                    t.append(t[j] + dt)
                    j += 1
                    t6 = t[j] * yr26yr      # 10^6 yr
                    t7 = t6 * s6yr27yr      # 10^7 yr

                        # Computation of the density of particles (cm^-3 GeV^-1)
                    #Rsb, Vsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)   # radius and velocity of the SB
                    #dr = Rsb/(number_bin_r - 1)                                             # space interval (pc)
                    Nr, r = diffusion_spherical(t[j], Rmax, t0, N_E[i], D[i], dr)          # density of particles (cm^-3 GeV^-1) and r array (pc)

                        # Computation of the particles distribution (GeV^-1)
                    N_part = shell_particles(Nr, r)
                    #N_part = N_part.tolist()
                    Nt.append(N_part)               # recording

                    """
                    if i == 0:
                            # In the ISM
                        pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
                        C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

                            # In the SB
                        Msb, Mswept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu)   # mass in the SB and swept-up mass (solar masses)

                            # Computation of the density of gas and recording
                        Nrgas = profile_gas_density(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu, Vsb, C02, Mswept, Msb, r)
                        Ngas.append(Nrgas)
                    """
                    """
                    if ((t[j] - t0 == 1000*dt) and (i == len(E)-1)):
                        popt, pcov = curve_fit(gauss, r, Nr)
                        plt.figure(figure_number)
                        plt.plot(r, Nr, '+', label = 'Simulation')
                        plt.plot(r, gauss(r, *popt), label='Fit')
                        plt.title('Density of CR in the SB at t=%.2e yr and at %.2e GeV' %(t[j], E[i]))
                        plt.xlabel('radius (pc)')
                        plt.ylabel(u'N(r) 'r'($GeV^{-1}$)')
                        plt.legend(loc = 'best')
                        figure_number += 1

                            # Verification
                        #delta_t = t[j]-t0
                        #print('The diffusion distance of the CR at dt = %.2e yr and E = %.2e GeV:' %(delta_t, E[i]))
                        #print(sqrt(6 * D[i] * t[j] * yr2s)/pc2cm)

                        #Dopt = popt[1]* (pc2cm)**2
                        #print('The simulated standard deviation of the CR at dt = %.2e and E = %.2e GeV:' %(delta_t, E[i]))
                        #print(sqrt(6 * Dopt * yr2s)/pc2cm)
                    """

                Ne.append(Nt)
                i += 1

            Ne = numpy.asarray(Ne)
            #print(2*Rsb)
            #Ngas = numpy.asarray(Ngas)
            t = numpy.asarray(t)

            my_CR_write = pickle.Pickler(CR_write)
            my_CR_write.dump(Ne)

            #my_gas_write = pickle.Pickler(gas_write)
            #my_gas_write.dump(Ngas)


    # number of particles at one time for each energy and each radius
ind = numpy.where(t == t0 + 1000*dt)[0]       # choose one time
n = len(N_part)
m = len(E)

y = numpy.zeros((n, m))
#label_name = []

for i in range (n):
    #label_name.append(r'$r_{in}$ = %.2f pc'%r[i])
    for j in range (m):
        y[i, j] = Ne[j, ind, i]

sum = numpy.zeros(m)
for k in range (m):
    sum[k] = numpy.sum(y[:,k])

label_name = 'none'
log_plot(figure_number, n, E, y, label_name, 'Number of CR in the SB at %.2e yr after the SN explosion' %t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '+')
log_plot(figure_number, 1, E, sum, 'Sum', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', 'x')
log_plot(figure_number, 1, E, N_E, 'Injected', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '-')
figure_number +=1

    # Verification
deltat = []             # time interval after the SN explosion (yr)
Ntot = numpy.sum(N_E)   # total number of particles injected in the SB (particles)
Nratio = []             # ratio of remained particles in time

for i in range (len(t)):

    if t[i] <= t0:
        continue

    else:
        deltat.append(t[i]-t0)
        z = numpy.zeros((m, n))
        for j in range (m):
            for k in range (n):
                z[j, k] = Ne[j, i, k]

        Ntot_t = numpy.sum(z)
        Nratio.append(Ntot_t/Ntot)

plot(figure_number, 1, deltat, Nratio, 'none', 'Ratio of remained particles in the considered volume', 'Time after the explosion (yr)', r'$N_{tot, t}/N_{tot, 0}$', '-')


plt.show()

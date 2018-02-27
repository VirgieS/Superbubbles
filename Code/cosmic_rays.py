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
#from Parameters_SB import *

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

                    # diffusion coefficient (cm^2 s^-1)
                        # D(p) = D0 * (p/p0)^(-delta)
                        # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
            delta = 1.0/3
            D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1
            D = D0 * (numpy.sqrt(E**2 + 2*mpgev*E)/p0)**delta

                # Computation of N(E,t,r)

                    # SN explosion: time (yr)
            t0 = 0      # in yr

                    # Initialization
            i = 0
            Ne = []             # matrix (len(E)xlen(t)xlen(r)) to register the particles distribution
            #Ngas = []           # matrix (len(t)xlen(r)) to register the gas density for each time and each radius
            figure_number = 0

            R_max = 1000

                    # density of population (GeV^-1)
            while i < len(E):
                t = []          # in yr

                    # Initialization
                t.append(0)
                j = 0
                t6 = t[j] * yr26yr      # 10^6 yr
                t7 = t6 * s6yr27yr      # 10^7 yr
                Nr, r = diffusion_spherical(t[j], R_max, t0, N_E[i], D[i]) # vector (len(r))

                Nt = []         # matrix (len(t)xlen(r))
                N_part = shell_particles(Nr, r)
                Nt.append(N_part)

                #if i == 0:
                #    Ngas = []           # matrix (len(t)xlen(r)) to register the gas density for each time and each radius

                while t[j] < lifetime:
                        # time (yr)
                    t.append(t[j] + dt)
                    j += 1
                    t6 = t[j] * yr26yr      # 10^6 yr
                    t7 = t6 * s6yr27yr      # 10^7 yr

                        # Density of population as function of the radius (cm^-3)
                    #diff = sqrt(6 * D[i] * t[j] * yr2s)/pc2cm                               # distance of diffusion (pc)
                    Nr, r = diffusion_spherical(t[j], R_max, t0, N_E[i], D[i])               # density (cm^-3) and vector r

                    #if i == 0:
                    #    Nrgas = profile_gas_density(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu, Vsb, C02, Mswept, Msb, r)
                    #    Ngas.append(Nrgas)

                    N_part = shell_particles(Nr, r)
                    Nt.append(N_part)

                    if ((t[j] - t0 == 1000*dt) and (i == len(E)-1)):
                        popt, pcov = curve_fit(gauss, r, Nr)
                        figure_number += 1
                        plt.figure(figure_number)
                        plt.plot(r, Nr, '+', label = 'Simulation')
                        plt.plot(r, gauss(r, *popt), label='Fit')
                        plt.title('Density of CR in the SB at t=%.2e yr and at %.2e GeV' %(t[j], E[i]))
                        plt.xlabel('radius (pc)')
                        plt.ylabel(u'N(r) 'r'($GeV^{-1}$)')
                        plt.legend(loc = 'best')

                            # Verification
                        #delta_t = t[j]-t0
                        #print('The diffusion distance of the CR at dt = %.2e yr and E = %.2e GeV:' %(delta_t, E[i]))
                        #print(sqrt(6 * D[i] * t[j] * yr2s)/pc2cm)

                        #Dopt = popt[1]* (pc2cm)**2
                        #print('The simulated standard deviation of the CR at dt = %.2e and E = %.2e GeV:' %(delta_t, E[i]))
                        #print(sqrt(6 * Dopt * yr2s)/pc2cm)

                Ne.append(Nt)
                i += 1

            Ne = numpy.asarray(Ne)
            t = numpy.asarray(t)

            my_CR_write = pickle.Pickler(CR_write)
            my_CR_write.dump(Ne)


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
log_plot(2, n, E, y, label_name, 'Number of CR in the SB at %.2e yr after the SN explosion' %t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '+')
log_plot(2, 1, E, sum, 'Sum', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', 'x')
log_plot(2, 1, E, N_E, 'Injected', 'Number of CR in the SB at %.2e yr after the SN explosion'%t[ind] , 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '-')

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

plot(3, 1, deltat, Nratio, 'none', 'Ratio of remained particles in the considered volume', 'Time after the explosion (yr)', r'$N_{tot, t}/N_{tot, 0}$', '-')


plt.show()

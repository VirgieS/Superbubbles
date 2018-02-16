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

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##--------------##
# Path for files #
##--------------##

os.chdir('/Users/stage/Documents/Virginie/Superbubbles/Files')

##-----------##
# Computation #
##-----------##

# Parameters for the system

    # Initialization
lifetime = 30e6         # average lifetime of the lowest mass B star (yr)
dt = 1e3                # time interval (yr)
Rsb = []                # radius of the SB (pc)
Rsb.append(0)

    # Computation of the diffusion coefficient and the initial density profile

        # Energies
Emin = 1            # minimum kinetic energy: Emin = 1GeV
Emax = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV
number_bin_E = 10.0
E = numpy.logspace(log10(Emin), log10(Emax), number_bin_E)  # GeV

        # power-law distribution of the cosmic rays (GeV^-1 cm^-3)
            # N(p) = N0 * (p/p0)^(-alpha)
            # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mp*E0)^(-alpha/2))
eta = 0.1           # efficiency of the cosmic rays acceleration
Esn = 1e51          # total kinetic energy released from the SN explosion (erg)
Esng = Esn*erg2GeV  # in GeV
mpgev = mp*MeV2GeV  # mass of the proton in GeV
p0 = 10             # normalization constant (GeV/c)
#E0 = sqrt(p0**2 + mpg**2) - mpgev
alpha = 2.0
integral_E = integrate.quad(lambda E: (E**2 + 2*mpgev*E)**(-(1 + alpha)/2.0) * (E + mpgev) * E, Emin, Emax)[0]
N0 = eta * Esng * cl**(1-alpha) * p0**(-alpha) * 1.0/integral_E        # normalization constant to have 0.1*Esn (GeV^-1)
N_E = N0/cl**(1-alpha) * (E**2 + 2*mpgev*E)**(-(1+alpha)/2.0) * (E + mpgev)/p0**(-alpha)

        # diffusion coefficient (cm^2 s^-1)
            # D(p) = D0 * (p/p0)^(-delta)
            # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
delta = 1.0/3
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1
D = D0 * (numpy.sqrt(E**2 + 2*mpgev*E)/p0)**delta

    # Computation of N(E,t,r)

        # SN explosion: position (pc) and time (yr)
t0 = 0      # in yr
r0 = 0      # in pc

        # R(t) = a * n0^alphar * L36^betar * t6^gammar                  (pc)
ar = 27.0
alphar = -1.0/5
betar = 1.0/5
gammar = 3.0/5
        # parameters of the OB association
n0 = 1              # mean density of the interstellar medium (particle/cm^3)
Nob = 1             # number of OB-stars in the association
Lob = 1e36          # mean luminosity of an OB-star (erg/s)
Pob = Nob*Lob       # mean power of the association (erg/s)
L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s

        # Initialization
i = 0
Ne = []             # matrix (len(E)xlen(t)xlen(r))
figure_number = 0

        # density of population (GeV^-1)
while i < len(E):
    t = []          # in yr

        # Initialization
    t.append(0)
    j = 0
    t6 = t[j] * yr26yr     # in 10^6 yr
    Rsb = 0
    Nr, r = diffusion_spherical(t[j], Rsb, t0, N_E[i], r0, D[i]) # vector (len(r))

    Nt = []         # matrix (len(t)xlen(r))
    N_part = shell_particles(Nr, r)
    Nt.append(N_part)

    while t[j] < lifetime:
        t.append(t[j] + dt)
        j += 1
        t6 = t[j] * yr26yr     # in 10^6 yr
        Rsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)[0]
        Nr, r = diffusion_spherical(t[j], Rsb, t0, N_E[i], r0, D[i]) # vector (len(r))
        N_part = shell_particles(Nr, r)
        Nt.append(N_part)

        if ((t[j-1] - t0 == 2e6) and (i == 0)):
            popt, pcov = curve_fit(gauss, r, Nr)
            print(popt)
            print(N_E[i])
            print(D[i])
            figure_number += 1
            plt.figure(figure_number)
            plt.plot(r, Nr, '+')
            plt.plot(r, gauss(r, *popt))
            plt.title('Density of CR in the SB at t=%.2e yr and at %.2e GeV' %(t[j], E[i]))
            plt.xlabel('radius (pc)')
            plt.ylabel(u'N(r) 'r'($GeV^{-1}$)')

    Ne.append(Nt)
    i += 1
Ne = numpy.asarray(Ne)
t = numpy.asarray(t)
print(Ne.shape)

    # number of particles at one time for each energy and each radius
ind = numpy.where(t == t0 + 1000*dt)[0]

n = len(N_part)
m = len(E)

y = numpy.zeros((n, m))
label_name = []

for i in range (n):
    label_name.append(r'$r_{in}$ = %.2f pc'%r[i])
    for j in range (m):
        y[i, j] = Ne[j, ind, i]

sum = numpy.zeros(m)
for k in range (m):
    sum[k] = numpy.sum(y[:,k])


log_plot(2, n, E, y, label_name, 'Number of CR in the SB', 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '+')
log_plot(2, 1, E, sum, 'Sum', 'Number of CR in the SB', 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', 'x')
log_plot(2, 1, E, N_E, 'injected', 'Number of CR in the SB', 'E (GeV)', u'N(E) 'r'($GeV^{-1}$)', '-')

plt.show()
"""
with open('data', 'wb') as data_write:
    my_data_write = pickle.Pickler(data_write)
    my_data_write.dump(Ne)

with open('data', 'rb') as data_load:
    my_data_load = pickle.Unpickler(data_load)
    N_read = my_data_load.load()
"""

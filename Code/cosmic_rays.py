##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from math import *
import scipy.integrate as integrate
from Functions import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##-----------##
# Computation #
##-----------##

# Parameters for the system

    # Initialization
lifetime = 1e6          # average lifetime of the lowest mass B star (yr)
dt = 1e3                # time interval (yr)
Rsb = []                # radius of the SB (pc)
Rsb.append(0)

    # Computation of the diffusion coefficient and the initial density profile

        # Energies
Emin = 1            # minimum kinetic energy: Emin = 1GeV
Emax = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV
number_bin_E = 10.0
E = numpy.logspace(log10(Emin), log10(Emax), number_bin_E)  # GeV

        # power-law distribution of the cosmic rays
            # N(p) = N0 * (p/p0)^(-alpha)
            # N(E) = N0/c * ((E^2 + 2*mp*c^2*E)^(-(1+alpha)/2) * (E + mp*c^2)/((E0^2 + 2*mpg*E0)^(-alpha/2))
eta = 0.1           # efficiency of the cosmic rays acceleration
Esn = 1e51          # total kinetic energy released from the SN explosion (erg)
Esng = Esn*erg2GeV  # in GeV
mpg = mp*MeV2GeV    # mass of the proton in GeV
p0 = 10             # normalization constant (GeV/c)
#E0 = sqrt(p0**2 + mpg**2) - mpg
alpha = 2.0
integral_E = integrate.quad(lambda E: (E**2 + 2*mpg*E)**(-(1 + alpha)/2.0) * (E + mpg) * E, Emin, Emax)[0]
N0 = eta * Esng * cl * p0**(-alpha) * 1.0/integral_E        # normalization constant to have 0.1*Esn (GeV^-1)
NE = N0/cl * (E**2 + 2*mpg*E)**(-(1+alpha)/2.0) * (E + mpg)/p0**(-alpha)

        # diffusion coefficient (cm^2 s^-1)
            # D(p) = D0 * (p/p0)^(-delta)
            # D(E) = D0 * (E^2 + 2*mpg*E)^(delta/2) * 1/p0^delta
delta = 1.0/3
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1
D = D0 * (numpy.sqrt(E**2 + 2*mpg*E)/p0)**delta

    # Computation of N(E,t,r)

        # SN explosion: position (pc) qnd time (yr)
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

        # density of population (GeV^-1)
while i < len(E):
    t = []          # in yr
    t.append(0)
    j = 0
    Nt = []         # matrix (len(t)xlen(r))
    while t[j] < lifetime:
        t.append(t[j] + dt)
        j += 1
        t6 = t[j] * yr26yr     # in 10^6 yr
        Rsb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)[0]
        Nr = diffusion_spherical(t[j], Rsb, t0, NE[i], r0, D[i]) # vector (len(r))
        Nt.append(Nr)
    Ne.append(Nt)
    i += 1
Ne = numpy.asarray(Ne)
print(Ne.shape)

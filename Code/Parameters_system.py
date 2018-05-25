"""
Here are all the free parameters for the system!
"""

import numpy

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##=============================##
# Superbubble and OB assocation #
##=============================##

    # OB asociation
Nob = 10           # number of OB-stars in the association
Lob = 1e36          # mean luminosity of an OB-star (erg/s)
Pob = Nob * Lob     # mean power of the association (erg/s)
L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s
lifetime = 4e6 #30e6     # average lifetime of the lowest B star (yr)

    # Fit parameters of the model
n0 = 10             # ambient density (cm^-3)
xH = 0.9            # mass fraction of hydrogen
xHe = 0.1           # mass fraction of helium
mu = xH + xHe * 4   # average molecular weight
percentage = 0.10   # percentage of the thickness of the shell in term of the outer radius of the SB
Ts = 1e2            # shell temperature (K)

    # R(t) = ar * n0^alphar * L36^betar * t6^gammar                  (pc)
    # Equation 51 of Weaver et al. (1977)
ar = 27.0
alphar = -1.0/5
betar = 1.0/5
gammar = 3.0/5

    # V(t) = av * n0^alphav * L36^betav * t6^gammav                  (km/s)
    # Equation 52 of Weaver et al. (1977)
av = 16.0
alphav = -1.0/5
betav = 1.0/5
gammav = -2.0/5

    # T(x) = at * n0^alphat * L38^betat * t7^gammat * (1-x)^deltat  (K)
    # Equation 4 of Mac Low and McCray (1987)
at = 7.4e6 #3.6e6
alphat = 2.0/35
betat = 8.0/35
gammat = -6/35
deltat = 2.0/5

    # n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan  (cm^-3)
    # Equation 5 of Mac Low and McCray (1987)
an = 1.4e-2 #4.0e-3
alphan = 19.0/35
betan = 6.0/35
gamman = -22.0/35
deltan = -2.0/5

    # p(t) = ap * n0^alphap * L_36^betap * t6^gammap               (dyne cm^-2)
    # Equation 22 of Weaver et al. (1977)
ap = 3.3e-11 #4.12e-12
alphap = 3.0/5
betap = 2.0/5
gammap = -4.0/5

    # LambdaT = al * T6^etal * zeta                                 (erg cm^3 s^-1)
    # Equation 6 of Mac Low and McCray (1987)
al = 1e-22
etal = -0.7
zeta = 1.0      # metallicity of the local ISM
epsilon = 1.2   # ratio ne/n

    # In the ISM
TISM = 88                           # temperature of the ambient gas (K)
C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)
C0 = numpy.sqrt(C02)
pISM = n0*kb*TISM                   # pressure in dyne cm^-2

##===========##
# Cosmic rays #
##===========##

    # Released energy by the acceleration of particles for one SN (erg)
eta = 0.1               # efficiency of the cosmic rays acceleration
Esn = 1e51              # total kinetic energy released from the SN explosion (erg)
Esng = Esn * erg2GeV    # GeV

    # Energy (GeV)
Emin_CR = 1            # minimum kinetic energy: Emin = 1GeV
Emax_CR = 1*PeV2GeV    # minimum kinetic energy: Emax = 1PeV in GeV

    # Power-law distribution of the cosmic rays (GeV^-1 cm^-3)
p0 = 10             # normalization constant (GeV/c)
alpha = 2.0         # exponent of the power-law distribution

    # Diffusion coefficient of the cosmic rays (cm^2 s^-1)
delta = 1.0/2       # exponent of the power-law of the diffusion coefficient
D0 = 1e28           # diffusion coefficient at 10 GeV/c in cm^2 s^-1  ==> prendre *10 et /10

##==============##
# Gamma emission #
##==============##

    # Energy of the gamma photons (GeV)
Emin_gamma = 100 * MeV2GeV      # 100 MeV (GeV)
Emax_gamma = 100 * TeV2GeV      # 100 TeV (GeV)

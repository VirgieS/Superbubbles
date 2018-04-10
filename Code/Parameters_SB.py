"""
Here are all the free parameters for the system!
"""

import numpy

# Parameters for the system

    # OB asociation
Nob = 88            # number of OB-stars in the association
Lob = 1e36          # mean luminosity of an OB-star (erg/s)
Pob = Nob * Lob     # mean power of the association (erg/s)
lifetime = 4e6 #30e6     # average lifetime of the lowest B star (yr)

    # Fit parameters of the model
n0 = numpy.array([1, 30, 100])              # mean density of the interstellar medium (particle/cm^3)
xH = 0.9            # mass fraction of hydrogen
xHe = 0.1           # mass fraction of helium
mu = xH + xHe * 4   # average molecular weight

    # R(t) = a * n0^alphar * L36^betar * t6^gammar                  (pc)
    # Equation 51 of Weaver et al. (1977)
ar = 27.0
alphar = -1.0/5
betar = 1.0/5
gammar = 3.0/5

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
TISM = 100                          # temperature of the ambient gas (K)

    # Released energy by the acceleration of particles for one SN (erg)
eta = 0.1           # efficiency of the cosmic rays acceleration
Esn = 1e51          # total kinetic energy released from the SN explosion (erg)

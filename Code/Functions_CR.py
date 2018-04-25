"""
Here are all functions needed to compute the CR
"""

##----------##
# Librairies #
##----------##
import matplotlib.pyplot as plt
import numpy
import scipy.integrate as integrate
from scipy.special import erf, erfc

##-----------------------------------------##
# Physical constants and Conversion factors #
##-----------------------------------------##
from Physical_constants import *
from Conversion_factors import *
from Parameters_system import *

##---------##
# Functions #
##---------##

def power_law_distribution(E):
    """
    Return the power-law distribution of the CR particles
    Inputs:
        E       :       energy array (GeV)
    Output:
        NE      :       intial particles distribution (GeV^-1)
    """
    mpgev = mp * MeV2GeV # mass of the proton in GeV

    integral_E = integrate.quad(lambda E: (E**2 + 2 * mpgev * E)**(-(1 + alpha)/2.0) * (E + mpgev) * E, Emin_CR, Emax_CR)[0]
    N0 = eta * Esng * cl**(1-alpha) * p0**(-alpha) * 1.0/integral_E         # normalization constant (GeV^-1 c)

    return N0/cl**(1 - alpha) * (E**2 + 2 * mpgev * E)**(-(1 + alpha)/2.0) * (E + mpgev)/p0**(-alpha)    # GeV^-1

def diffusion_coefficient(E):
    """
    Return the diffusion coefficient with a power-law
    Inputs:
        p0      :       normalization constant (GeV/c)
        D0      :       diffusion coefficient at p0 (cm^2 s^-1)
        E       :       energy array (GeV)
        delta   :       exponent of the power-law expression
    Output:
        D       :       diffusion coefficient (cm^2 s^-1)
    """
    mpgev = mp * MeV2GeV  # mass of the proton in GeV
    return D0 * (numpy.sqrt(E**2 + 2 * mpgev * E)/p0)**delta

def shell_particles(r_in, r_out, NE, D, deltat):
    """
    Return the number of particles in a shell of inner radius r_in and outer radius r_out
    Inputs:
        r_in        :       inner radius of the shell (pc)
        r_out       :       outer radisu of the shell (pc)
        NE          :       initial particles distribution (GeV^-1)
        D           :       diffusion coefficient (cm^2 s^-1)
        deltat      :       time after the SN explosion (yr)
    Output:
        N           :       number of particles in a shell (GeV^-1)
    """

    r_out = r_out * pc2cm   # in cm
    r_in = r_in * pc2cm     # in cm

    if deltat < 0:
        print('problem')

    if deltat < 1e-8 and r_in == 0:
        N = NE

    else:
        deltat = deltat * yr2s            # in s
        a = r_in/(numpy.sqrt(4 * D * deltat))
        b = r_out/(numpy.sqrt(4 * D * deltat))

        N = NE/(numpy.sqrt(numpy.pi)) * (numpy.sqrt(numpy.pi) * (erf(b) - erf(a)) + 2 * numpy.exp(-a**2) * a - 2 * numpy.exp(-b**2) * b)

    return N

def inf_particles(Rsb, NE, D, deltat):
    """
    Return the number of particles outside the superbubble.
    Inputs:
        Rsb         :       outer radius of the SB (pc)
        NE          :       initial particles distribution (GeV^-1)
        D           :       diffusion coefficient (cm^2 s^-1)
        deltat      :       time after the SN explosion (yr)
    Output:
        N           :       number of particles outside the SB (GeV^-1)
    """

    N = numpy.zeros_like(NE)
    if deltat > 1e-8:
        N = numpy.zeros_like(NE)
        a = (Rsb*pc2cm)/(numpy.sqrt(4 * D * deltat * yr2s))
        N = NE/(numpy.sqrt(numpy.pi)) * (numpy.sqrt(numpy.pi) * erfc(a) + 2 * numpy.exp(-a**2) * a)

    return N

def diffusion_spherical(delta_t, r, NE, D):
    """
    Return the density of the cosmic rays at each time and radius
        N(E, r, deltat) = N(E, 0, t0)/(4*pi*D(E)*deltat) * exp(-r^2/(4*D(E)*deltat))
    Inputs:
        delta_t :       time interval since the SN explosion (yr)
        r       :       distance (pc)
        NE      :       initial density of the population of CR (GeV^-1)
        D       :       diffusion coefficient (cm^2 s^-1)
    Output:
        N       :       density of CR (GeV^-1 cm^-3)
    """
    r = r * pc2cm               # cm
    delta_t = delta_t * yr2s    # s

    return NE/((4*numpy.pi * D * delta_t)**(3/2.0))*numpy.exp(-(r**2)/(4 * D * delta_t))

def gauss(r, A, Dt):
    """
    Fit a gaussian on the density profile as function of the radius.
    Input:
        r       :   vector radius (pc)
    Fit parameters:
        A       :   normalization of the gaussian
        Dt      :   factor related to the standard deviation (pc^2)
    Output:
        y = A/(4 * pi * Dt)**(3/2) * exp(- x**2/(4 * Dt)
    """
    return A/((4. * numpy.pi * Dt)**(3/2.0)) * numpy.exp(-(r)**2/(4. * Dt))

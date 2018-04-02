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

##---------##
# Functions #
##---------##

def power_law_distribution(Emin, Emax, E, alpha, eta, Esng, p0):
    """
    Function to compute the power-law distribution of the CR particles
    Inputs:
        Emin    :       minimum energy of the distribution (GeV)
        Emax    :       maximum energy of the distribution (GeV)
        E       :       energy array (GeV)
        alpha   :       exponent of the power-law distribution
        eta     :       efficiency of the cosmic rays acceleration
        Esng    :       energy released by the SN explosion (GeV)
        p0      :       normalization constant (GeV/c)
    """
    mpgev = mp*MeV2GeV  # mass of the proton in GeV

    integral_E = integrate.quad(lambda E: (E**2 + 2*mpgev*E)**(-(1 + alpha)/2.0) * (E + mpgev) * E, Emin, Emax)[0]
    N0 = eta * Esng * cl**(1-alpha) * p0**(-alpha) * 1.0/integral_E         # normalization constant (GeV^-1 c)

    return N0/cl**(1-alpha) * (E**2 + 2*mpgev*E)**(-(1+alpha)/2.0) * (E + mpgev)/p0**(-alpha)    # GeV^-1

def diffusion_coefficient(p0, D0, E, delta):
    """
    Function to compute the diffusion coefficient with a power-law
    Inputs:
        p0      :       normalization constant (GeV/c)
        D0      :       diffusion coefficient at p0 (cm^2 s^-1)
        E       :       energy array (GeV)
        delta   :       exponent of the power-law expression
    """
    mpgev = mp*MeV2GeV  # mass of the proton in GeV
    return D0 * (numpy.sqrt(E**2 + 2*mpgev*E)/p0)**delta

def diffusion_spherical(t, r, t0, NE, D):
    """
    Function to the density of the cosmic rays at each time and radius
        N(E, r, deltat) = N(E, 0, t0)/(4*pi*D(E)*deltat) * exp(-r^2/(4*D(E)*deltat))
    Inputs:
        t       :       time (yr)
        r       :       distance (pc)
        t0      :       time when the SN explode (yr) for the test if the SN has already exploded
        NE      :       initial density of the population of CR (GeV^-1)
        D       :       diffusion coefficient (cm^2 s^-1)
    """

    delta_t = t - t0            # time after the SN explosion (yr)

        # density of the particles in time and position (GeV^-1)
    N = numpy.zeros(len(r))

    if delta_t >= 0:            # if there is already an explosion
        if delta_t == 0:        # at the SN explosion
            print('explosion')
            N[0] = NE
        else:
            N = NE/((4*numpy.pi*D*(delta_t*yr2s))**(3/2.0))*numpy.exp(-(r*pc2cm)**2/(4*D*(delta_t*yr2s)))
    return N

def shell_particles(r_in, r_out, NE, D, deltat):
    """
    Function to compute the number of particles in a shell of inner radius r_in and outer radius r_out
    Inputs:
        r_in        :       inner radius of the shell (pc)
        r_out       :       outer radisu of the shell (pc)
        NE          :       initial particles distribution (GeV^-1)
        D           :       diffusion coefficient (cm^2 s^-1)
        deltat      :       time after the SN explosion (yr)
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
    Function to compute the number of particles outside the superbubble.
    Inputs:
        Rsb         :       outer radius of the SB (pc)
        NE          :       initial particles distribution (GeV^-1)
        D           :       diffusion coefficient (cm^2 s^-1)
        deltat      :       time after the SN explosion (yr)
    """

    N = numpy.zeros_like(NE)
    if deltat > 1e-8:
        N = numpy.zeros_like(NE)
        a = (Rsb*pc2cm)/(numpy.sqrt(4 * D * deltat * yr2s))
        N = NE/(numpy.sqrt(numpy.pi)) * (numpy.sqrt(numpy.pi) * erfc(a) + 2 * numpy.exp(-a**2) * a)

    return N

def gauss(x, A, Dt):
    """
    Function to fit a gaussian on the density profile as function of the radius.
    Input:
        r       :   vector radius (pc)
    Fit parameters:
        A       :   normalization of the gaussian (pc)
        Dt      :   factor related to the standard deviation (pc^2 s^-1 yr)
    """
    return A/((4.*numpy.pi*Dt)**(3/2.0))*numpy.exp(-(x)**2/(4.*Dt))

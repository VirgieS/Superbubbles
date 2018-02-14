"""
Here are all functions needed for the superbubble model
"""

##----------##
# Librairies #
##----------##
import matplotlib.pyplot as plt
import numpy
from math import *
import scipy.integrate as integrate

##-----------------------------------------##
# Physical constants and Conversion factors #
##-----------------------------------------##
from Physical_constants import *
from Conversion_factors import *

##---------##
# Functions #
##---------##

def log_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel):
    """
    Function to plot a log-log graphic
    Inputs:
        figure_number:      define the number of the figure
        number_of_plot:     define how many plot do you want on one figure (with the same axis)
        x:                  x-vector
        y:                  y-array (line = one y-array and row = each different y-plot)
        label_name:         legend of one plots
        title:              title of the plot
        xlabel:             label of the x-axis
        ylabel:             label of the y-axis
    """
    plt.figure(figure_number)
    if number_of_plot > 1:
        for i in range (number_of_plot):
            y_plot = y[i,:]
            plt.plot(x, y_plot, '+', label=label_name[i])
        plt.legend(loc='best')
    else:
        plt.plot(x, y, '+')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    return

def radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6):
    """
    Function to compute the radius and the velocity of the forward shock of a superbubble
    From the equation 52 and 53 of the article Weaver et al. 1977
        Rsb = ar * n0^alphar * L36^betar * t6^gammar                    (pc)
        Vsb = dRsb/dt = ar * gammar * n0^alphar  L36^betar * t6^gammar  (km/s)
    Inputs:
        ar      :       numerical factor
        alphar  :       exponent of the density
        betar   :       exponent of the luminosity
        gammar  :       exponent of the time
        n0      :       atomic density expressed in cm^-3
        L36     :       total luminosity injected by the stars and SN expressed in 10^36 erg/s
        t6      :       age of the system expressed in 10^6 yr
    """

        # Outer radius of the SB (pc)
            # R_sb = ar * n0^alphar * L36^betar * t6^gammar
    Rsb = ar * n0**alphar * L36**betar * t6**gammar

        # Velocity of the forward shock
            # Vsb = dRsb/dt = gammar * ar * n0^alphar * L36^betar * t6^(gammar-1)
    Vsb = ar*(gammar) * n0**alphar * L36**betar * t6**(gammar-1)

    return Rsb, Vsb

def masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb, mu):
    """
    Function to compute the mass in the superbubble and the swept-up mass by the superbubble
    From the equation 5 of the article Mac Low and McCray (1987)
        x = r/Rsb where r is the distance to the center of the superbubble (pc) and Rsb is the outer radius (pc)
        n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan        (cm^-3)
    Inputs:
        an      :       numerical factor
        alphan  :       exponent of the density
        betan   :       exponent of the luminosity
        gamman  :       exponent of the time
        deltan  :       exponent of the radius
        n0      :       atomic density expressed in cm^-3
        L38     :       total luminosity injected by the stars and the SNRs expressed in 10^38 erg/s
        t7      :       age of the system expressed in 10^8 yr
        Rsb     :       outer radius of the superbubble expressed in pc
        mu      :       average molecular weight
    """

        # mass within the superbubble (solar mass)
    integral_msb = integrate.quad(lambda x: (1-x)**deltan * x**2, 0, 1)[0]
    Msb = 4*numpy.pi * (Rsb*pc2cm)**3 * an * n0**alphan * L38**betan * t7**gamman * integral_msb * mu * mpg/Msun2g

        # swept-up mass (solar mass)
    Mswept = 4*numpy.pi/3.0 * (Rsb*pc2cm)**3 * n0 * mu * mpg/Msun2g

    return Msb, Mswept

def density_thickness_shell(mu, n0, Vsb, C02, Mswept, Msb, Rsb):
    """
    Function to compute the density and the thickness of the shell
    From the equation 67 of the article Weaver et al. (1977)
        ns = n0 * (Vsb^2 + C0^2)/Cs^2
        where Cs^2 = kb*Ts/(mu*mpg)
    Inputs:
        mu      :       average molecular weight
        n0      :       atomic density (cm^-3)
        Vsb     :       velocity of the forward shock expressed in km/s
        C02     :       isothermal sound speed in the ambient medium in km/s
        Mswept  :       swept-up mass by the suerbubble expressed in solar masses
        Msb     :       mass in the superbubble expressed in solar masses
        Rsb     :       outer radius of the superbubble expressed in pc
    """
        # density in the shell (cm^-3)
    Ts = 1e4                                    # temperature of the shell (K)
    Cs2 = kb*Ts/(mu*mpg)/(km2cm)**2             # isothermal sound speed in the shell (km/s)
    ns = n0 * (Vsb**2 + C02)/Cs2                # density in the shell (cm^-3)

        # thickness of the shell (pc)
    Ms = (Mswept-Msb)*Msun2g                                    # mass in the shell (g)
    Vs = Ms/(ns * mu * mpg)                                     # volume of the shell (cm^3)
    Rc = exp(log(Rsb**3 - 3/(4*numpy.pi)*Vs/(pc2cm**3))*1/3.0)  # radius of the continuity contact (pc)
    hs = Rsb-Rc                                                 # thickness of the shell (pc)

    return ns, hs

def pressure_SB(ap, alphap, betap, gammap, n0, L36, t6):
    """
    Function to compute the pressure in the superbubble
    From the equation 22 of the article Weaver et al. (1977)
        p = ap * n0^alphap * L36^betap * t6^gammap
    Inputs:
        ap      :       numerical factor
        alphap  :       exponent of the density
        betap   :       exponent of the luminosity
        gammap  :       exponent of the time
        n0      :       atomic density expressed in cm^-3
        L36     :       total luminosity injected by the stars and SN expressed in 10^36 erg/s
        t6      :       age of the system expressed in 10^6 yr
    """
    return ap * n0**(alphap) * L36**(betap) * t6**(gammap)

def luminosity_SB(al, etal, zeta, at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, epsilon, Rsb):
    """
    Function to compute the pressure in the superbubble
    From the equation 6 of the article Mac Low and McCray (1987)
        C = ne * n * LambdaT                (erg cm^-3 s^-1)
        where LambdaT = al * T6^etal * zeta (erg cm^3 s^-1)
              ne = epsilon * n (cm^-3)
        Lsb = int_Volumesb (C)          (erg/s)
    Inputs:
        # expression of LambdaT = al * T6^etal * zeta
        al      :       numerical factor in the expression of LambdaT
        etal    :       exponent of the temperature in the expression of LambdaT
        zeta    :       metallicity of the local ISM (in the expression of LambdaT)
        # expression of T = at * n0^alphat * L38^betat * t7^deltat * (1-x)^deltat
        at      :       numerical factor in the expression of T
        alphat  :       exponent of the density in the expression of T
        betat   :       exponent of the luminosity in the expression of T
        gammat  :       exponent of the time in the expression of T
        deltat  :       exponent of the radius in the expression of T
        # expression of n = an * n0^alphan * L38^betan * t7^deltan * (1-x)^deltan
        an      :       numerical factor in the expression of n
        alphan  :       exponent of the density in the expression of n
        betan   :       exponent of the luminosity in the expression of n
        gamman  :       exponent of the time in the expression of n
        deltan  :       exponent of the radius in the expression of n
        # others parameters
        n0      :       atomic density expressed in cm^-3
        L38     :       total luminosity injected by the stars and SN expressed in 10^38 erg/s
        t7      :       age of the system expressed in 10^7 yr
        epsilon :       ratio ne/n
        Rsb     :       outer radius of the superbubble (pc)
    """
    at6 = at * K26K
    deltax = 2 * deltan + etal * deltat
    integral_lsb = integrate.quad(lambda x: (1-x)**deltax * x**2, 1, 0)[0]
    return al * zeta * (at6 * n0**alphat * L38**betat * t7**gammat)**(etal) * epsilon * (an * n0**alphan * L38**betan * t7**gamman)**2 * integral_lsb * (Rsb*pc2cm)**3 # in erg/s

def profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb):
    """
    Function to compute the density and the temperature profiles in the superbubble
    From equations 4 and 5 of the article of Mac Low and McCray (1987)
        x = r/Rsb where r is the distance to the center of the SB (pc) and Rsb is the outer radius (pc)
        T(x) = at * n0^alphat * L38^betat * t7^gammat * (1-x)^deltat    (K)
        n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan    (cm^-3)
    Inputs:
        # expression of T = at * n0^alphat * L38^betat * t7^deltat * (1-x)^deltat
        at      :       numerical factor in the expression of T
        alphat  :       exponent of the density in the expression of T
        betat   :       exponent of the luminosity in the expression of T
        gammat  :       exponent of the time in the expression of T
        deltat  :       exponent of the radius in the expression of T
        # expression of n = an * n0^alphan * L38^betan * t7^deltan * (1-x)^deltan
        an      :       numerical factor in the expression of n
        alphan  :       exponent of the density in the expression of n
        betan   :       exponent of the luminosity in the expression of n
        gamman  :       exponent of the time in the expression of n
        deltan  :       exponent of the radius in the expression of n
        # others parameters
        n0      :       atomic density expressed in cm^-3
        L38     :       total luminosity injected by the stars and SN expressed in 10^38 erg/s
        t7      :       age of the system expressed in 10^7 yr
        Rsb     :       outer radius of the SB (pc)
    """
    r = numpy.linspace(0, Rsb, 100)
    x = r/Rsb
    Tsb = at * n0**alphat * L38**betat * t7**gammat * (1-x)**deltat
    nsb = an * n0**alphan * L38**betan * t7**gamman * (1-x)**deltan

    return Tsb, nsb

def diffusion_spherical(t, Rsb, t0, NE, r0, D):
    """
    This function compute the intensity density of the cosmic rays at each time and radius
    Inputs:
        t       :       time (yr)
        Rsb     :       radius of the SB (pc)
        t0      :       time when the SN explode (yr) for the test if the SN has already exploded
        NE      :       initial density of the population of CR (GeV^-1)
        r0      :       position of the SN explosion (pc)
        D       :       diffusion coefficient (cm^2 s^-1)
    """
    rmin = 0              # minimum radius (pc)
    rmax = Rsb            # maximum radius (pc)
    number_bin_r = 20.0
    r = numpy.linspace(rmin, rmax, number_bin_r)    # position in pc
    dr = (rmax - rmin)/number_bin_r     # in pc

        # density of the particles in time and position (GeV^-1)
    N = numpy.zeros(len(r))

    if t >= t0:
        if t == t0:
            ind = numpy.where((r < r0 + dr/2.0) & (r > r0 - dr/2.0))[0]
            N[ind] = NE
        else:
            N = NE/((4*numpy.pi*D*(t-t0)*yr2s)**(3/2.0))*numpy.exp(-((r-r0)*pc2cm)**2/(4*D*(t-t0)*yr2s))
    return N, r

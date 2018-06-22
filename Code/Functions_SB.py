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

## ------- ##
# Functions #
## ------- ##

def radius_velocity_SB(t6):
    """
    Returns the radius and the velocity of the forward shock of a superbubble
    From the equation 52 and 53 of the article Weaver et al. 1977
        Rsb = ar * n0^alphar * L36^betar * t6^gammar                    (pc)
        Vsb = dRsb/dt = av * n0^alphav * L36^betav * t6^gammav          (km/s)
    Inputs:
        t6      :       age of the system expressed in 10^6 yr
    Outputs:
        Rsb     :       outer radius of the SB (pc)
        Vsb     :       velocity of the forward shock (km/s)
    """

        # Outer radius of the SB (pc)
            # R_sb = ar * n0^alphar * L36^betar * t6^gammar
    Rsb = ar * n0**alphar * L36**betar * t6**gammar

        # Velocity of the forward shock
            # Vsb = dRsb/dt = av * n0^alphav * L36^betav * t6^gammav
    Vsb = av * n0**alphav * L36**betav * t6**gammav

    return Rsb, Vsb

def masses(t7, Rsb):
    """
    Returns the mass in the superbubble and the swept-up mass by the superbubble
    From the equation 5 of the article Mac Low and McCray (1987)
        x = r/Rsb where r is the distance to the center of the superbubble (pc) and Rsb is the outer radius (pc)
        n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan        (cm^-3)
    Inputs:
        t7      :       age of the system expressed in 10^8 yr
        Rsb     :       outer radius of the superbubble expressed in pc
    Outputs:
        Msb     :       mass inside the superbubble (solar mass)
        Mswept  :       swept-up mass by the growth of the SB (solar mass)
    """

        # mass within the superbubble (solar mass)
        # integral over x of the density in the SB
    integral_msb = integrate.quad(lambda x: (1-x)**deltan * x**2, 0, 1)[0]
    Msb = 4*numpy.pi * (Rsb*pc2cm)**3 * an * n0**alphan * L38**betan * t7**gamman * integral_msb * mu * mpg/Msun2g

        # swept-up mass (solar mass)
        # Mswept = Volume_sb * n0 * mu * mpg
    Mswept = 4*numpy.pi/3.0 * (Rsb*pc2cm)**3 * n0 * mu * mpg/Msun2g

    return Msb, Mswept

def density_thickness_shell_percentage(percentage, Rsb, Mswept, Msb):
    """
    Returns the density and the thickness of the shell with the correction of the Weaver's model

    Inputs:
        percentage      :   thickness of the shell is a percentage of the outer radius
        Rsb             :   outer radius of the SB (pc)
        Mswept          :   swept-up mass by the SB (solar masses)
        Msb             :   mass in the SB (solar masses)
    Outputs:
        ns              :   density in the shell (cm^-3)
        hs              :   thickness of the shell (pc)
    """
        # thickness (pc)
        # percentage of the outer radius
    hs = percentage * Rsb
    Rc = Rsb - hs

        # density (cm^-3)
        # from the corrected volume of the supershell
    Vs = 4 * numpy.pi/3.0 * ((Rsb * pc2cm)**3 - (Rc * pc2cm)**3)
    Ms = (Mswept - Msb) * Msun2g
    ns = Ms/(Vs * mu * mpg)

    return ns, hs

def density_thickness_shell(Vsb, Mswept, Msb, Rsb):
    """
    Returns the density and the thickness of the shell
    From the equation 67 of the article Weaver et al. (1977)
        ns = n0 * (Vsb^2 + C0^2)/Cs^2
        where Cs^2 = kb*Ts/(mu*mpg)
    Inputs:
        Vsb     :       velocity of the forward shock (km/s)
        C02     :       isothermal sound speed in the ambient medium (km/s)
        Mswept  :       swept-up mass by the suerbubble (solar masses)
        Msb     :       mass in the superbubble (solar masses)
        Rsb     :       outer radius of the superbubble (pc)
    Outputs:
        ns      :       density in the shell (cm^-3)
        hs      :       thickness in the shell (cm^-3)
    """
        # density in the shell (cm^-3)
    Cs2 = kb*Ts/(mu*mpg)/(km2cm)**2             # isothermal sound speed in the shell (km2/s2)
    ns = n0 * (Vsb**2 + C02)/Cs2                # density in the shell (cm^-3)

        # thickness of the shell (pc)
    Ms = (Mswept-Msb)*Msun2g                    # mass in the shell (g)
    Vs = Ms/(ns * mu * mpg)*1.0/(pc2cm)**3      # volume of the shell (pc^3)
    Rc = (Rsb**3 - 3.0/(4*numpy.pi)*Vs)**(1.0/3)# radius of the continuity contact (pc)
    hs = Rsb-Rc                                 # thickness of the shell (pc)

    return ns, hs

def pressure_SB(t6):
    """
    Returns the pressure in the superbubble
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
    Output:
        p       :       pressure in the SB (dyne cm^-2)
    """
    return ap * n0**(alphap) * L36**(betap) * t6**(gammap)

def luminosity_SB(t7, Rsb):
    """
    Returns the luminosity in the superbubble from the cooling rate.
    From the equation 6 of the article Mac Low and McCray (1987)
        C = ne * n * LambdaT                :   cooling rate per unit volume (erg cm^-3 s^-1)
        where LambdaT = al * T6^etal * zeta (erg cm^3 s^-1)
              ne = epsilon * n              : density of electrons (cm^-3)
        Lsb = int_Volumesb (C)          (erg/s)
    Inputs:
        t7      :       age of the system expressed in 10^7 yr
        Rsb     :       outer radius of the superbubble (pc)
    Output:
        Lsb     :       luminosity of the SB (erg s^-1)
    """
    at6 = at * K26K
    deltax = 2 * deltan + etal * deltat
    integral_lsb = integrate.quad(lambda x: (1-x)**deltax * x**2, 1, 0)[0]
    return al * zeta * (at6 * n0**alphat * L38**betat * t7**gammat)**(etal) * epsilon * (an * n0**alphan * L38**betan * t7**gamman)**2 * integral_lsb * (Rsb*pc2cm)**3 # in erg/s

def profile_density_temperature(t7, rsb, Rsb):
    """
    Returns the density and the temperature profiles in the superbubble
    From equations 4 and 5 of the article of Mac Low and McCray (1987)
        x = r/Rsb where r is the distance to the center of the SB (pc) and Rsb is the outer radius (pc)
        T(x) = at * n0^alphat * L38^betat * t7^gammat * (1-x)^deltat    (K)
        n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan    (cm^-3)
    Inputs:
        t7      :       age of the system expressed in 10^7 yr
        rsb     :       array of distance (pc)
        Rsb     :       outer radius of the SB (pc)
    Outputs:
        Tsb     :       temperature profile in the SB (K)
        nsb     :       density profile in the SB (cm^-3)
    """
    x = rsb/Rsb

        # temperature profile (K)
    Tsb = at * n0**alphat * L38**betat * t7**gammat * (1-x)**deltat

        # density profile (cm^-3)
    nsb = an * n0**alphan * L38**betan * t7**gamman * (1-x)**deltan

    return Tsb, nsb

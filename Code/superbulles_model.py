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

    # OB asociation
Nob = 100           # number of OB-stars in the association
Lob = 1e36          # mean luminosity of an OB-star (erg/s)
Pob = Nob*Lob       # mean power of the association (erg/s)

    # Fit parameters of the model
n0 = 1              # mean density of the interstellar medium (particle/cm^3)
xH = 0.9            # mass fraction of hydrogen
xHe = 0.1           # mass fraction of helium
mu = xH + xHe * 4   # average molecular weight

    # R(t) = a * n0^alphar * L36^betar * t6^gammar [pc]
    # L36 expressed in 10^36 erg/s, t6 in 10^6 yr and n0 in cm^-3
ar = 27.0
alphar = -1.0/5
betar = 1.0/5
gammar = 3.0/5

    # T(x) = at * n0^alphat * L38^betat * t7^gammat * (1-x)^deltat
    # L38 expressed in 10^38 erg/s, t7 in 10^7 yr and n0 in cm^-3
at = 3.6e6
alphat = 2.0/35
betat = 8.0/35
gammat = -6/35
deltat = 2.0/5

    # n(x) = an * n0^alphan * L38^betan * t7^gamman * (1-x)^deltan
    # L38 expressed in 10^38 erg/s, t7 in 10^7 yr and n0 in cm^-3
an = 4.0e-3
alphan = 19.0/35
betan = 6.0/35
gamman = -22.0/35
deltan = -2.0/5

    #  p(t) = ap * n0^alphap * L_36^betap * t6^gammap
    # L36 expressed in 10^36 erg/s, t6 in 10^6 yr and n0 in cm^-3
ap = 4.12e-12
alphap = 3.0/5
betap = 2.0/5
gammap = -4.0/5

# Computation of the radius of the superbubble during the time

    # time vector (must be given in yr)
t = []
t.append(0.0)   # birth of the association (yr)
dt = 1e3        # time interval (yr)

    # Important quantities
L36 = Pob/erg236erg     # mechanical energy expressed in 10^36 erg/s
L38 = L36/t36erg238erg  # mechanical energy expressed in 10^38 erg/s

    # Quantities that we want to compute during the time
Rsb = []        # outer radius of the superbubble (pc)
Vsb = []        # velocity of the forward shock (km/s)
# nsb = []        # density profile in the interior (cm^-3)
# Tsb = []        # temperature profile in the interior (K)
Msb = []        # mass within the superbubble (solar mass) fromt the integral of the density profile
Msbth = []      # mass within the superbubble (solar mass) from the mass loss rate
Mswept = []     # swept mass (solar mass)
Mratio = []     # ratio between the mass within the SB and the swept mass
h = []          # thickness of the shell (pc)
ns = []         # density in the shell (cm^-3)
hratio = []     # ratio between the thickness of the shell and the radius of the SB
psb = []        # pressure into the superbubble (dyne cm^-2)
pratio = []     # ratio between the pressure in the SB and the pressure in the ISM
Lsb = []        # luminosity of the superbubble (erg/s)

    # In the ISM
TISM = 100                          # temperature of the ambient gas (K)
pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

    # radiative cooling time (yr)
tc = 2.3e4 * n0**(-0.71) * L38**0.29

        # Initialization
Rsb.append(0)
Vsb.append(0)
#r = 0           # position from the center (pc)
#x = 0           # reduction variable: x = r/Rsb
# Tsb.append(0)
# nsb.append(0)
Msb.append(0)
Mswept.append(0)
Mratio.append(0)
h.append(0)
ns.append(0)
hratio.append(0)
psb.append(0)
pratio.append(0)
Lsb.append(0)

# Stop criterion
i = 0
lifetime = 30e6 # average lifetime of the lowest mass B star (yr)
#hgalac = 150 # tickness of the galactic plane
#while Rsb[i] < hgalac:
while t[i] < lifetime:

    t.append(t[i]+dt)
    i +=1
    t6 = t[i]/yr26yr
    t7 = t6/s6yr27yr

        # outer radius (pc) and velocity (km/s) of the forward shock
    Rsb.append(ar * n0**alphar * L36**betar * t6**gammar)               # in pc
    Vsb.append(ar*(gammar) * n0**alphar * L36**betar * t6**(gammar-1))  # in km/s

        # temperature and density profiles within the superbubble
    #r = numpy.linspace(0, Rsb[i], 100)
    #x = r/Rsb[i]

    #Tsb = at * n0**alphat * L38**betat * t7**gammat * (1-x)**deltat
    #nsb = an * n0**alphan * L38**betan * t7**gamman * (1-x)**deltan

        # mass within the superbubble (solar mass)
    integral_msb = integrate.quad(lambda x: (1-x)**deltan * x**2, 0, 1)[0]
    Msb.append(4*numpy.pi * (Rsb[i]*pc2cm)**3 * an * n0**alphan * L38**betan * t7**gamman * integral_msb * mu * mpg/Msun2g) # in solar mass

        # swept mass (solar mass)
    Mswept.append(4*numpy.pi/3.0 * (Rsb[i]*pc2cm)**3 * n0 * mu * mpg/Msun2g)    # in solar mass

        # ratio of mass
    Mratio.append(Msb[i]/Mswept[i])

        # density in the shell (cm^-3)
    Ts = 3*mu*mpg*(Vsb[i]*km2cm)**2/(16.0*kb)   # temperature of the ambiant gas (K)
    Ts = 1e4                                    # temperature of the shell (K)
    Cs2 = kb*Ts/(mu*mpg)/(km2cm)**2             # isothermal sound speed in the shell (km/s)
    ns.append(n0 * (Vsb[i]**2 + C02)/Cs2)       # density in the shell (cm^-3)

        # Volume of the shell (cm^3)
    Ms = (Mswept[i]-Msb[i])*Msun2g              # mass in the shell (g)
    Vs = Ms/(ns[i] * mu * mpg)                  # volume of the shell (cm^3)

        # Solution
    Rc = exp(log(Rsb[i]**3 - 3/(4*numpy.pi)*Vs/(pc2cm**3))/3)
    h.append(Rsb[i]-Rc)
    hratio.append(h[i]/Rsb[i])

        # pressure in the superbubble (dyne cm^-2)
    psb.append(ap * n0**(alphap) * L36**(betap) * t6**(gammap))
    pratio.append(psb[i]/pISM)

        # luminosity of the SB (erg/s)
    deltax = 2 * deltan -0.7 * deltat
    integral_lsb = integrate.quad(lambda x: (1-x)**deltax * x**2, 1, 0)[0]
    Lsb.append(1e-22 * (at/1e6 * n0**alphat * L38**betat * t7**gammat)**(-0.7) * 1.2 * (an * n0**alphan * L38**betan * t7**gamman)**2 * integral_lsb * (Rsb[i]*pc2cm)**3) # in erg/s

    if ((t[i]<tc + dt) and (t[i] > tc - dt)):
        print(Lsb[i]/(Pob))

pISM = np.ones_like(psb)*pISM
Pob = np.ones_like(Lsb)*Pob
tc = 2.3e4 * n0**(-0.71) * L38**0.29 # the radiative cooling time scale (yr)

# Plots
log_plot(1, 1, t, Rsb, 'none', 'Time evolution of the radius', 'Time (yr)', 'Radius (pc)')
log_plot(2, 1, t, Vsb, 'none', 'Time evolution of the velocity', 'Time (yr)', 'Velocity (km/s)')
log_plot(3, 2, t, np.array([Msb, Mswept]), ['Mass within the SB from the density profile', 'Swept-up mass'], 'Time evolution of the masses', 'Time (yr)', u'Mass'r'($M_\odot$)')
log_plot(4, 1, t, Mratio, 'none', 'Time evolution of the ratio between the two masses', 'Time (yr)', r'$M_{SB}/M_{swept-up}$')
log_plot(5, 1, t, ns, 'none', 'Time evolution of the density in the shell', 'Time (yr)', r'$n_s$($cm^{-3}$)')
log_plot(6, 1, t, h, 'none', 'Time evolution of the thickness of the shell', 'Time (yr)', 'h (pc)')
log_plot(7, 1, t, hratio, 'none', 'Time evolution of the ratio between the thickness of the shell and the radius of the SB', 'Time (yr)', r'$h_s/R_{SB}$')
log_plot(8, 2, t, np.array([psb, pISM]), ['within the SB', 'in the ISM'], 'Time evolution of the pressures', 'Time (yr)', u'Pressure 'r'(dyne $cm^{-2}$)')
log_plot(9, 1, t, pratio, 'none', 'Time evolution of the ratio between the two pressures', 'Time (yr)', r'$p_{SB}/p_{ISM}$')
log_plot(10, 2, t, np.array([Lsb, Pob]), ['of the SB', 'of the stars'], 'Time evolution of the luminosities', 'Time (yr)', u'L 'r'(erg $s^{-1}$)')
plt.show()

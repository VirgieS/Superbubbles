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
Lob = 1E36          # mean luminosity of an OB-star (erg/s)
Pob = Nob*Lob       # mean power of the association (erg/s)

    # Fit parameters of the model
n0 = 1              # mean density of the interstellar medium (particle/cm^3)
xH = 0.9            # mass fraction of hydrogen
xHe = 0.1           # mass fraction of helium
mu = xH + xHe * 4   # average molecular weight

    # R(t) = a * n0^alphar * L36^betar * t6^gammar [pc]
    # L36 expressed in 10^36 erg/s, t6 in 10^6 yr and n0 in cm^-3
ar = 27
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

    # n(x) = an * n_0^alphan * L38^betan * t7^gamman * (1-x)^deltan
    # L38 expressed in 10^38 erg/s, t7 in 10^7 yr and n0 in cm^-3
an = 4.0e-3
alphan = 19.0/35
betan = 6.0/35
gamman = -22.0/35
deltan = -2.0/5

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
nsb = []        # density profile in the interior (cm^-3)
Tsb = []        # temperature profile in the interior (K)
Msb = []        # mass within the superbubble (solar mass) fromt the integral of the density profile
Msbth = []      # mass within the superbubble (solar mass) from the mass loss rate
Mswept = []     # swept mass (solar mass)
Mratio = []     # ratio between the mass within the SB and the swept mass
h = []          # thickness of the shell (pc)
ns = []         # density in the shell (cm^-3)
hratio = []     # ratio between the thickness of the shell and the radius of the SB
psb = []        # pressure into the superbubble (dyne cm^-2)

    # Int the ISM
TISM = 100              # temperature of the ambient gas (K)
pISM = n0 * kb * TISM   # pressure in the ISM (dyne cm^-2)

        # Initialization
Rsb.append(0)
Vsb.append(0)
r = 0           # position from the center (pc)
x = 0           # reduction variable: x = r/Rsb
Tsb.append(0)
nsb.append(0)
Msb.append(0)
Mswept.append(0)
Mratio.append(0)
h.append(0)
ns.append(0)
hratio.append(0)

# Stop criterion
i = 0
lifetime = 30e6 # average lifetime of the lowest mass B star (yr)
#hgalac = 150 # tickness of the galactic plane
#while Rsb[i] < hgalac:
while t[i] < lifetime:

    t.append(t[i]+dt)
    i +=1
    t6 = t[i]/yr26yr

        # outer radius and velocity of the forward shock
    Rsb.append(ar * n0**alphar * L36**betar * t6**gammar)               # in pc
    Vsb.append(ar*(gammar) * n0**alphar * L36**betar * t6**(gammar-1))  # in km/s

        # temperature and density profiles within the superbubble
    #r = numpy.linspace(0, Rsb[i], 100)
    #x = r/Rsb[i]

    #Tsb = at * n0**alphat * L38**betat * (t6/s6yr27yr)**gammat * (1-x)**deltat
    #nsb = an * n0**alphan * L38**betan * (t6/s6yr27yr)**gamman * (1-x)**deltan

        # mass within the superbubble
    integral_sb = integrate.quad(lambda x: (1-x)**deltan * x**2, 0, 1)[0]
    Msb.append(4*numpy.pi * (Rsb[i]*pc2cm)**3 * an * n0**alphan * L38**betan * (t6/s6yr27yr)**gamman * integral_sb * mu * mpg/Msun2g) # in solar mass

        # swept mass
    Mswept.append(4*numpy.pi/3.0 * (Rsb[i]*pc2cm)**3 * n0 * mu * mpg/Msun2g)    # in solar mass

        # ratio of mass
    Mratio.append(Msb[i]/Mswept[i])

    # Computation of the thickness of the shell

        # density in the shell
    Ts = 3*mu*mpg*(Vsb[i]*km2cm)**2/(16.0*kb)   # temperature of the ambiant gas (K)
    C02 = kb*TISM/(mu*mpg)/(km2cm)**2           # isothermal sound speed in the ambiant gas (km/s)
    Ts = 1e4                                    # temperature of the shell (K)
    Cs2 = kb*Ts/(mu*mpg)/(km2cm)**2             # isothermal sound speed in the shell (km/s)
    ns.append(n0 * (Vsb[i]**2 + C02)/Cs2)       # density in the shell (cm^-3)

        # Volume of the shell
    Ms = (Mswept[i]-Msb[i])*Msun2g              # mass in the shell (g)
    Vs = Ms/(ns[i] * mu * mpg)                  # volume of the shell (cm^3)

        # Solution
    Rc = exp(log(Rsb[i]**3 - 3/(4*numpy.pi)*Vs/(pc2cm**3))/3)
    h.append(Rsb[i]-Rc)
    hratio.append(h[i]/Rsb[i])

        # pressure in the superbubble
    p.append(7/(3850 * numpy.pi)**(2.0/5) * Pob**(2/5) * (n0 * mu * mpg)**(3/5) * t[i]**(-4/5))

# Plots
log_plot(1, 1, t, Rsb, 'none', 'Time evolution of the radius', 'Time (yr)', 'Radius (pc)')
log_plot(2, 1, t, Vsb, 'none', 'Time evolution of the velocity', 'Time (yr)', 'Velocity (km/s)')
log_plot(3, 2, t, np.array([Msb, Mswept]), ['Mass within the SB from the density profile', 'Swept-up mass'], 'Time evolution of the masses', 'Time (yr)', u'Mass'r'($M_\odot$)')
log_plot(4, 1, t, Mratio, 'none', 'Time evolution of the ratio between the two masses', 'Time (yr)', r'$M_{SB}/M_{swept-up}$')
log_plot(5, 1, t, ns, 'none', 'Time evolution of the density in the shell', 'Time (yr)', r'$n_s$($cm^{-3}$)')
log_plot(6, 1, t, h, 'none', 'Time evolution of the thickness of the shell', 'Time (yr)', 'h (pc)')
log_plot(7, 1, t, hratio, 'none', 'Time evolution of the ratio between the thickness of the shell and the radius of the SB', 'Time (yr)', r'$h_s/R_{SB}$')
plt.show()

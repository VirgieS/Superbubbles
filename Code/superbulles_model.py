##------------------------##
# Librairies and functions #
##------------------------##
import matplotlib.pyplot as plt
import numpy
from math import *
import scipy.integrate as integrate

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *

##-----------##
# Computation #
##-----------##

# Parameters for the system

    # OB asociation
Nob = 10            # number of OB-stars in the association
Lob = 1E37          # mean luminosity of an OB-star (erg/s)
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
Msb = []        # mass within the superbubble (solar mass)
Mswept = []     # swept mass (solar mass)
Mratio = []     # ratio between the mass within the SB and the swept mass

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
    Rsb.append(27 * n0**alphar * L36**betar * t6**gammar)
    Vsb.append(27*(gammar) * n0**alphar * L36**betar * t6**(gammar-1))

        # temperature and density profiles within the superbubble
    #r = numpy.linspace(0, Rsb[i], 100)
    #x = r/Rsb[i]

    #Tsb = at * n0**alphat * L38**betat * (t6/s6yr27yr)**gammat * (1-x)**deltat
    #nsb = an * n0**alphan * L38**betan * (t6/s6yr27yr)**gamman * (1-x)**deltan

        # mass within the superbubble
    integral = integrate.quad(lambda x: (1-x)**deltan * x**2, 0, 1)[0]
    Msb.append(4*numpy.pi * (Rsb[i]*pc2cm)**3 * an * n0**alphan * L38**betan * (t6/s6yr27yr)**gamman * integral * mu * mpg/Msun2g)

        # swept mass
    Mswept.append(4*numpy.pi/3.0 * (Rsb[i]*pc2cm)**3 * n0 * mu * mpg/Msun2g)
    if t[i] == 1e3:
        M3 = Msb[i]

    if t[i] == 30e6:
        M6 = Msb[i]

        # ratio of mass
    Mratio.append(Msb[i]/Mswept[i])

# Plots

plt.figure(1)
plt.plot(t, Rsb, '+')
plt.title('Evolution of the superbubble radius')
plt.xlabel('Time after the formation of the association (yr)')
plt.ylabel('Radius (pc)')
plt.xscale('log')
plt.yscale('log')

plt.figure(2)
plt.plot(t, Vsb, '+')
plt.title('Evolution of the superbubble velocity')
plt.xlabel('Time after the formation of the association (yr)')
plt.ylabel('Velocity (km/s)')
plt.xscale('log')
plt.yscale('log')

plt.figure(3)
plt.plot(t, Msb, '+', label='Mass within the SB')
plt.plot(t, Mswept, '+', label='Swept mass')
plt.title('Evolution of the mass')
plt.xlabel('Time after the formation of the association (yr)')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(u'Mass'r'($M_\odot$)')
plt.legend(loc='best')

plt.figure(4)
plt.plot(t, Mratio, '+')
plt.title('Ratio between the swept mass and the mass within the SB')
plt.xlabel('Time after the formation of the association (yr)')
plt.ylabel(r'$M_{SB}/M_{swept}$')
plt.xscale('log')
plt.yscale('log')

plt.show()

print(M6/M3)

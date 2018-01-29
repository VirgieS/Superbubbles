# Librairies and functions
import matplotlib.pyplot as plt
import numpy
from math import *
from Physical_constants import *
from Conversion_factors import *

# Parameters for the system

    # OB asociation
Nob = 100            # number of OB-stars in the association
Lob = 1E37          # mean luminosity of an OB-star (erg/s)
Pob = Nob*Lob       # mean power of the association (erg/s)

    # Fit parameters of the model
    # L must be expressed in 10^36 erg, t in 10^6 yr and n0 in cm^-3
n0 = 1e2            # mean density of the interstellar medium (particle/cm^3)

    # R(t) = a * n0^alpha * L^beta * t^gamma [pc]
ar = 27
alphar = -1.0/5
betar = 1.0/5
gammar = 3.0/5

# Computation of the radius of the superbubble during the time
# example
t = []
t.append(0.0)       # birth of the association (yr)
deltat = 1e3        # time interval (yr)
L36 = Pob/1E36
t6 = t[0]/1E6
Rsb = []
Vsb = []

Rsb.append(0) # radius of the superbubble (pc)
Vsb.append(0) # velocity of the forward shock (km/s)


# Stop criterion
i = 0
hgalac = 150 # tickness of the galactic plane
while Rsb[i] < hgalac:
    t.append(t[i]+deltat)
    i +=1
    t6 = t[i]/1E6
    Rsb.append(27 * n0**alphar * L36**betar * t6**gammar)
    Vsb.append(27*(gammar) * n0**alphar * L36**betar * t6**(gammar-1))

# Plots

plt.figure(1)
plt.plot(t, Rsb, '+')
plt.title('Evolution of the superbubble radius')
plt.xlabel('Time after the formation of the association (yr)')
plt.ylabel('Radius (pc)')
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xscale('log')

plt.figure(2)
plt.plot(t, Vsb, '+')
plt.title('Evolution of the superbubble velocity')
plt.xlabel('Time after the formation of the association (yr)')
plt.ylabel('Velocity (km/s)')
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xscale('log')

plt.show()

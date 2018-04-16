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
D = 1       # difusion coefficient
N0 = 5
z0 = 2
t0 = 0

    # time (s)
tmin = 0.1            # intial time (s)
tmax = 1e2            # final time (s)
number_bin_t = 15.0
dt = (log10(tmax)-log10(tmin))/number_bin_t
t = numpy.logspace(log10(tmin), log10(tmax), number_bin_t))  # s
print(t)

    # position
number_bin_z = 200.0
zmin = 0            # initial position
zmax = 100          # final position
z = numpy.linspace(zmin, zmax, number_bin_z)
N = numpy.zeros((len(t), len(z)))

for i in range(len(t)):
    if t[i] < t0:
        continue
    elif ((t[i]>t0 - dt/2.0) and (t[i] < t0 +dt/2.0)):
        N[i,z0] = N0
        plt.plot(z, N[i, :], '+', label="at t=%.2f s" %t[i])
    else:
        N[i,:] = N0/(numpy.sqrt(4*numpy.pi*(t[i]-t0)))*numpy.exp(-(z[:]-z0)**2/(4*D*(t[i]-t0)))
        plt.plot(z, N[i, :], '+', label="at t=%.2f s" %t[i])

plt.title('Diffusion')
plt.xlabel('position')
plt.ylabel('N')
plt.legend(loc='best')
plt.show()

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
from Parameters_SB import *

##-----------##
# Computation #
##-----------##

# Computation of the radius of the superbubble during the time

    # Important quantities
L36 = Pob * erg236erg     # mechanical energy expressed in 10^36 erg/s
L38 = L36 * t36erg238erg  # mechanical energy expressed in 10^38 erg/s

    # Quantities that we want to compute during the time
Rsb = []        # outer radius of the superbubble (pc)
Vsb = []        # velocity of the forward shock (km/s)
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
#TISM = 100                          # temperature of the ambient gas (K)
pISM = n0 * kb * TISM               # pressure in the ISM (dyne cm^-2)
C02 = kb*TISM/(mu*mpg)/(km2cm)**2   # isothermal sound speed in the ambiant gas (km/s)

    # radiative cooling time (yr)
#tc = 2.3e4 * n0**(-0.71) * L38**0.29

        # Initialization
Rsb.append(0)
Vsb.append(0)
Msb.append(0)
Mswept.append(0)
Mratio.append(0)
h.append(0)
ns.append(0)
hratio.append(0)
psb.append(0)
pratio.append(0)
Lsb.append(0)

# Initialization
t0 = 0.01          # creation of the SB (yr)
tmin = t0
tmax = lifetime
number_bin_t = 30
t = numpy.logspace(log10(tmin), log10(tmax), number_bin_t)

figure_number = 1
#hgalac = 150 # tickness of the galactic plane
#while Rsb[i] < hgalac:
for i in range (1,len(t)):

    t6 = t[i] * yr26yr
    t7 = t6 * s6yr27yr

        # outer radius (pc) and velocity (km/s) of the forward shock
    R_sb, V_sb = radius_velocity_SB(ar, alphar, betar, gammar, n0, L36, t6)
    Rsb.append(R_sb)  # in pc
    Vsb.append(V_sb)  # in km/s

        # mass within the superbubble and swept-up mass (solar mass)
    M_sb, M_swept = masses(an, alphan, betan, gamman, deltan, n0, L38, t7, Rsb[i], mu)
    Msb.append(M_sb)
    Mswept.append(M_swept)

        # ratio of masses
    Mratio.append(Msb[i]/Mswept[i])

        # density (cm^-3) and thickness (pc) of the shell
    n_s, h_s = density_thickness_shell(mu, n0, Vsb[i], C02, Mswept[i], Msb[i], Rsb[i])
    ns.append(n_s)  # in cm^-3
    h.append(h_s)   # in pc

        # ratio of the thickness of the shell and the outer radius
    hratio.append(h[i]/Rsb[i])

        # pressure in the superbubble (dyne cm^-2)
    psb.append(pressure_SB(ap, alphap, betap, gammap, n0, L36, t6))
    pratio.append(psb[i]/pISM)

        # luminosity of the SB (erg/s)
    Lsb.append(luminosity_SB(al, etal, zeta, at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, epsilon, Rsb[i])) # erg/s
    """
    number_bin_r = 20
    if i == 25:
        rmin = 0                    # minimum radius (pc)
        rmax = Rsb[i]-h[i]          # maximum radius (pc)
        r = numpy.linspace(rmin, rmax, number_bin_r)    # position in pc
        Tsb, nsb = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, r)
        plot(1, 1, r/Rsb[i], nsb, 'at %.2f yr'%t[i], 'Density profile within the SB', 'radius (pc)', u'n(r) 'r'($cm^{-3}$)', '+')
        figure_number +=1

    if i == 5:
        rmin = 0                    # minimum radius (pc)
        rmax = Rsb[i]-h[i]          # maximum radius (pc)
        r = numpy.linspace(rmin, rmax, number_bin_r)    # position in pc
        Tsb, nsb, r = profile_density_temperature(at, alphat, betat, gammat, deltat, an, alphan, betan, gamman, deltan, n0, L38, t7, r)
        plot(2, 1, r/Rsb[i], nsb, 'at %.2f yr'%t[i], 'Density profile within the SB', 'radius (pc)', u'n(r) 'r'($cm^{-3}$)', '+')
        figure_number +=1
    """
pISM = numpy.ones_like(psb)*pISM
Pob = numpy.ones_like(Lsb)*Pob

# Plots
log_plot(figure_number, 1, t, Rsb, 'none', 'Time evolution of the radius', 'Time (yr)', 'Radius (pc)', '+')
figure_number +=1
log_plot(figure_number, 1, t, Vsb, 'none', 'Time evolution of the velocity', 'Time (yr)', 'Velocity (km/s)', '+')
figure_number +=1
log_plot(figure_number, 2, t, numpy.array([Msb, Mswept]), ['Mass within the SB from the density profile', 'Swept-up mass'], 'Time evolution of the masses', 'Time (yr)', u'Mass'r'($M_\odot$)', ['+', 'x'])
figure_number +=1
log_plot(figure_number, 1, t, Mratio, 'none', 'Time evolution of the ratio between the two masses', 'Time (yr)', r'$M_{SB}/M_{swept-up}$', '+')
figure_number +=1
log_plot(figure_number, 1, t, ns, 'none', 'Time evolution of the density in the shell', 'Time (yr)', r'$n_s$($cm^{-3}$)', '+')
figure_number +=1
log_plot(figure_number, 1, t, h, 'none', 'Time evolution of the thickness of the shell', 'Time (yr)', 'h (pc)', '+')
figure_number +=1
log_plot(figure_number, 1, t, hratio, 'none', 'Time evolution of the ratio between the thickness of the shell and the radius of the SB', 'Time (yr)', r'$h_s/R_{SB}$', '+')
figure_number +=1
log_plot(figure_number, 2, t, numpy.array([psb, pISM]), ['within the SB', 'in the ISM'], 'Time evolution of the pressures', 'Time (yr)', u'Pressure 'r'(dyne $cm^{-2}$)', ['+', 'x'])
figure_number +=1
log_plot(figure_number, 1, t, pratio, 'none', 'Time evolution of the ratio between the two pressures', 'Time (yr)', r'$p_{SB}/p_{ISM}$', '+')
figure_number +=1
log_plot(figure_number, 2, t, numpy.array([Lsb, Pob]), ['of the SB', 'of the stars'], 'Time evolution of the luminosities', 'Time (yr)', u'L 'r'(erg $s^{-1}$)', ['+', 'x'])
figure_number +=1

plt.show()

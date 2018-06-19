"""
It compute all the parameters of the superbubble from the Weaver's model.

All the parameters must be given in the Parameters_system.
"""

##------------------------##
# Librairies and functions #
##------------------------##

import matplotlib.pyplot as plt
import numpy
from Functions import *
from Functions_SB import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_system import *

##====##
# Path #
##====##

    # you need to change it
pathfigure_SB = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Weaver/'

##===========##
# Computation #
##===========##

    # selected time (yr)
t = 4.0e6           # in yr
t6 = t * yr26yr     # in Myr
t7 = t6 * s6yr27yr  # in 10 Myr

    # Correction factor
Rsb = 47.0                 # observed radius (pc)                              #you need to change it for your simulations
Rw, Vw = radius_velocity_SB(t6) # from Weaver's model (pc and km/s)
correction_factor = Rsb/Rw

    # Mass inside the SB and the sweptup mass (solar masses)
Msb, Mswept = masses(t7, Rsb)

    # thickness (pc), the density (cm^-3) of the shell
ns, hs = density_thickness_shell_percentage(percentage, Rsb, Mswept, Msb)

    # distance array of the shell (pc)
rmin = Rsb - hs   # in pc
rmax = Rsb        # in pc
number_bin_rs = 5
rshell = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_rs)

nshell = ns * numpy.ones(number_bin_rs)
Tshell = Ts * numpy.ones(number_bin_rs)

    # Density (cm^-3) and temperature (K) profiles of the Superbubble (cm^-3)
rsbmin = 1       # in pc
rsbmax = rmin       # in pc
number_bin_rsb = 20
rsb = numpy.logspace(numpy.log10(rsbmin), numpy.log10(rsbmax), number_bin_rsb)
Tsb, nsb = profile_density_temperature(t7, rsb, Rsb)

    # Plot
figure_number = 1
xlabel = '$r/R_{sb}$'
color = ['cornflowerblue', 'orange']
text = ''

r = []  # distance array (pc)
n = []  # density array (cm^-3)
T = []  # temperature array (K)

for i in range (number_bin_rsb):

    r.append(rsb[i])
    n.append(nsb[i])
    T.append(Tsb[i])

for j in range (number_bin_rs):

    r.append(rshell[j])
    n.append(nshell[j])
    T.append(Tshell[j])

n = numpy.asarray(n)/ns

        # Density/temperature profiles
title = 'none'
ylabel = ['$n(r)/n_{shell}$', '$T(r)$ [K]']
symbol = ['-.', ':']
y = [n, T]
label = ['Density', 'Temperature']
r = numpy.asarray(r)/Rsb
log_plot_multi(figure_number, r, y, label, title, xlabel, ylabel, symbol)
plt.savefig(pathfigure_SB+'Density_temperature_profiles.pdf')

figure_number += 1

        # Density
ylabel = '$n(r)/n_{shell}$'
symbol = ''
linestyle = '-.'
label = 'none'
colord = 'cornflowerblue'

log_plot(figure_number, 1, r, n, label, title, xlabel, ylabel, symbol, linestyle, colord, text)
plt.savefig(pathfigure_SB+'Density_profiles.pdf')

figure_number += 1

        # Temperature
ylabel = 'Temperature [K]'
symbol = ''

log_plot(figure_number, 1, r, T, label, title, xlabel, ylabel, symbol, linestyle, colord, text)
plt.savefig(pathfigure_SB+'Temperature_profiles.pdf')

figure_number += 1

    # time array (yr)
tmin = 0.01                    # Myr
tmax = 10                   # Myr
tmin = tmin/yr26yr          # yr
tmax = tmax/yr26yr          # yr
number_bin_t = 30
t_fix = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr
t7 = t6 * s6yr27yr                                  # 10 Myr

    # Evolution of the superbubble
        # Computation of the radius and the velocity of the forward shock
Rsb, Vsb = radius_velocity_SB(t6)
Rsb = correction_factor * Rsb

        # Computation of mass inside the SB and the sweptup mass
Msb, Mswept = masses(t7, Rsb)

        # Computation of the thickness and the density of the shell
ns, hs = density_thickness_shell(Vsb, Mswept, Msb, Rsb)

        # Computation of the pressure inside the SB
psb = pressure_SB(t6)

        # Computation of the luminosity of the SB due to the radiatiative loss
Lsb = luminosity_SB(t7, Rsb)

    # Plots
xlabel = 'Time [Myr]'
text = ''
symbol = ['', '']
linestyle = ['dashed', ':']

        # Radius vs velocity
ylabel = ['Radius [pc]', 'Velocity [km s$^{-1}$]']
symbol = ['-.', ':']
label = ['Radius', 'Velocity']
y = [Rsb, Vsb]
log_plot_multi(figure_number, t6, y, label, title, xlabel, ylabel, symbol)
plt.savefig(pathfigure_SB+'Radius_velocity.pdf')

figure_number += 1

        # Velocities comparison
ylabel = 'Velocity [km s$^{-1}$]'
label = ['V$_{sb}$', 'C$_0$']
y = [Vsb, C0*numpy.ones_like(t6)]

color = ['cornflowerblue', 'orange']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Velocities.pdf')

figure_number += 1

        # Masses inside the SB and swept-up
label = ['$M_{sb}$', '$M_{swept-up}$']
ylabel = 'Mass [$M_{\odot}$]'
y = [Msb, Mswept]

log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Masses.pdf')

figure_number += 1

        # Percentage of mass in the SB
label = 'none'
ylabel = '$M_{sb}/M_{su}$'
y = Msb/Mswept

log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, '', '-', 'cornflowerblue', text)
plt.savefig(pathfigure_SB+'Mass_percentage.pdf')

figure_number += 1
        # Thickness vs radius
label = ['$h_s$', '$R_{sb}$']
ylabel = 'Lenght [pc]'
y = [hs, Rsb]

log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Thickness_radius.pdf')

figure_number += 1

        # Percentage of thickness
label = 'none'
ylabel = '$h_s/R_{sb}$'
y = hs/Rsb

log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, '', '-.', 'cornflowerblue', text)
plt.savefig(pathfigure_SB+'Thickness_percentage.pdf')

figure_number += 1

        # Pressure inside the SB
symbol = ['+', 'x']
linestyle = ['', '']

label = ['$p_{sb}$','$p_{ism}$']
ylabel = 'p [dyne cm$^{-2}$]'
y = [psb, pISM*numpy.ones_like(t6)]

color = ['cornflowerblue', 'orange']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Pressures.pdf')

figure_number += 1

        # Percentage of pressure
label = 'none'
ylabel = '$p_{sb}/P_{ISM}$'
y = psb/pISM

log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, '+', '', 'cornflowerblue', text)
plt.savefig(pathfigure_SB+'Pressure_percentage.pdf')

figure_number += 1

        # Luminosity of the SB
ylabel = '$L_{sb}$ [erg s$^{-1}$]'
y = Lsb

log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, '+', '', 'cornflowerblue', text)
plt.savefig(pathfigure_SB+'Luminosity.pdf')

plt.show()

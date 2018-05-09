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

pathfigure_SB = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/Weaver/'

##===========##
# Computation #
##===========##

    # selected time (yr)
t = 4.5e6           # in yr
t6 = t * yr26yr     # in Myr
t7 = t6 * s6yr27yr  # in 10 Myr

    # Correction factor
Robs = 47.0                 # observed radius (pc)                              #you need to change it for your simulations
Rsb = radius_velocity_SB(t6)[0] # from Weaver's model (pc and km/s)
correction_factor = Robs/Rsb

    # Mass inside the SB and the sweptup mass (solar masses)
Msb, Mswept = masses(t7, Rsb)

    # mass inside the supershell (solar masses)
#Ms = Mswept - Msb

    # thickness (pc), the density (cm^-3) of the shell
ns, hs = density_thickness_shell_percentage(percentage, Rsb, Mswept, Msb)

    # distance array of the shell (pc)
rmin = Robs - hs   # in pc
rmax = Robs        # in pc
number_bin_rs = 5
rs = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_rs)

ns = ns * numpy.ones(number_bin_rs)
Ts = Ts * numpy.ones(number_bin_rs)

    # Density (cm^-3) and temperature (K) profiles of the Superbubble (cm^-3)
rsbmin = 1e-2       # in pc
rsbmax = rmin       # in pc
number_bin_rsb = 20
rsb = numpy.logspace(numpy.log10(rsbmin), numpy.log10(rsbmax), number_bin_rsb)
Tsb, nsb = profile_density_temperature(t7, rsb, Robs)

    # Plot
figure_number = 1
xlabel = 'Radius [pc]'
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

    r.append(rs[j])
    n.append(ns[j])
    T.append(Ts[j])

        # Density/temperature profiles
title = ''
ylabel = ['$n(r)$ [cm$^{-3}$]', '$T(r)$ [K]']
symbol = ['-.', ':']
y = [n, T]
label = ['Density', 'Temperature']
log_plot_multi(figure_number, r, y, label, title, xlabel, ylabel, symbol)
plt.savefig(pathfigure_SB+'Density_profiles.pdf')

figure_number += 1
"""
        # density
title = ''
ylabel = 'Density [cm$^{-3}$]'
symbol = ''
linestyle = '-.'
label = 'none'
print(percentage)

log_plot(figure_number, 1, r, n, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Density_profiles.pdf')

figure_number += 1

        # temperature
title = ''
ylabel = 'Temperature [K]'
symbol = ''
linestyle = ':'
label = 'none'

log_plot(figure_number, 1, r, T, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Temperature_profiles.pdf')

figure_number += 1
"""

    # time array (yr)
tmin = 0.01                    # Myr
tmax = 37                   # Myr
tmin = tmin/yr26yr          # yr
tmax = tmax/yr26yr          # yr
number_bin_t = 300
t_fix = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)    # yr
t6 = t_fix * yr26yr                                 # Myr
t7 = t6 * s6yr27yr                                  # 10 Myr

    # Evolution of the superbubble
        # Computation of the radius and the velocity of the forward shock
Rsb, Vsb = radius_velocity_SB(t6)

        # Computation of mass inside the SB and the sweptup mass
Msb, Mswept = masses(t7, Rsb)

        # Computation of the thickness and the density of the shell
ns, hs = density_thickness_shell(Vsb, Mswept, Msb, Rsb)

        # Computation of the pressure inside the SB
psb = pressure_SB(t6)

        # Computation of the luminosity of the SB due to the radiatiative loss
Lsb = luminosity_SB(t7, Rsb)

    # Plots
#figure_number = 1
xlabel = 'Time [Myr]'
text = ''
#text = '$n_{ob}$ = %d\n$n_c$ = %.2e cm$^{-3}$, $T_c$ = %.2e K, $p_c$ = %.2e dynes cm$^{-2}$\n$T_s$ = %.2e K, $n_0$ = %.2e cm$^{-3}$'%(Nob, an, at, ap, Ts, n0)

        # Radius vs velocity
#title = 'Time evolution of the radius and the velocity'
title = ''
ylabel = ['Radius [pc]', 'Velocity [km s$^{-1}$]']
symbol = ['x', '+']
label = ['Radius', 'Velocity']
y = [Rsb, Vsb]
log_plot_multi(figure_number, t6, y, label, title, xlabel, ylabel, symbol)
plt.savefig(pathfigure_SB+'Radius_velocity.pdf')

figure_number += 1

        # Velocities comparison
#title = 'Comparison of the velocities'
ylabel = 'Velocity [km s$^{-1}$]'
symbol = ['x', '+']
label = ['V$_{sb}$', 'C$_0$']
y = [Vsb, C0*numpy.ones_like(t6)]
linestyle = ['', '']
color = ['cornflowerblue', 'orange']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Velocities.pdf')

figure_number += 1


        # Masses inside the SB and swept-up
#title = 'Time evolution of the masses'
label = ['$M_{sb}$', '$M_{swept-up}$']
ylabel = 'Mass [$M_{\odot}$]'
y = [Msb, Mswept]
symbol = ['+', 'x']
linestyle = ['', '']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Masses.pdf')

figure_number += 1

        # Thickness vs radius
#title = 'Time evolution of the thickness of the shell compare to the radius'
label = ['$h_s$', '$R_{sb}$']
ylabel = 'Lenght [pc]'
y = [hs, Rsb]
symbol = ['+', 'x']
linestyle = ['', '']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Thickness_radius.pdf')

figure_number += 1

        # Percentage of thickness
#title = 'Time evolution of the ratio of the thickness'
label = 'none'
ylabel = '$h_s/R_{sb}$'
y = hs/Rsb
symbol = '+'
linestyle = ''
color = 'cornflowerblue'
log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Thickness_percentage.pdf')

figure_number += 1

        # Pressure inside the SB
#title = 'Time evolution of the pressure inside the SB \ncompare to the pressure in the ambient medium'
label = 'none'
ylabel = 'p [dyne cm$^{-2}$]'
y = [psb, pISM*numpy.ones_like(t6)]
symbol = ['+', 'x']
linestyle = ['', '']
color = ['cornflowerblue', 'orange']
log_plot(figure_number, 2, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Pressures.pdf')

figure_number += 1

        # Percentage of pressure
#title = 'Time evolution of the ratio of the pressure'
label = 'none'
ylabel = '$p_{sb}/P_{ISM}$'
y = psb/pISM
symbol = '+'
linestyle = ''
color = 'cornflowerblue'
log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Pressure_percentage.pdf')

figure_number += 1

        # Luminosity of the SB
#title = 'Time evolution of the luminosity of the SB by cooling rate'
label = 'none'
ylabel = '$L_{sb}$ [erg s$^{-1}$]'
y = Lsb
symbol = '+'
linestyle = ''
log_plot(figure_number, 1, t6, y, label, title, xlabel, ylabel, symbol, linestyle, color, text)
plt.savefig(pathfigure_SB+'Luminosity.pdf')

plt.show()

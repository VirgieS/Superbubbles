"""
It computes the CR distribution to check the shape
"""

##------------------------##
# Librairies and functions #
##------------------------##

import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit
from Functions import *
from Functions_SB import radius_velocity_SB
from Functions_CR import *

# Physical constants and conversion factors
from Physical_constants import *
from Conversion_factors import *
from Parameters_system import *

##====##
# Path #
##====##

pathfigure_CR = '/Users/stage/Documents/Virginie/Superbubbles/figures/30_Dor_C/CR/'

##===========##
# Computation #
##===========##

    # SN explosion time
t0min = 3       # Myr
t0max = 37      # Myr
t0 = random_SN(t0min, t0max)/yr26yr   # SN explosion time in yr

    # time array (yr) ONLY AFTER THE SN EXPLOSION
tmin = t0min                # Myr
tmax = t0max + 1            # Myr
tmin = tmin/yr26yr          # yr
tmax = tmax/yr26yr          # yr
number_bin_t = 300
t = numpy.logspace(numpy.log10(tmin), numpy.log10(tmax), number_bin_t)    # yr
t6 = t * yr26yr         # Myr
t7 = t6 * s6yr27yr          # 10 Myr

    # Energy of the cosmic rays (GeV)
number_bin_E = 10
ECR = numpy.logspace(numpy.log10(Emin_CR), numpy.log10(Emax_CR), number_bin_E)

    # Cosmic rays distribution (GeV^-1) for each time step

        # Initialization
number_bin_r = 2000
NCR = numpy.zeros((number_bin_t, number_bin_r, number_bin_E))
figure_number = 1

indt = numpy.where(t > t0)[0]

for i in (indt):

    delta_time = t[i] - t0

        # Initial cosmic rays distribution (GeV^-1)
    NE = power_law_distribution(ECR)

        # Diffusion coefficient (cm^2 s^-1)
    D = diffusion_coefficient(ECR)

    for j in range (number_bin_E):

            # distance from the center (pc)
        rmin = 0.1      # pc
        rmax = 3 * numpy.sqrt(6 * D[j] * delta_time * yr2s)/pc2cm  # pc
        r = numpy.logspace(numpy.log10(rmin), numpy.log10(rmax), number_bin_r)  # pc

        # Density of cosmic rays (GeV^-1 cm^-3)
        NCR[i, :, j] = diffusion_spherical(delta_time, r, NE[j], D[j])

            # Fit
        if j == 0:

            popt, pcov = curve_fit(gauss, r, NCR[i, :, j])
            fit = gauss(r, *popt)
            d_diff_s = numpy.sqrt(6 * D[j] * delta_time * yr2s)/pc2cm
            Dopt = popt[1]
            d_diff_f = numpy.sqrt(6 * Dopt)

            Title = 'Density of cosmic rays at %.2e GeV and %.2e yr after the SN explosion\n' %(ECR[j], delta_time)
            y = [NCR[i, :, j], fit]
            label = ['Simulation', 'Fit']
            xlabel = 'Distance [pc]'
            ylabel = 'N(r, t, E) [GeV$^{-1}$ cm$^{-3}$]'
            symbol = ['+', '']
            linestyle = ['', '-']
            text = '$d_{diff}$ = %.2f pc (simulations)\n$d_{diff}$ = %.2f pc (fit)'%(d_diff_s, d_diff_f)
            k = i - indt[0]
            plot(figure_number, 2, r, y, label, Title, xlabel, ylabel, symbol, linestyle, text)
            plt.savefig(pathfigure_CR+'CR_E%dgev_t%dyr.pdf'%(ECR[j], k))
            figure_number += 1

                # Verification of the fit
            print('dt = %.2e yr' %delta_time)
            print('The diffusion distance of the CR at dt = %.2e yr and E = %.2e GeV:' %(delta_time, ECR[j]))
            print(d_diff_s)
            print('The simulated standard deviation of the CR at dt = %.2e and E = %.2e GeV:' %(delta_time, ECR[j]))
            print(d_diff_f)

plt.show()

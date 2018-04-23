How to use the different codes.

##===========##
# Funcions.py #
##===========##

There are all the basics functions.
- log_plot      :   returns the plot in log-log scale for both axis
- plot          :   returns the plot in linear scale for both axis
- histogramme   :   returns the histogramme of data
- random_PL     :   returns a random number with size=size (default is 1) for a power-law distribution
- random_SN     :   returns a random number with size=size (default is 1) for a uniform distribution from xmin to xmax
- interpolation :   returns the linear interpolation of a function from a specific data set

##==============##
# Funcions_SB.py #
##==============##

There are all the functions to compute the parameters of the superbubble from the Weaver's model.
- radius_velocity_SB                  :   returns the outer radius (pc) and the velocity (km/s) of the forward shock of the superbubble
- masses                              :   returns the mass inside the superbubble (solar masses) and the swept-up mass by the growth of the superbubble (solar masses)
- density_thickness_shell_percentage  :   returns the density (cm^-3) and the thickness (pc) of the shell if you consider the correction of the Weaver's model
- density_thickness_shell             :   returns the density (cm^-3) and the thickness (pc) of the shell
- pressure_SB                         :   returns the pressure inside the superbubble (dyne cm^-2)
- luminosity_SB                       :   returns the luminosity of the superbubble from the cooling rate (erg s^-1)
- profile_density_temperature         :   returns the temperature (K) and density (cm^-3) profiles inside the superbubble

##==============##
# Funcions_CR.py #
##==============##

Tere are all the functions to compute the CR particles distribution.
- power_law_distribution  :     returns the particles distribution for a power-law
- diffusion_coefficient   :     returns the diffusion coefficient for a power-law
- shell_particles         :     returns the number of particles in a shell from r_in to r_out
- inf_particles           :     returns the number of particles outside the superbubble from r to infinity
- diffusion_spherical     :     returns the density of CR at each time and distance step (from the solution of the diffusion equation without energy loss)
- gauss                   :     returns the gaussian fit from the solution of the diffusion equation without energy loss

##=================##
# Funcions_gamma.py #
##=================##

There are all the functions to compute the gamma emission of the superbubble.
- luminosity      :   returns the gamma luminosity in a specific range of energy (erg s^-1)
- spectral_index  :   returns the photon spectral index for a specific range of energy
- data            :   returns the gamma luminosities in specific energy range and all the energy range of the SB, the spectral index, the luminosity unit, the number of pwne, the number of OB-stars


##=====================##
# Parameters_systems.py #
##=====================##

There are all the parmeters of the system.
- SB parameters   :   free parameters to compute all the parameters of the SB following the Weaver's model
- CR parameters   :   free parameters to compute the cosmic rays production of the SB
- Gamma emission  :   free parameters to compute the gamma emission of the SB

##==================##
# Conversion_factors #
##==================##

All the conversion factors for all units.

##==================##
# Physical_constants #
##==================##

All the physical constants.

##============##
# Statistic.py #
##============##

The main program to make statistic about the gamma emission of one specific SB.
All you need to give in this program are (look for '%#you need to change it for your simulations'):
BEFORE THE ITERATIONS
- nit       :   the number of iterations that you want to do
- t_end     :   the estimated age of the SB (yr) (for 30 Dor C it is 4.5 Myr)
- Robs      :   the size of the SB that you observe to compute the correction factor from the Weaver's model (pc)
- zone      :   which zone you want to compute (1: inside the superbubble, 2: in the shell, 3: outside the superbubble)
- t0min     :   the minimum of the SN explosions time (yr) (default is 3 Myr)
IN THE ITERATIONS
- n         :   the number of SN explosions that you want

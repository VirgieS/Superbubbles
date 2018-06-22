How to use the different codes.

##===========##
# Funcions.py #
##===========##

There are all the basics functions.
- log_plot            :   returns the plot in log-log scale for both axis
- plot                :   returns the plot in linear scale for both axis
- semilog_plot        :   returns the plot in y-log scale and x-linear scale
- log_plot_multi      :   returns a plot in log-log scale with multiple y-axes and the same x-axis
- plot_multi          :   returns a plot in linear scale with multiple y-axes and the same x-axis
- histogramme         :   returns the histogramme of data
- random_PL           :   returns a random number with size=size (default is 1) for a power-law distribution
- interpolation1d     :   returns the linear interpolation of a 1d-function from a specific data set
- interpolation2d     :   returns the linear interpolation of a 2d-function from a specific data set
- loglog_interpolation:   returns the log-log interpolation of a 1d-function from a specific data set
- probability         :   returns the different probabilities for the gamma-ray emission of the superbubble

##==============##
# Funcions_SB.py #
##==============##

There are all the functions to compute the parameters of the superbubble from the Weaver's model.
- radius_velocity_SB                  :   returns the radius (pc) of the superbubble and the velocity (km/s) of the forward shock
- masses                              :   returns the mass inside the superbubble (solar masses) and the swept-up mass by the growth of the superbubble (solar masses)
- density_thickness_shell_percentage  :   returns the density (cm^-3) and the thickness (pc) of the shell if you consider the correction of the Weaver's model
- density_thickness_shell             :   returns the density (cm^-3) and the thickness (pc) of the shell from the Weaver's model
- pressure_SB                         :   returns the pressure inside the superbubble (dyne cm^-2)
- luminosity_SB                       :   returns the luminosity of the superbubble from the cooling rate (erg s^-1)
- profile_density_temperature         :   returns the temperature (K) and density (cm^-3) profiles inside the superbubble

##==============##
# Funcions_CR.py #
##==============##

There are all the functions to compute the CR particles distribution.
- power_law_distribution  :   returns the injected particles distribution for a power-law (GeV^-1)
- diffusion_coefficient   :   returns the diffusion coefficient for a power-law (cm^2 s^-1)
- diffusion_time          :   returns the typical diffusion time scale (yr)
- shell_particles         :   returns the number of particles in a shell from r_in to r_out (GeV^-1)
- inf_particles           :   returns the number of particles outside the superbubble from r to infinity (GeV^-1)
- diffusion_spherical     :   returns the density of CR at each time and distance step (from the solution of the diffusion equation for a homogeneous and isotropic diffusion) (GeV^-1 cm^ -3)
- gauss                   :   returns the gaussian fit from the solution of the diffusion equation for a homogeneous and isotropic diffusion

##=================##
# Funcions_gamma.py #
##=================##

There are all the functions to compute the gamma emission of the superbubble.
- cosmicray_lis   :   returns the cosmic rays spectrum from a table of kinetic energy (GeV^-1 cm^-3)
- pwn_emission    :   returns the TeV emission of a pulsar wind nebula (erg s^-1)
- psr_emission    :   returns the GeV emission of a pulsar (erg s^-1)
- luminosity      :   returns the gamma luminosity in a specific range of energy (erg s^-1)
- spectral_index  :   returns the photon spectral index for a specific range of energy
- data            :   returns the gamma-rays luminosity, the differential gamma-ray luminosity in the whole energy range, the TeV and GeV emission of PWN and pulsar, the number of remained OB-stars and the parameters of the superbubble to check the values
- energy_gamma    :   returns the energy radiation by gamma photons (erg)

##=====================##
# Parameters_systems.py #
##=====================##

There are all the parmeters of the system.
- SB parameters   :   free parameters to compute all the parameters of the SB following the Weaver's model and beyond this model
- CR parameters   :   free parameters to compute the cosmic rays production of the SB
- Gamma emission  :   free parameters to compute the gamma emission of the SB
- SN and time     :   time array of the computation and tsnmin and tsnmax

##==================##
# Conversion_factors #
##==================##

All the conversion factors for all units.

##==================##
# Physical_constants #
##==================##

All the physical constants.

## ============ ##
# Cosmic_rays.py #
## ============ ##

This program compares the gaussian fit and the computatio of the CR distribution by diffusion.

## ============= ##
# Weaver_model.py #
## ============= ##

This program computes the parameters of the system from the Weaver model.

All you need to give is the time at which you want to compute the density and temperature profiles and the time array to compute the evolution of the parameters of the superbubble.

## =========== ##
# Iterations.py #
## =========== ##

The main program to compute each sampling of the system for future statistic analyzes.
All you need to give in this program are (look for '#you need to change it for your simulations'):

BEFORE THE SAMPLINGS
- nit             :   the number of sampling that you want to do
- zone            :   which zone you want to compute (1: inside the superbubble, 2: in the shell, 3: outside the superbubble)
- need_correction :   if you correct the outer radius from Weaver's model by the observed radius
- t_end           :   if you correct the outer radius, then you need to give the estimated age of the SB (yr) (for 30 Dor C it is 4.5 Myr)
- Rsb             :   if you correct the outer radius, then you need to give the size of the SB that you observe to compute the correction factor from the Weaver's model (pc)

IN THE SAMPLING
- Emin            :   minimum energy to compute the spectral index and the associated gamma-ray luminosity (GeV)
- Emax            :   maximum energy to compute the spectral index and the associated gamma-ray luminosity (GeV)

It returns a file  with:

- Lum_HESS_it     :   gamma-ray luminosity in the H.E.S.S. energy range (erg s^-1)
- Lum_Fermi_it    :   gamma-ray luminosity in the Fermi energy range (erg s^-1)
- Lum_it          :   gamma-ray luminosity in the whole energy range (erg s^-1)
- Gamma_HESS_it   :   photon spectral index in the H.E.S.S. energy range
- Gamma_GeV_it    :   photon spectral index from 1 GeV to 10 GeV
- Gamma_MeV_it    :   photon spectral index from 100 MeV to 1 GeV
- Lum_pwn_it      :   TeV emission of PWN (erg s^-1)
- Lum_psr_it      :   GeV emission of pulsar (erg s^-1)
- nob_it          :   number of remained massive stars
- Rsb             :   outer radius of the superbubble (pc)
- Vsb             :   velocity of the forward shock (km/s)
- Ms              :   mass in the shell (solar masses)
- ns              :   density in the shell (cm^-3)


## ========= ##
# Plotting.py #
## ========= ##

The main program to compute the statistic analyzes of the samplings.
Make sure that you first run Iterations.py

Load the data from files to built one single file with:
TOTAL
- Lum_HESS_it     :   gamma-ray luminosity in the H.E.S.S. energy range (erg s^-1)
- Lum_Fermi_it    :   gamma-ray luminosity in the Fermi energy range (erg s^-1)
- Lum_it          :   gamma-ray luminosity in the whole energy range (erg s^-1)
- Gamma_HESS_it   :   photon spectral index in the H.E.S.S. energy range
- Gamma_GeV_it    :   photon spectral index from 1 GeV to 10 GeV
- Gamma_MeV_it    :   photon spectral index from 100 MeV to 1 GeV
- Lum_pwn_it      :   TeV emission of PWN (erg s^-1)
- Lum_psr_it      :   GeV emission of pulsar (erg s^-1)

Except for 30 Dor C

The program builds two files for 30 Dor C
GENERAL:
- Lum_HESS_it     :   gamma-ray luminosity in the H.E.S.S. energy range (erg s^-1)
- Lum_Fermi_it    :   gamma-ray luminosity in the Fermi energy range (erg s^-1)
- Lum_it          :   gamma-ray luminosity in the whole energy range (erg s^-1)
- Gamma_HESS_it   :   photon spectral index in the H.E.S.S. energy range
- Gamma_GeV_it    :   photon spectral index from 1 GeV to 10 GeV
- Gamma_MeV_it    :   photon spectral index from 100 MeV to 1 GeV
- Lum_pwn_it      :   TeV emission of PWN (erg s^-1)
- Lum_psr_it      :   GeV emission of pulsar (erg s^-1)
SB:
- tsn_it          :   supernova explosion times (yr)
- nsn_it          :   number of happened SN

The program plot the graphics.

## ============= ##
# Plotting_tot.py #
## ============= ##

When the previous is already run and the file TOTAL is written, this program can make the statistical analyzes of the samplings.

Plot the mean and the standard deviation from the statistical analyzes and compute the different probabilities of the superbubble.

## ================= ##
# Plotting_one_run.py #
## ================= ##

When Plotting.py is already run and the file TOTAL is written, this program can make a comparison between samplings in both energy ranges to show the necessities to make statistical analyzes.

## ==================== ##
# Plotting_comparison.py #
## ==================== ##

When Plotting.py is already run and the file TOTAL is written for all different parameters that we want to compare, this program plots the mean value of each gamma-ray luminosities and spectral indices to show a comparison between the infleunce of parameters.

## ================== ##
# Plotting_30_Dor_C.py #
## ================== ##

When Plotting.py is run for the case 30 Dor C and GENEREAL and SB are written for all samplings, the program computes the statistical analyzes of the system and compare it to the H.E.S.S. obersvations.

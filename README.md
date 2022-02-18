# EBTEL3
Modified version of the Enthalpy-Based Thermal Evolution of Loops model accounting for kinetic energy as described in Abhishek Rajhans et al. 2022 ApJ 924 13

For more information regarding the EBTEL model see:

+ <a href="http://adsabs.harvard.edu/abs/2008ApJ...682.1351K">Klimchuk et al. 2008, ApJ, 682:1351-1362 (Paper 1)</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...752..161C">Cargill et al. 2012A, ApJ, 752:161 (Paper 2)</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...758....5C">Cargill et al. 2012B, ApJ, 758:5 (Paper 3)</a>
+ <a href="https://ui.adsabs.harvard.edu/abs/2016ApJ...829...31B">Barnes et al. 2016, ApJ, 829:31 (Paper 4)</a>
+ <a href="https://ui.adsabs.harvard.edu/abs/2022ApJ...924...13R">Rajhans et al. 2022, ApJ, 924:13 (Paper 5)</a>

Papers 1-3 developed the basic framework of EBTEL which is based on 0D hydrodynamical description of coronal loops. It models single fluid plasma and assumed the flows to be subsonic at all stages and the default timestep was 1 second. See also [ebtel++](https://github.com/rice-solar-physics/ebtelPlusPlus), a C++ implementation of the EBTEL model (Paper 4) that treats the electron and ion populations separately. Note that while the IDL version of EBTEL only solves the single fluid equations, there is little or no difference between the two-fluid and single-fluid models below roughly 5 MK. An adaptive time grid was used in ebtel++ instead of default time grid of 1 second.  

EBTEL 3 (Paper 5) relaxes the assumption of subsonic flows in papers (1-4) and uses an adaptive time grid ensuring proper numerical resolution of the loop's evolution in the impulsive phase. It is coded in IDL. Additionally the Mach numbers and velocities produced by EBTEL3 are in better agreement with field aligned 1D simulations. EBTEL3 models single fluid plasma. Users interested in two fluid picture will find it useful to compare results obtained from ebtel++ with ebtel3 with the latter giving more reliable estimates of flows in loops. 

## Terms of Use

This code is the authorized version of EBTEL3 dated 2022. Updates will be made as and when necessary.This code is distributed as is. Modifications by the user are unsupported. If you believe you have encountered an error in the original (i.e. unaltered) code, create an issue or submit a pull request with the requested changes.

Use of EBTEL3 should be acknowledged by referencing Papers 1,2,3,4, and 5 as listed above.  

## Purpose

Compute the evolution of spatially-averaged (along the field) loop quantities using simplified equations.  The instantaneous differential emission measure of the transition region is also computed. This version incorporates all the modifications from Paper 2 and is written in modular form for clarity. DEM parts unchanged except for TR pressure correction (see Paper 2).

## Inputs

+ ttime  = time array (s) (The user defined time array should be preferably given at intervals of 1 sec. The code will however return results at nonuniform timesteps depending on instantaneous timescales)
+ heat   = heating rate array (erg cm^-3 s^-1)   (direct heating only; the first element, heat(0), determines the initial static equilibrium)
+ length = loop half length (top of chromosphere to apex) (cm)

## Optional Keyword Inputs

+ classical = set to use the Unsaturated classical heat flux
+ dem_old   = set to use old technique of computing DEM(T) in the trans. reg. (weighted average of demev, demcon, and demeq)
+ flux_nt   = energy flux array for nonthermal electrons impinging the chromosphere (input as a positive quantity; erg cm^-2 s^-1)
+ energy_nt = mean energy of the nonthermal electrons in keV (default is 50 keV)
+ rtv       = set to use Rosner, Tucker, Vaiana radiative loss function (Winebarger form)

## Outputs (ALL outputs listed below are returned at entries in non-uniform timeout array and NOT the input uniform time array)
+ timeout = (non-uniform timesteps at which results are returned. It is different from the input time array given.)
+ heatout = (heating function corresponding to timeout.)
+ t (t_a) = temperature array corresponding (avg. over coronal section of loop / apex)
+ n (n_a) = electron number density array (cm^-3)(coronal avg. / apex) 
+ p (p_a) = pressure array (dyn cm^-2) (coronal avg. /apex)
+ v = velocity array (cm s^-1) (r4 * velocity at base of corona) 
+ c11 = C1 (or r3 in this code) 
+ dem_tr = differential emission measure of transition region, dem(time,T), both legs (dem = n^2 * ds/dT  cm^-5 K^-1) (Note:  dem_tr is not reliable when a nonthermal electron flux is used.)
+ dem_cor = differential emission measure of corona, dem(time,T), both legs(Int{dem_cor+dem_tr dT} gives the total emission measure of a loop strand having a cross sectional area of 1 cm^2)
+ logtdem = logT array corresponding to dem_tr and dem_cor
+ f_ratio = ratio of heat flux to equilibrium heat flux (ratio of heat flux to tr. reg. radiative loss rate)
+ rad_ratio = ratio of tr. reg. radiative loss rate from dem_tr and from r3*(coronal rate)
+ cond = conductive loss from corona
+ rad_cor =  coronal radiative loss

## Correspondence with Variables in ApJ Articles

+ Paper 1
 + r1 = c_3
 + r2 = c_2
 + r3 = c_1
+ Papers 1-3
 + f, ff = F_0
 + f_eq, ff_eq = - R_tr
 + dem_eq = DEM_se

## Usage

+ To include the transition region DEM:
  + `IDL> ebtel2, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout,logtdem`
+ To include a nonthermal electron energy flux:
  + `IDL> ebtel2, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout, logtdem`, flux_nt=flux_nt`
+ To compute rad_ratio (Requires 25% more computing time):
  + `IDL> ebtel2, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout, logtdem, f_ratio, rad_ratio`

### Example Run

+ Define the necessary parameters
    + `IDL> time = findgen(10000)` 										Define input time array
    + `IDL> heat = fltarr(10000)` 										Define corresponding heating array
    + `IDL> heat0 = 0.01` 												Amplitude of (nano)flare (erg cm^-3 s^-1)
    + `IDL> for i = 0, 250 do heat(i) = heat0*time(i)/250.`  				Triangular profile rise
    + `IDL> for i = 251, 500 do heat(i) = heat0*(500. - time(i))/250.` 	Decay
    + `IDL> heat_bkg = 1.e-6`     										Low level constant background heating
    + `IDL> heat = heat + heat_bkg`
    + `IDL> length = 7.5e9`           									Loop half length (cm)
    + `IDL> .compile ebtel3`                Compile

+ Hydro simulation
    + `IDL> ebtel3, time, heat, length, t, n, p, v, ta, na, pa, c11, timeout, heatout, dem_tr, dem_cor, logtdem, /classical`

+ Other examples listed below are exactly the same as those provided in [ebtel](https://github.com/rice-solar-physics/EBTEL) and are repeated here for completeness. 

+ Plot temperature evolution
    + `IDL> plot, time, t, xtit='Time (s)', ytit='Temperature (K)'`

+ Total differential emission measure (corona plus footpoint)
    + `IDL> dem_tot = dem_cor + dem_tr`

+ FeXII (195) G(T) function from CHIANTI
    + `IDL> gofnt, 'fe_12', 190, 200, t_array, g_array, density=1.e9`

+ Compute FeXII intensity
    + `IDL> intensity_ebtel, dem_tot, logtdem, g_array, t_array, int, t, n, length, int_avg`

+ Plot intensity evolution
    + `IDL> plot, time, int, xtit='Time (s)', ytit='FeXII (195) Intensity'`

+ Integrate DEM(T) over 60 s interval starting at t = 1000 s
    + `IDL> dem60 = total(dem_tot(1000:1059,*),1)`

+ Plot DEM(T) for 60 s integration
    + `IDL> plot, logtdem, alog10(dem60), xtit='log T (K)',ytit='log DEM (cm!U-5!N K!U-1!N)', tit='1000-1059 s Integration', xran=[5.,7.], /ynoz`

## Intensities

For observations in which temperature response function, G(T), has units of DN s^-1 pix^-1 cm^5 and the loop diameter, d, is larger than the pixel dimension, l_pix:

+ I_cor_perp = d/(2L) * Int{G(T) * dem_cor(T) * dT}
+ I_tr_perp = d/l_pix * Int{G(T) * dem_tr(T) * dT}
+ I_tr_parallel = Int{G(T) * dem_tr(T) * dT}

for lines-of-sight perpendicular and parallel to the loop axis. I_tr_perp assumes that the transition region is thinner than l_pix. See `intensity_ebtel.pro` for more information.


## CHANGELOG
+ January 2022, Modified to relax the assumption of subsonic flows and incorporated adaptive timegrd in IDL version of EBTEL.
+ May 2012. PC version. Modular form.
+ 2013 Jan 15, JAK, Fixed a bug in the compution of the radiation function at 10^4 K; important for computing radiation losses based on dem_tr; ge vs. gt in computing rad;  lt vs. le in computing rad_dem

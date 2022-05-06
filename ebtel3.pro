;This version of EBTEL solves the 0D equations for a single fluid plasma and does not assume it to be subsonic.
   pro ebtel2, ttime, heat, length, t, n, p, v, ta, na, pa, c11,timeout,heatout, dem_tr, dem_cor,  $
       logtdem, f_ratio, rad_ratio, cond, rad_cor,$
       classical=classical, dynamic=dynamic, dem_old=dem_old, $
       flux_nt=flux_nt, energy_nt=energy_nt, rtv=rtv
     ;
     ; NAME:  Enthalpy-Based Thermal Evolution of Loops (EBTEL)
     ;
     ; PURPOSE:
     ;   Compute the evolution of spatially-averaged (along the field) loop quantities
     ;   using simplified equations.  The instantaneous differential emission measure of
     ;   the transition region is also computed. This version incorporates all the modifications 
     ;   from Cargill et al (2012) and is written in modular form for clarity. 
     ;   DEM parts unchanged except for TR pressure correction (see
     ;   Paper 2). Original list of JK modifications removed from preamble.
     ;
     ; INPUTS:
     ;   ttime  = time array (s)
     ;   heat   = heating rate array (erg cm^-3 s^-1)   (direct heating only)
     ;              (the first element, heat(0), determines the initial static equilibrium)
     ;   length = loop half length (top of chromosphere to apex) (cm)
     ;
     ; OPTIONAL KEYWORD INPUTS:
     ;   classical = set to use the UNsaturated classical heat flux
     ;   dynamic   = set to use dynamical r1 and r2 (NOT recommended, especially when T > 10 MK). Now redundant.
     ;   dem_old   = set to use old technique of computing DEM(T) in the trans. reg.
     ;               (weighted average of demev, demcon, and demeq)
     ;   flux_nt   = energy flux array for nonthermal electrons impinging the chromosphere
     ;               (input as a positive quantity; erg cm^-2 s^-1)
     ;   energy_nt = mean energy of the nonthermal electrons in keV (default is 50 keV)
     ;   rtv       = set to use Rosner, Tucker, Vaiana radiative loss function (Winebarger form)
     ;
     ; OUTPUTS:
     ;   t (t_a) = temperature array corresponding to time (avg. over coronal section of loop / apex)
     ;   n (n_a) = electron number density array (cm^-3) (coronal avg. / apex)
     ;   p (p_a) = pressure array (dyn cm^-2) (coronal avg. /apex)
     ;   v = velocity array (cm s^-1) (r4 * velocity at base of corona)
     ;   c11 = C1 (or r3 in this code)
     ;   dem_tr = differential emission measure of transition region, dem(time,T), both legs
     ;             (dem = n^2 * ds/dT  cm^-5 K^-1)
     ;             (Note:  dem_tr is not reliable when a nonthermal electron flux is used.)
     ;   dem_cor = differential emission measure of corona, dem(time,T), both legs
     ;             (Int{dem_cor+dem_tr dT} gives the total emission measure of a
     ;             loop strand having a cross sectional area of 1 cm^2)
     ;   logtdem = logT array corresponding to dem_tr and dem_cor
     ;   f_ratio = ratio of heat flux to equilibrium heat flux
     ;             (ratio of heat flux to tr. reg. radiative loss rate)
     ;   rad_ratio = ratio of tr. reg. radiative loss rate from dem_tr and from r3*(coronal rate)
     ;   cond = conductive loss from corona
     ;   rad_cor =  coronal radiative loss
     ;
     ; CORRESPONDENCE WITH VARIABLES IN ASTROPHYSICAL JOURNAL ARTICLES:
     ;          (Klimchuk et al., 2008, ApJ, 682, 1351; Cargill et al., ApJ 2012)
     ;   r1 = c_3
     ;   r2 = c_2
     ;   r3 = c_1
     ;   f, ff = F_0
     ;   f_eq, ff_eq = - R_tr
     ;   dem_eq = DEM_se
     ;
     ; USAGE:
     ;   To include the transition region DEM:
     ;      IDL> ebtel3, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout, logtdem
     ;   To include a nonthermal electron energy flux:
     ;      IDL> ebtel3, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout, logtdem, flux_nt=flux_nt
     ;   To compute rad_ratio:
     ;      IDL> ebtel3, ttime, heat, t, n, p, v, ta, na, pa, dem_tr, dem_cor, c11, timeout, heatout, logtdem, f_ratio, rad_ratio
     ;      (Takes 25% more computing time.)
     ;
     ; INTENSITIES:
     ;   For observations in which temperature response function, G(T), has units of
     ;   DN s^-1 pix^-1 cm^5 and the loop diameter, d, is larger than the pixel dimension, l_pix:
     ;
     ;      I_cor_perp = d/(2L)* Int{G(T)*dem_cor(T)*dT}
     ;      I_tr_perp = d/l_pix * Int{G(T)*dem_tr(T)*dT}
     ;      I_tr_parallel = Int{G(T)*dem_tr(T)*dT} ,
     ;
     ;   for lines-of-sight perpendicular and parallel to the loop axis.  I_tr_perp assumes that
     ;   the transition region is thinner than l_pix.
     ;
     ; MISCELLANEOUS COMMENTS:
	 ;   A 1 sec time is generally adequate, but in applications where exceptionally strong conductive 
	 ;      cooling is expected (e.g., intense short duration heating events, especially in short loops), 
	 ;      a shorter time step may be necessary. If there is a question, users should compare runs with 
	 ;      different time steps and verify that there are no significant differences in the results.
     ;   Runs much more quickly if the transition region DEM is not computed.
     ;   Speed can be increased by increasing the minimum DEM temperature from 10^4 to, say, 10^5 K
     ;      or by decreasing the maximum DEM temperature from 10^8.5 to, say, 10^7.5
     ;      (search on 450 and 451).
     ;   The equilibrium base heat flux coefficient of 2/7 is appropriate for uniform heating;
     ;      a coefficient of 4/7 is more appropriate for apex heating.
     ;   To have equal amounts of thermal and nonthermal heating:  flux_nt = heat*length.
     ;   It is desirable to have a low-level background heating during the cooling phase so that the
     ;      coronal temperature does not drop below values at which the corona DEM is invalid.
     ;   r1 = c_3 = 0.7 gives more accurate coronal evolution than the original 0.5, especially in
     ;      the late phase of cooling.  However, it produces excess DEM at the very hottest temperatures
     ;      during impulsive events, and the transition region DEM is somewhat elevated.  We have
     ;      therefore introduced r1_tr = 0.5, which provides a more accurate transition region DEM at
     ;      the same time that r1 = 0.7 provides a more accurate radiative cooling.
     ;   v = (c_2/c_3)*(t_tr/t)*v_0 = (r2/r1)*(t_tr/t)*(v/r4) at temperature t_tr in the transition
     ;      region, where t is the average coronal temperature.
     ;
       
     ; HISTORY:
     ; May 2012. PC version. Modular form. 
     ; See original ebtel.pro for many additional comments.
     ; 2013 Jan 15, JAK, Fixed a bug in the compution of the radiation function at 10^4 K; 
     ;      important for computing radiation losses based on dem_tr;
     ;      ge vs. gt in computing rad;  lt vs. le in computing rad_dem
     ; ---------------------------------------------------------------------
     common params, k_b, mp, kappa_0   
     ntot = n_elements(ttime)
     ; Physical constants Can comment out Hydrad lines if needed.
     k_b = 1.38e-16
     mp = 1.67e-24
     n_he_n_p = 0.075;   He/p abundance.
     z_avg = (1 + 2*n_he_n_p)/(1 + n_he_n_p); Include Helium
     z_avg = 1.; For Hydrad comparison.
     kb_fact = 0.5*(1.+1./z_avg)
     k_b = k_b*kb_fact; Modify equation of state for non-e-p plasma
     m_fact = (1 + n_he_n_p*4.)/(2 + 3.*n_he_n_p); Include Helium
     m_fact = (1 + n_he_n_p*4.)/2.; For Hydrad comparison
     mp = mp*m_fact*(1 + z_avg)/z_avg; Average ion mass
     ;print,k_b/mp
     kappa_0 = 1.e-6 
     kappa_0 = 8.12e-7; Hydrad value
     ; Set ksX in loss function (T-breaks).
     radloss, rad, 1.e6, 0, rtv=rtv
     ; Calculate initial Cs (rs in this code).
     calc_c3, r1
     r1=0.6
     ; Ratio of temperature at top of transition region to apex temperature,
     ;      used only for calculation of transition region DEM, dem_tr
     ;      See miscellaneous comments above.
     ; r1_tr = 0.5
     r1_tr = r1
     ; Ratio of average to apex temperature (c_2 in ApJ paper)
     calc_c2, r2
     r2 = 0.9
     ; Ratio of conductive to radiative losses in equilibrium (c_1 in ApJ paper)
     ; Initial value of C1. Fix value now then later iteration on initial state sets it properly
     r3 = 2
     ; Ratio of average to base velocity
     r4 = 1.0
     ; Invalid heat array
     if (heat(0) eq 0.) then begin
       print, '* * * * No initial loop heating * * * *'
       goto, jump99
     endif
     
     identnum1 = round(1000000000.*randomu(seed,1)) 
     identnum2 = round(1000000000.*randomu(seed,1)) 
     identnum3 = round(1000000000.*randomu(seed,1)) 
     ;Only for the purpose of creating temporary files so that multiple ebtel codes can run in same directory
     q = heat ; user defined heating array
     qtemporary = heat ; temporary value of heating taken while converting uniform grid into adaptive grid
     timetemporary = ttime ; temporary value of time taken while converting uniform grid into adaptive grid
     ntot = n_elements(ttime)
     nntot=ntot 
     t = fltarr(nntot) ; creates temperature array
     n = fltarr(nntot) ; creates density array
     ta = fltarr(nntot) ; creates apex temperature array
     na = fltarr(nntot) ; creates apex density array
     p = fltarr(nntot)  ; creates pressure array
     pa = fltarr(nntot) ; creates apex pressure array
     v = fltarr(nntot)  ; creates velocity array
     ve = fltarr(nntot) ; creates array of velocity solved by analytical solution to cubic equation
     c11 = fltarr(nntot); creates array of ratio of total radiation loss from TR and corona
     cond = fltarr(nntot)
     rad_cor = fltarr(nntot)
     dz = fltarr(nntot)
     dp = fltarr(nntot-1)
     dn = fltarr(nntot-1)
     timeout = fltarr(nntot) ; creates array for timesteps but on adaptive grid
     heatout = fltarr(nntot) ; creates array for heating but on adaptive grid
     flux_nt_out  = fltarr(nntot)
     j_nt_out  = fltarr(nntot)
     heatout(0) = q(0) ; fixes initial value of heating
     timeout(0) = ttime(0) ; fixes initial value of time		
     c1 = -(2./7.)*kappa_0 ; Set up coefficient of Spitzer thermal conduction flux
     c_sat = -1.5*k_b^1.5/9.1e-28^0.5 ; Set up coefficient of saturation thermal conduction flux
     sat_limit = 1./6. ;saturation limit used in saturated thermal conduction flux
     if not keyword_set(flux_nt) then flux_nt = timeout*0.     
     ; Set up flux of nonthermal electrons as 0 if keyword flux_nt is not used
     flux_nt = -flux_nt 
     ; converts user defined positive nonthermal flux into negative showing it is a loss from corona
     maxfluxnt = max(abs(flux_nt))
     if not keyword_set(energy_nt) then energy_nt = 50.
     ; if energy of non thermal particles not specified by keyword energy_nt, then default value is 50 keV
     j_nt = 6.241e8*flux_nt/energy_nt 
     ; converts energy flux of non thermal particles into number flux
     dlogt_cor = -alog10(r2)
     dj = fix(100*dlogt_cor)
     nj = 2*dj + 1
     ; Set up DEM in transition region
     if n_params() gt 12 then begin
       logtdem = findgen(451)/100. + 4.
       tdem = 10.^logtdem
       root_tdem = tdem^0.5
       fourth_tdem = tdem^0.25
       dem_tr = fltarr(nntot,451)
       dem_cor = fltarr(nntot,451)
       rad_dem = fltarr(451)
       demev = fltarr(nntot,451)
       demcon = demev
       demeq = demev
       rad_ratio = fltarr(nntot)
       f_ratio = fltarr(nntot)
       root_c2 = (kappa_0/(20.*k_b))^0.5/k_b    ; root_c2 to avoid overflow in dem_ev
       c3 = -5.*k_b
       c4 = (kappa_0/14.)^0.5/k_b
       ;    Radiation in transition region
       for i = 0, 450 do begin
         radloss, rad,tdem(i),1, rtv=rtv
         rad_dem(i)=rad
;         if tdem(i) le 1.e4 then begin
         if tdem(i) lt 1.e4 then begin
           rad_dem(i) = 1.
         endif
       endfor
       root_rad_dem = rad_dem^0.5
     endif
     ; ---------------------------
     ; Initial static equilibrium
     ; 2 methods. (a) Use EBTEL eqm. (b) Use scalings laws. (a) recommended.
     ; ---------------------------
     ; Set up trial values for C1 = 2
     tt_old = r2*(3.5*r3/(1. + r3)*length*length*q(0)/kappa_0)^(2./7.)
     radloss, rad, tt_old, 1, rtv=rtv
     nn = (q(0)/((1. + r3)*rad))^0.5
     nn_old = nn
     ; Iterate on TT and r3
     for i=1,100 do begin
       calc_c1, tt_old, nn, length, rad, r3
       tt_new = r2*(3.5*r3/(1. + r3)*length*length*q(0)/kappa_0)^(2./7.)
       radloss, rad, tt_new, 1, rtv=rtv
       nn = (q(0)/((1. + r3)*rad))^0.5
       err = tt_new - tt_old
       err_n = nn_old - nn
       if abs(err) lt 1e3 then  begin
         i=100
       endif
       print,r3,tt_new,tt_old,err_n,i
       tt_old = tt_new
       nn_old = nn
     endfor
     tt = tt_old
     nn = (q(0)/((1. + r3)*rad))^0.5
     ;   If want to fix out of eqm start, e.g. cooling flare. Section 4.2, Paper 3.
;        tt=1.e7*r2
;        nn = 1e9/r2   
     print, ' '
     print, 'Model parameters'
     print, '  r1 = ', r1
     print, '  r2 = ', r2
     print, '  R3 = ', r3
     print, ' '  
     t(0)=tt ; initial temperature
     n(0)=nn ; initial density
     p(0) = 2.*k_b*n(0)*t(0) ; inital pressure
     v(0) = 0. ; initial velocity = 0
     ta(0) = t(0)/r2 ; initial apex temperature
     calc_lambda, t(0), sc ; calculates initial equilibrium scale height
     na(0) = n(0)*r2*exp(-2*length*(1.-sin(3.14159/5.))/3.14159/sc); calculates initial apex density
     pa(0) = 2*k_b*na(0)*ta(0) ; calculates initial apex pressure
     print, 'Initial static equilibrium'
     print, '  L = ', length
     print, '  Q = ', q(0)
     print, '  T, T_a = ', tt,ta(0)
     print, '  n, n_a = ', nn, na(0)
     ; Scaling law values
     lambda_0 = 1.78e-19    
     lambda_0 = 1.95e-18  ;  lambda = lambda_0*T^bb
     bb = -0.5
     bb = -2./3.
     q_0 = heat(0)       ;  heating rate (erg cm^-3 s^-1)
     t_0 = r2*(3.5/kappa_0*q_0)^(2./7.)*length^(4./7.) ;  temperature (K) (avg., not apex)
     p_0 = r2^(-7./4.)*(8./7.*kappa_0/lambda_0)^0.5*k_b       $
       *t_0^((11.-2.*bb)/4.)/length  ;  total pressure (dyn cm^-2)
     n_0 = 0.5*p_0/(t_0*k_b)  ;  electron number density (cm^-3)
     v_0 = 0.  ;  velocity
     print, ' '
     print, 'Scaling law values'
     print, '  T = ', t_0
     print, '  P = ', p_0
     print, '  n = ', n_0
     ; ----------------------
     ; Time-dependent heating
     ; ----------------------
     print,'MAX TIMETEMPORARY',max(timetemporary)
     h_tot = 0.
     mmm=0
     ; Loop over time steps
     timeout(0) =0.
     compquant = 0.
     while compquant le max(timetemporary) do begin
     for i = 0, nntot-2 do begin ;ntot-3 taken to make interpolation used for conversion to adaptive grid easier
      if  timeout(i) gt max(timetemporary)  then break  
       f_cl = c1*((t(i)/r2)^3.5)/length ; Spitzer conduction flux
       if keyword_set(classical) then begin
         f = f_cl ; if keyword classical is used then only spitzer conduction flux is considered
       endif else begin
         f_sat = sat_limit*c_sat*n(i)*t(i)^1.5 ; else saturation flux is also taken
         f = -f_cl*f_sat/(f_cl*f_cl + f_sat*f_sat)^0.5 ; effective conduction flux
       endelse
       radloss, rad, t(i), 1, rtv=rtv ; calculates power loss function at current temperature
       condtime = 4.e-10*n(i)*length^2./t(i)^2.5 ; calculates conduction cooling timescale
       heattime = 1.e5
       if i eq 0 then begin
         if qtemporary(i+1) ne heatout(i) then heattime = heatout(i)/abs(heatout(i)-qtemporary(i+1))
       endif  
        if i ne 0 and (heatout(i) ne heatout(i-1)) then heattime = heatout(i)/abs(heatout(i)-heatout(i-1)) 
       denominator =10.
       timescalearray = [condtime,heattime] 
       dt = min(timescalearray)/denominator
        ; ensures that timestep at each step is significantly lower than the lowest timescale of system
       calc_c1, t(i), n(i), length, rad, r3 
       ; calculates ratio of total radiation loss from TR and corona at current time
       calc_c2, r2
       calc_c3, r1
       r12 = (r1/r2)
       r12_tr = r1_tr/r2
       c11(i) = r3; Store for output
       ; Equilibrium thermal conduction flux at base (-R_tr in ApJ paper)
       machn = abs(v(i)*1.2216)/(1.5e4*t(i)^0.5) ; mach number
       mcorr1 = exp(machn*machn*0.5*1.*3.)
       ; factor accounting for departure of average quantities from hydrostatic equilibrium
       f_eq = -r3*n(i)*n(i)*rad*length  ;radiation loss from transition region  
       xxx0 =  exp((7.4839e-5)*length/t(i)) ; relates average pressure to base pressure
       pfac2 = v(i)*v(i)*mp/2. ; kinetic energy of a proton
       pfac3 = n(i)*v(i)/length ; flux at base /Length  
       k1pp =  2./3.*(heatout(i) + ((1. + 1./(r3))*f_eq/(length))-xxx0*(pfac3*pfac2*2.*mcorr1)/(r12)) 
       ;  derivative of pressure
       k1p = k1pp*dt ; integrated pressure equuation
       k1nn = (0.5*p(i)*v(i)*xxx0*mcorr1/(r12*k_b*t(i)*length)) 
       ;derivative of number density in absence of non-thermal particles
       k1n =k1nn*dt ; integrated density equuation
       p(i+1) = p(i)+k1p ; pressure at new timestamp
       n(i+1) = n(i)+k1n ; density at new timestamp
       t(i+1) = p(i+1)/(2.*n(i+1)*k_b) ; temperature at new timestamp
       coff2= (2.*(5.0*k_b*t(i+1)*r12/mp)+0.*2.74956e11) ;defined as b in cubic eqn
       coff3_1 = (mp*(n(i+1))*mcorr1*xxx0)
       coff3_2 = (2.*((f)-f_eq)*r12) 
       ;---------------------FLUX NON-THERMAL -----------------------------------------
       if maxfluxnt ne 0. then begin
         fnt = dspline(timetemporary,flux_nt,timeout(i))
         jnt = dspline(timetemporary,j_nt,timeout(i))
         flux_nt_out(i) = fnt
         j_nt_out(i) = jnt
         k1pp = k1pp-(2./3.)*(1. - (0.5*mp*v(i)*v(i)+1.5*k_b*t(i))/energy_nt)*fnt/length
         k1nn = k1nn+jnt/length  
         coff3_2 = (2.*((f+fnt)-f_eq)*r12) 
       endif
  ;---------------------FLUX NON-THERMAL -----------------------------------------
       coff3 = coff3_2/coff3_1  ; defined as c in cubic eqn  
       facp = abs(coff2)
       facq= abs(coff3)
       pow_phi = 3.*alog10((facp))-2.*alog10((facq))-alog10(6.75) ; discriminant in log
       pos_phi = (1.+signum(coff2)*10.^pow_phi)^0.5 ; value of discriminant if positive
       if (coff3 gt 0.) and ((-1.+pos_phi) ge 0.) then sugar=((-1.+pos_phi)^0.333)-((1.+pos_phi)^0.333)
       ; q > 0 and -1 + discrim > 0
       if (coff3 gt 0.) and ((-1.+pos_phi) lt 0.) then sugar=(-(abs(-1.+pos_phi))^0.333)-((1.+pos_phi)^0.333)
       ; q > 0 and -1 + discrim < 0
       if (coff3 lt 0.) and ((1.-pos_phi) ge 0.) then sugar= ((1.+pos_phi)^0.333)+((1.-pos_phi)^0.333)
       ; q < 0 and 1 - discrim > 0
       if (coff3 lt 0.) and ((1.-pos_phi) lt 0.) then sugar= ((1.+pos_phi)^0.333)-((abs(1.-pos_phi))^0.333)
       ; q < 0 and 1 - discrim < 0           
       vbeforeis= ((0.5*facq)^0.333)*sugar ; velocity in absence of  non-thermal particles 
       if i gt 0 then begin
        if abs(vbeforeis-v(i)) ge abs(v(i)*0.1) then vbeforeis  = (v(i)+vbeforeis)/2.
       endif  
       v(i+1) =vbeforeis
       varia = [abs(n(i+1)-n(i))/n(i),abs(p(i+1)-p(i))/p(i),abs(t(i+1)-t(i))/t(i)] 
       ; array of fractional changes in T, n and P between two timesteps
       epsilon = max(varia) ; maximum fractional change
       timeout(i+1) = timeout(i)+dt ; new time stamp
       temptime = timeout(i+1) ; stores it for temporary usage
       heatout(i+1) = interpol(qtemporary,timetemporary,temptime)
       ; interpolates heating on new timestamp, which is different from user defined time stamps
       tempq = heatout(i+1)
       if  epsilon ge 0.1 then begin 
       ; ensures fractional change in T, n and P should not be more than 10 % between two consecutive timesteps
	gre = 0 ; keeps an account of iterations required to do the above stated task
	while epsilon ge 0.1 do begin ; unless and until fractional change is below 10 % keeps iterating
	 gre = gre+1
	 timeout(i+1) = timeout(i)+dt/2. ; divides original time difference into half
	 temptime = timeout(i+1) ; new time stamp
         heatout(i+1) = interpol(qtemporary,timetemporary,temptime) ; new interpolated heating
	 tempq = heatout(i+1) ; store it for temporary use
	 p(i+1) = p(i)+k1pp*dt/2. ;  pressure at new timestamp 
	 tempp = p(i+1)
	 n(i+1) = n(i)+k1nn*dt/2. ;  density at new timestamp 
	 tempn = n(i+1)
	 t(i+1) = p(i+1)/(2.*n(i+1)*k_b) ;  temperature at new timestamp 
	 tempt = t(i+1)
	 coff2= (2.*(5.0*k_b*t(i+1)*r12/mp)+0.*2.74956e11) ;defined as b in cubic eqn
	 coff3_1 = (mp*(n(i+1))*mcorr1*xxx0)
	 coff3 = (coff3_2/coff3_1) ; defined as c in cubic eqn
	 facp = abs(coff2)
	 facq= abs(coff3)
	 pow_phi = 3.*alog10((facp))-2.*alog10((facq))-alog10(6.75) ; discriminant in log
	 pos_phi = (1.+signum(coff2)*10.^pow_phi)^0.5 ; value of discriminant if positive
         if (coff3 gt 0.) and ((-1.+pos_phi) ge 0.) then sugar=((-1.+pos_phi)^0.333)-((1.+pos_phi)^0.333) 
         ; q > 0 and -1 + discrim > 0
         if (coff3 gt 0.) and ((-1.+pos_phi) lt 0.) then sugar=(-(abs(-1.+pos_phi))^0.333)-((1.+pos_phi)^0.333) 
         ; q > 0 and -1 + discrim < 0
         if (coff3 lt 0.) and ((1.-pos_phi) ge 0.) then sugar= ((1.+pos_phi)^0.333)+((1.-pos_phi)^0.333) 
         ; q < 0 and 1 - discrim > 0
         if (coff3 lt 0.) and ((1.-pos_phi) lt 0.) then sugar= ((1.+pos_phi)^0.333)-((abs(1.-pos_phi))^0.333) 
         ; q < 0 and 1 - discrim < 0           
         vbeforeisol= ((0.5*facq)^0.333)*sugar
         if i gt 0  then begin
          if abs(vbeforeisol-v(i)) ge abs(v(i)*0.1) then vbeforeisol  = (v(i)+vbeforeisol)/2.
       endif  
       v(i+1) = vbeforeisol
       varia = [abs(n(i+1)-n(i))/n(i),abs(p(i+1)-p(i))/p(i),abs(t(i+1)-t(i))/t(i)]
       ; array of fractional changes in T, n and P between two timesteps
       epsilon = max(varia) ; maximum fractional change
       dt = dt/2. ; if fractional change of max 10 % not met again divide time interval into 2
       if gre ge 20 then break 
       ; allows maximum 20 iterations , useful in pointing out any arithmetic error			
	endwhile
	p(i+1) = tempp ; values stored
	n(i+1) = tempn
	t(i+1) = tempt
	endif
       h_tot = h_tot+heatout(i)
       ; Calculate scale height
       calc_lambda, t(i+1), sc
       ; calculate apex quantities 
       ta(i+1) = t(i+1)/r2;
       na(i+1) = n(i+1)*r2*exp(-2.*length*(1.-sin(3.14159/5.))/3.14159/sc);
       pa(i+1) = 2*k_b*na(i+1)*ta(i+1)
       ; Differential emission measure
       if n_params() gt 12 then begin
         ;   Transition Region
         if (r12_tr*t(i) gt tdem(450)) then begin
           print, 'Transition region temperature outside DEM range'
           return
         endif
         if (f ne f_eq) then        $
           cf = f*f_eq/(f - f_eq)   $
         else                       $
           cf = 1.e10*f 
         for j = 0, 450 do begin
           if (tdem(j) lt r12_tr*t(i)) then begin
             if not keyword_set(dem_old) then begin
               ;    New method
               aaa = kappa_0*tdem(j)^1.5
               bbb = -5.*k_b*n(i)*v(i)
               p2kt2 = (p(i)/(2.*k_b*tdem(j)))^2
               ; Next line has TR fix
               p2kt2 = (p(i)*exp(2*sin(3.14159/5.)*length/3.14159/sc)/(2.*k_b*tdem(j)))^2
               ccc = -p2kt2*rad_dem(j)
               dtds1 = (-bbb + (bbb*bbb - 4.*aaa*ccc)^0.5)/(2.*aaa)
               dtds2 = (-bbb - (bbb*bbb - 4.*aaa*ccc)^0.5)/(2.*aaa)
               dtds = max(dtds1, dtds2)
               dem_tr(i,j) = 2.*p2kt2/dtds  ; factor 2 for both legs
             endif else begin
               ;    Old method
               ;           approximation to tr. reg. dem when evaporation dominates
               dem_ev = (root_c2/n(i))*(root_c2*p(i)*p(i)  $
                 /root_tdem(j))/v(i)
               dem_con = c3*n(i)*v(i)/rad_dem(j) 
               ;           approximation to tr. reg. dem under equilibrium conditions
               dem_eq = c4*p(i)/(root_rad_dem(j)*fourth_tdem(j))
               
               dem_tr(i,j) = 2.*(f*dem_ev + cf*dem_eq - f_eq*dem_con)   $
                 /(f + cf - f_eq)              ; factor 2 for both legs
               demev(i,j) = 2.*dem_ev
               demcon(i,j) = 2.*dem_con
               demeq(i,j) = 2.*dem_eq
             endelse
             f_ratio(i) = f/f_eq
             if keyword_set(classical) then begin
               f_ratio(i) = f/f_eq
             endif else begin
               f_ratio(i) = f_cl/f_sat
             endelse
           endif 
         endfor
         ss = where(dem_tr(i,*) lt 0., nneg)
         if (nneg gt 0) then begin
           print, ' '
           print, '***** Negative DEM;   i = ', i
           print, ' '
           if (i ne 0) then begin
             for j = 0, 450 do dem_tr(i,j) = dem_tr(i-1,j)
           endif else begin
             dem_tr(0,*) = 0.  ; to avoid problems at start of run with saturated heat flux
           endelse
         endif
         ;   Corona   (EM distributed uniformly over temperture interval [t_min, t_max])
         t_max = max([t(i)/r2, 1.1e4])
         t_min = max([t(i)*(2. - 1/r2), 1.e4])
         j_max = fix((alog10(t_max) - 4.0)*100)
         j_min = fix((alog10(t_min) - 4.0)*100)
         em = 2.*n(i)*n(i)*length            ; factor of 2 for both legs
         ;         dem0 = em/(t_max - t_min)
         delta_t = 10.^4*(10.^((j_max+0.5)/100.)    $
           - 10.^((j_min-0.5)/100.))
         dem0 = em/delta_t
         for j = j_min, j_max do   $
           dem_cor(i,j) = dem0  
         ;   Transition region radiation losses based on DEM
         if n_params() gt 11 then begin
           rad_loss = 0. 
           for j = 0, 450 do begin
             if (tdem(j) lt r12_tr*t(i)) then $
               rad_loss = rad_loss + dem_tr(i,j)*rad_dem(j)*tdem(j)*0.01*2.3  ; 2.3=ln(10) 
           endfor
           rad_ratio(i) = -rad_loss/f_eq
         endif 
       endif
       cond(i)=f
       rad_cor(i)=f_eq/r3
       ;print,timeout(i),heatout(i),t(i),n(i),v(i),c11(i)
     endfor
     dem_tr(nntot-1,*) = dem_tr(nntot-2,*)
     dem_cor(nntot-1,*) = dem_cor(nntot-2,*)
     c11(nntot-1)=c11(nntot-2)
     compquant = max(timeout,xyz)
     timestamppp = fltarr(2001)
     for j =0,2000 do timestamppp(j) = timeout(0) + (max(timeout)-timeout(0))*j/2000.  
     ;print,minmax(timestamppp)
     hheat = interpol(heatout[0:xyz],timeout[0:xyz],timestamppp) 
     tt = interpol(t[0:xyz],timeout[0:xyz],timestamppp)
     nn = interpol(n[0:xyz],timeout[0:xyz],timestamppp)
     vv = interpol(v[0:xyz],timeout[0:xyz],timestamppp)
     tta = interpol(ta[0:xyz],timeout[0:xyz],timestamppp)
     nna = interpol(na[0:xyz],timeout[0:xyz],timestamppp)
     pp = interpol(p[0:xyz],timeout[0:xyz],timestamppp)
     ppa = interpol(pa[0:xyz],timeout[0:xyz],timestamppp)
     cc11 = interpol(c11[0:xyz],timeout[0:xyz],timestamppp)
     ddem_cor = fltarr(2001,451)
     for ii = 0,450 do ddem_cor(*,ii) = abs(interpol(reform(dem_cor(0:xyz,ii)),timeout(0:xyz),timestamppp))
     ddem_tr = fltarr(2001,451)
     for ii = 0,450 do ddem_tr(*,ii) = abs(interpol(reform(dem_tr(0:xyz,ii)),timeout(0:xyz),timestamppp))
     save,filename = 'temp'+STRTRIM(identnum1,2)+STRTRIM(identnum2,2)+STRTRIM(identnum3,2)+STRTRIM(mmm,2)+'.sav',tt,nn,pp,vv,tta,nna,ppa,cc11,timestamppp,hheat,ddem_tr,ddem_cor 
     t = fltarr(nntot) ; creates temperature array
     n = fltarr(nntot) ; creates density array
     ta = fltarr(nntot) ; creates apex temperature array
     na = fltarr(nntot) ; creates apex density array
     p = fltarr(nntot)  ; creates pressure array
     pa = fltarr(nntot) ; creates apex pressure array
     v = fltarr(nntot)  ; creates velocity array
     ve = fltarr(nntot) ; creates array of velocity solved by analytical solution to cubic equation
     c11 = fltarr(nntot); creates array of ratio of total radiation loss from TR and corona
     cond = fltarr(nntot)
     rad_cor = fltarr(nntot)
     dz = fltarr(nntot)
     dp = fltarr(nntot-1)
     dn = fltarr(nntot-1)
     timeout = fltarr(nntot) ; creates array for timesteps but on adaptive grid
     heatout = fltarr(nntot) ; creates array for heating but on adaptive grid
     flux_nt_out  = fltarr(nntot)
     j_nt_out  = fltarr(nntot)
     dem_cor = fltarr(nntot,451)
     dem_tr = fltarr(nntot,451)
     t(0) = tt(i-1)
     n(0) = nn(i-1)
     p(0) = pp(i-1)
     v(0) = vv(i-1)
     ta(0) = tta(i-1)
     na(0) = nna(i-1)
     pa(0) = ppa(i-1)
     c11(0) = cc11(i-1)
     dem_cor(0,*) = ddem_cor(i-1,*)
     dem_tr(0,*) = ddem_tr(i-1,*)     
     timeout(0) = timestamppp(i-1)
     heatout(0) = hheat(i-1)
     mmm = mmm+1
    endwhile
    print,'mm',mmm,identnum
     heatout = fltarr((2001.*mmm)-mmm+1) 
     timeout = fltarr((2001.*mmm)-mmm+1) 
     t = fltarr((2001.*mmm)-mmm+1)
     p = fltarr((2001.*mmm)-mmm+1)     
     n = fltarr((2001.*mmm)-mmm+1)
     ta = fltarr((2001.*mmm)-mmm+1)
     pa = fltarr((2001.*mmm)-mmm+1)     
     na = fltarr((2001.*mmm)-mmm+1)
     v =  fltarr((2001.*mmm)-mmm+1)
     c11 = fltarr((2001.*mmm)-mmm+1)     
     dem_cor = fltarr((2001.*mmm)-mmm+1,451)
     dem_tr = fltarr((2001.*mmm)-mmm+1,451)
     for xxx = 0,mmm-1 do begin
     restore,'temp'+STRTRIM(identnum1,2)+STRTRIM(identnum2,2)+STRTRIM(identnum3,2)+STRTRIM(xxx,2)+'.sav',/v   
     if xxx eq 0 then yyy = 0 else yyy = 1
     t(2000.*xxx+yyy:2000.*(xxx+1.))   =   tt(yyy:2000)
     n(2000.*xxx+yyy:2000.*(xxx+1.))   =   nn(yyy:2000)
     p(2000.*xxx+yyy:2000.*(xxx+1.))   =   pp(yyy:2000)
     v(2000.*xxx+yyy:2000.*(xxx+1.))   =   vv(yyy:2000)
     ta(2000.*xxx+yyy:2000.*(xxx+1.))   =   tta(yyy:2000)
     pa(2000.*xxx+yyy:2000.*(xxx+1.))   =   ppa(yyy:2000)
     na(2000.*xxx+yyy:2000.*(xxx+1.))   =   nna(yyy:2000)
     c11(2000.*xxx+yyy:2000.*(xxx+1.))   =   cc11(yyy:2000)
     dem_cor(2000.*xxx+yyy:2000.*(xxx+1.),*) = ddem_cor(yyy:2000,*)
     dem_tr(2000.*xxx+yyy:2000.*(xxx+1.),*) = ddem_tr(yyy:2000,*)
     heatout(2000.*xxx+yyy:2000.*(xxx+1.)) = hheat(yyy:2000)
     timeout(2000.*xxx+yyy:2000.*(xxx+1.)) = timestamppp(yyy:2000)
     endfor
     afiles = findfile('temp'+STRTRIM(identnum1,2)+STRTRIM(identnum2,2)+STRTRIM(identnum3,2)+'*.sav')
     file_delete,afiles[*]
     jump99:
     return
   end

pro arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
; Save Planet information Here
; Select planet in ait_retrieval.pro INPUT variable 'planet'
; New Planat? Copy and edit to make new entry if needed.

 Rpl=0
 Rjup  = 71492D5         ; cm

 ;================================================================
 ;Hat32 System Parameters
 if (planet eq 'Hat32') then begin
  Rp    = 1.789D0     ; Radius of Planet (minimum of transit spectra)
  Teff  = 1786D0          ; Original Temperature
  g     = 6.6988D2           ; cm/s^2 gravity
  rstar = 1.219D0           ; Rsun
  rorbit= 0.03426        ;AU
  stellarfile='stellar_spectra/HAT-P-32_lte062-4.5-0.0a+0.0.BT-Settl.7.ncdf'
  Tstar = 6269D0          ;K
 endif
 ;
 ;================================================================
 ;WASP 107 System Parameters
 if (planet eq 'Wasp107') then begin
  Rp    = 0.94D0     ; Radius of Planet (minimum of transit spectra)
  Teff  = 770D0          ; Original Temperature
  g     = 3.090295D2           ; cm/s^2 gravity
  rstar = 0.66           ; Rsun
  rorbit= 0.055         ;AU
  stellarfile='stellar_spectra/WASP-107_lte044-4.5-0.0a+0.0.BT-Settl.7.ncdf'
  Tstar = 4430D0          ;K
 endif
 ;
 ;================================================================
 ;WASP 43 System Parameters
 if (planet eq 'Wasp43') then begin
  Rp    = 0.93D0     ; Radius of Planet (minimum of transit spectra)
  Teff  = 1440D0          ; Original Temperature
  g     = 51.168D2           ; cm/s^2 gravity
  rstar = 0.598           ; Rsun
  rorbit= 0.01424         ;AU
  stellarfile='stellar_spectra/Wasp43A.Phoenix.ncdf'
  Tstar = 4400D0          ;K
 endif
 ;
 ;================================================================
 ;WASP 39 System Parameters
 if (planet eq 'Wasp39') then begin
  Rp    = 1.27D0     ; Radius of Planet (minimum of transit spectra)
  Teff  = 1120D0          ; Original Temperature
  g     = 4.3D2           ; cm/s^2 gravity
  rstar = 0.895           ; Rsun
  rorbit= 0.0486         ;AU
  stellarfile='stellar_spectra/Wasp39A.Phoenix.ncdf'
  Tstar = 5400D0          ;K
 endif
 ;
 ;================================================================
 ;WASP 121 System Parameters
 if (planet eq 'Wasp121') then begin
  Rp    = 1.6938182D0     ; Radius of Planet (minimum of transit spectra)
  Teff  = 1700D0          ; Original Temperature
  g     = 8.43D2          ; cm/s^2 gravity
  rstar = 1.458           ; Rsun
  rorbit= 0.02544         ;AU
  stellarfile='stellar_spectra/WASP121A.Phoenix.ncdf'
  Tstar = 6460D0          ;K
 endif
 ;
 ;================================================================
 ;HAT26 System Parameters
 if (planet eq 'Hat26') then begin
  Rp    = 0.565D0         ; Radius of Planet (minimum of transit spectra)
  Teff  = 1100D0          ; Original Temperature
  g     = 4.58D2          ; cm/s^2 gravity
  rstar = 0.788           ; Rsun
  rorbit= 0.04787         ;AU
  stellarfile='stellar_spectra/HAT26A.Phoenix.ncdf'
  Tstar = 5079D0          ;K
 endif
 ;
 ;================================================================
 ;HAT18 System Parameters
 if (planet eq 'Hat18') then begin
  Rp    = 0.947D0         ; Radius of Planet (minimum of transit spectra)
  Teff  = 841D0           ; Original Temperature
  g     = 5.42D2          ; cm/s^2 gravity
  rstar = 0.717           ; Rsun
  rorbit= 0.0559          ;AU
  stellarfile='stellar_spectra/HAT18A.Phoenix.ncdf'
  Tstar = 4870D0          ;K
 endif
 ;================================================================
 ;================================================================

  Rpl   = Rp*Rjup         ; Rjup (cm)      
  logg  = alog10(g)       ; logg
  Rs    = rstar*6.96D10   ;cm
 ;================================================================

 if (Rpl eq 0) then print,'ERROR: Exoplanet Parameters Not Loaded, see ait_load_system_parameters.pro'

end

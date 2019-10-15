FUNCTION atmospheric_density, temperature, pressure_pa, humidity, xc, $
                              force_xw=force_xw, water_vapor=water_vapor, $
                              dry_air=dry_air, verbose=verbose, $
                              atmosphere_values=atmosphere_values

;+
; NAME:
;	  atmospheric_density
;
; PURPOSE:
;	  Calculate the atmospheric density of air at the supplied conditions
; EXPLANATION:
;	  This procedure will take as inputs a set of physical conditions 
;     (temperature, density, humidity, and CO2 concentration) and 
;     calculate the refractivity (the index of refraction minus 1)
;     of air at a specific wavelength.
; CALLING SEQUENCE:
;	  r = refractivity(wavelength, temp, pressure, humidity, xc)
;
; INPUTS:
;	  temperature = temperature of air, in Celsius (C)
;                   (default = 20)
;	  pressure_pa = local air pressure, in pascals (P)
;                   (default = 100000)
;	  humidity    = relative humidity, in percent (%)
;                   (default = 75)
;	  xc          = concentration of CO2, in parts per million (ppm)
;                   (default = 380)
;
; OPTIONAL INPUT KEYWORDS:
;	  force_xw    = use the value of this keyword as the molar fraction
;                   of water vapor instead of calculating x_w from the 
;                   supplied conditions
;	  water_vapor = if supplied, then calculated the partial pressure
;                   of pure water vapor
;	  dry_air     = if supplied, then calculated the partial pressure
;                   of completely dry air
;	  verbose     = indicates level of output information
;                   0 = no output
;                   1 = standard output
;                   2 = extended output
;
; OPTIONAL OUTPUT KEYWORD:
;	  atmosphere_values = structure containing calculated atmospheric
;                         parameters for given conditions
;     R           = the gas constant
;     Z           = compressibility of moist air
;     Ma          = molar mass of dry air (at a given xc)
;     Mw          = molar mass of water vapor
;     svp         = saturation vapor pressure of water vapor in air 
;                   at temperature TT_K
;     f           = enhancement factor of water vapor in air
;     density     = calculated density
;     TT_C        = temperature of air, in Celsius (C)
;     TT_K        = temperature of air, in Kelvin (K) = TT_C + 273.15
;     pressure_Pa = local air pressure, in pascals (P)
;     humidity    = relative humidity, in percent (%)
;     xw          = concentration of CO2, in parts per million (ppm)
;	
; RESULTS:
;	  density = calculated density
;
; EXAMPLE:
;	  density_a    = atmospheric_density(12, 70000, 35, 400, /dry_air)
; COMMON BLOCKS:
;     none	
;
; PROCEDURE:
;	  Ciddor, P.E., 1996, "Refractive index of air: new equations for 
;         the visible and near infrared", Applied Optics LP, vol. 35, 
;         Issue 9, p.1566
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, INAF/Osservatorio Astrofisico di Arcetri, 2005
;-

IF NOT KEYWORD_SET(temperature) THEN temperature = 20
TT_C           = temperature
TT_K           = TT_C + 273.15

IF N_ELEMENTS(humidity) NE 1 THEN    humidity       = 75
humidity_partial = humidity / 100.

IF N_ELEMENTS(pressure_pa) NE 1 THEN pressure_pa = 100000
IF N_ELEMENTS(xc) NE 1 THEN          xc          = 450
IF NOT KEYWORD_SET(verbose) THEN     verbose     = 0

IF verbose GE 2 THEN BEGIN 
    PRINT
    PRINT, FORMAT='(%"Temp: %6.2f C  <>  Pressure - %8.1f (Pa)  <>  ' + $
                  'Humidity - %6.1f (\%)  <>  CO2 - %3.3d (ppm)")',$
           tt_C,  pressure_Pa, humidity, xc
ENDIF

; *****************  Constants  *******************
; Ciddor, 1996, Appendix A
AA =  1.2378847d * 1e-5   ; K^(-2)
BB = -1.9121316d * 1e-2   ; K^(-1)
CC = 33.93711047d         ;
DD = -6.3431645d * 1e3    ; K

alpha = 1.00062d
beta  = 3.14d * 1e-8      ; Pa^(-1)
gamma = 5.6d * 1e-7       ; deg C^(-2)

a0 =  1.58123d * 1e-6     ; K Pa^(-1)
a1 = -2.9331d * 1e-8      ; Pa^(-1)
a2 =  1.1043d * 1e-10     ; K^(-1) Pa^(-1)
b0 =  5.707d * 1e-6       ; K Pa^(-1)
b1 = -2.051d * 1e-8       ; Pa^(-1)
c0 =  1.9898d * 1e-4      ; K Pa^(-1)
c1 = -2.376d * 1e-6       ; Pa^(-1)
d  =  1.83d * 1e-11       ; K^2 Pa^(-2)
e  = -0.765d * 1e-8       ; K^2 Pa^(-2)

; Ciddor, 1996, Section 3
; gas constant
R = 8.314510d            ; J mol^(-1) K^(-1) - gas constant
; molar mass of water vapor
Mw = 0.018015d           ; kg/mol - molar mass of water vapor
; molar mass of dry air containing a concentration of CO2 of Xc ppm
Malpha = 1e-3 * (28.9635 + 12.011 * 1e-6 * (xc - 400))

; *****************             *******************

; saturation vapor pressure of water vapor in air at temperature T
; Ciddor, 1996, Section 3
svp = EXP( AA * TT_K^2 + BB * TT_K + CC + DD/TT_K )

; enhancement factor of water vapor in air
f   = (alpha + beta * pressure_Pa + gamma * TT_C^2)

IF N_ELEMENTS(force_xw) LT 1 THEN BEGIN
    xw  = f * humidity_partial * svp / pressure_Pa
ENDIF ELSE BEGIN
    xw  = force_xw
ENDELSE
 
; Ciddor, 1996, Appendix
IF N_ELEMENTS(Z) LT 1 THEN BEGIN
    Z = 1 - (pressure_Pa/TT_K) *   (a0 + a1 * TT_C + a2 * TT_C^2  + $
                                   (b0 + b1 * TT_C) * xw          + $
                                   (c0 + c1 * TT_C) * xw^2 )      + $
            (pressure_Pa/TT_K)^2 * (d + e * xw^2)
ENDIF

IF KEYWORD_SET(water_vapor) THEN BEGIN
   density = (pressure_Pa * Mw * xw / (Z * R * TT_K))
   ;PRINT,pressure_Pa,Mw,xw,Z,R,TT_K
ENDIF ELSE IF KEYWORD_SET(dry_air) THEN BEGIN
   density = (pressure_Pa * Malpha / (Z * R * TT_K)) * (1 - xw)
ENDIF ELSE BEGIN
   density = (pressure_Pa * Malpha / (Z * R * TT_K)) * (1 - xw * (1 - Mw/Malpha))
ENDELSE

IF verbose GE 2 THEN BEGIN 
    PRINT, FORMAT='(%"svp: %8.2f  <>  f: %9.7f  <>  xw: %12.4f  <>  ' + $
                  'Z: %14.6f  <>  density: %8.4f")',$
           svp, f, xw, Z, density
ENDIF

atmosphere_values = CREATE_STRUCT('R', R, 'Z', Z, 'Ma', Malpha, 'Mw', Mw, $
                                  'svp', svp, 'f', f, 'density', density, $
                                  'TT_C', TT_C, 'TT_K', TT_K, 'pressure_Pa', pressure_Pa, $
                                  'humidity', humidity, 'xw', xw)

RETURN, density

END

; ***************************************

FUNCTION refractivity, wavelength, temperature, pressure_pa, humidity, xc, verbose=verbose
;+
; NAME:
;	  refractivity
;
; PURPOSE:
;	  Calculate the refractivity of air at the supplied conditions
; EXPLANATION:
;	  This procedure will take as inputs a set of physical conditions 
;     (temperature, density, humidity, and CO2 concentration) and 
;     calculate the refractivity (the index of refraction minus 1)
;     of air at a specific wavelength.
; CALLING SEQUENCE:
;	  r = refractivity(wavelength, temp, pressure, humidity, xc)
;
; INPUTS:
;	  wavelength  = wavelength at which to calculate refractivity,
;                   in nanometers (nm)
;	  temperature = temperature of air, in Celsius (C)
;                   (default = 20)
;	  pressure_pa = local air pressure, in pascals (P)
;                   (default = 100000)
;	  humidity    = relative humidity, in percent (%)
;                   (default = 75)
;	  xc          = concentration of CO2, in parts per million (ppm)
;                   (default = 380)
;
; OPTIONAL INPUT KEYWORDS:
;	  verbose     = indicates level of output information
;                   0 = no output
;                   1 = standard output
;                   2 = extended output
;
; OPTIONAL OUTPUT KEYWORD:
;	
; RESULTS:
;	  nprop = calculated refractivity
;
; EXAMPLE:
;	  r = refractivity(630, 12, 70000, 35, 400)
; COMMON BLOCKS:
;     none	
;
; PROCEDURE:
;	  Ciddor, P.E.: 1996, "Refractive index of air: new equations for 
;         the visible and near infrared", Applied Optics LP vol. 35, 
;         Issue 9, p.1566
;
;	  Ciddor, P.E. and Hill, R.J.: 1999, Applied Optics vol. 38, p. 1663.  
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, INAF/Osservatorio Astrofisico di Arcetri, 2005
;-

wavelength_nm = 633.0d     ; nm
wavelength_nm = wavelength     ; nm
wavelength_mic = wavelength_nm / 1000.
; convert wavelengths in air to vacuum wavelengths [ lamda(air) = lamba(vacuum) / n(air) ]
; using mean index of refraction of air = 1.00027
wavelength_vac = wavelength_mic * 1.00027
wavenumber     = 1 / wavelength_vac

IF N_ELEMENTS(temperature) NE 1  THEN temperature = 20
IF N_ELEMENTS(humidity) NE 1 THEN     humidity    = 75
IF N_ELEMENTS(pressure_pa) NE 1 THEN  pressure_pa = 100000
IF N_ELEMENTS(xc) NE 1 THEN           xc          = 380
IF NOT KEYWORD_SET(verbose) THEN      verbose     = 0

; *****************  Constants  *******************
; Ciddor, 1996, Appendix A
; from Peck and Reeder (1962)
k0 = 238.0185d    ; microns^(-2)
k1 = 5792105.0d   ; microns^(-2)
k2 = 57.362d      ; microns^(-2)
k3 = 167917.0d    ; microns^(-2)

; from Owens (1967)
w0 = 295.235d     ; microns^(-2)
w1 =   2.6422d    ; microns^(-2)
w2 =  -0.032380d  ; microns^(-4)
w3 =   0.004028d  ; microns^(-6)
; *****************             *******************

; Refractivity of Air at 15 deg C, 101325 Pa, 0% humidity, and
; a fixed concentration of 450 ppm of CO2
; Ciddor, 1996, Eq. (1)
nas  = ( 1e-8 * ( k1 / (k0 - wavenumber^2) + k3 / (k2 - wavenumber^2) ) ) + 1

; Refractivity of Air at 15 deg C, 101325 Pa, 0% humidity, and 
; a variable concentration of CO2 of Xc ppm
; Ciddor, 1996, Eq. (2)
naxs = ((nas - 1) * ( 1 + 0.534e-6 * (xc - 450)) ) + 1

; Refractivity of water vapor at 20 deg C, 1333 Pa, 0% humidity
; correction factor derived by Ciddor, 1996 by fitting to measurements
cf = 1.022
; Ciddor, 1996, Eq. (3)
nws = ( 1e-8 * cf * ( w0 + (w1 * wavenumber^2) + (w2 * wavenumber^4) + $
                         (w3 * wavenumber^6) ) ) + 1

; density of dry air at standard conditions
density_axs  = atmospheric_density( 15, 101325, 0, xc, force_xw=0, /dry_air)
; density of water vapor at standard conditions
density_ws   = atmospheric_density( 20, 1333, 100, xc, force_xw=1)

; density of dry air at input/actual conditions
density_a    = atmospheric_density( temperature , pressure_pa, humidity, xc, /dry_air)
; density of water vapor at input/actual conditions
density_w    = atmospheric_density( temperature , pressure_pa, humidity, xc, /water)

density      = atmospheric_density (temperature, pressure_pa, humidity, xc)

; Ciddor, 1996, Eq. (5)

IF (density_ws NE 0) THEN nprop_a = (density_a/density_axs) * (naxs - 1) ELSE nprop_a = 0
IF (density_ws NE 0) THEN nprop_w = (density_w/density_ws) * (nws - 1)   ELSE nprop_w = 0
nprop   = nprop_a + nprop_w


;PRINT, density_a, density_axs, density_w, density_ws
IF (verbose GE 1) THEN PRINT, FORMAT='(%"n(axs): %8.1f  <>  n(ws): %8.1f  <>  rho(a/axs): %15.6f  <>  rho(w/ws):: %15.6f  <>  n(prop): %8.1f")',$
       (naxs-1)*1e8, (nws-1)*1e8, (density_a/density_axs), (density_w/density_ws), nprop * 1e8
IF (verbose GE 2) THEN PRINT,FORMAT='(%"n(air): %8.1f <>  n(water): %8.1f")', $
       (density_a/density_axs) * (naxs - 1) * 1e8, (density_w/density_ws) * (nws - 1) * 1e8

RETURN, nprop

END

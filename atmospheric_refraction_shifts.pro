
FUNCTION atmospheric_refraction_shifts, times_in, wavelengths, Obs_Coord=Obs_Coord, Meteo_Conditions=Meteo_Conditions

;+
; NAME:
;	  atmospheric_reftraction_shifts
;
; PURPOSE:
;	  Calculate the atmospheric refraction for the Sun at the input times and wavelengths
; EXPLANATION:
;	  This procedure will calculate the magnitude of the atmospheric
;     refraction at a series of times and wavelengths. 
;
;	  This procedure will also compute the offsets in heliocentric coordinates
;     which vary in time due to changes in both the magnitude of the refraction 
;     and the orientation of the parallactic angle relative to the solar rotation axis.
;
;     The default
;     location is the Sacramento Peak Observatory, and typical 
;     atmospheric conditions (temperature, density, humidity, 
;     and CO2 concentration) are assumed.
; CALLING SEQUENCE:
;	  atm_refract = atmospheric_reftraction_shifts(times_in, wavelengths)
;
; INPUTS:
;	  times_in    = array of times for the calculation, in Julian dates
;	  wavelengths = aaray of wavelengths for which to perform the calculation, in nanometers
;                   (default = [600, 850])
;
; OPTIONAL INPUT KEYWORDS:
;	  Obs_Coord    = coordinates of observing locations, as an array:
;                    [longitude (deg, E+,W-), latitude (deg, N+,S-), altitude (m) ]
;                    default = [ -105.8, 32.8, 2900] = Sac Peak
;
;     Meteo_Conditions = local atmospheric conditions
;                    [temperature (deg C), pressure (Pa), humidity (%), co2 (ppm) ]
;                    default = [ 14, 73000, 35, 380] = typical Sac Peak
;
; OPTIONAL OUTPUT KEYWORD:
;	
; RESULTS:
;	  atm_refract_str = a structure containing information on the calculated refraction and dispersion.
;         refraction - the downward deflection of the Sun in the sky at the given wavelength [arcsecods].
;                       FLTARR(num_times, num_wavelengths)
;         parallactic_ang - the angle between the dispersion direction (perpendicular to the horizon) and
;                           the direction of celestial (or terrestrial) north [degrees].           
;                           FLTARR(num_times)
;         sfts_heliocent_ew - the shift of the image in heliocentric East-West direction, with an eastward shift
;                             given a negative value [arcsecods].
;                             FLTARR(num_times, num_wavelengths)
;         sfts_heliocent_ns - the shift of the image in heliocentric Nort-South direction, with a southward shift
;                             given a negative value [arcsecods].          
;                             FLTARR(num_times)
;         times_jd - the input times, in Julian days         
;                      DBLARR(num_times)
;         solar_elevation - the elevation of the Sun at the input times [degrees]           
;                      FLTARR(num_times)
;         solar_azimuth - the azimuth of the Sun at the input times [degrees] 
;                      FLTARR(num_times)
;         wavelengths - the array of input wavelengths [nanometers]
;                      FLTARR(num_wavelengths) 
; EXAMPLE:
;	  times_24Sep2011 = JULDAY(09,24,2011,0.,0.,0.) + FINDGEN(1440)/1440.
;     atm_sft = atmospheric_reftraction_shifts( times_24Sep2011 )
;     dispersion_ew = atm_sft.SFTS_HELIOCENT_EW(800:*,0)-atm_sft.SFTS_HELIOCENT_EW(800:*,1)
;     dispersion_ns = atm_sft.SFTS_HELIOCENT_NS(800:*,0)-atm_sft.SFTS_HELIOCENT_NS(800:*,1)
;     plot,dispersion_ew,dispersion_ns
; COMMON BLOCKS:
;     none	
;
; PROCEDURE:
;
; NOTES:
;    Calls refractivity.pro, which is described in 
;        Reardon, K.P. 2006, Solar Physics, v239, pg 503, DOI: 10.1007/s11207-006-0283-2
;
;    The total angle between the solar north pole and the parallactic angle is calculated as:
;        Parallactic Angle - Solar P-anlge
;    I think that's right...
;
;    Also calls eq2hor.pro and sun.pro from the IDL Astronomy Library 
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, INAF/Osservatorio Astrofisico di Arcetri, 2005
;	Written, Kevin Reardon, National Solar Observatory, 2014
;-

hatm         = 11
earth_radius = 6375.

IF N_ELEMENTS(wavelengths) LT 1 THEN wavelengths = 600 
num_waves    = N_ELEMENTS(wavelengths)
num_times    = N_ELEMENTS(times_in)

sunpos, times_in, ra_sun, dec_sun

IF N_ELEMENTS(obs_coord) LT 3 THEN BEGIN 
    ;Sac Peak Coordinates
    longitude =     -105.8
    latitude  =     32.8
    altitude  =     2900
ENDIF ELSE BEGIN
    ;Sac Peak Coordinates
    longitude =     obs_coord[0]       ; longitude, east positive
    latitude  =     obs_coord[1]       ; latitude, north positive
    altitude  =     obs_coord[2]       ; altitude, meters
ENDELSE

IF N_ELEMENTS(Meteo_Conditions) LT 4 THEN BEGIN 
    ; Default meteorological conditions 
	temp         = 14.                 ; temperature (deg C)
	pressure     = 73400               ; pressure (Pa)
	humidity     = 35.                 ; humidity (%)
	co2          = 380.                ; CO2 concentration (ppm)
	; This would be STP, but not very useful for high-altitude observatories
	;temp         = 0.
	;pressure     = 100000
	;humidity     = 0.
	;co2          = 380.
ENDIF ELSE BEGIN
    ;Sac Peak Coordinates
	temp         = Meteo_Conditions[0]
	pressure     = Meteo_Conditions[1]
	humidity     = Meteo_Conditions[2]
	co2          = Meteo_Conditions[3]
ENDELSE

refractivity = refractivity(wavelengths,temp,pressure,humidity,co2,/VERBOSE)

; calculate solar position in sky of observer
; refract = 0 forces eq2hor to calculate true celestial position, not refracted
eq2hor, ra_sun, dec_sun, times_in, elevation_sun, azimuth_sun, hourangle_sun, $
    refract=0, lat=latitude, lon=longitude, alt=altitude

parang_sun_sine    =  SIN(hourangle_sun * !DTOR) / SIN((90-elevation_sun)*!DTOR) * SIN((90-latitude)*!DTOR)
parang_sun         = ASIN( parang_sun_sine >(-1.0d) <1.0d ) * !RADEG

h0 = (hatm - altitude/1000.) / (earth_radius + altitude/1000.)
coeff_a = refractivity * (1 - h0)
coeff_b = refractivity * (h0 - refractivity / 2.)

refraction_atm = FLTARR(N_ELEMENTS(times_in), num_waves)
FOR wv = 0,num_waves - 1 DO BEGIN
    dispersion_1wv = coeff_a(wv) * TAN((90-elevation_sun)*!DTOR) - coeff_b(wv) * (TAN((90-elevation_sun)*!DTOR)^3)
    dispersion_1wv = dispersion_1wv * !RADEG * 3600.
    refraction_atm(*,wv) = dispersion_1wv
ENDFOR

; ------------------------------------------------------------------------------------
;     calculate the offsets in heliocentric coordinates due to atmospheric refraction
; ------------------------------------------------------------------------------------

caldat, times_in, mon, day, year, time

p_angle_sun = FLTARR(num_times)

FOR tt=0,num_times-1 DO BEGIN
    sun, year[tt], mon[tt], day[tt], time[tt], PA=p_angle_temp
    p_angle_sun[tt] = p_angle_temp
ENDFOR
PRINT,p_angle_temp
parallactic_to_solar = parang_sun - p_angle_sun

num_waves = N_ELEMENTS(wavelengths)
sfts_heliocent_ew = refraction_atm * 0.0
sfts_heliocent_ns = refraction_atm * 0.0

FOR wv=0,num_waves-1 DO BEGIN
    sfts_heliocent_ew(*,wv) = SIN((180 - parallactic_to_solar) * !DTOR) * refraction_atm(*,wv)
    sfts_heliocent_ns(*,wv) = COS((180 - parallactic_to_solar) * !DTOR) * refraction_atm(*,wv)
ENDFOR


; ------------------------------------------------------------------------------------
;     calculate the image scale change in direction of atmospheric refraction
; ------------------------------------------------------------------------------------

image_fov     = 60.      ; arcseconds
image_fov_deg = image_fov / 3600.

alt_range   = [MIN(elevation_sun)>5,MAX(elevation_sun)<90]    ; range of degrees present on selected day
alt_scale   = 60.                       ; number of divisions per degree on interpolated altitude scale
alt_stepn   = ROUND (((alt_range[1] - alt_range[0]) ) * alt_scale)  ; number of steps for altitude scale
alt_scl_new = dindgen(alt_stepn)/alt_scale + alt_range(0)
alt_scl_ave = (alt_scl_new[0:-2] + alt_scl_new[1:-1]) / 2.
alt_scl_deg = alt_scl_new[1:-1] - alt_scl_new[0:-2]

image_scale_factor = FLTARR(N_ELEMENTS(elevation_sun),num_waves)

FOR wv=0,num_waves-1 DO BEGIN
    sun_is_up                = WHERE(elevation_sun GE 1)
    HELP,sun_is_up
    ; calculate atmospheric refraction on a grid evenly spaced in elevation
    refraction_mag           = INTERPOL(refraction_atm[sun_is_up,wv], elevation_sun[sun_is_up], alt_scl_new,/SPLINE)
    ; calculate the difference in refraction between adjacent elevation steps
    refraction_diff          = MEDIAN(refraction_mag[0:-2] - refraction_mag[1:*],15)
    ; normalize the refraction differences to the step size of the degree scale
    refraction_diff_perdeg   = refraction_diff / alt_scl_deg
    ; multiply the refraction difference (per degree) by the size of the field of view
    refraction_diff_fov      = refraction_diff_perdeg * image_fov_deg
    image_scale_cor          = 1 - refraction_diff_fov / image_fov
    image_scale_cor_time     = INTERPOL(image_scale_cor,alt_scl_ave, elevation_sun)
    image_scale_factor[*,wv] = image_scale_cor_time
ENDFOR

atm_refract_str = CREATE_STRUCT( 'refraction',        refraction_atm, $
                                 'parallactic_ang',   parang_sun, $
                                 'sfts_heliocent_ew', sfts_heliocent_ew, $
                                 'sfts_heliocent_ns', sfts_heliocent_ns, $
                                 'image_scl_factor',  image_scale_factor, $
                                 'solar_elevation',   elevation_sun, $
                                 'solar_azimuth',     azimuth_sun, $
                                 'times_jd',          times_in, $
                                 'wavelengths',       wavelengths )
                           
RETURN, atm_refract_str

END
                           

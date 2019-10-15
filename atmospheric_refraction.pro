FUNCTION atmospheric_refraction, wavelengths, input_times, $
                                 latitude=latitude, longitude=longitude, $
                                 altitude=altitude, air_temp=air_temp, $
                                 air_pressure=air_pressure, humidity=humidity, $
                                 co2_conc=co2_conc, $
                                 verbose=verbose,alt_sun=alt_all,az_sun=az_all

;+
; NAME:
;	  atmospheric_refraction
;
; PURPOSE:
;	  Calculate the refraction of the Sun at different wavelengths
; EXPLANATION:
;	  This program calculates the refractivity of the air for a given set 
;	  of atmospheric conditions and the position of the Sun for an array
;	  of input times. 
;
; CALLING SEQUENCE:
;	  atm_ref_struct = atmospheric_refraction(wavelengths, input_times)
;
; INPUTS:
;	  wavelength   = wavelengths at which to calculate refractivity,
;                    in nanometers (nm)
;                    (default = [400, 500, 600, 700, 800])
;	  input_times  = times for which to calculate the position of the Sun
;                    and the related atmospheric refraction (in Julian days)
;                    (default = current time)
;
; OPTIONAL INPUT KEYWORDS
;	  latitude     = latitude of observation location, in decimal degrees,
;                    positive for north latitudes     (default = 20.71)
;	  longitude    = longitude of observation location, in decimal degrees,
;                    positive for east longitudes     (default = -156.25)
;	  altitude     = altitude of observation location, in meters
;                    (default = 3055)
;	  air_temp     = temperature of air, in Celsius (C)
;                    (default = 20)
;	  air_pressure = local air pressure, in pascals (P)
;                    (default = 100000)
;	  humidity     = relative humidity, in percent (%)
;                    (default = 75)
;	  co2_conc     = concentration of CO2, in parts per million (ppm)
;                    (default = 380)
;	  verbose     = indicates level of output information
;                   0 = no output
;                   1 = standard output
;                   2 = extended output
;
; OPTIONAL OUTPUT KEYWORD:
;	
; RESULTS:
;     Program returns a structure containing refraction magnitude for each
;         input time and wavelength, and the parallactic angle for each 
;         input time. The structure also contains descriptors of the units
;         used for each array.
;
;     ** Structure <574550>, 4 tags, length=52, data length=52, refs=1:
;	  REFRACTION_MAG             FLOAT     Array[1, 5]
;     PARALLACTIC_ANGLE          DOUBLE          -13.435534
;     REFRACTION_MAG_UNIT        STRING    'arcsec'
;     PARALLACTIC_ANGLE_UNIT     STRING    'degrees'
;
; EXAMPLE:
;	  
; COMMON BLOCKS:
;     none	
;
; REQUIRES:
;     Requires the following program from the IDL Standard Library. 
;         poly.pro
;     Requires the following programs from the IDL Astronomy User's Library.
;         (http://idlastro.gsfc.nasa.gov/homepage.html)
;         Note: there is a precess.pro also in the JHU/APL IDL Library, but
;               this version from the Astronomy User's Library should be used.
;         cirrange.pro
;         co_aberration.pro
;         co_nutate.pro
;         ct2lst.pro
;         eq2hor.pro
;         hadec2altaz.pro
;         nutate.pro
;         precess.pro
;         sunpos.pro
;         premat.pro
;     Requires the following program from the JHU/APL IDL Library (also included
;         in the IDL Astronomy User's Library). 
;         (http://fermi.jhuapl.edu/s1r/idl/s1rlib/local_idl.html)
;         isarray.pro
;     Requires the following program that should accompany this file. 
;         refractivity.pro
;
; NOTES:
;
;     Example: 
;         Calculates the atmospheric refraction at five different wavelengths  
;         for 22 UTC (near local noon at Haleakala) for each day during the year.
;
;         times_noon_1year = FINDGEN(365) + JULDAY(01,01,2006,22.,0.,0.)
;         wavelengths      = [400,525,630,850,1600]
;         atm_ref_struct   = atmospheric_refraction(wavelengths, times_noon_1year)
;         help,atm_ref_struct,/st
;            ** Structure <560540>, 4 tags, length=10244, data length=10244, refs=1:
;            REFRACTION_MAG           FLOAT     Array[365, 5]
;            PARALLACTIC_ANGLE        DOUBLE    Array[365]
;            REFRACTION_MAG_UNIT      STRING    'arcsec'
;            PARALLACTIC_ANGLE_UNIT   STRING    'degrees'
;
;         Dispersion is calculated by simply taking the difference between 
;             the refraction at two different wavelengths
;         PRINT,MIN(atm_ref_struct.REFRACTION_MAG(*,0) - atm_ref_struct.REFRACTION_MAG(*,3),MAX=max),max
;            0.125757      1.37185
;
;         Calculates the atmospheric refraction at five different wavelengths  
;         for an entire day, 21 March, 2006
;
;         times_1day       = FINDGEN(1440)/1440. + JULDAY(03,21,2006,0.,0.,0.)
;         wavelengths      = [400,525,630,850,1600]
;         atm_ref_struct   = atmospheric_refraction(wavelengths, times_1day)
;         PLOT,(atm_ref_struct.REFRACTION_MAG(*,0) - atm_ref_struct.REFRACTION_MAG(*,3))>(-20)<20
;
;	  This program takes a selection of wavelengths and the four parameters 
;     describing the atmospheric conditions and calculates the refractivity
;     of the air at each wavelength.
;     The position of the Sun is calculated for each input time. These two
;     values are then combined to provide the magnitude and direction of the
;     atmospheric refraction for the Sun.
;
;     The atmospheric conditions are those at the telescope aperture. The 
;     calculated absolute refraction is probably accurate only for zenith 
;     distances less than 75 degrees. The atmospheric dispersion may be 
;     accurate up to zenith distances of 80 degrees.
;
;     The refraction is calculated only for the Sun, but it would be relatively
;     trivial to change the program to allow calculations for any position
;     on the celestial sphere.
;
;     Haleakala Coordinates: -156.25, 20.71, 3055.
;         longitude = -156.25, latitude = 20.71, altitude = 3055.
;
;     La Palma Coordinates: -17.88, 28.76, 2350.
;         longitude = -17.88, latitude = 28.76, altitude = 2350.
;
;     Sac Peak Coordinates: -105.8, 32.8, 2900.
;         longitude = -105.8, latitude = 32.8, altitude = 2900.
;
;     Default location is Haleakala, Maui, Hawai`i
;
;     see "Reardon, K.P.: 2006, Solar Physics, in press" for more details
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, INAF/Osservatorio Astrofisico di Arcetri, 2006
;-

arcsec_conversion = !RADEG * 3600.
IF NOT KEYWORD_SET(verbose) THEN      verbose     = 1

wavelengths  = FLOAT(wavelengths)
num_waves    = N_ELEMENTS(wavelengths)
IF N_ELEMENTS(num_waves) LT 1    THEN wavelengths = FLOAT([400, 500, 600, 700, 800])
num_waves    = N_ELEMENTS(wavelengths)
;print,wavelengths

IF N_ELEMENTS(input_times) LT 1    THEN input_times = systime(/JULIAN,/UTC)
input_times  = DOUBLE(input_times)
num_times    = N_ELEMENTS(input_times)

IF N_ELEMENTS(latitude) LT 1     THEN latitude  = 20.71
IF N_ELEMENTS(longitude) LT 1    THEN longitude = -156.25
IF N_ELEMENTS(altitude) LT 1     THEN altitude  = 50.
latitude     = FLOAT(latitude(0))
longitude    = FLOAT(longitude(0))
altitude     = FLOAT(altitude(0))

IF N_ELEMENTS(air_temp) LT 1     THEN air_temp = 20.
IF N_ELEMENTS(air_pressure) LT 1 THEN air_pressure = 100000.
IF N_ELEMENTS(humidity) LT 1     THEN humidity = 75.
IF N_ELEMENTS(co2_conc) LT 1     THEN co2_conc = 380.
air_temp     = FLOAT(air_temp(0))
air_pressure = FLOAT(air_pressure(0))
humidity     = FLOAT(humidity(0))
co2_conc     = FLOAT(co2_conc(0))

refractivity = refractivity(wavelengths, air_temp, air_pressure, humidity, $
               co2_conc, VERBOSE=verbose - 1)

sunpos,input_times,ra_all,dec_all 

eq2hor,ra_all,dec_all,input_times,alt_all,az_all,ha_all,$
    refract=0,lat=latitude,lon=longitude,alt=altitude
;print,latitude,longitude,altitude
;hlp,ra_all,dec_all,input_times,alt_all,az_all,ha_all
;plot,alt_all

beta = 0.001254 * (273.15 + air_temp) / 273.15
coeff_a = refractivity * (1 - beta)
coeff_b = refractivity * (beta - refractivity / 2.)         

refraction_calc = FLTARR(num_times, num_waves)
FOR wv = 0,num_waves - 1 DO BEGIN
    refraction_wv = coeff_a(wv) * TAN((90-alt_all)*!DTOR) - coeff_b(wv) * (TAN((90-alt_all)*!DTOR)^3)
    refraction_wv = refraction_wv * arcsec_conversion
    refraction_calc(*,wv) = refraction_wv
ENDFOR 

parallactic_angle_sin =  SIN(ha_all * !DTOR) / SIN((90-alt_all)*!DTOR) * SIN((90-latitude)*!DTOR)
parallactic_angle     = ASIN( parallactic_angle_sin >(-1.0d) <1.0d ) * !RADEG

atmospheric_refraction = CREATE_STRUCT('refraction_mag', refraction_calc, $
                                       'parallactic_angle', parallactic_angle, $
                                       'refraction_mag_unit', 'arcsec', $
                                       'parallactic_angle_unit', 'degrees')
STOP

RETURN,atmospheric_refraction

END

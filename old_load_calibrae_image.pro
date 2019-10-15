FUNCTION load_calibrate_image, filename, extension, channel, wavelength_nb=wavelength_nb, dark_file=dark_file, gain_file=gain_file, $
                               rotate_array=rotate_array, keep_size=keep_size, calibration_location=calibration_location, $
                               target_scale=target_scale, verbose=verbose, header_out=header_out

;+
; NAME:
;	  load_calibrate_image
;
; PURPOSE:
;	  
; EXPLANATION:
;	  
; CALLING SEQUENCE:
;	  calibrated_image = load_calibrate_image(filename, extension, channel)
;
; INPUTS:
;	  filename     = name of FITS file containing raw image of interest
;	  extension    = extension in FITS file for raw image of interest
;         channel      = instrument channel suitable for raw image of interest
; 
; OPTIONAL INPUT KEYWORDS:
;     wavelength_nb - wavelength of loaded array, for narrowband data; if not provided, it is read from the FITS header
;     dark_file     - name of dark file that can be used for calibration; if not provided it will be 
;                         determined from load_calibration_info parameters
;     gain_file     - name of gain file that can be used for calibration; if not provided it will be 
;                         determined from load_calibration_info parameters
;     rotate_array  - three options for how to rotate array:
;                    '*solar*' - rotate so solar north is aligned with vertical axis of output array
;                    '*grid*'  - rotate so dot grid would be aligned with vertical axis of output array
;                     number   - rotate array by input value [degrees]
;                    [default = "solar"]
;     target_scale  - desired output spatial scale, in arcsec/pixel; identical for both axes
;                    [default = 0.096]
;     keep_size     - make output array have same spatial dimensions as original input array, even
;                        after rescaling [default = 1]
;     calibration_location - directory path to calibration files, absolute or relative (e.g. 'calibration_files')
;     verbose       - determines if any information is printed out during execution
;
; OPTIONAL OUTPUT KEYWORD:
;     header_out    - header corresponding to loaded image	
; RESULTS:
;
; EXAMPLE:
;	  
; COMMON BLOCKS:
;     none	
;
; PROCEDURE:
;
; NOTES:
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, National Solar Observatory, 2018-2019
;-


IF N_ELEMENTS(verbose)   LE 0 THEN verbose=0
IF verbose LT 2 THEN restore_verbose = 0 ELSE restore_verbose = 1

IF N_ELEMENTS(keep_size) LE 0 THEN keep_size=1

; set up desired orientation and scale 
IF NOT KEYWORD_SET(rotate_array) THEN rotate_array='solar_north'
IF NOT KEYWORD_SET(target_scale) THEN target_scale=0.096

; determine whether input image is narrowband or whitelight
IF StrMatch(channel,'*nb*') THEN channel_id = 'nb' ELSE channel_id = 'wl'
; check input channel specification against typical file directory scheme
channel_guess = channel
IF StrMatch(filename, '*spectral*') THEN channel_guess='nb'
IF StrMatch(filename, '*whitelight*') THEN channel_guess='wl'
IF channel_guess NE channel THEN PRINT,"Warning: Input channel type doesn't seem to match filename path!"

; read in primary header 
noimage         = readfits(filename, exten=0, main_hdr,/Silent)
; read in image to be loaded and extension header
image_array     = readfits(filename, exten=extension, image_hdr, /Silent)
header_out      = [main_hdr, image_hdr]
; determine array dimensions of input image
image_array_sz  = SIZE(image_array,/str)
image_array_szx = image_array_sz.dimensions[0]
image_array_szy = image_array_sz.dimensions[1]

; get time - in string and Julian Day - of the image
image_dateobs    = sxpar(image_hdr, 'DATE-OBS')
image_dateobs_jd = fits_date_convert(image_dateobs)
date_str         = STRMID(image_dateobs,0,10)

; check input channel specification against information in FITS header
image_channel    = sxpar(header_out, 'CHANNEL')
IF StrMatch(image_channel,'*whitelight*',/Fold_Case) THEN channel_guess='wl' ELSE IF StrMatch(image_channel,'*narrowband*',/Fold_Case) THEN channel_guess='nb'
IF channel_guess NE channel THEN PRINT,"Warning: Input channel type doesn't seem to match FITS header information!"

; if wavelength isn't set, determine it from the FITS header information
IF NOT keyword_set(wavelength_nb) THEN wavelength_nb = ROUND(sxpar(image_hdr, 'WAVELNTH'))

; get calibration information from external function
; date_str is actually ignored, so it doesn't matter
cal_params = load_calibration_info(date_str, 'ibis_' + channel_id)

IF NOT Keyword_Set(calibration_location) THEN BEGIN
    repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
    calibration_location = repository_location + 'calibration_files/'
ENDIF

; determine if calibration files are valid and accessible
IF N_ELEMENTS(dark_file) EQ 1 THEN dark_file_use = dark_file ELSE dark_file_use = cal_params.dark_file
dark_name       = cal_params.dark_name
IF FILE_TEST(calibration_location) THEN dark_file_use = calibration_location + '/' + dark_file_use
dark_file_valid = FILE_TEST(dark_file_use)

; identify proper gain file and array for narrowband image of interest
IF N_ELEMENTS(gain_file) EQ 1 THEN gain_file_use = gain_file ELSE gain_file_use = cal_params.gain_file
gain_name = cal_params.gain_name
IF channel_id EQ 'nb' THEN BEGIN
    gain_waves     = FIX(cal_params.gain_file[1,*])
    gain_idx       = get_closest(gain_waves, wavelength_nb)
    gain_file_use  = gain_file_use[0, gain_idx]
    
    gain_waves     = FIX(cal_params.gain_name[1,*])
    gain_idx       = get_closest(gain_waves, wavelength_nb)
    gain_name      = cal_params.gain_name[0, gain_idx]
ENDIF
IF FILE_TEST(calibration_location) THEN gain_file_use = calibration_location + '/' + gain_file_use
gain_file_valid = FILE_TEST(gain_file_use)

IF dark_file_valid NE 1 THEN BEGIN
    dark_files = File_Search(calibration_location,'*dark*' + channel_id + '*', count=dark_count)
    dark_file_use  = dark_files[0]
ENDIF

IF gain_file_valid NE 1 THEN BEGIN
    gain_files = File_Search(calibration_location,'*gain*' + channel_id + '*', count=dark_count)
    gain_file_use  = gain_files[0]
ENDIF

; restore identified calibration files
RESTORE,Verbose=restore_verbose,dark_file_use
RESTORE,Verbose=restore_verbose,gain_file_use

; set defined calibration array to a common name
res = EXECUTE('dark_cal = ' + dark_name[0])
res = EXECUTE('gain_cal = ' + gain_name[0])

; apply dark and flat correction to data
image_array = (image_array - dark_cal) / gain_cal

; rotate image to roughly match solar Cartesian coordinates (North up, East left)
image_array = ROTATE(image_array,cal_params.transpose)

; Determine rotation angle to be applied to image
; we now allow the rotation angle to specified as a time-dependent polynomial
; the time scale is defined as fractional days after midnight
; the suitable rotation angle will be calculated based on the time of the input image
; if only a single term is given, the rotation angle is taken to be constant
; note that time-dependent rotation angle is only applied if rotation to solar north is desired
rot_fit_params     = cal_params.rot_to_solnorth
rot_fit_params_num = N_ELEMENTS(rot_fit_params)
image_time_scale   = image_dateobs_jd - cal_params.rot_to_sol_reftime
image_rot_angle    = 0.0

FOR nn=0,rot_fit_params_num-1 DO BEGIN
    image_rot_angle += rot_fit_params[nn] * (image_time_scale)^nn
ENDFOR

IF STRMATCH(rotate_array,'*solar*') THEN BEGIN
    image_array = ROT(image_array,image_rot_angle,CUBIC=-0.5)
    IF verbose GE 2 THEN PRINT,'Rotating: ', image_rot_angle
ENDIF ELSE IF STRMATCH(rotate_array,'*grid*') THEN BEGIN
    image_array = ROT(image_array,cal_params.rot_to_grid[0],CUBIC=-0.5)
    IF verbose GE 2 THEN PRINT,'Rotating: ', cal_params.rot_to_grid[0]
ENDIF ELSE IF (size(rotate_array,/str)).TYPE_NAME NE 'STRING' THEN BEGIN
    image_array = ROT(image_array,rotate_array,CUBIC=-0.5)
    IF verbose GE 2 THEN PRINT,'Rotating: ', rotate_array
ENDIF

; rescale image to desired, equal plate scale
; retrieve plate scale of loaded image (based on filter wavelength for narrowband data)
IF channel_id EQ 'wl' THEN BEGIN
    plate_scale_im = cal_params.plate_scale[0:1]
ENDIF ELSE BEGIN
    ; filter_idx_select = where(fix(cal_params.plate_scale[2,*])  eq wavelength)
    filter_idx_select = get_closest(cal_params.plate_scale[2,*],wavelength_nb)
    plate_scale_im = cal_params.plate_scale[0:1, filter_idx_select]
ENDELSE

; determine appropriate (pixel) size of output image to achieve the desired spatial scale
target_scale_ratio = plate_scale_im / target_scale
target_scale_pixel = ROUND([image_array_szx,image_array_szy] * target_scale_ratio)
target_scale_diff  = target_scale_pixel - [image_array_szx,image_array_szy]
image_array        = CONGRID(image_array, target_scale_pixel[0], target_scale_pixel[1], CUBIC=-0.5)

IF channel_id EQ 'wl' THEN BEGIN
    image_shift = cal_params.optical_shift
ENDIF ELSE BEGIN
    ; filter_idx_select = where(fix(cal_params.plate_scale[2,*])  eq wavelength)
    filter_idx_select = get_closest(cal_params.optical_shift[2,*],wavelength_nb)
    image_shift       = cal_params.optical_shift[0:1, filter_idx_select]
ENDELSE

; image size will typically change following rescaling
; place the rescaled image in an array the same size as the input array,
; either through trimming or adding a border.
; this seems to be a complicated approach, but it appears to work
IF KEYWORD_SET(keep_size) THEN BEGIN

    pix_rangex_arr = [0,target_scale_pixel[0] - 1]
    pix_rangex_out = [0,image_array_szx-1]
    IF target_scale_diff[0] GT 0 THEN pix_rangex_arr = [ROUND( target_scale_diff[0]/2.), ROUND( target_scale_diff[0]/2.) + image_array_szx-1]
    IF target_scale_diff[0] LT 0 THEN pix_rangex_out = [ROUND(-target_scale_diff[0]/2.), ROUND(-target_scale_diff[0]/2.) + target_scale_pixel[0] -1]

    pix_rangey_arr = [0,target_scale_pixel[1] - 1]
    pix_rangey_out = [0,image_array_szy-1]
    IF target_scale_diff[1] GT 0 THEN pix_rangey_arr = [ROUND( target_scale_diff[1]/2.), ROUND( target_scale_diff[1]/2.) + image_array_szy-1]
    IF target_scale_diff[1] LT 0 THEN pix_rangey_out = [ROUND(-target_scale_diff[1]/2.), ROUND(-target_scale_diff[1]/2.) + target_scale_pixel[1] -1]

    image_array_keep = FLTARR(image_array_szx, image_array_szy) + MEDIAN(image_array)
    image_array_keep[pix_rangex_out[0]:pix_rangex_out[1],pix_rangey_out[0]:pix_rangey_out[1]] = $
         image_array[pix_rangex_arr[0]:pix_rangex_arr[1],pix_rangey_arr[0]:pix_rangey_arr[1]]
    image_array = image_array_keep
ENDIF

; apply shifts to images so the whitelight and different narrowband images are approximately co-aligned
IF channel_id EQ 'wl' THEN BEGIN
    image_shift = cal_params.optical_shift
ENDIF ELSE BEGIN
    ; filter_idx_select = where(fix(cal_params.plate_scale[2,*])  eq wavelength)
    filter_idx_select = get_closest(cal_params.optical_shift[2,*],wavelength_nb)
    image_shift       = cal_params.optical_shift[0:1, filter_idx_select]
ENDELSE

IF verbose GE 1 THEN PRINT,"Applying a shift of : ",image_shift," pixels"
; use ROT to apply the shift, because it doesn't cause wrapping at the edges of the array, which can be annoying (and wrong)
image_array = ROT(image_array, 0.0, 1.0, (image_array_szx-1)/2. - image_shift[0], (image_array_szy-1)/2. - image_shift[1], /Interp, Missing = MEDIAN(image_array))

RETURN,image_array

END

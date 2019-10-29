FUNCTION load_calibration_info, date_str, instrument_channel, cal_directory = cal_directory

; sample usage:
; nb_21apr_info = load_alma_calibration_info('21Apr2017','ibis_nb') 

; these orientations should be pretty standard for IBIS 
; they shouldn't need to be changed for most observations.
ibis_wl_transpose    =  3
ibis_nb_transpose    =  7
; these values for ROSA may change depending on how the feed optics were setup
;rosa_gband_transpose =  6
;rosa_cak_transpose   =  1

; this is a "standard" value which was adopted but could be changed
target_plate_scale   = 0.096

; look for a "calibration_files" directory as a subdirectory of the directory
; where this program is found

IF NOT KEYWORD_SET(cal_directory) THEN BEGIN
calibration_location = repository_location + 'calibration_files/'
repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
calibration_location = repository_location + 'calibration_files/'
ENDIF ELSE BEGIN
calibration_location = cal_directory
ENDELSE

    ; name of IDL save file containing calibration frames and
    ; name of variable in save file for desired calibration type
    ; the user may save dark and gain files with whatever names desired and update information here

    wl_dark_cal_file 	  = 'avg.whitelight.darks.20Dec2018.sav'
    wl_dark_name          = 'wl_dark_ave'
    wl_gain_cal_file      = 'avg.whitelight.gains.20Dec2018.sav'
    wl_gain_name          = 'wl_gain_ave'

    nb_dark_cal_file      = 'nb.darks.20181220.sav'
    nb_dark_name          = 'nb_darks'

    nb_6173_gain_cal_file = ['Fe6173.gain.blueshift.info.20181220_165244.sav', '6173']
    nb_6173_gain_cal_name = ['nb_gain_1652', '6173', 'nb_gain_info_1652']

    nb_7699_gain_cal_file = ['K7699.gain.blueshift.info.20181220_194457.sav', '7699']
    nb_7699_gain_cal_name = ['nb_gain_1944', '7699', 'nb_gain_info_1944']

    nb_6563_gain_cal_file = ['Halpha.gain.blueshift.info.20181220_194457.sav', '6563']
    nb_6563_gain_cal_name = ['nb_gain_1944', '6563', 'nb_gain_info_1944']

    nb_8542_gain_cal_file = ['Ca8542.gain.blueshift.info.20181220_161321.sav', '8542']
    nb_8542_gain_cal_name = ['nb_gain_1613', '8542', 'nb_gain_info_1613']

    nb_gain_cal_file      = [[nb_6173_gain_cal_file], [nb_6563_gain_cal_file], [nb_7699_gain_cal_file], [nb_8542_gain_cal_file]]
    nb_gain_name          = [[nb_6173_gain_cal_name], [nb_6563_gain_cal_name], [nb_7699_gain_cal_name], [nb_8542_gain_cal_name]]

    nb_bad_pixel_file      = 'narrowband.bad.pixel.index.txt'
    wl_bad_pixel_file      = 'whitelight.bad.pixel.index.txt'

    ;pick a reference time for some time-dependent values
    ; typically pick the time of midnight before the start of the observations
    time_ref                = (fits_date_convert('2018-12-20T14:00:00'))[0]
    num_minutes = 500
    
    ; this is how a fixed value for rotation would be defined
    ;rot_wl_to_sol_north     = [0.32]    ; determined from fitting to HMI
    ; this is a quadratic fit to the observed rotation
    ; where the independent variable is fractional days since the reference time
    rot_wl_to_sol_north     = [-1.28] ; determined from fitting to HMI

    ; reference time is midnight in UT
    rot_wl_to_sol_reftime   = time_ref

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    ;scale_ibis_wl          = [ 0.0975688, 0.0977281]  ; determined from dot grid

    ; IBIS narrowband plate scales - determined from fitting to HMI using stx_findlocation.pro
    ; this generally requires running alignment over many whitelight vs. HMI combinations to beat down the noise
    ; but in this case we found a plate scale that was about 2% larger than from the dot grid
    ; Don't know which is correct in an absolute sense, but the key is to match the whitelight image
    ; to the HMI plate scale, so we'll go with HMI-derived value
    scale_ibis_wl           = [ 0.0976, 0.0976 ]   ; mean value, determined from fitting to HMI, about 0.2% larger
    scale_ibis_wl           = [ 0.09545, 0.09765] ; asymmetric values, combining HMI average with ratio derived from dot grid

    ; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    rot_ibis_wl             = [ -0.26 ]

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    
    scale_ibis_8542   = [0.09535, 0.09756, 8542]
    scale_ibis_6563   = [0.09522, 0.09747, 6563]
    scale_ibis_7699   = [0.09533, 0.09742, 7699]
    scale_ibis_6173   = [0.09518, 0.09739, 6173]
    scale_ibis_nb     = [[scale_ibis_6173],[scale_ibis_6563],[scale_ibis_7699],[scale_ibis_8542]]

; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    ; the rotation angle is assumed to be the same for all filters/wavelengths
    rot_ibis_nb             = [-0.10]

    ; now we can take the difference of the grid angles to determine the relative rotation between the
    ; narrowband and whitelight channels. 
    ; this value will be applied to the science observations as well
    rot_nb_to_wl            = [[time_ref],[rot_ibis_nb - rot_ibis_wl]]
    
    ; values of any integer shifts that should be applied to the narrowband or whitelight data
    ; not used by load_calibrated_image.pro (as of July, 2019)
    shift_ibis_wl_even      = [0, 0]
    shift_ibis_nb_even      = [0, 0]

    ; optical shifts between images through different filters
    ; these are determined by computing the cross correlation between the Air Force target images taken through different filters
    ; one filter is choses as the "reference wavelength," or some average position could be chosen
    ; in this case, the H-alpha filter was chosen, since it's position was most "central"
    ; (i.e. the other filters had shifts on either side of the H-alpha position
    ; these shifts are due to optical effects from the (tilted) prefilter or the automatic repositioning of the camera 
    ; to adjust focus for each filter.

    shift_target_8542 = [0.42895508,-1.45053, 8542 ]
    shift_target_6563 = [1.19836, -0.641968, 6563 ]
    shift_target_7699 = [2.142761, 1.323975, 7699 ]
    shift_target_6173 = [0.812685, 0.454041, 6173 ]
    shift_target_nb   = [ [shift_target_6173], [shift_target_6563], [shift_target_7699], [shift_target_8542] ]
    
    ; optical shifts between narrowband reference wavelength and whitelight channel
    ; this value is also calculated from the cross-correlation between the whitelight target image and the 
    ; target image at the reference wavelength
    shift_wl_to_nb          = [ 0.0, 0.0 ]
    

  ; load_destretch vectors


IF FILE_TEST(calibration_location + 'destr.components.nb2wl.new.v2.txt',/Read) THEN BEGIN
	vects_in    = read_ascii(calibration_location + '/destr.components.nb2wl.new.v2.txt',type='float',record_start=0,data_start=1)
	disp_nb2wl  = reform(vects_in.field1[2:3,*],2,49,49)
	rdisp_nb2wl = reform(vects_in.field1[0:1,*],2,49,49)
    ENDIF ELSE BEGIN
        rdisp_nb2wl = -1
        disp_nb2wl  = -1
    ENDELSE



    ; define which filters were in the wheel for a given observing day, and which filters were used
    filter_ids              = ['', '8542', '7773', '6563', '', '7699', '7090', '5876']
    ;filters_used           = [1,3,6,7]   ; filter wheel order
    filters_used            = [6,3,1,7]   ; wavelength sampling order  @ TODO: i might have to change this.
    num_filters             = N_ELEMENTS(filters_used)
; this is a placeholder for the future definition of the time-dependent offsets between wavelengths due to
    ; atmospheric dispersion or time-dependent changes in the offset between whitelight and narrowband channels
    
    num_timesteps           = 300
    wl_to_nb_drift          = DBLARR(4,num_timesteps,num_filters) + 1
    wl_optical_drifts       = FLTARR(2,num_timesteps) + 1
    

    ; this is a placeholder for the future definition of the time-dependent offsets between wavelengths due to
    ; atmospheric dispersion or time-dependent changes in the offset between whitelight and narrowband channels
    ;wl_to_nb_drift          = DBLARR(4,num_timesteps,num_filters) + 1
    wl_to_nb_drift          = [0,0]
    wl_to_nb_drift          = REFORM(wl_to_nb_drift, 2, 1)

    found_atm_disp_file = FILE_TEST(FILE_SEARCH(calibration_location, 'atmospheric.dispersion.calc.20Dec2018.sav'),/Read)		
    IF STRLOWCASE(instrument_channel) EQ 'ibis_nb' THEN BEGIN
        atm_dispersion_file = (FILE_SEARCH(calibration_location, 'atmospheric.dispersion.calc.20Dec2018.sav'))[0]
        IF FILE_TEST(atm_dispersion_file,/Read) THEN BEGIN
            RESTORE,Verbose=0,atm_dispersion_file
            ; provides ATM_DISP_CALC and SEQUENCE_TIMES

            ; pull out refraction values decomposed into x- and y-shifts (in solar heliocent ric coordinates)
            dispcalc_sfts       = [REFORM(atm_disp_calc.SFTS_HELIOCENT_EW,1,1080,5),REFORM(atm_disp_calc.SFTS_HELIOCENT_NS,1,1080,5)]
            dispcalc_sfts_size  = SIZE(dispcalc_sfts)
            num_timesteps       = dispcalc_sfts_size[2]
            ; ATM_DISP_CALC.times_jd also stores the times for the calculated refraction
            atm_dispersion_times = sequence_times

            ; take difference between the refraction at the wavelengths of the whitelight im age
            ; and the refraction at each narrowband wavelength to get theoretical dispersion offsets 
            ; between each of the channels/filters
            ; the atmospheric dispersion values are stored in increasing wavelength order - [5434,7090,7200,7699]
            ; 7200 (index 2) is the wavelength of the whitelight images, so all the shifts 
            ; are differenced with respect to that wavelength
            atm_dispersion_nb            = FLTARR(3,num_timesteps,num_filters)
            atm_dispersion_nb[2,*,0]     = atm_disp_calc.wavelengths[1]*10
            atm_dispersion_nb[2,*,1]     = atm_disp_calc.wavelengths[3]*10
            atm_dispersion_nb[2,*,2]     = atm_disp_calc.wavelengths[0]*10
            atm_dispersion_nb[0:1,*,0,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,1]
            atm_dispersion_nb[0:1,*,1,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,3]
            atm_dispersion_nb[0:1,*,2,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,0]
        ENDIF ELSE BEGIN
            num_timesteps              = 300
            atm_dispersion_nb          = FLTARR(3,2,num_filters)
            FOR filtn = 0,num_filters-1 DO BEGIN & 
                atm_dispersion_nb[2,*,filtn] = FLOAT(filter_ids[filters_used[filtn]]) & 
                ENDFOR
            atm_dispersion_times       = [time_ref, time_ref+1]
        ENDELSE
        wl_optical_drifts       = FLTARR(2,num_timesteps)
    ENDIF 


    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout        = FLTARR(1,5)
    hmi_pixel_cutout[0,*]   = [31, 1030, 232, 1231, time_ref] 
    ;hmi_pixel_cutout[0,*]    = [87, 246, 77, 236, time_ref]
    
    ; ---------- ROSA Corrections --------
    ; values are saved here for future use but are not currently returned
    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband  = [0.0, 0.0]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.0, 0.0]
    rot_rosa_cak      = [+0.41]
    scale_rosa_gband  = [0.07571, 0.07748]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.1543,  0.1542]
    rot_rosa_cak      = [+0.41]


; package all the above information into a structure to be returned to the user

rot_to_solnorth_combo = rot_wl_to_sol_north

CASE STRLOWCASE(instrument_channel) OF
    'ibis_wl' : BEGIN
        output_info = CREATE_STRUCT('transpose',          ibis_wl_transpose, $
                                    'rot_to_solnorth',    rot_to_solnorth_combo, $
                                    'rot_to_sol_reftime', rot_wl_to_sol_reftime, $
                                    'rot_to_grid',        rot_ibis_wl, $
                                    'plate_scale',        scale_ibis_wl, $
                                    'uniform_plate_scale', target_plate_scale, $
                                    'shift_even_scale',   shift_ibis_wl_even, $
                                    'optical_shift',      shift_wl_to_nb, $
                                    'hmi_pixel_cutout',   hmi_pixel_cutout, $
                                    'dark_file',          wl_dark_cal_file, $
                                    'dark_name',          wl_dark_name, $
                                    'gain_file',          wl_gain_cal_file, $
                                    'gain_name',          wl_gain_name, $
                                    'bad_pixel_file',     wl_bad_pixel_file)
    END
    'ibis_nb' : BEGIN
        ; add the relative rotation between the narrowband and whitelight channels to 
        ; the rotation-to-solar-north value
        rot_to_solnorth_combo[0]    += rot_nb_to_wl[1]
        output_info = CREATE_STRUCT('transpose',          ibis_nb_transpose, $
                                    'rot_to_solnorth',    rot_to_solnorth_combo, $
                                    'rot_to_sol_reftime', rot_wl_to_sol_reftime, $
                                    'rot_to_grid',        rot_ibis_nb, $
                                    'plate_scale',        scale_ibis_nb, $
                                    'uniform_plate_scale', target_plate_scale, $
                                    'shift_even_scale',   shift_ibis_nb_even,$
                                    'optical_shift',      shift_target_nb, $
                                    'hmi_pixel_cutout',   hmi_pixel_cutout, $
                                    'nb_to_wl_destr_ref', rdisp_nb2wl, $
                                    'nb_to_wl_destr_sft', disp_nb2wl, $
                                    'wl_to_nb_drift',     wl_to_nb_drift, $
                                    'wl_drifts',          wl_optical_drifts, $
                                    'atm_dispersion',     atm_dispersion_nb, $
                                    'atm_dispersion_times',atm_dispersion_times, $
                                    'filter_ids',         filter_ids, $
                                    'dark_file',          nb_dark_cal_file, $
                                    'dark_name',          nb_dark_name, $
                                    'gain_file',          nb_gain_cal_file, $
                                    'gain_name',          nb_gain_name, $
				    'bad_pixel_file',     nb_bad_pixel_file)
    END
ENDCASE

RETURN,output_info

END

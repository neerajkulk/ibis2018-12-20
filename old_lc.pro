FUNCTION load_calibration_info, date_str, instrument_channel

; sample usage:
; nb_20dec_info = load_calibration_info('20Dec2018','ibis_nb') 

ibis_wl_transpose    =  3
ibis_nb_transpose    =  7
;rosa_gband_transpose =  6
;rosa_cak_transpose   =  1

target_plate_scale   = 0.096

repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
calibration_location = repository_location + 'calibration_files/'

   
;'20dec2018' 

    
    wl_dark_cal_file = 'whitelight.darks.20Dec2018.sav'
    wl_dark_name     = 'wl_dark_ave'
    wl_gain_cal_file = 'whitelight.gains.20Dec2018.sav'
    wl_gain_name     = 'wl_gain'

; @   nb_dark_cal_file = 'nb.dark.30Sep2014.ave.sav'
; @   nb_dark_name     = 'nb_darks'
    nb_6173_gain_cal_file = ['Fe6173.gain.blueshift.info.20181220_165244.sav', '6173']
    nb_6173_gain_cal_name = ['nb_gain_1652', '6173']
    nb_7699_gain_cal_file = ['K7699.gain.blueshift.info.20181220_194457.sav', '7699']
    nb_7699_gain_cal_name = ['nb_gain_1944', '7699']
    nb_6563_gain_cal_file = ['Halpha.gain.blueshift.info.20181220_194457.sav', '6563']
    nb_6563_gain_cal_name = ['nb_gain_1944', '6563']
    nb_8542_gain_cal_file = ['Ca8542.gain.blueshift.info.20181220_161321.sav', '8542']
    nb_8542_gain_cal_name = ['nb_gain_1613', '8542']

    nb_gain_cal_file  = [[nb_6173_gain_cal_file], [nb_7699_gain_cal_file], [nb_6563_gain_cal_file], [nb_8542_gain_cal_file]]
    nb_gain_name      = [[nb_6173_gain_cal_name], [nb_7699_gain_cal_name], [nb_6563_gain_cal_name], [nb_8542_gain_cal_name]]


    num_minutes   = 500
    time_start    = 14.5
    time_axis_min = (DINDGEN(num_minutes)/60.) + time_start
    time_ref      = (fits_date_convert('2018-12-20T14:00:00'))[0]

    rot_grid_to_sol_north = [[time_ref],[-0.249]]
    rot_wl_to_sol_north   = [[time_ref],[0.321]]    ; determined from fitting to HMI

    ;scale_ibis_wl     = [ 0.09668, 0.09680 ]
    ;scale_ibis_wl     = [ 0.09761, 0.09548 ]  ; determined from dot grid
    scale_ibis_wl     = [ 0.09099, 0.09312 ]   ; determined from fitting to HMI, about 0.2% larger
    rot_ibis_wl       = [ -0.26 ]

    scale_ibis_8542   = [0.09535, 0.09756, 8542]
    scale_ibis_6563   = [0.09522, 0.09747, 6563]
    scale_ibis_7699   = [0.09533, 0.09674, 7699]
    scale_ibis_6173   = [0.09518, 0.09739, 6173]
    scale_ibis_nb     = [[scale_ibis_6173],[scale_ibis_6563],[scale_ibis_7699],[scale_ibis_8542]]
    rot_ibis_nb       = [ 0.069 ]
    
    shift_ibis_wl_even = [0, 0]
    shift_ibis_nb_even = [0, 0]  

    shift_target_8542 = [ -0.0,  -0.0, 8542 ] 
    shift_target_6563 = [ -0.0,   0.0, 6563 ] 
    shift_target_7699 = [ -0.0,  -0.0, 7699 ] 
    shift_target_6173 = [ -0.0,  -0.0, 6173 ] 
    shift_target_nb   = [ [shift_target_6173], [shift_target_6563], [shift_target_7699], [shift_target_8542] ]
    
    filter_ids        = [6173, 6563, 7699, 8542]
     ;filters_used     = [1,3,5,7]   ; filter wheel order
    filters_used      = [1,3,5,7]   ; wavelength sampling order
    num_filters       = N_ELEMENTS(filters_used)

    wl_to_nb_drift    =  DBLARR(4,num_minutes,num_filters)
    wl_optical_drifts = FLTARR(2,num_minutes)
    atm_dispersion_nb = FLTARR(2,num_minutes,num_filters)

    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout_target_1 = [346,1345,359,1358] 
    hmi_pixel_cutout_target_2 = [375,1374,324,1323] 
    hmi_pixel_cutout          = FLTARR(2,4)
    hmi_pixel_cutout[0,*]     = hmi_pixel_cutout_target_1
    hmi_pixel_cutout[1,*]     = hmi_pixel_cutout_target_2

    ; ---------- ROSA Corrections --------
    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband  = [0.0, 0.0]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.0, 0.0]
    rot_rosa_cak      = [+0.41]
    scale_rosa_gband  = [0.07571, 0.07748]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.1543,  0.1542]
    rot_rosa_cak      = [+0.41]


rot_to_solnorth_combo = rot_grid_to_sol_north

CASE STRLOWCASE(instrument_channel) OF
    'ibis_wl' : BEGIN
        rot_to_solnorth_combo[*,1] += rot_ibis_wl[0]
        output_info = CREATE_STRUCT('transpose',        ibis_wl_transpose, $
                                    'rot_to_solnorth',  rot_to_solnorth_combo, $
                                    'rot_to_grid',      rot_ibis_wl, $
                                    'plate_scale',      scale_ibis_wl, $
                                    'shift_even_scale', shift_ibis_wl_even, $
                                    'hmi_pixel_cutout', hmi_pixel_cutout)
    END
    'ibis_nb' : BEGIN
        rot_to_solnorth_combo[*,1] += rot_ibis_nb[0]
        output_info = CREATE_STRUCT('transpose',       ibis_nb_transpose, $
                                    'rot_to_solnorth', rot_to_solnorth_combo, $
                                    'rot_to_grid',     rot_ibis_nb, $
                                    'plate_scale',     scale_ibis_nb, $
                                    'shift_even_scale',shift_ibis_nb_even,$
                                    'optical_shift',   shift_target_nb, $
                                    'wl_to_nb_drift',  wl_to_nb_drift, $
                                    'wl_drifts',       wl_optical_drifts, $
                                    'atm_dispersion',  atm_dispersion_nb, $
                                    'filter_ids',      filter_ids)
    END
ENDCASE

RETURN,output_info

END)

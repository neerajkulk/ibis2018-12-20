restore_verbose   = 0

; ------------------------------------------------------------
;      Define Dataset and Paths
; ------------------------------------------------------------

day_dir      = '20Dec2018'
day_id       = '20181220'
day_tag      = STRMID(day_dir,0,5)

base_dir     = '/SMdata4/IBIS/Dec2018/'    ;i changed the directory, let's see if it works...

ave_ser_dir  = base_dir + day_dir + '/ibis/averaged_series/'

; ------------------------------------------------------------
;      Load Dark Calibration Images
; ------------------------------------------------------------
dark_wl_files      = file_search(ave_ser_dir,'DarkCalibration.whitelight.combineall.' + day_id + '*.series.ave.sav', count=num_darks)     
wl_darks_all       = FLTARR(1000,1000,num_darks)
wl_darks_exptime   = FLTARR(num_darks)
wl_darks_imcount   = INTARR(num_darks)
wl_darks_timerange = DBLARR(2,num_darks)

FOR nn=0,num_darks-1 DO BEGIN
    RESTORE,Verbose=0,dark_wl_files[nn]
    wl_darks_all[*,*,nn]      = Series_Ave
    wl_darks_exptime[nn]      = MEDIAN(FLOAT(series_ave_input[6,*,0]))
    wl_darks_imcount[nn]      = FIX(MEDIAN(Series_Cnt,/Even))
    validp                    = (WHERE(strlen(Series_Ave_Input[5,*,*]) GE 1))
    wl_darks_timerange[*,nn]  = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
ENDFOR

output_darkvar = 'wl_darks_' + day_tag
result = execute(output_darkvar + '= wl_darks_all')

save_variables = output_darkvar + ', dark_wl_files, wl_darks_timerange'
result = execute('SAVE,' + save_variables + ', FILENAME=''whitelight.darks.' + day_dir + '.sav'',/COMP')

; ------------------------------------------------------------

; ------------------------------------------------------------
;      Load Flat Calibration Images
; ------------------------------------------------------------
flat_wl_files      = file_search(ave_ser_dir, 'FlatFieldCalibration.whitelight.combineall.' + day_id + '*.series.ave.sav', count=num_flats)     
wl_flats_all       = FLTARR(1000,1000,num_flats)
wl_flats_exptime   = FLTARR(num_flats)
wl_flats_imcount   = INTARR(num_flats)
wl_flats_timerange = DBLARR(2,num_flats)

FOR nn=0,num_flats-1 DO BEGIN
    RESTORE,Verbose=0,flat_wl_files[nn]
    TVSCL,Series_Ave
    wl_flats_all[*,*,nn]      = Series_Ave
    wl_flats_exptime[nn]      = MEDIAN(FLOAT(Series_Ave_Input[6,*,0]))
    wl_flats_imcount[nn]      = FIX(MEDIAN(Series_Cnt,/Even))
    validp                    = (WHERE(strlen(Series_Ave_Input[5,0,*]) GE 1))
    wl_flats_timerange[*,nn]  = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
    dark_match                = get_closest(REBIN(wl_darks_timerange,1,num_darks),MEAN(wl_flats_timerange[*,nn]))
    wl_flats_all[*,*,nn]     -= wl_darks_all[*,*,dark_match]
    wl_flats_all[*,*,nn]     /= MEAN(wl_flats_all[*,*,nn]) 
ENDFOR

output_flatvar = 'wl_flats_' + day_tag
result = execute(output_flatvar + '= wl_flats_all')

save_variables = output_flatvar + ', flat_wl_files, wl_flats_timerange'
result = execute('SAVE,' + save_variables + ', FILENAME=''whitelight.gains.' + day_dir + '.sav'',/COMP')

END

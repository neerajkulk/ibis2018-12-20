; ----------------------------------------------------------------------
; get information on observing time and conditions                                                                                  
; ----------------------------------------------------------------------                                                            
; we need some information on the acquisition time of each image file                                                               
; we'll first see if this has already been parsed out of the FITS headers                                                           
; and stored in a IDL .sav file                                                                                                     
basedir                = '/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/'                                                           
dst_header_info_file   = 'wl.20Dec2018.dst.seeing.header.info.new.sav'  ; looks for this file in working directory                                                  
IF FILE_TEST(dst_header_info_file) NE 1 THEN BEGIN                                                                                 
    ; if that file is not available, we will walk through the FITS files and extensions                                             
    ;     and parse the information from the headers                                                                                
    ; this procedure does that and creates the missing save file, which will then be restored below                                 
    extract_dst_header_info, wl_data_files, output_save_file='wl.20Dec2018.dst.seeing.header.info.new.sav'                         
ENDIF                         
                                                                                                                                    
restore,/verbose,'wl.20Dec2018.dst.seeing.header.info.new.sav'                                                                                             
; provides - date_obs_jd,date_end_jd                                                                                                
sequence_times         = MEAN((date_obs_jd + date_end_jd)/2.,dim=1)                                                                
sequence_reftime       = LONG(MIN(sequence_times)) - 0.5d
sequence_times_fracday = sequence_times - ; @@@ is this specific to observations?                                                                               

; calculate atmospheric dispersion in different wavelength
; 680 nm is the central wavelength of filter used in the whitelight channel

;atm_disp_times = findgen(1080)*30/86400. + (julday(12,20,2018,14.d,00,00)
atm_disp_calc     = atmospheric_refraction_shifts(sequence_times,[617.3,656.3,680.0,769.9,854.2])

; instead, one can simply restore this SAVE file, which contains the precalculated dispersion information
;RESTORE,/Verbose,'calibrations/atmospheric.dispersion.calc.20Dec2018.sav'

; pull out refraction values decomposed into x- and y-shifts (in solar heliocentric coordinates)
dispcalc_sfts     = [REFORM(atm_disp_calc.SFTS_HELIOCENT_EW,1,10,5),REFORM(atm_disp_calc.SFTS_HELIOCENT_NS,1,10,5)]

; take difference to get theoretical shifts between different IBIS filter bands
; this provides the dispersion between pairs of wavelengths

sftcalc_6173_6563 = dispcalc_sfts[*,*,0] - dispcalc_sfts[*,*,1]
sftcalc_6173_7699 = dispcalc_sfts[*,*,0] - dispcalc_sfts[*,*,3]
sftcalc_6173_8542 = dispcalc_sfts[*,*,0] - dispcalc_sfts[*,*,4]
sftcalc_6563_7699 = dispcalc_sfts[*,*,1] - dispcalc_sfts[*,*,3]
sftcalc_6563_8542 = dispcalc_sfts[*,*,1] - dispcalc_sfts[*,*,4]
sftcalc_7699_8542 = dispcalc_sfts[*,*,3] - dispcalc_sfts[*,*,4]

; take difference to get theoretical shifts between wavelength of whitelight 
; images and the wavelength of different IBIS filters

sftcalc_wl_6173 = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,0]
sftcalc_wl_6563 = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,1]
sftcalc_wl_7699 = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,3]
sftcalc_wl_8542 = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,4]

; ----------------------------------------------------------------------
; identify all available files
; ----------------------------------------------------------------------

ibisdata_basedir = '/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/'
wl_data_files = file_search(ibisdata_basedir + 'whitelight/ScienceObservation/20181220_160454*/s*fits')
num_files     = N_ELEMENTS(wl_data_files)
nb_data_files = wl_data_files
startcut      = STRPOS(wl_data_files[0],'ScienceObservation')
for nn=0,num_files-1 do begin & nb_data_files[nn] = ibisdata_basedir + 'spectral/' + STRMID(nb_data_files[nn],startcut,1000) & endfor

nb_data_files = file_search(ibisdata_basedir + 'spectral/ScienceObservation/20181220_160454*/s*fits')
   
; ----------------------------------------------------------------------
; now calculate observed offsets
; ----------------------------------------------------------------------

; Choose narrowband wavelengths where image is mostly smooth continuum and 
; a relaible comparison to the whitelight image can be made
;ca8542_cont = [] cant really see the photospheric granulation at 8542 and halpha
;h6563_cont = []
k7699_cont = indgen(21,start = 44, increment = 1)
fe6173_cont = indgen(5, start = 71, increment = 1)

; wWe can use this definition to calculate the image shifts only for those narrowband images
; taken at near-continuum wavelengths that can be most reliably compared to the whitelight images.
; The advantage is that it will be faster if it doesn't load and calculate offsets for
; unnecessary images.
align_extn_all = [k7699_cont, fe6173_cont]
; Instead, this will calculate the image shifts for all wavelength positions
;align_extn_all = indgen(56)

num_extn_tot  = 76
num_extn_calc = N_ELEMENTS(align_extn_all)

; run through each image and calculate shifts with respect to whitelight
; -- this might take a while ...
shifts_wl2nb = FLTARR(2,num_extn_tot,num_files)
for filn=0,num_files-1 do begin & $
    for extn=0,num_extn_calc-1 do begin & $ 
        file_extension = align_extn_all[extn] & $
        wl_extension=file_extension & $
	wlim = load_calibrate_image(wl_data_files[filn],wl_extension+1,'wl',cal_directory='calibrations',/keep,verbose=0) & $                         
	nbim = load_calibrate_image(nb_data_files[filn],file_extension+1,'nb',cal_directory='calibrations',/keep,verbose=0, apply_dist=1) & $
	nbim_sfit = sfit(nbim,3) & wlim_sfit = sfit(wlim,3) & $
	shifts_wl2nb[*,file_extension,filn] = xyoff(nbim/nbim_sfit,wlim/wlim_sfit,768,768,/quiet) & $
    endfor & $
    IF (filn MOD 20) EQ 0 THEN print,systime(),filn & $
endfor 

; extract between whitelight and IBIS narrowband filters
offsets_wl_7699   = (MEDIAN(shifts_wl2nb[*,k7699_cont,*],dim=2))*0.096 
offsets_wl_6173   = (MEDIAN(shifts_wl2nb[*,fe6173_cont,*],dim=2))*0.096


; calculate shifts between different combinations of IBIS narrowband filters

 offsets_7699_6173 = (offsets_wl_7699 - offsets_wl_6173)

; subtract off the theoretical atmospheric dispersion values from the observed shifts,
; in order to isolate the (optical) shifts due to the drift between the whitelight
; and narrowband channel
drifts_wl2nb = FLTARR(2,1048,2)

drifts_wl2nb[*,*,0] = offsets_wl_6173 - sftcalc_wl_6173
drifts_wl2nb[*,*,1] = offsets_wl_7699 - sftcalc_wl_7699

; apply a median smoothing filter and correct some time steps with spurious values
for wv = 0,1 do begin & for dir = 0,1 do begin & $
    drifts_wl2nb[dir,*,wv] = MEDIAN(REFORM(drifts_wl2nb[dir,*,wv]),3) & $
endfor & endfor

;drifts_wl2nb[1,1034:1037,0] = MEDIAN(drifts_wl2nb[1,1024:1033,0])
;drifts_wl2nb[1,826,2]       = MEDIAN(drifts_wl2nb[1,816:825,0])
;drifts_wl2nb[1,835,2]       = MEDIAN(drifts_wl2nb[1,827:834,0])
;drifts_wl2nb[1,1035,2]       = MEDIAN(drifts_wl2nb[1,1025:1034,0])

;drifts_wl2nb[0,1034:1038,1] = MEDIAN(drifts_wl2nb[0,1024:1033,1])
;drifts_wl2nb[0,1036,2]      = MEDIAN(drifts_wl2nb[0,1031:1041,2])
;drifts_wl2nb[1,1033:1038,1] = MEDIAN(drifts_wl2nb[1,1024:1033,1])

;drift_fits_param = FLTARR(2,6)
;drift_fits       = FLTARR(2,1048)
;drift_fits_param[0,*] = poly_fit(sequence_times_fracday,rebin(drifts_wl2nb[0,*,*],1,1048,1),5,yfit=fit_xdrift)
;drift_fits_param[1,*] = poly_fit(sequence_times_fracday,rebin(drifts_wl2nb[1,*,*],1,1048,1),5,yfit=fit_ydrift)
;drift_fits[0,*]       = fit_xdrift
;drift_fits[1,*]       = fit_ydrift

;print,"The following lines can be copied into load_calibration_info" 
;print,"to provide the definition of the optical drift:"
;print,drift_fits_param,format='("    wl_to_nb_drift = [", 11( F10.2, ", "),F10.2,"]")' 
;print,"    wl_to_nb_drift = REFORM(wl_to_nb_drift,2,6)"


END

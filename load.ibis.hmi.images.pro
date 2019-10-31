
; before starting, make sure the proper calibration info routine is loaded!
; if this isn't in the IDL library path (or current directory), one may need to specify 
; the path to the program
;.r load_calibration_info.25Apr2019
; also, compile the destretch code before executing this the first time
; .r reg

IF N_ELEMENTS(define_variables) EQ 0 THEN define_variables = 1

; what level of verbosity to use when loading the IBIS images
load_verbosity = 0
; sometime it may not be necessary to load the corresponding narrowband file
; and since it can take some time, we can skip it by setting this flag to zero
load_nb        = 0

; ------------------------------------------------------------------
; determine which data is available and define what ranges to process
; ------------------------------------------------------------------
wl_files     = file_search('/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/whitelight/ScienceObservation/20181220_145901/','*fits')
; there are extra spectral files, so we'll construct an array of matching narrowband
; files by simply changing the name of the leading directory
dir_name_len = strlen('whitelight')
nb_files     = 'spectral' + strmid(wl_files,dir_name_len,999)
num_files    = N_ELEMENTS(wl_files)

wlszx = 1000 & wlmidx = wlszx/2.
wlszy = 1000 & wlmidy = wlszy/2.
; this is the average time (in seconds) that it takes to go through a full scan of all the spectral lines
seq_cadence            = 18.0

; define which images to load
seqstart = 0
;seqend   = 0
; uncomment or change this to do more than one sequence
seqend   = num_files - 1


; how many extensions for each file to process
; start at the first extensions (extensions are 1-indexed, not 0)
exten_start      = 1
;num_exten        = 1
; see how many extensions there are in the files
fits_open,wl_files[1],fcb
max_exten        = fcb.nextend     ; should be 56 for the 25 April 2019 data
fits_close,fcb
; uncomment this (and change the value) if you want to do more than just the first extension

num_exten       = max_exten

wl_dstr_all = FLTARR(1000,1000,num_exten,num_files) ; array to store the destretched images

; ------------------------------------------------------------------
; define output variables, if necessary
; ------------------------------------------------------------------
define_variables = 1; set this so that arrays and kernels are initialized 

IF define_variables THEN BEGIN
    IF N_ELEMENTS(wl2hmi_sfts) LE 2 * num_files THEN BEGIN
        wl2hmi_sfts = FLTARR(2,num_exten,num_files)
    ENDIF

    ; how big of a kernel to use for destretching
    destretch_kernel = BYTARR(16,16)
    ; do a test destretch so we can determine how big the grid of destretch vectors will be
    test_destr = reg(fltarr(wlszx,wlszy),fltarr(wlszx,wlszy),destretch_kernel, rdisp=rdispout, disp=dispout)
    disp_sz    = SIZE(dispout)
    disp_all   = FLTARR(disp_sz[1],disp_sz[2],disp_sz[3],num_exten, num_files)
    rdisp_all  = disp_all
    
    ; shouldn't have to revisit this block after the first execution
    ; but if you do need to, just set define_variables = 1 from the command line 
    ; before executing
    define_variables = 0
ENDIF

; ------------------------------------------------------------------
; load HMI header and determine some timing information about the cube
; ------------------------------------------------------------------
; again, assuming the (relative) path to the HMI data directory
;hmi_cont_cube_file = '../sdo/target/cubes/hmicont.fits'
hmi_cont_cube_file  = '/SMdata1/nkulkarni/sdo_data/middle/target/cubes/hmicont.fits'
hmiim_cube          = readfits(hmi_cont_cube_file, hmihdr, /Silent)

; the HMI data is stored as a 3-D data cube
; we need to find out which image in that cube corresponds to the image/sequence loaded above
; first we find the time for the first image in the cube
hmi_starttime_string = sxpar(hmihdr,'STARTTIM')
; need to trim the final character off of the string
hmi_starttime_string = strmid(hmi_starttime_string,0,STRLEN(hmi_starttime_string)-1)
hmi_starttime_jd     = fits_date_convert(hmi_starttime_string)
; then we determine how many seconds between images in the cube
hmi_time_step        = sxpar(hmihdr,'CADENCE')
hmi_image_scale      = sxpar(hmihdr,'CDELT1')

; ------------------------------------------------------------------
; start main loop through sequences and images
; ------------------------------------------------------------------

FOR seqnum = seqstart,seqend DO BEGIN
    IF num_exten EQ max_exten THEN PRINT,systime(),seqnum,' - ',wl_files[seqnum] $ 
    ELSE IF (seqnum MOD 50) EQ 0 THEN PRINT,systime(),seqnum
    ;seqnum = 100
    
    FOR exten_cnt = 0,num_exten-1 DO BEGIN
        exten = exten_cnt + exten_start
	; assumes WL and NB calibrations are in a subdirectory 'calibrations/' of the current directory 
	wlim     = load_calibrate_image(wl_files[seqnum],exten,'wl',header=wl_header, cal_info=cal_info_wl)
	IF load_nb THEN BEGIN
	    nbim7090 = load_calibrate_image(nb_files[seqnum],exten,'nb',header=nb_header, cal_info=cal_info_nb)
	ENDIF

	extension_header_start = where(STRMATCH(wl_header,'*XTENSION=*'))
	; this the acquisition time of the loaded image
	wl_image_time          = fits_date_convert(sxpar(wl_header[extension_header_start:*],'DATE-OBS',/silent))
	; taken from the primary header, this is the time when the sequence started
	wl_seq_time            = fits_date_convert(sxpar(wl_header[0:extension_header_start],'DATE-OBS',/silent))
	; this is the middle time of the scan
	wl_seq_mid_time        = wl_seq_time + (seq_cadence/2.)/86400.

	; again, assuming the (relative) path to the HMI data directory
	;hmi_cont_cube_file = '../sdo/target/cubes/hmicont.fits'
	hmi_cont_cube_file = '/SMdata1/nkulkarni/sdo_data/disk_center/target/level1/target/cubes/hmicont.fits'
	hmiim = readfits(hmi_cont_cube_file,nslice=0,hmihdr,/Silent)

	; the HMI data is stored as a 3-D data cube
	; we need to find out which image in that cube corresponds to the image/sequence loaded above
	; first we find the time for the first image in the cube
	hmi_starttime_string = sxpar(hmihdr,'STARTTIM')
	; need to trim the final character off of the string
	hmi_starttime_string = strmid(hmi_starttime_string,0,STRLEN(hmi_starttime_string)-1)
	hmi_starttime_jd     = fits_date_convert(hmi_starttime_string)
	; then we determine how many seconds between images in the cube
	hmi_time_step        = sxpar(hmihdr,'CADENCE')

	; here we determine the number of seconds between the whitelight image and the started
	; of the HMI cube
	; we use the sequence middle time, rather than the time of the specific image because that
	; way a constant HMI image will be found for all images in a sequence, avoiding any potential
	; changes in the HMI reference within a single scan (which could cause discontinuities in
	; the alignment).
	; Note that an IBIS scan takes ~12 seconds, while the HMI data were originally sampled at
	; a 45 second cadence (upsampled to 12 seconds to match AIA in the data cubes)
	wl2_hmi_time_diff    = (wl_seq_mid_time - hmi_starttime_jd)*86400
	; finally, we divide the number of seconds by the cadence to find the appropriate index
	; for the HMI image corresponding to the whitelight image
	hmi_image_num        = ROUND(wl2_hmi_time_diff / hmi_time_step)

	; now read in the matching HMI image
	; this can be slow if the image number is at the end of the file
	;hmiim = readfits(hmi_cont_cube_file,nslice=hmi_image_num,hmihdr,/Silent)
	; so we read in the whole cube before the loop and just extract the appropriate image
	hmiim = hmiim_cube[*,*,hmi_image_num]

	; scale the HMI image to the same plate scale as the IBIS data
	hmiim_size = size(hmiim)
	; determine number of pixels for scaled up image
	; need to know target IBIS plate scale = 0.096 arcsec/pixel
	hmiim_ibis_pix = ROUND(hmiim_size[1:2] * hmi_image_scale / 0.096)
	hmi_ibisscl     = CONGRID(hmiim,hmiim_ibis_pix[0],hmiim_ibis_pix[1], CUBIC=-0.5)
	
	; now cut down HMI image to the same number of pixels as the IBIS data, 
	; with a window that is approximately centered on the average IBIS field of view
	hmicut          = cal_info_wl.hmi_pixel_cutout
	hmi_ibisscl     = hmi_ibisscl[hmicut[0]:hmicut[1],hmicut[2]:hmicut[3]]

	; calculate rigid shifts between IBIS and HMI images
	wl_sft = xyoff(hmi_ibisscl, wlim, 768,768, /quiet)
	wl2hmi_sfts[*,exten_cnt,seqnum] = wl_sft

	; calculate destretch vectors
        ; before destretching, take out the bulk shifts since they might confuse the destretching if
        ; the a comparable to or greater than the size of the subfield/kernel
        ; two options to shift images 
        ;     - the first (SHIFT) causes the image to wrap around, which may be bad when destretching
        ;     - the second (ROT) properly fills the unknown borders with the array mean
	;wlim_sft = SHIFT(wlim,ROUND(wl_sft[0]), ROUND(wl_sft[1]))
	wlim_sft = ROT(wlim,0,1.0,wlmidx-ROUND(wl_sft[0]),wlmidy-ROUND(wl_sft[1]),missing=MEAN(wlim))
	wlim_destr = reg(wlim_sft, hmi_ibisscl, bytarr(16,16), rdisp=rdispout, disp=dispout)
	disp_all[*,*,*,exten_cnt,seqnum]  = dispout
	rdisp_all[*,*,*,exten_cnt,seqnum] = rdispout
	wl_dstr_all[*,*,exten_cnt,seqnum] = wlim_destr

    ENDFOR
    	SAVE,wl2hmi_sfts,disp_all,rdisp_all,wl_dstr_all,filename='145901.wl.to.hmi.alignment.params.20Dec2018.sav'

ENDFOR

END

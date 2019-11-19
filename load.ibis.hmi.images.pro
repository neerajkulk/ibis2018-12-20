
; before starting, make sure the proper calibration info routine is loaded!
; if this isn't in the IDL library path (or current directory), one may need to specify 
; the path to the program
;.r load_calibration_info.25Apr2019
; also, compile the destretch code before executing this the first time
; .r reg
; .r reg_loop

IF N_ELEMENTS(define_variables) EQ 0 THEN define_variables = 1

verbose        = 1
; what level of verbosity to use when loading the IBIS images
load_verbosity = 0

; sometime it may not be necessary to load the corresponding narrowband file
; and since it can take some time, we can skip it by setting this flag to zero
load_nb        = 1

; The total number of destretched images may be too large to hold in memory
; since we could have tens of thousands of images at 4 MB each.
; So in order to save the destretched arrays, we can write them out to a file. 
; This involves opening a file, putting on a rudimentary FITS header, and then
; dumping each destretched array into that file as it is processed.
write_destr_arrays = 1

; ------------------------------------------------------------------
; determine which data is available and define what ranges to process
; ------------------------------------------------------------------
wl_files     = file_search('/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/whitelight/ScienceObservation/20181220_145901/','*fits')
; there are extra spectral files, so we'll construct an array of matching narrowband
; files by simply changing the name of the leading directory
dir_name_len = strlen('whitelight')
nb_files = file_search('/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/spectral/ScienceObservation/20181220_145901/','*fits')
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
fits_read,fcb,data,hdrorig_wl,exten=1,/header_only
fits_close,fcb
; uncomment this (and change the value) if you want to do more than just the first extension
num_exten        = max_exten
wlpos            = 0

; if we want to write out the destretched arrays, we will prepare the suitable FITS header and
;     open the output files for writing.
IF write_destr_arrays GE 1 THEN BEGIN
    ; define the name of the output file(s)
    filename_out_base = '.destr.all.20Dec2018.v1.fits'
    wl_filename_out = '/SMdata1/nkulkarni/destretching/whitelight' + filename_out_base
    
    ; This is a list of the FITS keywords we will copy from the original data files.
    copykeywd = ['ORIGIN','TELESCOP','INSTRUME','CHANNEL','DETECTOR','FILEORIG','GAIN_PRE','WL_PRFLT','EXPTIME']
    numkeywd  = N_ELEMENTS(copykeywd)

    ; make a new FITS header of the size appropriate for the expected data array output
    mkhdr,wlhdr_out,4,[wlszx,wlszy,max_exten,seqend-seqstart+1]
    hdrlen = N_ELEMENTS(wlhdr)
    FOR nn=0,numkeywd-1 DO sxaddpar,wlhdr_out,copykeywd[nn],sxpar(hdrorig_wl,copykeywd[nn]),before='COMMENT'
    sxaddpar,wlhdr_out,'HISTORY','Destretched using load.ibis.hmi.images.pro on ' + SYSTIME(),before='COMMENT'
    writefits,wl_filename_out,nulldata,wlhdr_out
    OPENU,wllun,wl_filename_out,/swap_if_little_endian,/get_lun
    point_lun,wllun,2880
    point_lun,-wllun,wlpos

    IF load_nb THEN BEGIN
        nb_filename_out = '/SMdata1/nkulkarni/destretching/narrowband' + filename_out_base

        fits_open,nb_files[1],fcbnb
        fits_read,fcbnb,data,hdrorig_nb,exten=1,/header_only
        fits_close,fcbnb
        mkhdr,nbhdr_out,4,[wlszx,wlszy,max_exten,seqend-seqstart+1]
        FOR nn=0,numkeywd-1 DO sxaddpar,nbhdr_out,copykeywd[nn],sxpar(hdrorig_nb,copykeywd[nn]),before='COMMENT'
        sxaddpar,nbhdr_out,'HISTORY','Destretched using load.ibis.hmi.images.pro on ' + SYSTIME(),before='COMMENT'
        writefits,nb_filename_out,nulldata,nbhdr_out
        OPENU,nblun,nb_filename_out,/swap_if_little_endian,/get_lun
        point_lun,nblun,2880
    ENDIF

ENDIF
        
; ------------------------------------------------------------------
; define output variables, if necessary
; ------------------------------------------------------------------
IF define_variables THEN BEGIN
    IF verbose GE 1 THEN PRINT,'Defining Variables'
    IF N_ELEMENTS(wl2hmi_sfts) LE 2 * num_files THEN BEGIN
        wl2hmi_sfts = FLTARR(2,num_exten,num_files)
    ENDIF

    ; how big of a kernel to use for destretching
    destretch_kernel = BYTARR(32,32)
    destretch_kernel = [64,32,16]
    ; do a test destretch so we can determine how big the grid of destretch vectors will be
    ;test_destr = reg(fltarr(wlszx,wlszy),  fltarr(wlszx,wlszy),  destretch_kernel, rdisp=rdispout, disp=dispout)
    test_destr = reg_loop(fltarr(wlszx,wlszy),fltarr(wlszx,wlszy), destretch_kernel, rdisp=rdispout, disp=dispout)
    disp_sz    = SIZE(dispout)
    disp_all   = FLTARR(disp_sz[1],disp_sz[2],disp_sz[3],num_exten, num_files)
    rdisp_all  = disp_all
    
    wlim_destr = FLTARR(1000,1000,num_exten)
    IF load_nb THEN nbim_destr = FLTARR(1000,1000,num_exten)
    
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
    IF num_exten EQ max_exten THEN PRINT,systime(),seqnum,' - ',wl_files[seqnum],' - ',wlpos $ 
    ELSE IF (seqnum MOD 50) EQ 0 THEN PRINT,systime(),seqnum
    ;seqnum = 100
    
    FOR exten_cnt = 0,num_exten-1 DO BEGIN
        exten = exten_cnt + exten_start
	; assumes WL and NB calibrations are in a subdirectory 'calibrations/' of the current directory 
	wlim     = load_calibrate_image(wl_files[seqnum],exten,'wl',header=wl_header,verbose=load_verbosity, cal_info=cal_info_wl)
	IF load_nb THEN BEGIN
	    nbim = load_calibrate_image(nb_files[seqnum],exten,'nb',header=nb_header,verbose=load_verbosity, cal_info=cal_info_nb, apply_dist=1)
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
        hmi_cont_cube_file = '/SMdata1/nkulkarni/sdo_data/middle/target/cubes/hmicont.fits'
        
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
	hmi_ibisscl     = CONGRID(hmiim, hmiim_ibis_pix[0], hmiim_ibis_pix[1], CUBIC=-0.5)
	
	; now cut down HMI image to the same number of pixels as the IBIS data, 
	; with a window that is approximately centered on the average IBIS field of view
	hmicut          = cal_info_wl.hmi_pixel_cutout
	hmi_ibisscl     = hmi_ibisscl[hmicut[0]:hmicut[1],hmicut[2]:hmicut[3]]

	; calculate rigid shifts between IBIS and HMI images
	wl_sft = xyoff(hmi_ibisscl, wlim, 768,768, /quiet)
	wl2hmi_sfts[*,exten_cnt,seqnum] = wl_sft

	; calculate destretch vectors
        ; before destretching, take out the bulk shifts since they might confuse the destretching if
        ; they are comparable to or greater than the size of the subfield/kernel.
        ; There are two options to shift images 
        ;     - the first (SHIFT) is simple but causes the image to wrap around, which may be bad when destretching
        ;     - the second (ROT) instead properly fills the unknown borders with the array mean
	;wlim_sft = SHIFT(wlim,ROUND(wl_sft[0]), ROUND(wl_sft[1]))
	wlim_sft                          = ROT(wlim,0,1.0,wlmidx-ROUND(wl_sft[0]),wlmidy-ROUND(wl_sft[1]),missing=MEAN(wlim))
	;wlim_destr[*,*,exten_cnt]         = reg(wlim_sft, hmi_ibisscl, bytarr(32,32), rdisp=rdispout, disp=dispout)
	wlim_destr[*,*,exten_cnt]         = reg_loop(wlim_sft, hmi_ibisscl, destretch_kernel, rdisp=rdispout, disp=dispout)
	
	disp_all[*,*,*,exten_cnt,seqnum]  = dispout
	rdisp_all[*,*,*,exten_cnt,seqnum] = rdispout


	IF load_nb THEN BEGIN
            nbim_sft                          = ROT(nbim,0,1.0,wlmidx-ROUND(wl_sft[0]),wlmidy-ROUND(wl_sft[1]),missing=MEAN(nbim))
            nbim_destr[*,*,exten_cnt]         = doreg(nbim_sft, rdispout, dispout)
        ENDIF

    ENDFOR
    
    ; if requested dump the destretched array(s) into the FITS file as raw values
    IF write_destr_arrays GE 1 THEN BEGIN
        WRITEU,wllun,wlim_destr
        IF load_nb THEN WRITEU,nblun,nbim_destr
    ENDIF
    
    ;STOP
    FLUSH,wllun,nblun
    point_lun,-wllun,wlpos
    point_lun,-nblun,nbpos
    ;print,wlpos,nbpos
    
    IF ((seqnum+1) MOD 50) EQ 0 THEN BEGIN
        SAVE,wl2hmi_sfts,disp_all,rdisp_all,filename='wl.to.hmi.alignment.params.20Dec2018.level2.sav'
    ENDIF

ENDFOR

; we are done writing so close any of the destretched array output files that were open
IF write_destr_arrays GE 1 THEN FREE_LUN,wllun,nblun

END

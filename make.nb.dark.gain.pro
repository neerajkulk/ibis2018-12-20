; script to generate an array of dark and flat calibrations for the narrowband
; channel on a specific day.
;
; Assumes: * that average_series_images has already been run for 
;            the whitelight Dark and Flat series one the given day
;          * at least one Dark and Flat exists for the given day
;          * that the data is all 1000 x 1000 pixels
;          * that the darks have the same exposure time as the flat frames,
;            or that if the exposure times are different that that doesn't
;            make a difference in the measured signal.
;          * dark frame closest-in-time to each flat field image is the 
;            most appropriate correction to apply to said flat

; do you want the RESTORE commands to list all the extracted variables?
restore_verbose   = 0
; do you want to automatically save the generated darks and flats?
do_save_data      = 1
; should the resulting gain arrays be duplicated into an array
;    that contains the series id? This allows it to be compared 
;    to the gain table generated from other series, but increases
;    the amount of data in memory
do_label_data   = 1

; ------------------------------------------------------------
; set day, filter, and series index to process
; ------------------------------------------------------------
day_dir      = '20Dec2018/'
day_id       = '20181220'

base_dir     = '/SMdata4/kreardon/IBIS/Dec2018/'
ave_ser_dir  = base_dir + day_dir + '/ibis/averaged_series/'

; select which filter to process
;flat_type    = 'halpha' & filter_name = 'Halpha'
;flat_type    = 'ca8542'  & filter_name = 'Ca8542'
;flat_type    = 'na5896' & filter_name = 'Na5896'
;flat_type    = 'fe6173' & filter_name = 'Fe6173'
flat_type    = 'k7699'  & filter_name = 'K7699'
;flat_type    = 'o7774'  & filter_name = 'O7774'
;flat_type    = 'fe5434' & filter_name = 'Fe5434'


; define list of series which contain flats for each filter
; this might vary if dedicated flats were done for each filter,
; instead of doing all the filters in the same flat.
; Right now this needs to be determined by user inspection of
; the logs or file contents, but yes, this determination 
; should be automated...

k7699_flat_series    = ['20181220_194457']
;ca8542_flat_series   = ['20181220_161321'] 

;na5896_flat_series   = ['20170421_192733','20170421_180800']
;k7699_flat_series    = ca8542_flat_series
;fe5434_flat_series   = ca8542_flat_series

;fe6173_flat_series    = ['20181220_165244']

; define the index of the flat series list to process
flat_index  = 0

; not currently used, but it could be used to keep all calibration
; files for a given dataset in the same location
calibration_dir = File_Dirname(Routine_Filepath(/Either),/Mark_Directory) + 'calibration_files/'
FORWARD_FUNCTION filter_fringe_fft

; if the the average files are generated at the DST (using the data
; on /dstdata) that are labeled "andor1" and "andor2" instead of 
; "spectral" and "whitelight" but this isn't usually the case
narrowband_tag  = 'spectral'

; ------------------------------------------------------------

; waveref - a nominal central wavelength for each filter
; flat_linepos_wavestep - the size of the spectral step for the averaged profile
; prefilter_profile_cor - prefilter profile correction type - there may be 
;                         different options to refine prefilter correction

CASE flat_type OF 
    'halpha': BEGIN
        waveref = 6562.8
        prefilter_profile_cor = '2D'
        flat_linepos_wavestep  = 0.04
    END
    'fe6173': BEGIN
	waveref = 6173.34
	prefilter_prof_cor = '1D'
	flat_linepos_wavestep = 0.01
    END
    'ca8542': BEGIN
        waveref = 8542.14
        prefilter_profile_cor = '1D'
        flat_linepos_wavestep = 0.01
    END
    'k7699': BEGIN
	waveref = 7698.96
	prefilter_profile_cor = '1D'
	flat_linepos_wavestep = 0.01
    END
    'na5896': BEGIN
        waveref = 5896.0
        prefilter_profile_cor = '1D'
        flat_linepos_wavestep  = 0.005
    END
ENDCASE

; ------------------------------------------------------------
;      Load Dark Calibration Images
; ------------------------------------------------------------
; find all averaged dark files
dark_nb_files       = file_search(ave_ser_dir,'DarkCalibration.spectral.combineall.' + day_id + '*.series.ave.sav', count=num_darks)     

;hard-coded HACK to remove arrays that are not 1000x1000

dark_nb_files = dark_nb_files[0:5]

num_darks = n_elements(dark_nb_files)

Print, 'num darks is: '
Print, num_darks
nb_darks            = FLTARR(1000,1000,num_darks)
nb_darks_exptime    = FLTARR(num_darks)
nb_darks_imcount    = INTARR(num_darks)
nb_darks_timerange  = DBLARR(2,num_darks)

FOR nn=0,num_darks-1 DO BEGIN
    RESTORE,Verbose=0,dark_nb_files[nn]
    nb_darks[*,*,nn]        = Series_Ave
    nb_darks_exptime[nn]    = MEDIAN(FLOAT(series_ave_input[6,*,0]))
    nb_darks_imcount[nn]    = FIX(MEDIAN(Series_Cnt,/Even))
    validp                  = (WHERE(strlen(Series_Ave_Input[5,*,*]) GE 1))
    nb_darks_timerange[*,nn]  = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
ENDFOR

; ------------------------------------------------------------
; ------------------------------------------------------------

; select appropriate series list for chosen filter
flat_series_res = EXECUTE('flat_series_use = ' + flat_type + '_flat_series')
flat_series     = flat_series_use[flat_index]
series_id       = STRMID(flat_series,STRPOS(flat_series,day_id)+9,4)
PRINT,'Processing Series ID ' + series_id

; a list of bad pixels
; this worked okay when the bad pixels were few and static.
; that might not be the case for newer data...
bad_pix         = [[158LL,156,157,158,156,157,257,134,621,893],$
                   [289,  290,290,290,291,291,207,341,100,842]]
; convert 2-D pixel coordinates to 1D coordinates
bad_pix_idx     = bad_pix[*,1]*1000LL + bad_pix[*,0] 

flat_file_nb = File_Search(ave_ser_dir,'FlatFieldCalibration.' + narrowband_tag + '*' + flat_series + '*.sav',count=num_flat_nb)

        ; Restore averaged flat series, get it's time, and find the closest dark calibration
        PRINT,flat_file_nb
        RESTORE,flat_file_nb, Verbose=restore_verbose
        nb_flat_timerange = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
        dark_match = get_closest(REBIN(nb_darks_timerange,1,num_darks),MEAN(nb_flat_timerange))
        
        ; read in some information from the saved headers
        numhdr         = N_ELEMENTS(HEADERS_ALL[0,*,0])
        wavelnth_all   = FLTARR(numhdr)
        filter_num_all = BYTARR(numhdr)
        fp1_volt_all   = FLTARR(numhdr)
        FOR nn=0,numhdr-1 DO BEGIN
            wavelnth_all[nn]     = FLOAT(sxpar(HEADERS_ALL(*,nn,0),'WAVELNTH'))
            fp1_volt_all[nn]     = FLOAT(sxpar(HEADERS_ALL(*,nn,0),'FP1_VOLT'))
            filter_num_all[nn]   = BYTE(sxpar(HEADERS_ALL(*,nn,0),'FILTER'))
        ENDFOR
        ; find all the wavelengths associated with the wavelength of the filter of interest
        ; i.e. within 5 Å of the central wavelength
        wavematch = WHERE(ABS(wavelnth_all - waveref) LE 5)
        filter_id = MEDIAN(filter_num_all[wavematch])
        
        filter_match = WHERE(FIX(series_ave_input[2,*,0]) EQ filter_id)
    
        flat_seq         = REFORM(series_ave[*,*,filter_match])
        flat_seq_size    = SIZE(flat_seq)
        xsz              = flat_seq_size[1]
        ysz              = flat_seq_size[2]
        num_wave         = flat_seq_size[3]
        
        exptime                     = MEAN(FLOAT(SERIES_AVE_INPUT[6,filter_match,*])) * 1000
        flat_waves                  = REFORM(FLOAT(SERIES_AVE_INPUT[3,filter_match,0]))
        flat_ave_mask               = BYTARR(xsz,ysz)
        flat_ave_mask(0:xsz-1,0:ysz-1) += 1
        flat_ave_mask_r200          = RADIAL_DISTANCES([1,xsz,ysz],[500,500]) LE 600
        flat_ave_mask_r600          = RADIAL_DISTANCES([1,xsz,ysz],[500,500]) LE 600
        flat_ave_mask_full          = RADIAL_DISTANCES([1,xsz,ysz],[500,500]) LE 1000
        
        flat_power_sum              = FLTARR(xsz,ysz)
        flat_seq_nosurf_ave         = FLTARR(xsz,ysz)
        flat_seq_surf               = FLTARR(xsz,ysz,num_wave)

        ; variables to define which filter should be used for the FFT filtering
        fft_filt_halpha = 0
        fft_filt_ca8542 = 0
        fft_filt_na5896 = 0
        fft_filt_k7699  = 0
        fft_filt_fe5434 = 0
        ; set the flag for the chosen prefilter
        res = EXECUTE('fft_filt_' + filter_name + ' = 1')
        
        ; load the reference prefilter profile
        prefilter_file = 'prefilter.' + STRING(ROUND(waveref),FORMAT='(I4.4)') + '.reference.profile.*.sav'
        prefilter_file_match = File_Search('~kreardon/git/IBIS/prefilter_reference/',prefilter_file)
        RESTORE,Verbose=restore_verbose,prefilter_file_match[0]
        
        ; calculate the prefilter intensities at the specific wavelengths of the flat field series
        ; by interpolating the reference prefilter profile onto the observed wavelengths
        res = EXECUTE('prefilter_wave  = prefilt' + STRING(ROUND(waveref),FORMAT='(I4.4)') + '_ref_wvscl')
        res = EXECUTE('prefilter_trans = prefilt' + STRING(ROUND(waveref),FORMAT='(I4.4)') + '_ref_main')
        prefilter = INTERPOL(prefilter_trans, prefilter_wave, flat_waves)
        
        PRINT,'Calculating flats for filter ' + filter_name
        
        ; cycle through each wavelength position of the flat field series                             
        FOR imn=0,num_wave - 1 DO BEGIN
            ; remove dark + bias
            flat_seq[*,*,imn]      -= nb_darks[*,*,dark_match]
            
            ; Correct smearing due to transfer of array without shuttering
            ; maybe not really necessary unless exposures are very short
            ; or data was taken at the limb, for example.
            ; flat_seq[*,*,imn]       = frame_transfer_cor(flat_seq[*,*,imn], exptime, 0.001, binning=1, ratio=0.5, /readout)
            
            ; Apply correction for non-linear response of detector
            ; assumes narrowband flat field data was taken with Andor Camera ID 2748
            flat_seq[*,*,imn]       = ibis_linearity_correction_andor(flat_seq[*,*,imn], camera_id='2748', $
                                                                        data_includes_bias=0, verbose=0,lin_ptr=lin_ptr)
            
            ; replace bad (hot, cold) pixels with local median value
            temp_im                     = flat_seq[*,*,imn]
            temp_im_med                 = MEDIAN(temp_im,5)
            temp_im(bad_pix_idx)        = temp_im_med(bad_pix_idx)
            flat_seq[*,*,imn]           = temp_im

            ; remove high-frequency fringes
            flat_seq[*,*,imn]           = filter_fringe_fft(flat_seq[*,*,imn], halpha=fft_filt_halpha, ca8542=fft_filt_ca8542, $
                                                            na5896=fft_filt_na5896, fft_power=imtest_pow, fft_mask=imtest_mask)
            ; for reference, save average power for all flat field images
            flat_power_sum             += ABS(imtest_pow)
            
            ; determine and remove 2-D polynomial surface from each flat image
            ; this surface will be dominated by prefilter and blueshift variations
            ; this will aid in determining small-scale, static (i.e. not wavelength dependent) errors
            flat_seq_sfit               = SFIT(flat_seq[*,*,imn],2)
            flat_seq_nosurf_ave        += flat_seq[*,*,imn] / (flat_seq_sfit>100)
            flat_seq_surf[*,*,imn]      = flat_seq_sfit
                        
            ; divide out prefilter transmission at wavelength of each flat-field image
            flat_seq[*,*,imn]          /= prefilter[imn]

            tvscl,flat_seq[*,*,imn] 
        ENDFOR
    
        flat_seq_nosurf_ave     /= num_wave
        ; the following would be another way to get the "static" errors, but some residual parabolic errors remain 
        ; flat_seq_static_errors             = REBIN(flat_seq,xsz,ysz,1)
        flat_seq_static_errors   = flat_seq_nosurf_ave
        
        ; divide each flat field image by the map of static errors
        FOR imn=0,num_wave - 1 DO flat_seq[*,*,imn] /= flat_seq_static_errors
        
        ; fit line core of each profile in the map to get distribution of shifts
        ;   flat_wave_interp = output wavelength scale
        ;   line_start,line_ends = bounds on where to search for line minimum
        ;   flat_linepos_wavestep = wavelength sampling of evenly spaced profile
        ;   flat_ave_mask_full = mask defining region over which to perform line fitting
        ;   fit_points = number of wavelength points to fit around line center
        print,'flat seq is:'
	print,flat_seq
	print,'flat_waves is :'
        print,flat_waves

	flat_linepos             = find_line_shifts(flat_seq, flat_waves, wavelength_scale_interp=flat_wave_interp,$
                                       interp=flat_linepos_wavestep,radial=flat_ave_mask_full,fit_points=31)

        ; No fit a 2-D polynomial to the measured blueshifts to get the blueshift map
        ; the function returns the pixel position of the optical center of the beam (i.e. the
        ; minimum of the polynomial).
        ;   flat_ave_mask_r600 = mask for area over which to fit polynomial
        ;   flat_blueshift_coeff = polynomial coefficients
        ;   flat_bluefit = array of fitted blueshifts
        flat_optcent             = ibis_optical_center(flat_linepos(*,*,0)*flat_linepos_wavestep,$
                                       mask=flat_ave_mask_r600,coeff=flat_blueshift_coeff,$
                                       blueshift_fit_matrix=flat_bluefit)

        ; fix minimum position of blueshift map to zero
        flat_bluefit_use         = flat_bluefit - MIN(flat_bluefit(xsz/2.-100:xsz/2.+100,ysz/2.-100:ysz/2.+100))
        
        ; align and sum profiles (within masked area) acccording to blueshifts determined above
        ;   flat_waves = observed wavelengths
        ;   flat_ave_mask_r200 = mask for area over which to sum profiles
        ;   interp = wavelength step for interpolated profile
        flat_alignprof           = align_sum_profiles(flat_seq,flat_bluefit_use,flat_waves,$
                                       radial=flat_ave_mask_r200,interp=0.01,/PLOT)
               
        ; take average profile and reapply blueshifts to construct array of expected
        ; intensities at each [x,y,lambda] coordinate in the observed array.
        ;   flat_alignprof.alignedprof = aligned and averaged reference profile
        ;   flat_alignprof.wave_scale_orig = wavelength scale of flat field observations
        ;   flat_alignprof.wave_scale_interp = interpolated wavelength scale of reference profile
        ;   flat_ave_mask_full = area over which to calculated shifted profiles
        flat_flatalign           = construct_shifted_profile_array(flat_alignprof.alignedprof,$
                                       flat_bluefit_use,flat_alignprof.wave_scale_orig,$
                                       flat_alignprof.wave_scale_interp,radial=flat_ave_mask_full)

        ; construct flat field by dividing observed by expected profiles
        nb_gain                  = flat_seq /  flat_flatalign
        
        ; reapply the flat_seq_static_errors (which should just be multaplicative with
        ; the other system response variations
        FOR imn=0,num_wave - 1 DO BEGIN
            nb_gain[*,*,imn] *= flat_seq_static_errors 
        ENDFOR
        
        ; normalize the whole gain array to a mean of unity
        nb_gain /= MEAN(nb_gain)
        
        ; create a structure defining how gain array was constructed (i.e. what corrections
        ; were applied
        nb_gain_info = CREATE_STRUCT(  'flat_file',       flat_file_nb, $
                                        'wavelengths',    flat_waves, $
                                        'bad_pixels',     bad_pix, $
                                        'prefilter_profile', prefilter, $
                                        'frame_transfer', 'not corrected',$
                                        'linearity',      'corrected',$
                                        'fringe_fft',     'corrected',$
                                        'prefilter',      'not included' ) 
        
        IF do_label_data THEN BEGIN
            res = EXECUTE('nb_flat_seq_' + series_id     + ' = flat_seq')
            res = EXECUTE('nb_flat_static_' + series_id  + ' = flat_seq_static_errors')
            res = EXECUTE('nb_gain_' + series_id         + ' = nb_gain')
            res = EXECUTE('nb_gain_info_' + series_id    + ' = nb_gain_info')
            res = EXECUTE('flat_linepos_' + series_id      + ' = flat_linepos')
            res = EXECUTE('flat_bluefit_' + series_id    + ' = flat_bluefit')
            res = EXECUTE('flat_blueshift_coeff_' + series_id    + ' = flat_blueshift_coeff')
            res = EXECUTE('flat_alignprof_' + series_id      + ' = flat_alignprof')
            res = EXECUTE('prefilter_' + series_id    + ' = prefilter')
        ENDIF

        IF do_save_data THEN BEGIN
            variables = 'nb_darks, nb_darks_exptime, nb_darks_timerange, dark_nb_files'
            output_dark_filename = 'nb.darks.' + day_id + '.sav'
   	    saveoutput = EXECUTE('SAVE, /Compress, Filename='''' + output_dark_filename + '''',' + variables)

            variables = 'nb_gain_' + series_id + ', nb_gain_info_' + series_id + ', flat_linepos_' + series_id
            variables = variables + ', flat_bluefit_' + series_id + ', flat_alignprof_' + series_id
            variables = variables + ', nb_flat_static_' + series_id + ', flat_blueshift_coeff_' + series_id
            variables = variables + ', prefilter_' + series_id

            output_filename = filter_name + '.gain.blueshift.info.' + flat_series + '.sav'

            saveoutput = EXECUTE('SAVE, /Compress, Filename='''' + output_filename + '''',' + variables)
        ENDIF

END



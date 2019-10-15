;cd,'/SMdata4/IBIS/Dec2018/20Dec2018/ibis/'

xsz = 1000
ysz = 1000

; ********************** Load Darks *******************************
RESTORE,/VE,'DarkCalibration.spectral.combineall.20181220_163252.series.ave.sav'
dark_cal1 = SERIES_AVE

; ********************** Load Flats *******************************
RESTORE,/VE,'FlatFieldCalibration.spectral.byfilter+wave.20181220_161321.series.ave.sav'
flat8542_seq1_wv = reform(float(SERIES_AVE_INPUT[3,*,1]))

flat_8542_ser1613    = series_ave
for nn=0,20 do flat_8542_ser1613[*,*,nn] -= dark_cal1 

; ********************** CCD Level Corrections to Flats *******************************
flat_8542_ser1613cor = flat_8542_ser1613
for nn=0,20 do begin & flat_8542_ser1613cor[*,*,nn] = frame_transfer_cor(flat_8542_ser1613cor[*,*,nn], 50, 0.001, binning=1, ratio=0.5, /readout) & $
    flat_8542_ser1613cor[*,*,nn] = ibis_linearity_correction_andor(flat_8542_ser1613cor[*,*,nn], camera_id='2748', data_includes_bias=0, verbose=0,lin_ptr=lin_ptr) & $
    flat_8542_ser1613cor[*,*,nn] = filter_fringe_fft(flat_8542_ser1613cor[*,*,nn],ca8542=1,fft_power=imtest_pow, fft_mask=imtest_mask) & $ 
endfor

; ********************** Determine Static Errors in Flats *******************************
flat_8542_ser1613_fitres = FLTARR(xsz,ysz)
flat_8542_ser1613_fitave = FLTARR(xsz,ysz)
for nn=0,20 do begin & flatfit = sfit(flat_8542_ser1613cor[*,*,nn],2) & $
    flat_8542_ser1613_fitres += flat_8542_ser1613cor[*,*,nn]/flatfit & flat_8542_ser1613_fitave += flatfit & endfor
flat_8542_ser1613_fitave /= 21.
flat_8542_ser1613_fitres /= 21.
flat_8542_ser1613corfix = flat_8542_ser1613cor
for nn=0,20 do begin & flat_8542_ser1613corfix[*,*,nn] = flat_8542_ser1613cor[*,*,nn]/flat_8542_ser1613_fitres/flat_8542_ser1613_fitave & endfor

; ********************** Load Spectral Scan *******************************

RESTORE,/VE,'ScienceObservation.spectral.byfilter+wave.20181220_160919.series.ave.sav'
scan8542_full1 = series_ave
scan8542_full1_wv = reform(float(SERIES_AVE_INPUT[3,*,1]))
for nn=0,100 do scan8542_full1[*,*,nn] -= dark_cal1

; ********************** CCD Level Corrections to Spectral Scan *******************************
scan8542_full1_cor = scan8542_full1
for nn=0,100 do begin & scan8542_full1_cor[*,*,nn] = frame_transfer_cor(scan8542_full1_cor[*,*,nn], 50, 0.001, binning=1, ratio=0.5, /readout) & $ 
   scan8542_full1_cor[*,*,nn] = ibis_linearity_correction_andor(scan8542_full1_cor[*,*,nn], camera_id='2748', data_includes_bias=0, verbose=0,lin_ptr=lin_ptr) & $
   scan8542_full1_cor[*,*,nn] = filter_fringe_fft(scan8542_full1_cor[*,*,nn],ca8542=1,fft_power=imtest_pow, fft_mask=imtest_mask) & $
endfor


; ********************** Determine Static Errors in Spectral Scans *******************************
scan8542_full1_fitave = FLTARR(xsz,ysz)
scan8542_full1_fitres = FLTARR(xsz,ysz)

for nn=0,100 do begin & flatfit = sfit(scan8542_full1_cor[*,*,nn],2) & $
    scan8542_full1_fitres += scan8542_full1_cor[*,*,nn]/flatfit & scan8542_full1_fitave += flatfit & endfor
scan8542_full1_fitave /= 101
scan8542_full1_fitres /= 101

; ********************** Remove Static Errors in Spectral Scans *******************************
scan8542_full1_corfix = scan8542_full1_cor
for nn=0,100 do begin & scan8542_full1_corfix[*,*,nn] = scan8542_full1_cor[*,*,nn]/flat_8542_ser1613_fitres/scan8542_full1_fitave & endfor

; ********************** Determine Blueshift from Flats *******************************

flat_linepos_wavestep = 0.02
fitmask_all = bytarr(xsz,ysz)+1
fitmask_rad = radial_distances([1,xsz,ysz],[xsz/2,ysz/2]) LE 580

flat_linepos          = find_line_shifts(flat_8542_ser1613corfix, flat8542_seq1_wv, wavelength_scale_interp=flat_wave_interp,$
                            line_start = 3, line_end=17, interp=flat_linepos_wavestep,radial=fitmask_all,fit_points=11)


flat_optcent          = ibis_optical_center(flat_linepos(*,*,0)*flat_linepos_wavestep,$
                                mask=fitmask_rad,coeff=flat_blueshift_coeff,$
                                blueshift_fit_matrix=flat_bluefit)

flat_bluefit_zero = flat_bluefit - MIN(flat_bluefit(xsz/2.-100:xsz/2.+100,ysz/2.-100:ysz/2.+100))

END
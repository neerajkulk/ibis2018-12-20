; ---------------------------------------------------------------
;   enter observation/data specific details here
; ---------------------------------------------------------------

; find all the relevant grid files
your_data_directory = '/SMdata4/kreardon/IBIS/Dec2018/20Dec2018/ibis/'
nb_grid_file = your_data_directory + 'spectral/GridImages/20181220_163909/s000.GridImages.fits'
wl_grid_file = your_data_directory + 'whitelight/GridImages/20181220_163909/s000.GridImages.fits'                            

im_per_file = 760
; number of grid points in the destretch array 
;     49 x 49 is what results for a 20 x 20 kernel and a 1000 x 1000 pixel image
ngrdx = 49LL
ngrdy = 49LL

num_prefilt = 4
nbgrid_byfilt = FLTARR(2,ngrdx,ngrdy,num_prefilt)
prefilt_index = FLTARR(2, num_prefilt)
prefilt_index[*,0] = [0,209]
prefilt_index[*,1] = [210,439]
prefilt_index[*,2] = [440,649]
prefilt_index[*,3] = [650,759]

output_distortion_file = 'calibrations/destr.components.nb2wl.new.v2.txt'

; ---------------------------------------------------------------
; ---------------------------------------------------------------

; load each corresponding grid image from the narrowband and whitelight channels
; and then destretch the two grid images two find the subarray offsets

disp_nb_grid   = FLTARR(2,ngrdx,ngrdy,im_per_file)
rdisp_nb_grid  = FLTARR(2,ngrdx,ngrdy,im_per_file)

for nn=0,im_per_file-1 do begin & $ 
    PRINT,nn & $
    ; we use load_calibrate_image to make sure the grid data are treated the same way as the science data
    nbgrid = load_calibrate_image(nb_grid_file,nn+1,'nb',cal_dir='calibrations',cal_info=cal_param_nb,header=header_nb,rot='solar',verbose=0) & $
    wlgrid = load_calibrate_image(wl_grid_file,nn+1,'wl',cal_dir='calibrations',cal_info=cal_param_nb,header=header_nb,rot='solar',verbose=0) & $
    ; divide out a fitted surface (sfit) to remove any intensity gradients due to wavelength shifts that might disturb the destretching
    nbgrid_reg = reg((nbgrid/sfit(nbgrid,4)),wlgrid/sfit(wlgrid,4),bytarr(20,20),disp=disp_out,rdisp=rdisp_out) & $                               
    disp_nb_grid[*,*,*,nn] = disp_out & rdisp_nb_grid[*,*,*,nn] = rdisp_out & $
endfor

; now average up the destretch vectors grouped by prefilter
FOR filtnum = 0,num_prefilt - 1 DO BEGIN
    nbgrid_byfilt[*,*,*,filtnum] = MEDIAN((disp_nb_grid - rdisp_nb_grid)[*,*,*,prefilt_index[0,filtnum]:prefilt_index[1,filtnum]],dim=4)
ENDFOR

; for each prefilter, subtract off the mean value of the destretch vectirs
; this corresponds to any uncorrected bulk shift that is specific to a individual filter
nbgrid_byfilt_zeroed = nbgrid_byfilt
for filtnum=0,num_prefilt - 1 do begin & $
    for dir=0,1 do begin & $
        nbgrid_byfilt_zeroed[dir,*,*,filtnum] -= MEAN(nbgrid_byfilt[dir,8:40,8:40,filtnum]) & $
endfor & endfor

; now average up the bulk-shift-corrected destretch vectors to get a single distortion map for all the prefilters
; we assume (after checking) the optical distortions are roughly the same in all wavelengths
disp_nb2wl        = REFORM(MEAN(nbgrid_byfilt_zeroed[*,*,*,*],dim=4)) 
disp_nb2wl[1,*,0] = disp_nb2wl[1,*,1]
disp_nb2wl        = SMOOTH(disp_nb2wl,[1,3,3],/edge_trunc)
disp_nb2wl        = disp_nb2wl + rdisp_nb_grid[*,*,*,0]
rdisp_nb2wl       = REFORM(rdisp_nb_grid[*,*,*,0])

; write the map of destretch vectors (and reference positions) out to a text file, which may be the fastest way to 
; load in the values from load_calibration_info.pro
destr_nb2wl = reform([rdisp_nb2wl,disp_nb2wl],4,ngrdx*ngrdy)
openw,11,output_distortion_file
for nn=0,ngrdx*ngrdy-1 do begin & printf,11,STRING(destr_nb2wl[*,nn],FORMAT='(F8.4)') & endfor
close,11

END

FUNCTION reg_loop, data, reference, kernels, rdisp=rdisp, disp=disp, data_loop=data_loop

kernels_sz = SIZE(kernels, /ST)

num_kernels = kernels_sz.DIMENSIONS[0]
max_size    = 0
data_loop = data

FOR knum=0,num_kernels-1 DO BEGIN

    kernel_in = BYTARR(kernels[knum], kernels[knum])
    data_in = data_loop
    data_loop = reg(data_in, reference, kernel_in, rdisp=rdisp_out, disp=disp_out)

    return_val = EXECUTE('rdisp_out_' + STRING(knum,FORMAT='(I2.2)') + ' = rdisp_out')
    return_val = EXECUTE('disp_out_'  + STRING(knum,FORMAT='(I2.2)') + ' = disp_out')
    
    IF N_ELEMENTS(rdisp_out) GT max_size THEN target_kernel = knum

ENDFOR

;PRINT,target_kernel

IF knum GT 1 THEN BEGIN
    return_val = EXECUTE('rdisp_target = rdisp_out_' + STRING(target_kernel,FORMAT='(I2.2)') )
    
    disp_out_sum = rdisp_target * 0.0
    
    FOR knum=0,num_kernels-1 DO BEGIN
        return_val = EXECUTE('rdisp_sm = rdisp_out_' + STRING(knum,FORMAT='(I2.2)') )
        return_val = EXECUTE('disp_sm  =  disp_out_' + STRING(knum,FORMAT='(I2.2)') )

        IF knum EQ target_kernel THEN BEGIN
            disp_out_sum += (disp_sm - rdisp_sm)
        ENDIF ELSE BEGIN
            ;PRINT,knum
            disp_sm_int_x = min_curve_surf(REFORM(disp_sm(0,*,*)),REFORM(rdisp_sm(0,*,*)),REFORM(rdisp_sm(1,*,*)),$
                xpout=REFORM(rdisp_target(0,*,*)),ypout=REFORM(rdisp_target(1,*,*)))
            disp_sm_int_y = min_curve_surf(REFORM(disp_sm(1,*,*)),REFORM(rdisp_sm(0,*,*)),REFORM(rdisp_sm(1,*,*)),$
                xpout=REFORM(rdisp_target(0,*,*)),ypout=REFORM(rdisp_target(1,*,*)))
            
            disp_out_sum(0,*,*) += disp_sm_int_x - rdisp_target(0,*,*)
            disp_out_sum(1,*,*) += disp_sm_int_y - rdisp_target(1,*,*)

        ENDELSE            

    ENDFOR
    
    disp_out_sum += rdisp_target

    rdisp = rdisp_target
    disp  = disp_out_sum
    
ENDIF ELSE BEGIN
    rdisp = rdisp_out
    disp  = disp_out
ENDELSE

data_out = doreg(data, rdisp, disp)

RETURN, data_out

END

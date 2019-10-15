FUNCTION make_arr, fname
  num_frames = 76
  xpix = 1000
  ypix = 1000
  
  obs = MAKE_ARRAY(xpix, ypix, num_frames, /INTEGER, VALUE = 0) 
  
  FOR frame=1,num_frames DO BEGIN
     PRINT, 'reading frame number:', frame
     current_frame = readfits(fname, exten_no = frame) 
     FOR i = 0, xpix-1 DO BEGIN 
        FOR j = 0, ypix-1 DO BEGIN 
           obs[i,j,frame-1] = current_frame[i,j] 
        ENDFOR 
     ENDFOR 
  ENDFOR
  RETURN, obs
END



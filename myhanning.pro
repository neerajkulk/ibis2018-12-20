FUNCTION myhanning, nxx, nyy, per

wind = FLTARR(nxx,nyy)

xx = FLTARR(nxx)
xx(*) = 1.0
nx = ROUND(nxx*per/100.)
hax = HANNING(2*nx)
xx(0:nx-1) = hax(0:nx-1)
xx(nxx-nx:*) = hax(nx:*)

yy = FLTARR(nyy)
yy(*) = 1.0
ny = ROUND(nyy*per/100.)
hay = HANNING(2*ny)
yy(0:ny-1) = hay(0:ny-1)
yy(nyy-ny:*) = hay(ny:*)

wind = xx # yy

RETURN, wind

END

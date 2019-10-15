@initrh
@files.common
@geometry.common
@atmos.common
@opacity.common
@spectrum.common

   metalFile      = "metals.out"
   moleculeFile   = "molecules.out"
   opacFile       = 'opacity.out'
   backgroundFile = 'background.dat'

   result = readInput('input.out')
   result = readgeometry('geometry.out')
   result = readatmos('atmos.out')
   result = readspectrum('spectrum.out')
   result = openOpacity(opacFile)

   WAVENO = 0  & RAYNO = 0  ;; Assuming 500 nm is the first wavelength.

   readOpacity, WAVENO, RAYNO

   tau500 = fltarr(geometry.Nx, geometry.Ny, geometry.Nz)

   FOR m=0, geometry.Ny-1 DO BEGIN & FOR n=0, geometry.Nx-1 DO BEGIN & tau500[n, m, *] = gettau(geometry.z, chi_c[n, m, *]) & ENDFOR & ENDFOR

   
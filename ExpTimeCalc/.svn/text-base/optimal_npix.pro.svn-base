function optimal_npix, imag
  
  log_npix =        2.9906515     -0.21442814*imag
  npix = 10.^log_npix
  npix = floor(npix) + 1
  too_few = where(npix lt 2)
  if (too_few[0] ne -1) then npix[too_few] = 2

  return, npix

end

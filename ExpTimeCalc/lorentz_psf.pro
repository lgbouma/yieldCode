function lorentz_psf, x, y
  r2 = x^2+y^2
  r = sqrt(r2)
  fwhm = 1.2
  ;return, (fwhm/(2.0*!PI)) / ((fwhm/2.0)^2 + r2)
  return, (fwhm/(!DPI)^2)*((fwhm/2.0)^2 + r2)^(-2)
end

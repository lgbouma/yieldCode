function calc_shz, lum, teff
; lum is bolometric luminosity in solar units
; teff is effective temperature in Kelvins
;
; calculate inner and outer S values corresponding to the HZ
; using Rappapoopoo Erratum Table 3

  c_in = [1.0146, 8.1884D-5, 1.9394D-9, -4.3618D-12, -6.8260D-16] ; moist greenhouse
  c_out= [0.3507, 5.9578D-5, 1.6707D-9, -3.0058D-12, -5.1925D-16] ; max greenhouse

  shz1 = 1d0*poly(teff, c_in)
  shz2 = 1d0*poly(teff, c_out)

  ahz1 = sqrt(lum/shz1)
  ahz1 = sqrt(lum/shz2)

  return, [ahz1, ahz2]

end

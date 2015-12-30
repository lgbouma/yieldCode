function calc_s_hz, teff

; calculate inner and outer S values corresponding to the HZ
; using Rappapoopoo Erratum Table 3

  c_in = [1.0146, 8.1884D-5, 1.9394D-9, -4.3618D-12, -6.8260D-16] ; moist greenhouse
  c_out= [0.3507, 5.9578D-5, 1.6707D-9, -3.0058D-12, -5.1925D-16] ; max greenhouse

  s_hz_in = 1d0*poly(teff, c_in)
  s_hz_out = 1d0*poly(teff, c_out)

  return, [s_hz_in, s_hz_out]

end

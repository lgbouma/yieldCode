function number_of_stars_expz, d, n0, h
;
; returns total number of stars within a maximum distance d in
; assuming the number density of stars varies as
; n = n0 exp(-|z|/h)
;
; d should be in parsecs
; n0 should be in stars per cubic parsec
; h should be in parsecs
;

  e = d/h
  near = where(e lt 1.0D-3, complement=far)
; this is to avoid numerical instability. Should be good to within 0.05%

  n_tot = (4./3.)*!PI*d^3. * n0
  n_tot[far] = n_tot[far] * 3.0 * (0.5*(e[far])^2. + exp(-e[far])*(1.+e[far]) - 1.0)/(e[far])^3.

  return, n_tot

end

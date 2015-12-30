function jlf_model, mj, p=p

  if (keyword_set(p)) then begin
     p=p
  endif else begin
     p=[81.5356,4.97171,0.577618,3.06956,0.411137]
  endelse

  dmj = 0.*mj
  mj_peak = 7.5
  hi = where(mj ge mj_peak, complement=lo)

;  if (hi[0] ne -1) then dmj[hi] = (abs(mj[hi]-p[1])/p[2])^p[3]
;  if (lo[0] ne -1) then dmj[lo] = (abs(mj[lo]-p[1])/p[4])^p[5]
;  jlf_calc = p[0] * exp(-dmj)

  jlf_calc = 0.0*mj

  if (hi[0] ne -1) then begin
     dmj = (abs(mj[hi]-mj_peak)/p[1])
     jlf_calc[hi] = p[0] * (1.0 + dmj^2.)^(-p[3])
  endif

  if (lo[0] ne -1) then begin
     dmj = (abs(mj[lo]-mj_peak)/p[2])
     jlf_calc[lo] = p[0] * (1.0 + dmj^2.)^(-p[4])
  endif
  
  return, jlf_calc

end

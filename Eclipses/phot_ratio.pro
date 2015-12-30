function phot_ratio, teff1, teff2, tmag1, tmag2, ph_p

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  nstars = n_elements(teff1)

  ph_1_filt = dblarr(nfilt, nstars)
  ph_2_filt = dblarr(nfilt, nstars)
  
  recipteff1 = 4000./teff1
  recipteff2 = 4000./teff2
  for j=0, nfilt-1 do begin
    ph_1_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff1 + $
         ph_p[j,2]*recipteff1^2. + ph_p[j,3]*recipteff1^3.
    ph_2_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff2 + $
         ph_p[j,2]*recipteff2^2. + ph_p[j,3]*recipteff2^3.
  endfor
  
  ph_1_filt[where(ph_1_filt lt 0)] = 0.0
  ph_2_filt[where(ph_2_filt lt 0)] = 0.0
  
  ph_1 = total(ph_1_filt, 1)*10.^(-0.4*(tmag1-10.))
  ph_2 = total(ph_2_filt, 1)*10.^(-0.4*(tmag2-10.))
  
  phr1 = ph_1/(ph_1+ph_2)
  return, phr1
  ;phr2 = 1.0-phr1

END

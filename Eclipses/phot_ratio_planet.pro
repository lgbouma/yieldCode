function phot_ratio_planet, teff1, teff2, tmag1, dm, r2, ph_p, tband

  kb = 1.38e-16
  c = 2.998e10
  h = 6.626e-27

  ; 10 pc in Rsun
  d = 206265.*215.*10.

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  nstars = n_elements(teff1)

  ph_1_filt = dblarr(nfilt, nstars)
  
  recipteff1 = 4000./teff1
  for j=0, nfilt-1 do begin
    ph_1_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff1 + $
         ph_p[j,2]*recipteff1^2. + ph_p[j,3]*recipteff1^3.
  endfor
  
  ph_1_filt[where(ph_1_filt lt 0)] = 0.0
  
  ph_1 = total(ph_1_filt, 1)*10.^(-0.4*(tmag1-10.-dm))
 
;  readcol, 'tband.csv', lam, tband ; for lambda in angstroms
  lam = tband[*,0]
  tran = tband[*,1]
  lam_cm = lam*1.E-8
  dlam = lam_cm[1]-lam_cm[0]
  
  tones = fltarr(nstars)+1.
  lones = fltarr(n_elements(lam_cm))+1.
  bph = 2.*c/((tones#(lam_cm^4.))*(exp(h*c/(kb*(teff2#lam_cm)))-1.))*!dpi*((r2#lones)/d)^2.*dlam
  filt = tones#tran
  ph_2 = total(bph*filt, 2)
 
  phr2 = ph_2/(ph_1+ph_2)
  return, phr2
  ;phr2 = 1.0-phr1
END

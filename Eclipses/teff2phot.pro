PRO teff2phot, teff, tmag, ph_p, ph

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  nstars = n_elements(teff)

  ph_filt = dblarr(nfilt, nstars)
  
  recipteff = 4000./teff
  for j=0, nfilt-1 do begin
    ph_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff + $
         ph_p[j,2]*recipteff^2. + ph_p[j,3]*recipteff^3.
  endfor
  
  ph_filt[where(ph_filt lt 0)] = 0.0

  ph = total(ph_1_filt, 1)*10.^(-0.4*(tmag-10.))

END

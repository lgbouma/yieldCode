PRO stack_prf_eclip, imag, teff, ph_p, prf, ph_star, dx=dx, dy=dy, fov_ind=fov_ind, mask=mask, sind=sind
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_prf = size(prf)
  npix_tot = sz_prf[4]
  nstars = n_elements(imag)
  if (keyword_set(mask)) then mask=mask else $
	mask = intarr(npix_tot)+1 

  ph_10_filt = dblarr(nfilt, nstars)
  ph_star = dblarr(total(mask), nstars)
  sind = intarr(total(mask), nstars)
 
  minds = where(mask)
  
  recipteff = 4000./teff
  for j=0, nfilt-1 do begin
    ph_10_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff + $
         ph_p[j,2]*recipteff^2. + ph_p[j,3]*recipteff^3.
  endfor

  ph_10_filt[where(ph_10_filt lt 0.0)] = 0.0

  for i=0, nstars-1 do begin
    prf_this = reform(prf[dx[i],dy[i],fov_ind[i],*,0:(nfilt-1)])
    ph_star_full = prf_this#ph_10_filt[*,i]
    ph_star_sub = ph_star_full[minds] 
    sind[*,i] = reverse(sort(ph_star_sub))
   ; ph_star[*,i] = 10.0^(-0.4*(imag[i]-10.))*total(ph_star_sub[sind[*,i]], /cumulative)
    ph_star[*,i] = 10.0^(-0.4*(imag[i]-10.))*ph_star_sub; [sind[*,i]]
  end

END

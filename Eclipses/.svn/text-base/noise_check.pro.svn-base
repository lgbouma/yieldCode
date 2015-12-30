PRO noise_check, sys=sys 
  imag = findgen(181)/10.+3.
  print, imag
  frac_file = 'bigfrac24_105_f3p33.fits' ; prf file 
;  rad_file = 'bigrad24_105_f3p33.fits' ; radius file 
  ph_file = 'ph_T_filt.fits' ; photon fluxes for T=10 vs Teff
  cr_file = 'crnoise.fits' ; photon fluxes for T=10 vs Teff
  tic_file = 'tic_teff.fits'
  ; Open the fits files
  frac_fits = mrdfits(frac_file)
; rad_fits = mrdfits(rad_file)/60. ; put into pixels
  ph_fits = mrdfits(ph_file)
  cr_fits = fltarr(100,64)
  ;  cr_fits = mrdfits(cr_file)
  tic_fits = mrdfits(tic_file)

  saturation=150000.
  aspix = 21.1
  npix_max = 64
  npix_min = 3
  mask2d = intarr(16,16)
  mid = indgen(8)+4
  mask1 = mask2d
  mask2 = mask2d
  mask1[*,mid] = 1
  mask2[mid,*] = 1
  mask2d = mask1*mask2
  mask1d = reform(mask2d, 16*16)
 
  teff = 4500. + dblarr(n_elements(imag))
  dx = 4 + dblarr(n_elements(imag))
  dy = 6 + dblarr(n_elements(imag))
  fov_ind = 1 + dblarr(n_elements(imag))
  field_angle = 6. + dblarr(n_elements(imag))
  elat = 30. + dblarr(n_elements(imag))
  tmag = imag + interpol(tic_fits[*,1], tic_fits[*,0], teff)

  if (keyword_set(sys)) then sys=sys else sys=60.
  readnoise=10.
  noises = dblarr(n_elements(imag), npix_max)
  dilution = dblarr(n_elements(imag), npix_max)
  diln = dblarr(n_elements(imag))
  satn = dblarr(n_elements(imag))
  noise_star = dblarr(n_elements(imag))
  noise_ro  = dblarr(n_elements(imag))
  noise_sys = dblarr(n_elements(imag))
  noise_zodi = dblarr(n_elements(imag))
  exptime = dblarr(n_elements(imag)) + 3600. ; what's the noise per hour?
  
  ; ph_star is npix x nstar
  stack_prf_eclip, tmag, teff, ph_fits, frac_fits, star_ph, $
      dx=dx, dy=dy, fov_ind=fov_ind, mask=mask1d, sind=sind
  dil_ph = dblarr(total(mask1d), n_elements(imag))
  for jj=0,n_elements(imag)-1 do begin
      star_ph[*,jj] = total(star_ph[sind[*,jj],jj], /cumulative)
  end

  zodi_flux, elat, aspix, zodi_ph

  noises = dblarr(n_elements(imag), npix_max)
  dilution = dblarr(n_elements(imag), npix_max)
  ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  exptime = dblarr(n_elements(imag)) + 3600.
  for ii=0,(npix_max-1) do begin
;    thiscr = cr[*,ii]
    calc_noise_eclip, star_ph, dil_ph, exptime, $
                 readnoise, sys, noise, $
                 npix_aper=(ii+1), $
                 field_angle=field_angle, $
;                 cr_noise = thiscr[crind[obs]]/sqrt(60.0/ffi_len)*star[obsid].ffi, $
                 subexptime=2., $
                 geom_area = 69.1, $
                 aspix=aspix, $
                 zodi_ph=zodi_ph, $
                 dilution=dil, $
                 e_tot_sub=estar
                 ;noise_star=shot_noise, $
	         ;noise_sky=bknd_noise, $
 		 ;noise_ro=read_noise, $
		 ;noise_sys=sys_noise
    dilution[*,ii] = dil
    noises[*,ii] = dil*noise
    ;shot_noises[*,ii] = shot_noise*1d6
    ;print, median(shot_noise*1d6)
    if (ii eq 0) then satn = (estar gt SATURATION)
  end
  noises = noises[*,(npix_min-1):(npix_max-1)]
  dilution = dilution[*,(npix_min-1):(npix_max-1)]
  noise = min(noises, ind, dimension=2)
  
  npix = ind / n_elements(imag) + npix_min
 
  for ii=0,n_elements(imag)-1 do begin 
    calc_noise_eclip, star_ph[*,ii], dil_ph[*,ii], exptime[ii], $
                 readnoise, sys, n, $
                 npix_aper=npix[ii], $
                 field_angle=field_angle[ii], $
;                 cr_noise = thiscr[crind[obs]]/sqrt(60.0/ffi_len)*star[obsid].ffi, $
                 subexptime=2., $
                 geom_area = 69.1, $
                 aspix=aspix, $
                 zodi_ph=zodi_ph[ii], $
                 dilution=dil, $
                 e_tot_sub=estar, $
                 noise_star=sn, $
                 noise_dil=dn, $
	         noise_zodi=zn, $
 		 noise_ro=rn, $
		 noise_sys=syn
  
                 diln[ii] = dil
                 noise[ii] = n
                 noise_star[ii]=sn
	         noise_zodi[ii]=zn
 		 noise_ro[ii]=rn 
		 noise_sys[ii]=syn
  end
  fitsout = [[imag], [npix], [noise], [noise_star], [noise_zodi], [noise_ro], [noise_sys], [satn], [diln]]
  mwrfits, fitsout, 'noises_newphot.fits'
END



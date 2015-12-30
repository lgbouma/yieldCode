PRO noise_calc, sys=sys 
  imag = findgen(151)/10.+5.
  ;print, imag
  npix_min=1
  npix_max=49
  frac_fits = mrdfits('frac24_1p0.fits')
  bk_fits = mrdfits('fov24_1p0_bkgnd_param.fits')
  frac_3d = mean(frac_fits, dimension=1)
  frac_2d = mean(frac_3d, dimension=1)
  frac = frac_2d[1,*]
  elon=110.
  elat=30.
  if (keyword_set(sys)) then sys=sys else sys=60.
  noises = dblarr(n_elements(imag), npix_max-npix_min+1)
  dilution = dblarr(n_elements(imag), npix_max-npix_min+1)
  diln = dblarr(n_elements(imag))
  satn = dblarr(n_elements(imag))
  exptime = dblarr(n_elements(imag)) + 3600. ; what's the noise per hour?
  for ii=(npix_min-1),(npix_max-1) do begin
     calc_noise, imag, exptime, noise, $
                 npix_aper=(ii+1), $
                 frac_aper=frac[ii], $
                 fov_ind=1, $
                 teff=4500., $
                 e_pix_ro = 10.,$
                 subexptime=2., $
                 geom_area = 61.2, $
                 sys_lim = sys, $
                 pix_scale = 24.*3600./4096., $
                 elon=elon, $
                 elat=elat, $
                 dilution=dil, $
                 e_star_sub=estar
  dilution[*,ii] = dil
    noises[*,ii] = (1.0+dil)*noise
    if (ii eq 0) then begin
	satn = (estar gt 150000.)
        print, noise
    end
  end
  minnoise = min(noises, ind, dimension=2)
  npix = ind / n_elements(imag) + 1
  print, median(npix)
  calc_noise, imag, exptime, noise, $
                 npix_aper=npix, $
                 frac_aper=frac[npix-1], $
                 fov_ind=1, $
                 teff=4500, $
                 e_pix_ro = 10.,$
                 subexptime=2., $
                 geom_area = 61.2, $
                 sys_lim = sys, $
                 pix_scale = 24.*3600./4096., $
                 elon=elon, $
                 elat=elat, $
                 dilution=diln, $
                 e_star_sub=estar, $
                 noise_star=shot_noise, $
			  noise_sky=bknd_noise, $
 			  noise_ro=read_noise, $
			  noise_sys=sys_noise
 
  fitsout = [[imag], [npix], [noise], [shot_noise], [bknd_noise], [read_noise], [sys_noise], [satn], [diln]]
  mwrfits, fitsout, 'noises_newphot.fits'
  
END



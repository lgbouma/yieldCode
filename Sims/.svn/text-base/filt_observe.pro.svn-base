pro filt_observe, sstruct=sstruct, pstruct=pstruct, sfile=sfile, pfile=pfile, outfile=outfile, $
	fov=fov, geomarea=geomarea, readnoise=readnoise, tranmin=tranmin, thresh=thresh, $
	nodil=nodil, ffi_len=ffi_len, $
	sys_limit=sys_limit, keep_ntra=keep_ntra, duty_cycle=duty_cycle, $
        bk_file=bk_file, sp_file=sp_file, ph_file=ph_file, prf_file=prf_file
	

 REARTH_IN_RSUN = 0.0091705248
;;;;;; basic parameters here

  DWELL_TIME = 13.66d0 ;27.4d0             ; days per field
  DOWNLINK_TIME = 16.0d0/24.0d0 ; correct for downlink time
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tranmin)) then NTRA_OBS_MIN = tranmin else NTRA_OBS_MIN = 2
  if (keyword_set(geomarea)) then GEOM_AREA=geomarea else GEOM_AREA=61.2; cm^2
  if (keyword_set(fov)) then fov=fov else fov=24.0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  if (keyword_set(ffi_len)) then ffi_len=ffi_len else ffi_len=30.
  apo_blank = (DWELL_TIME-DOWNLINK_TIME)*(1.0-duty_cycle/100.0)
 ;SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0 
  SUB_EXP_TIME = 2.0
  SATURATION = 150000.
  npix_max = 49
  npix_min = 3
  ; photon counts
  if (keyword_set(ph_file)) then ph_fits=mrdfits(ph_file) else ph_fits=0
  ; big frac file
  if (keyword_set(prf_file)) then prf_file=prf_file else prf_file='../ExpTimeCalc/bigfrac24_105_4700k.fits'
  prf_fits = mrdfits(prf_file)
  ; background file
  if (keyword_set(bk_file)) then bk_fits = mrdfits(bk_file) else bk_fits=0
  ; spline fits to PSRR
  if (keyword_set(sp_file)) then begin
	sp_fits = mrdfits(sp_file) 
 	sp_fits[4,*] = sp_fits[4,*]-sp_fits[4,0]
 	sp_fits[5,*] = sp_fits[5,*]-sp_fits[5,0]
 	sp_fits[6,*] = sp_fits[6,*]-sp_fits[6,0]
 	sp_fits[7,*] = sp_fits[7,*]-sp_fits[7,0]
  endif
  ccd_pix = 4096.0
  PIX_SCALE = fov*3600./ccd_pix
  pix_scale_deg = fov/ccd_pix


  if keyword_set(sstruct) then star=sstruct  else restore, sfile
  if keyword_set(pstruct) then planet=pstruct else restore, pfile
  ;restore, infile
  nstars = n_elements(star)
  nplanets = n_elements(planet)
  print, 'Observing ', nplanets, ' planets around ', nstars,' stars.'
  
; for each transiting planet, how many transits were observed?
  tra = where(planet.tra gt 0)
  traid = planet[tra].hostid

   print, 'Calculating number of eclipses'
  if ~(keyword_set(keep_ntra)) then begin
    planet[tra].ntra_obs = $
      n_eclip_blank(planet[tra].p, DWELL_TIME, $
      2.0*double(star[traid].npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank)
  endif

  print, 'Diluting FFIs'
  tra_ps = where((planet.tra gt 0) and (star[planet.hostid].ffi lt 1))
  if (tra_ps[0] ne -1) then begin
    dur_min = planet[tra_ps].dur*24.0*60.0
    planet[tra_ps].dep_eff = planet[tra_ps].dep*dil_ffi(dur_min, 1.0, ffis=nps)
    planet[tra_ps].dur_eff = (nps*1.0)/(24.0*60.0)
  endif

  tra_ffi = where((planet.tra gt 0) and (star[planet.hostid].ffi gt 0))
  if (tra_ffi[0] ne -1) then begin
    dur_min = planet[tra_ffi].dur*24.0*60.0
    planet[tra_ffi].dep_eff = planet[tra_ffi].dep*dil_ffi(dur_min, float(ffi_len), ffis=ffis)
    planet[tra_ffi].dur_eff = (ffis*ffi_len)/(24.0*60.0)
  endif
; for each observed transiting planet, calculate snr
 
   obs = where(planet.tra and (planet.ntra_obs gt 0))
   obsid = planet[obs].hostid
   nobs = n_elements(obs)
   print, 'Observing ', nobs, ' transits'
;  obs = indgen(n_elements(star))

  fov_ind = intarr(n_elements(obs))
  fov_ind[where((star[obsid].coord.fov_r ge 0.104*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.365*CCD_PIX))] = 1
  fov_ind[where((star[obsid].coord.fov_r ge 0.365*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.592*CCD_PIX))] = 2
  fov_ind[where(star[obsid].coord.fov_r  ge 0.592*CCD_PIX)] = 3
  field_angle = star[obsid].coord.fov_r / CCD_PIX * FOV 
 
  print, 'Stacking and sorting PRFs'
  ; ph_star is npix x nstar
  dx = floor(10*randomu(seed, nobs))
  dy = floor(10*randomu(seed, nobs))
  stack_prf, star[obsid].mag.i, star[obsid].teff, ph_fits, prf_fits, ph_star, dx=dx, dy=dy, fov_ind=fov_ind

  bin_sys = (star[obsid].pri or star[obsid].sec)
  bin_sep = star[obsid].companion.sep
  bin_imag = star[star[obsid].companion.ind].mag.i
  bin_teff = star[star[obsid].companion.ind].teff
  bins = where(bin_sys)

  stack_prf, bin_imag, bin_teff, ph_fits, prf_fits, ph_bin, dx=dx, dy=dy, fov_ind=fov_ind
 
  pix_sep = bin_sep[bins]/pix_scale   ; from arcsec to pixels
  r = sp_fits[fov_ind[bins],*]       ; distance (in pixels)
  di = sp_fits[fov_ind[bins]+4,*]    ; imag attenuation
  dibin = interpol(di, r, pix_sep)    ; interpolate over spline fit
  ph_bin = 10.0^(-0.4*dibin) * ph_bin ; attenuate the count rate
  
  print, 'Calculating noise'
  noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  dilution = dblarr(n_elements(obs), npix_max-npix_min+1)
  shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  exptime = dblarr(n_elements(obs)) + 3600.
  for ii=0,(npix_max-npix_min) do begin
     calc_noise_filt, ph_star, exptime, noise, $
		 npix_aper=(ii+npix_min), $
                 field_angle=field_angle, $
  		 fov_ind=fov_ind, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
		 bk_p = bk_fits, $
                 elon=star[obsid].coord.elon, $
                 elat=star[obsid].coord.elat, $
		 dilution=dil, $
		 e_star_sub=estar, $
		 noise_star=shot_noise, $
		 bin_sys= bin_sys, $
 		 bin_ph = ph_bin, $
  		 bin_sep = bin_sep
    dilution[*,ii] = dil
    if (keyword_set(nodil)) then noises[*,ii] = noise $
	else noises[*,ii] = dil*noise
    shot_noises[*,ii] = shot_noise*1d6
    ;print, median(shot_noise*1d6)
    ;if (ii eq 0) then star[obsid].sat = (estar gt SATURATION)
  end
  minnoise = min(noises, ind, dimension=2)
  star[obsid].npix = ind / n_elements(obs) + npix_min
  star[obsid].snr = 1.0/minnoise
  star[obsid].dil = dilution[ind]
  ; Calculate SNR in phase-folded lightcurve
  et_folded = double(planet[obs].ntra_obs) * $
    planet[obs].dur_eff * 24.0 * 3600
  et_tran = planet[obs].dur_eff * 24.0 * 3600
  planet[obs].snr = planet[obs].dep_eff / (minnoise) * sqrt(et_folded/3600.)
  planet[obs].snrtran = planet[obs].dep_eff / (minnoise) * sqrt(et_tran/3600.)
  planet[obs].snrgress = planet[obs].snrtran * $
 		sqrt(2. * planet[obs].ntra_obs * $
		     REARTH_IN_RSUN * planet[obs].r / $
		    (6.0*star[obsid].r*(1.0+planet[obs].b^2.)))
;   decide if it is 'detected'.
  det = where((planet.tra gt 0) and $
  	      (planet.ntra_obs ge NTRA_OBS_MIN) and $
	      (planet.snr ge SNR_MIN))
  print, 'Detected ', n_elements(det), ' planets.'
  planet[det].det = 1
;  detected = where(star.planet_hz.tra gt 0 and star.planet_hz.ntra_obs ge NTRA_OBS_MIN and star.planet_hz.snr ge SNR_MIN)
;  star[detected].planet_hz.det = 1
  if keyword_set(sstruct) then sstruct=star 
  if keyword_set(pstruct) then pstruct=planet 
  if keyword_set(outfile) then save, filen=outfile, planet

end

pro mult_observe, sstruct=sstruct, pstruct=pstruct, sfile=sfile, pfile=pfile, outfile=outfile, $
	fov=fov, geomarea=geomarea, readnoise=readnoise, tranmin=tranmin, thresh=thresh, $
	frac_file=frac_file, nodil=nodil, red=red, al_bk=al_bk, al_phot=al_phot, al_aper=al_aper, $
	sys_limit=sys_limit, keep_ntra=keep_ntra, duty_cycle=duty_cycle, bk_file=bk_file

 REARTH_IN_RSUN = 0.0091705248
;;;;;; basic parameters here

  DWELL_TIME = 13.66d0 ;27.4d0             ; days per field
  DOWNLINK_TIME = 16.0d0/24.0d0 ; correct for downlink time
;  DWELL_TIME = 27.4d0/2.0                ; days per field
;  DWELL_TIME = DWELL_TIME - 8.0d0/24.0d0  ; correct for downlink time
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tranmin)) then NTRA_OBS_MIN = tranmin else NTRA_OBS_MIN = 2
  if (keyword_set(geomarea)) then GEOM_AREA=geomarea else GEOM_AREA=61.2; cm^2
  if (keyword_set(fov)) then fov=fov else fov=24.0
  if (keyword_set(red)) then red=red else red=0
  if (keyword_set(al_bk)) then al_bk=al_bk else al_bk=0
  if (keyword_set(al_phot)) then al_phot=al_phot else al_phot=0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  apo_blank = (DWELL_TIME-DOWNLINK_TIME)*(1.0-duty_cycle/100.0)
 ;SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0 
  SUB_EXP_TIME = 2.0
  ; PIX_SCALE = 20.0; arcsec per pixel side
  ;SURVEY_FILE = '../Survey/survey_fov23_seg13.sav'
  SATURATION = 150000.
  RMAX_HAB = 2.0 & SHZ_IN = 0.5 & SHZ_OUT = 1.5  ; radii, incident fluxes corresponding to inner/outer edge of HZ
  npix_max = 49
  npix_min = 3
  if (keyword_set(frac_file)) then frac_file=frac_file else frac_file='../ExpTimeCalc/frac24_1p0.fits'
  if (keyword_set(bk_file)) then bk_fits = mrdfits(bk_file) else bk_fits=0
  frac_fits = mrdfits(frac_file)
  ;fov = 24.0
  ;n_segs = 13
  ;n_cams = 4
  ccd_pix = 4096.0
;;;;;;; set up some plotting stuff
 
  PIX_SCALE = fov*3600./ccd_pix
  pix_scale_deg = fov/ccd_pix


;;;;;;;; let's get to work
  if keyword_set(sstruct) then star=sstruct  else restore, sfile
  if keyword_set(pstruct) then planet=pstruct else restore, pfile
  ;restore, infile
  nstars = n_elements(star)
  nplanets = n_elements(planet)
  print, 'Observing ', nplanets, ' planets around ', nstars,' stars.'
  
; for each transiting planet, how many transits were observed?
  tra = where(planet.tra gt 0)
  traid = planet[tra].hostid

  if ~(keyword_set(keep_ntra)) then begin
    planet[tra].ntra_obs = $
      n_eclip_blank(planet[tra].p, DWELL_TIME, $
      2.0*double(star[traid].npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank)
  endif

; for each observed transiting planet, calculate snr
 
   obs = where(planet.tra and (planet.ntra_obs gt 0))
   obsid = planet[obs].hostid

;  obs = indgen(n_elements(star))

  fov_ind = intarr(n_elements(obs))
  fov_ind[where((star[obsid].coord.fov_r ge 0.104*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.365*CCD_PIX))] = 1
  fov_ind[where((star[obsid].coord.fov_r ge 0.365*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.592*CCD_PIX))] = 2
  fov_ind[where(star[obsid].coord.fov_r  ge 0.592*CCD_PIX)] = 3
  field_angle = star[obsid].coord.fov_r / CCD_PIX * FOV 
  ;print, total(fov_ind eq 0)
  ;print, total(fov_ind eq 1)
  ;print, total(fov_ind eq 2)
  ;print, total(fov_ind eq 3)
 
  noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  dilution = dblarr(n_elements(obs), npix_max-npix_min+1)
  shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  ;exptime = double(star[obs].planet.ntra_obs) * star[obs].planet.dur * 24.0 * 3600.
  exptime = dblarr(n_elements(obs)) + 3600.
  for ii=0,(npix_max-npix_min) do begin
     ;print, 'npix = ', ii+1
     thisfrac = frac_fits[*,*,*,(ii+npix_min-1)]
     frac = thisfrac[star[obsid].dx, star[obsid].dy,fov_ind]
     calc_noise, star[obsid].mag.i, exptime, noise, $
		 npix_aper=(ii+npix_min), $
		 frac_aper=frac, $
                 field_angle=field_angle, $
  		 fov_ind=fov_ind, $
                 teff=star[obsid].teff, $
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
                 red=red, $
  		 al_bk=al_bk, $
  		 al_phot=al_phot, $
		 noise_star=shot_noise
    dilution[*,ii] = dil
    if (keyword_set(nodil)) then noises[*,ii] = noise $
	else noises[*,ii] = (1.0+dil)*noise
    shot_noises[*,ii] = shot_noise*1d6
    ;print, median(shot_noise*1d6)
    if (ii eq 0) then star[obsid].sat = (estar gt SATURATION)
  end
  minnoise = min(noises, ind, dimension=2)
  star[obsid].npix = ind / n_elements(obs) + npix_min
  star[obsid].snr = 1.0/minnoise
  star[obsid].dil = dilution[ind]
  
  ; Calculate SNR in phase-folded lightcurve
  exptime = double(planet[obs].ntra_obs) * $
    planet[obs].dur * 24.0 * 3600.
  if (keyword_set(al_aper)) then begin
	frac=1.0
	npix_use=9.
  endif  else begin 
	frac = frac_fits[star[obsid].dx, star[obsid].dy, fov_ind,star[obsid].npix-1]
	npix_use = star[obsid].npix
  endelse
  calc_noise, star[obsid].mag.i, exptime, noise, $
		 npix_aper = npix_use, $
		 frac_aper = frac, $
  		 field_angle = field_angle, $
		 fov_ind = fov_ind, $
                 teff=star[obsid].teff, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
		 bk_p = bk_fits, $
                 elon=star[obsid].coord.elon, $
                 elat=star[obsid].coord.elat, $
		 dilution=dil, $
                 red=red, $
   		 al_phot=al_phot, $
   		 al_bk=al_bk, $
		 e_star_sub=estar
  if (total(planet[obs].dep) gt 0) then begin
    print, 'Calculating SNR of pre-defined depth'
    if (keyword_set(nodil)) then planet[obs].snr = planet[obs].dep / noise $
    else planet[obs].snr = planet[obs].dep / ((1.0+dil)*noise)
  endif else begin
    if (keyword_set(nodil)) then planet[obs].dep = SNR_MIN * noise $
      else planet[obs].dep = SNR_MIN * noise * (1.0+dil)
    planet[obs].r = star[obsid].r * sqrt(planet[obs].dep) / REARTH_IN_RSUN
    planet[obs].snr = SNR_MIN
  endelse

  ; Calculate SNR per transit
  exptime =  planet[obs].dur * 24.0 * 3600.
  if (keyword_set(al_aper)) then begin
	frac=1.0
	npix_use=9.
  endif  else begin 
	frac = frac_fits[star[obsid].dx, star[obsid].dy, fov_ind,star[obsid].npix-1]
	npix_use = star[obsid].npix
  endelse
  calc_noise, star[obsid].mag.i, exptime, noise, $
		 npix_aper = npix_use, $
		 frac_aper = frac, $
 		 field_angle = field_angle, $
		 fov_ind=fov_ind, $
                 teff = star[obsid].teff, $
                 e_pix_ro = E_PIX_RO, $
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
		 bk_p = bk_fits, $
                 elon=star[obsid].coord.elon, $
                 elat=star[obsid].coord.elat, $
		 dilution=dil, $
                 red=red, $
 		 al_phot=al_phot, $
 		 al_bk=al_bk, $
		 e_star_sub=estar
     
  if (total(planet.dep) gt 0) then begin
    if (keyword_set(nodil)) then planet[obs].snrtran =  planet[obs].dep / noise $
    else planet[obs].snrtran =  planet[obs].dep / ((1.0+dil)*noise)
    planet[obs].snrgress = planet[obs].snrtran * $
 		sqrt(2. * planet[obs].ntra_obs * $
		     REARTH_IN_RSUN * planet[obs].r / $
		    (6.0*star[obsid].r*(1.0+planet[obs].b^2.)))
    endif else begin
      if (keyword_set(nodil)) then planet[obs].snrtran = planet[obs].dep / noise $
         else planet[obs].snrtran = planet[obs].dep / ((1.0+dil)*noise)
    endelse
  
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

pro rad_observe, struct=struct, infile=infile, outfile=outfile, filen=filen, $
	fov=fov, geomarea=geomarea, readnoise=readnoise, tranmin=tranmin, thresh=thresh, $
	frac_file=frac_file, nodil=nodil, red=red, al_bk=al_bk, al_phot=al_phot, sys_limit=sys_limit, $ 
	keep_ntra=keep_ntra, duty_cycle=duty_cycle

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
  if (keyword_set(duty_cyce)) then duty_cycle=duty_cycle else duty_cycle=100.0
  apo_blank = (DWELL_TIME-DOWNLINK_TIME)*duty_cycle/100.0
  ;SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0 
  SUB_EXP_TIME = 2.0
  ; PIX_SCALE = 20.0; arcsec per pixel side
  ;SURVEY_FILE = '../Survey/survey_fov23_seg13.sav'
  SATURATION = 150000.
  RMAX_HAB = 2.0 & SHZ_IN = 0.5 & SHZ_OUT = 1.5  ; radii, incident fluxes corresponding to inner/outer edge of HZ
  npix_max = 49
  npix_min = 1
  if (keyword_set(frac_file)) then frac_file=frac_file else frac_file='../ExpTimeCalc/frac24_1p0.fits'
  frac_fits = mrdfits(frac_file)
  ;fov = 24.0
  ;n_segs = 13
  ;n_cams = 4
  ccd_pix = 4096.0
;;;;;;; set up some plotting stuff
 
  PIX_SCALE = fov*3600./ccd_pix
  pix_scale_deg = fov/ccd_pix


;;;;;;;; let's get to work
  if keyword_set(struct) then star=struct  else restore, infile
  ;restore, infile
  nstars = n_elements(star)
  print, 'Observing ', nstars,' stars.'
  

; to plot a band representing the galactic plane
;
;  m = 100
;  glon = -180. + 360.0*dindgen(m)/double(m-1)
;  glat = 0.0*glon
;  euler, glon, glat, elon, elat, select=6
;  aitoff, elon, elat, x_gal, y_gal
;  oplot, x_gal, y_gal, psym=8, syms=0.5*syms_multiplier, color=fsc_color('Black') 

; for each transiting planet, how many transits were observed?

  n_per = n_elements(star[0].planet.p)
  print, 'periods: ', n_per
 
  if (n_per gt 1) then begin
    tra = where(total(star.planet.tra,1) gt 0)
  endif else begin
    tra = where(star.planet.tra gt 0)
  endelse

;  tra = indgen(n_elements(star))
  if ~(keyword_set(keep_ntra)) then begin
    for kk=0,n_per-1 do begin
       star[tra].planet[kk].ntra_obs = $
	  n_eclip_blank(star[tra].planet[kk].p, DWELL_TIME, $
	  2.0*double(star[tra].npointings), periblank=DOWNLINK_TIME)
    endfor
  endif
; for each observed transiting planet, calculate snr
  ;if (total(star.planet.dep) eq 0) then begin
 
  if (n_per gt 1) then begin
   obs = where((total(star.planet.tra,1) gt 0) and $
		(total(star.planet.ntra_obs,1) gt 0))
  endif else begin
    obs = where(star.planet.ntra_obs gt 0)
  endelse
;  obs = indgen(n_elements(star))
  fov_ind = intarr(n_elements(obs))
  fov_ind[where((star[obs].coord.fov_r ge 0.104*CCD_PIX) and $
		(star[obs].coord.fov_r lt 0.365*CCD_PIX))] = 1
  fov_ind[where((star[obs].coord.fov_r ge 0.365*CCD_PIX) and $
		(star[obs].coord.fov_r lt 0.592*CCD_PIX))] = 2
  fov_ind[where(star[obs].coord.fov_r  ge 0.592*CCD_PIX)] = 3
  field_angle = star[obs].coord.fov_r / CCD_PIX * FOV 
  ;print, total(fov_ind eq 0)
  ;print, total(fov_ind eq 1)
  ;print, total(fov_ind eq 2)
  ;print, total(fov_ind eq 3)
 
  ;endif else begin
  ;  obs = findgen(n_elements(star))
  ;  star.planet.ntra_obs[where(star.planet.ntra_obs lt 1)] = 1
  ;endelse
  noises = dblarr(n_elements(obs), npix_max)
  dilution = dblarr(n_elements(obs), npix_max)
  shot_noises = dblarr(n_elements(obs), npix_max)
  ;exptime = double(star[obs].planet.ntra_obs) * star[obs].planet.dur * 24.0 * 3600.
  exptime = dblarr(n_elements(obs)) + 3600.
  for ii=(npix_min-1),(npix_max-1) do begin
     ;print, 'npix = ', ii+1
     thisfrac = frac_fits[*,*,*,ii]
     frac = thisfrac[star[obs].dx, star[obs].dy,fov_ind]
     calc_noise, star[obs].mag.i, exptime, noise, $
		 npix_aper=(ii+1), $
		 frac_aper=frac, $
                 field_angle=field_angle, $
                 teff=star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat, $
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
    if (ii eq 0) then star[obs].sat = (estar gt SATURATION)
  end
  minnoise = min(noises, ind, dimension=2)
  star[obs].npix = ind / n_elements(obs) + 1
  star[obs].snr = 1.0/minnoise
  star[obs].dil = dilution[ind]
  ; Minimum-radius mode
;  if (n_per gt 1) then begin
;    star[obs].planet.dep = SNR_MIN * minnoise
;    star[obs].planet.r = star[obs].r * sqrt(star[obs].planet.dep) / REARTH_IN_RSUN
;    star[obs].planet.snr = SNR_MIN
;    obs = where(star.planet_hz.tra gt 0 and star.planet_hz.ntra_obs gt 0)
  for ii=0,n_per-1 do begin
  ; Calculate SNR in phase-folded lightcurve
    exptime = double(star[obs].planet[ii].ntra_obs) * $
       star[obs].planet[ii].dur * 24.0 * 3600.
    frac = frac_fits[star[obs].dx, star[obs].dy, fov_ind,star[obs].npix-1]
    calc_noise, star[obs].mag.i, exptime, noise, $
		 npix_aper = star[obs].npix, $
		 frac_aper = frac, $
  		 field_angle = field_angle, $
                 teff=star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat, $
		 dilution=dil, $
                 red=red, $
   		 al_phot=al_phot, $
   		 al_bk=al_bk, $
		 e_star_sub=estar
    if (total(star.planet[ii].dep) gt 0) then begin
      print, 'Calculating SNR of pre-defined depth'
      if (keyword_set(nodil)) then star[obs].planet[ii].snr = star[obs].planet[ii].dep / noise $
	else star[obs].planet[ii].snr = star[obs].planet[ii].dep / ((1.0+dil)*noise)
    endif else begin
      if (keyword_set(nodil)) then star[obs].planet[ii].dep = SNR_MIN * noise $
        else star[obs].planet[ii].dep = SNR_MIN * noise * (1.0+dil)
      star[obs].planet[ii].r = star[obs].r * sqrt(star[obs].planet[ii].dep) / REARTH_IN_RSUN
      star[obs].planet[ii].snr = SNR_MIN
    endelse

    ; Calculate SNR per transit
    exptime =  star[obs].planet[ii].dur * 24.0 * 3600.
    frac = frac_fits[star[obs].dx, star[obs].dy, fov_ind, star[obs].npix-1]
	calc_noise, star[obs].mag.i, exptime, noise, $
		 npix_aper = star[obs].npix, $
		 frac_aper = frac, $
 		 field_angle = field_angle, $
                 teff = star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat, $
		 dilution=dil, $
                 red=red, $
 		 al_phot=al_phot, $
 		 al_bk=al_bk, $
		 e_star_sub=estar
     
    if (total(star.planet[ii].dep) gt 0) then begin
      if (keyword_set(nodil)) then star[obs].planet[ii].snrtran =  star[obs].planet[ii].dep / noise $
         else star[obs].planet[ii].snrtran =  star[obs].planet[ii].dep / ((1.0+dil)*noise)
      star[obs].planet[ii].snrgress = star[obs].planet[ii].snrtran * $
		sqrt(2. * star[obs].planet.ntra_obs * $
		     REARTH_IN_RSUN * star[obs].planet[ii].r / $
		    (6.0*star[obs].r*(1.0+star[obs].planet.b^2.)))
    endif else begin
      if (keyword_set(nodil)) then star[obs].planet[ii].snrtran = star[obs].planet[ii].dep / noise $
         else star[obs].planet[ii].snrtran = star[obs].planet[ii].dep / ((1.0+dil)*noise)
    endelse
  
;   decide if it is 'detected'.

    detected = where((star.planet[ii].tra gt 0) and $
		(star.planet[ii].ntra_obs ge NTRA_OBS_MIN) and $
		(star.planet[ii].snr ge SNR_MIN))
    star[detected].planet[ii].det = 1
    print, 'Period ', ii, ' detected ', n_elements(detected), ' planets.'
  endfor
;  detected = where(star.planet_hz.tra gt 0 and star.planet_hz.ntra_obs ge NTRA_OBS_MIN and star.planet_hz.snr ge SNR_MIN)
;  star[detected].planet_hz.det = 1
  if keyword_set(struct) then struct=star 
  if keyword_set(outfile) then save, filen=outfile, star

; Report basic parameters


end

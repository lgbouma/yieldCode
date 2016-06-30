pro eclip_observe, eclipse, star, bk, deep, frac, ph_p, cr, var, $
	aspix=aspix, effarea=effarea, readnoise=readnoise, $
  tranmin=tranmin, thresh=thresh, $
	ps_len=ps_len, ffi_len=ffi_len, saturation=saturation, $
	sys_limit=sys_limit, duty_cycle=duty_cycle, $
  dwell_time=dwell_time, downlink=downlink, $
  subexptime=subexptime, extMission=extMission, randSeed=randSeed, $
  ps_only=ps_only, run=run
;+
; NAME: eclip_observe
; PURPOSE: take in an eclipStruct and find out which are detected.
; Dilution due to background stars and binary companions happens here. It gets called
; if stars on this healpix tile get eclipsing planets.
; IN:
;   1. eclipse: "eclip" object (eclipStruct) with transit parameters filled out
;   2. star: outStarStruct of targets
;   3. bk: background stars (to be uniformly sampled)
;   4. deep: deeper dilution
;   5. frac: PRF file (for us, as-built from deb woods). [10, 20, 4, 144, 9]:
;       4 field angles. 9 wavelengths. 12 Teffs (table 1). 144 total pixels in image.
;   6. ph_p: photon fluxes in TESS bandpass (vs t_eff) (array of n_pix * n_star)
;   7. cr: a 100*64 array of zeros, supposed to be cosmic rays (presumably was tried then dropped)
;   8. var: 100*4 array of numbers to account for stellar variability
;       sigma_v is computed for Teff<4500K, 4500-5000K, 5000-6000K, and T_eff>6000K. Each star
;       gets assigned the same varibility stat btwn runs
;   9. aspix: arcsec per pixel (21.1)
;   10. effarea: of CCDs (~69cm^2)
;   11. readnoise: 10 electrons per subexposure
;   12. tranmin: 1 (minimum # eclipses for a detection)
;   13. thresh: 5 (detection threshold in phase-folded light curve; can be refined in postProcessing)
;   14. ps_len: 2 (minutes)
;   15. ffi_len: 30 (minutes)
;   16. saturation: 150,000 (electrons)
;   17. sys_limit: 60 ppm/hr (presumably systematic noise floor)
;   18. duty_cycle: array number hp tiles long (2908) of 100 (time blanked at perigee in minutes)
;   19. dwell_time: set to orbit period (13.66days)
;   20. downlink: 16./24. (downlink time in days)
;   21. subexptime: seconds per subexposure (2 seconds)
; OUT:
;-
  REARTH_IN_RSUN = 0.0091705248
  if (keyword_set(ps_len)) then ps_len=ps_len else ps_len = 2.0
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tranmin)) then NTRA_OBS_MIN = tranmin else NTRA_OBS_MIN = 2
  if (keyword_set(effarea)) then effarea=effarea else effarea=61.2; cm^2
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  if (keyword_set(ffi_len)) then ffi_len=ffi_len else ffi_len=30.
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.
  if (keyword_set(saturation)) then saturation=saturation else saturation=150000.
  if (keyword_set(dwell_time)) then dwell_time=dwell_time else dwell_time=13.66d0 ; days per orbit
  if (keyword_set(downlink)) then downlink=downlink else downlink=16.0d/24.0d
  apo_blank = (dwell_time-downlink)*(1.0-duty_cycle/100.0)

  assert, min(cr eq make_array(100,64)) eq 1, 'Legacy cr holder should be zeros'
  sz_frac = size(frac)  
  npix_img = sqrt(sz_frac[4]) ;sz_frac[4] is total # of px in img. Sqrt is # px / img side = 12.
  npix_min = 3
  npix_max = (npix_img/2)^2 ; 36, cf. Fig 14 S15, corresponds to I_c = 6 and brighter.
  mask2d = intarr(npix_img,npix_img)
  mid = indgen(npix_img/2)+npix_img/4
  mask1 = mask2d
  mask2 = mask2d
  mask1[*,mid] = 1
  mask2[mid,*] = 1
  mask2d = mask1*mask2 ; 12*12 array of 1's and 0's, presumably masking for stack
  mask1d = reform(mask2d, npix_img*npix_img)

  cind = findgen(npix_img)
  ones = fltarr(npix_img)+1.
  xx = cind#ones
  yy = ones#cind
  
  nstars = n_elements(star)
  neclipses = n_elements(eclipse)
  print, 'Observing ', neclipses, ' eclipses around ', nstars,' stars.'
  dx = floor(10*randomu(randSeed, neclipses)) ; same tile gets same seed
  dy = floor(20*randomu(randSeed, neclipses))
  crind = dy*5 + dx

  vtf = [1D7, 6000., 5000., 4500., 1D0]
  
  ecid = eclipse.hostid
  print, 'Calculating number of eclipses'
  fnow = eclipse.f
  ecc = eclipse.ecc
  w = eclipse.w
  p = eclipse.p
  f1 = !dpi/2. - w
  f2 = -!dpi/2. - w
  ; Kepler problem in reverse
  e1 = 2.0*atan(sqrt(1.0-ecc)*sin(f1/2.0), sqrt(1.0+ecc)*cos(f1/2.0))
  e2 = 2.0*atan(sqrt(1.0-ecc)*sin(f2/2.0), sqrt(1.0+ecc)*cos(f2/2.0))
  enow = 2.0*atan(sqrt(1.0-ecc)*sin(fnow/2.0), sqrt(1.0+ecc)*cos(fnow/2.0))
  m1 = (e1-ecc*sin(e1))
  m2 = (e2-ecc*sin(e2))
  mnow = (enow-ecc*sin(enow))
  dr1 = (m1-mnow)/(2.*!dpi) + 4.0
  dr2 = (m2-mnow)/(2.*!dpi) + 4.0
  
  ; Assign variability to eclipses
  for ii=0, n_elements(vtf)-2 do begin
    thisteff = where(star[ecid].teff lt vtf[ii] and star[ecid].teff ge vtf[ii+1])
    if (thisteff[0] ne -1) then begin
      thisvar = var[*,ii]
      eclipse[thisteff].var = thisvar[crind[thisteff]]
    end
  end

  dayoff1 = (dr1 mod 1)*p 
  dayoff2 = (dr2 mod 1)*p
 
  ; Get number of *observed* eclipses
  if run eq 0 then begin
    eclipse.pri.neclip_obs1 = $
        n_eclip(eclipse.p, dwell_time, $
        double(eclipse.pri.npointings), dayoff1, periblank=downlink, apoblank=apo_blank)
    eclipse.pri.neclip_obs2 = $
        n_eclip(eclipse.p, DWELL_TIME, $
        double(eclipse.pri.npointings), dayoff2, periblank=downlink, apoblank=apo_blank)
  endif
  if run eq 1 then begin
    eclipse.ext.neclip_obs1 = $
        n_eclip(eclipse.p, dwell_time, $
        double(eclipse.ext.npointings), dayoff1, periblank=downlink, apoblank=apo_blank)
    eclipse.ext.neclip_obs2 = $
        n_eclip(eclipse.p, DWELL_TIME, $
        double(eclipse.ext.npointings), dayoff2, periblank=downlink, apoblank=apo_blank)
  endif

  ; Get eclipse durations and depths given whether they're obsd as PSs or pseudo-FFIs
  assert, extMission, 'Have not written otherwise'
  if extMission and ps_only then begin
    tra_ps = where(star[eclipse.hostid].ffi eq 0 or star[eclipse.hostid].nPntgs eq 1)
    tra_ffi = -1
  endif

  if extMission and ~ps_only then begin
    if run eq 0 then begin
      tra_ps = where(eclipse.ffiClass eq 8 or eclipse.ffiClass eq 5 or eclipse.ffiClass eq 2)
      tra_ffi = where(eclipse.ffiClass eq 3 or eclipse.ffiClass eq 6 or eclipse.ffiClass eq 9)
    endif
    if run eq 1 then begin
      tra_ps = where(eclipse.ffiClass eq 7 or eclipse.ffiClass eq 8 or eclipse.ffiClass eq 9)
      tra_ffi = where(eclipse.ffiClass eq 4 or eclipse.ffiClass eq 5 or eclipse.ffiClass eq 6)
    endif
  endif

  if (tra_ps[0] ne -1) then begin
    dur1_min = eclipse[tra_ps].dur1*24.0*60.0
    dur2_min = eclipse[tra_ps].dur2*24.0*60.0
    if run eq 0 then begin
      eclipse[tra_ps].pri.dep1_eff = eclipse[tra_ps].dep1*dil_ffi_eclip(dur1_min, float(ps_len), $
                                  ffis=nps, randSeed=randSeed)
      eclipse[tra_ps].pri.dur1_eff = (nps*ps_len)/(24.0*60.0)
      eclipse[tra_ps].pri.dep2_eff = eclipse[tra_ps].dep2*dil_ffi_eclip(dur2_min, float(ps_len), $
                                  ffis=nps, randSeed=randSeed+9001) ; diff seeds for each
      eclipse[tra_ps].pri.dur2_eff = (nps*ps_len)/(24.0*60.0)
    endif
    if run eq 1 then begin
      eclipse[tra_ps].ext.dep1_eff = eclipse[tra_ps].dep1*dil_ffi_eclip(dur1_min, float(ps_len), $
                                  ffis=nps, randSeed=randSeed)
      eclipse[tra_ps].ext.dur1_eff = (nps*ps_len)/(24.0*60.0)
      eclipse[tra_ps].ext.dep2_eff = eclipse[tra_ps].dep2*dil_ffi_eclip(dur2_min, float(ps_len), $
                                  ffis=nps, randSeed=randSeed+9001) ; diff seeds for each
      eclipse[tra_ps].ext.dur2_eff = (nps*ps_len)/(24.0*60.0)
    endif
  endif

  if (tra_ffi[0] ne -1) then begin
    dur1_min = eclipse[tra_ffi].dur1*24.0*60.0
    dur2_min = eclipse[tra_ffi].dur2*24.0*60.0
    if run eq 0 then begin
      eclipse[tra_ffi].pri.dep1_eff = eclipse[tra_ffi].dep1*dil_ffi_eclip(dur1_min, float(ffi_len), $
                                    ffis=ffis, randSeed=randSeed+9002)
      eclipse[tra_ffi].pri.dur1_eff = (ffis*ffi_len)/(24.0*60.0)
      eclipse[tra_ffi].pri.dep2_eff = eclipse[tra_ffi].dep2*dil_ffi_eclip(dur2_min, float(ffi_len), $
                                    ffis=ffis, randSeed=randSeed+9003)
      eclipse[tra_ffi].pri.dur2_eff = (ffis*ffi_len)/(24.0*60.0)
    endif
    if run eq 1 then begin
      eclipse[tra_ffi].ext.dep1_eff = eclipse[tra_ffi].dep1*dil_ffi_eclip(dur1_min, float(ffi_len), $
                                    ffis=ffis, randSeed=randSeed+9002)
      eclipse[tra_ffi].ext.dur1_eff = (ffis*ffi_len)/(24.0*60.0)
      eclipse[tra_ffi].ext.dep2_eff = eclipse[tra_ffi].dep2*dil_ffi_eclip(dur2_min, float(ffi_len), $
                                    ffis=ffis, randSeed=randSeed+9003)
      eclipse[tra_ffi].ext.dur2_eff = (ffis*ffi_len)/(24.0*60.0)
    endif
  endif

  ; For each observed transiting eclipse, calculate preliminary snr (pre-dilution)
  if run eq 0 then obs = where((eclipse.pri.neclip_obs1 + eclipse.pri.neclip_obs2) gt NTRA_OBS_MIN)
  if run eq 1 then obs = where((eclipse.ext.neclip_obs1 + eclipse.ext.neclip_obs2) gt NTRA_OBS_MIN)
  if (obs[0] ne -1) then begin
    obsid = eclipse[obs].hostid
    nobs = n_elements(obs)
    print, 'Observing ', nobs, ' transits'
    print, 'Stacking and sorting PRFs'
    ; ph_star is npix x nstar
    assert, max(eclipse[*].coord.fov_ind) eq 0, 'error: fov_ind zeroing'
    stack_prf_eclip, star[obsid].mag.t, star[obsid].teff, ph_p, frac, star_ph, $
        dx=dx[obs], dy=dy[obs], fov_ind=eclipse[obs].coord.fov_ind, mask=mask1d, sind=sind

    ; BEBs and HEBs
    bebdil = where(eclipse[obs].class eq 3 or eclipse[obs].class eq 4 or eclipse[obs].class eq 5)
    assert, bebdil eq -1
    beb_ph = dblarr(total(mask1d), nobs)
    for jj=0,nobs-1 do begin
      star_ph[*,jj] = total(star_ph[sind[*,jj],jj], /cumulative)
      beb_ph[*,jj] = total(beb_ph[sind[*,jj],jj], /cumulative)
    end      
    assert, max(beb_ph) eq 0, 'What are background eclipsing binaries doing in eclip_observe?'

    zodi_flux, eclipse[obs].coord.elat, aspix, zodi_ph
    eclipse[obs].zodi_ph = zodi_ph

    print, 'Calculating noise'
    noises = dblarr(n_elements(obs), npix_max)
    dilution = dblarr(n_elements(obs), npix_max)
    ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
    exptime = dblarr(n_elements(obs)) + 3600.
    for ii=0,(npix_max-1) do begin ; see p13 S+15. This gets noises for all apertures; takes min below
      thiscr = cr[*,ii]
      calc_noise_eclip, star_ph, beb_ph, exptime, $ ; computes sigma-1hr (then sqrt(t_int) below)
          readnoise, sys_limit, noise, $
          npix_aper=(ii+1), field_angle=eclipse[obs].coord.field_angle, $
          cr_noise = star[obsid].ffi*0., $ ; wonky call for array size compatibility
          subexptime=subexptime, $
          geom_area = effarea, $
          aspix=aspix, $
          zodi_ph=zodi_ph, $
          dilution=dil, $
          e_tot_sub=estar, $
          noise_star=shot_noise
      dilution[*,ii] = dil
      if (keyword_set(nodil)) then noises[*,ii] = noise $
      else noises[*,ii] = dil*noise
      if (ii eq 0) then eclipse[obs].sat = (estar gt SATURATION)
    end
    noises = noises[*,(npix_min-1):(npix_max-1)]
    dilution = dilution[*,(npix_min-1):(npix_max-1)]
    minnoise = min(noises, ind, dimension=2)
    eclipse[obs].npix = ind / n_elements(obs) + npix_min
    eclipse[obs].dil = dilution[ind]

    for ss=0,nobs-1 do eclipse[obs[ss]].star_ph = star_ph[eclipse[obs[ss]].npix-1,ss]
    ; Calculate SNR (pre-dilution) in "phase-folded" (in this idealized way) lightcurve
    if run eq 0 then begin
      et1_eclp = eclipse[obs].pri.dur1_eff * 24.0
      et2_eclp = eclipse[obs].pri.dur2_eff * 24.0
      eclipse[obs].pri.snreclp1 = eclipse[obs].pri.dep1_eff $
                                  / sqrt(minnoise^2./et1_eclp + eclipse[obs].var^2.)
      eclipse[obs].pri.snreclp2 = eclipse[obs].pri.dep2_eff $
                                  / sqrt(minnoise^2./et2_eclp + eclipse[obs].var^2.)
      eclipse[obs].pri.snr1 = sqrt(eclipse[obs].pri.snreclp1^2. * double(eclipse[obs].pri.neclip_obs1))
      eclipse[obs].pri.snr2 = sqrt(eclipse[obs].pri.snreclp2^2. * double(eclipse[obs].pri.neclip_obs2))
      eclipse[obs].pri.snr  = sqrt(eclipse[obs].pri.snr1^2. + eclipse[obs].pri.snr2^2.)
    endif
    if run eq 1 then begin
      et1_eclp = eclipse[obs].ext.dur1_eff * 24.0
      et2_eclp = eclipse[obs].ext.dur2_eff * 24.0
      eclipse[obs].ext.snreclp1 = eclipse[obs].ext.dep1_eff $
                                  / sqrt(minnoise^2./et1_eclp + eclipse[obs].var^2.)
      eclipse[obs].ext.snreclp2 = eclipse[obs].ext.dep2_eff $
                                  / sqrt(minnoise^2./et2_eclp + eclipse[obs].var^2.)
      eclipse[obs].ext.snr1 = sqrt(eclipse[obs].ext.snreclp1^2. * double(eclipse[obs].ext.neclip_obs1))
      eclipse[obs].ext.snr2 = sqrt(eclipse[obs].ext.snreclp2^2. * double(eclipse[obs].ext.neclip_obs2))
      eclipse[obs].ext.snr  = sqrt(eclipse[obs].ext.snr1^2. + eclipse[obs].ext.snr2^2.)
      eclipse[obs].snrf  = sqrt(eclipse[obs].pri.snr1^2. + eclipse[obs].pri.snr2^2. + $
                           eclipse[obs].ext.snr1^2. + eclipse[obs].ext.snr2^2.)
    endif

    ; decide if it is 'detected' before dilution.
    if run eq 0 then det = where(((eclipse.pri.neclip_obs1 + eclipse.pri.neclip_obs2) ge NTRA_OBS_MIN) $
                           and (eclipse.pri.snr ge SNR_MIN))
    if run eq 1 then det = where(((eclipse.ext.neclip_obs1 + eclipse.ext.neclip_obs2) ge NTRA_OBS_MIN) $
                           and (eclipse.snrf ge SNR_MIN))

    ;;; Dilute if there are planets detected (preliminarily) with SNR > 5 (nb dilution makes SNR decrease)
    ;;; This is "second half" towards most precise SNR we derive
    if (det[0] ne -1) then begin
      detid = eclipse[det].hostid
      ndet = n_elements(det)
      print, 'Re-observing ', ndet, ' transits'
      
      print, 'Stacking and sorting PRFs'
      ; ph_star is npix x nstar
      stack_prf_eclip, star[detid].mag.t, star[detid].teff, ph_p, frac, star_ph, $
          dx=dx[det], dy=dy[det], fov_ind=eclipse[det].coord.fov_ind, mask=mask1d, sind=sind
      ;bk_ph = eclipse[det].bk_ph

      zodi_flux, eclipse[det].coord.elat, aspix, zodi_ph
      eclipse[det].zodi_ph = zodi_ph

      dil_ph = dblarr(total(mask1d), ndet)
      beb_ph = dblarr(total(mask1d), ndet)
      ;tgt_ph = dblarr(total(mask1d), ndet)
      
      ; print, "Diluting with binary companions"
      ; Binary companions dilute planet detections. Not EBs or HEBs
      bindil = where(eclipse[det].class eq 1)
      if (bindil[0] ne -1) then begin
        nbindil = n_elements(bindil)
        dilute_binary, eclipse[det[bindil]], star, frac, ph_p, $
            dx[det[bindil]], dy[det[bindil]], dilvec, aspix=aspix, radmax=6.0, randSeed=randSeed
        for jj=0,nbindil-1 do dil_ph[*,bindil[jj]] = dil_ph[*,bindil[jj]] + dilvec[*,jj]
      end
    
      ; Everything is diluted by deep stars
      print, "Diluting with deep stars"
      dilute_eclipse_img, eclipse[det], deep, frac, ph_p, $
          dx[det], dy[det], dilvec, aspix=aspix, sq_deg=0.0134, radmax=2.0, randSeed=randSeed+9004
      for jj=0,ndet-1 do dil_ph[*,jj] = dil_ph[*,jj] + dilvec[*,jj]
      
      ; BEB and BTP hosts will be diluted by the background star. Don't add others
      print, "Diluting with background stars"
      bkdil = where(eclipse[det].class ne 3 and eclipse[det].class ne 5) 
      if (bkdil[0] ne -1) then begin
        nbkdil = n_elements(bkdil)
        dilute_eclipse_img, eclipse[det[bkdil]], bk, frac, ph_p, $
            dx[det[bkdil]], dy[det[bkdil]], dilvec, aspix=aspix, sq_deg=0.134, $
            radmax=4.0, randSeed=randSeed+9005
        for jj=0,nbkdil-1 do dil_ph[*,bkdil[jj]] = dil_ph[*,bkdil[jj]] + dilvec[*,jj]
      end
   
      ; Everything
      print, "Diluting with other target stars"
      dilute_eclipse_img, eclipse[det], star, frac, ph_p, $
          dx[det], dy[det], dilvec, aspix=aspix, sq_deg=13.4, radmax=6.0, randSeed=randSeed+9006
      for jj=0,ndet-1 do dil_ph[*,jj] += dilvec[*,jj] ; 160423: uh.. not +=?
   
      ; Sort into the same pixel order as target star flux
      for jj=0,ndet-1 do begin
        star_ph[*,jj] = total(star_ph[sind[*,jj],jj], /cumulative)
        dil_ph[*,jj] = dil_ph[*,jj] + beb_ph[*,jj]; + tgt_ph[*,jj]
        dil_ph[*,jj] = total(dil_ph[sind[*,jj],jj], /cumulative)
        beb_ph[*,jj] = total(beb_ph[sind[*,jj],jj], /cumulative)
      end      

      noises = dblarr(n_elements(det), npix_max)
      dilution = dblarr(n_elements(det), npix_max)
      exptime = dblarr(n_elements(det)) + 3600.
      for ii=0,(npix_max-1) do begin
        thiscr = cr[*,ii]
        calc_noise_eclip, star_ph, dil_ph, exptime, $
            readnoise, sys_limit, noise, $
            npix_aper=(ii+1), $
            field_angle=eclipse[det].coord.field_angle, $
            cr_noise=star[detid].ffi*0., $
            subexptime=subexptime, $
            geom_area=effarea, $
            aspix=aspix, $
            zodi_ph=zodi_ph, $
            dilution=dil, $
            e_tot_sub=estar, $
            noise_star=shot_noise
        dilution[*,ii] = dil
        if (keyword_set(nodil)) then noises[*,ii] = noise $
        else noises[*,ii] = dil*noise
        if (ii eq 0) then eclipse[det].sat = (estar gt SATURATION)
      end
      noises = noises[*,(npix_min-1):(npix_max-1)]
      dilution = dilution[*,(npix_min-1):(npix_max-1)]
      minnoise = min(noises, ind, dimension=2)
      eclipse[det].npix = ind / n_elements(det) + npix_min
      eclipse[det].dil = dilution[ind]
      if run eq 0 then eclipse[det].pri.snrhr = 1. / minnoise
      if run eq 1 then eclipse[det].ext.snrhr = 1. / minnoise
      for ss=0,ndet-1 do eclipse[det[ss]].star_ph = star_ph[eclipse[det[ss]].npix-1,ss]
      
      ; Calculate SNR in phase-folded lightcurve
      if run eq 0 then begin
        et1_eclp = eclipse[det].pri.dur1_eff * 24.0
        et2_eclp = eclipse[det].pri.dur2_eff * 24.0
        eclipse[det].pri.snreclp1 = eclipse[det].pri.dep1_eff $
                                    / sqrt(minnoise^2./et1_eclp + eclipse[det].var^2.)
        eclipse[det].pri.snreclp2 = eclipse[det].pri.dep2_eff $
                                    / sqrt(minnoise^2./et2_eclp + eclipse[det].var^2.)
        eclipse[det].pri.snr1 = sqrt(eclipse[det].pri.snreclp1^2 * double(eclipse[det].pri.neclip_obs1))
        eclipse[det].pri.snr2 = sqrt(eclipse[det].pri.snreclp2^2 * double(eclipse[det].pri.neclip_obs2))
        eclipse[det].pri.snr  = sqrt(eclipse[det].pri.snr1^2. + eclipse[det].pri.snr2^2.)
      endif
      if run eq 1 then begin
        et1_eclp = eclipse[det].ext.dur1_eff * 24.0
        et2_eclp = eclipse[det].ext.dur2_eff * 24.0
        eclipse[det].ext.snreclp1 = eclipse[det].ext.dep1_eff $
                                    / sqrt(minnoise^2./et1_eclp + eclipse[det].var^2.)
        eclipse[det].ext.snreclp2 = eclipse[det].ext.dep2_eff $
                                    / sqrt(minnoise^2./et2_eclp + eclipse[det].var^2.)
        eclipse[det].ext.snr1 = sqrt(eclipse[det].ext.snreclp1^2 * double(eclipse[det].ext.neclip_obs1))
        eclipse[det].ext.snr2 = sqrt(eclipse[det].ext.snreclp2^2 * double(eclipse[det].ext.neclip_obs2))
        eclipse[det].ext.snr  = sqrt(eclipse[det].ext.snr1^2. + eclipse[det].ext.snr2^2.)
        eclipse.snrf  = sqrt(eclipse.pri.snr1^2. + eclipse.pri.snr2^2. + $
                             eclipse.ext.snr1^2. + eclipse.ext.snr2^2.)
      endif

      if run eq 0 then begin
        det1 = where((eclipse.pri.neclip_obs1 ge NTRA_OBS_MIN) and $
             (eclipse.pri.snr1 ge SNR_MIN))
        det2 = where((eclipse.pri.neclip_obs2 ge NTRA_OBS_MIN) and $
             (eclipse.pri.snr2 ge SNR_MIN))
        det = where(((eclipse.pri.neclip_obs1 + eclipse.pri.neclip_obs2) ge NTRA_OBS_MIN) and $
            (eclipse.pri.snr ge SNR_MIN))
        if (det1[0] ne -1) then eclipse[det1].pri.det1 = 1
        if (det2[0] ne -1) then eclipse[det2].pri.det2 = 1
        if (det[0] ne -1) then begin
          eclipse[det].pri.det = 1
          print, 'Detected ', n_elements(det), ' eclipses in primary.'
        endif
      endif
      if run eq 1 then begin
        det1 = where((eclipse.ext.neclip_obs1 ge NTRA_OBS_MIN) and $
             (eclipse.ext.snr1 ge SNR_MIN))
        det2 = where((eclipse.ext.neclip_obs2 ge NTRA_OBS_MIN) and $
             (eclipse.ext.snr2 ge SNR_MIN))
        det = where(((eclipse.ext.neclip_obs1 + eclipse.ext.neclip_obs2) ge NTRA_OBS_MIN) and $
            (eclipse.ext.snr ge SNR_MIN))
        detf = where(eclipse.snrf ge SNR_MIN)
        if (det1[0] ne -1) then eclipse[det1].ext.det1 = 1
        if (det2[0] ne -1) then eclipse[det2].ext.det2 = 1
        if (det[0] ne -1) then begin
          eclipse[det].ext.det = 1
          print, 'Detected ', n_elements(det), ' eclipses in extended.'
        endif
        if (detf[0] ne -1) then begin
          eclipse[detf].detf = 1
          print, 'Detected ', n_elements(detf), ' eclipses over combined missions.'
        endif
      endif ; run (extended mission) if
    endif ;det if
  endif ; obs if
  assert, max(cr) eq 0, 'error: nonzero cosmic ray flux. means blowback to stack_prf'
end

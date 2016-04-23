pro eclip_observe, eclipse, star, bk, deep, frac, ph_p, cr, var, $
	aspix=aspix, effarea=effarea, readnoise=readnoise, $
  tranmin=tranmin, thresh=thresh, $
	ps_len=ps_len, ffi_len=ffi_len, saturation=saturation, $
	sys_limit=sys_limit, duty_cycle=duty_cycle, $
  dwell_time=dwell_time, downlink=downlink, $
  subexptime=subexptime, extMission=extMission, randSeed=randSeed, $
  ps_only=ps_only
;+
; NAME: eclip_observe
; PURPOSE: take in an eclipStruct and find out which are detected. I think 
; there is also "dilution" due to binaries going on in here. It gets called 
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
;   13. thresh: 5 (detection threshold in phase-folded light curve)
;   14. ps_len: 2 (minutes)
;   15. ffi_len: 30 (minutes)
;   16. saturation: 150,000 (electrons)
;   17. sys_limit: 60 ppm/hr (presumably systematic noise floor)
;   18. duty_cycle: array number hp tiles long (2908) of 100 (time blanked at perigee in minutes)
;   19. dwell_time: set to orbit period (13.66days)
;   20. downlink: 16./24. (downlink time in days)
;   21. seconds per subexposure (2 seconds)
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
  eclipse.neclip_obs1 = $
      n_eclip(eclipse.p, dwell_time, $
      2.0*double(eclipse.npointings), dayoff1, periblank=downlink, apoblank=apo_blank)
  eclipse.neclip_obs2 = $
      n_eclip(eclipse.p, DWELL_TIME, $
      2.0*double(eclipse.npointings), dayoff2, periblank=downlink, apoblank=apo_blank)

  print, 'Diluting FFIs'
  if ~extMission then tra_ps = where(star[eclipse.hostid].ffi lt 1)
  if extMission then tra_ps = where(star[eclipse.hostid].ffi eq 0 or $
                                    star[eclipse.hostid].nPntgs eq 1)
  if (tra_ps[0] ne -1) then begin
    dur1_min = eclipse[tra_ps].dur1*24.0*60.0
    dur2_min = eclipse[tra_ps].dur2*24.0*60.0
    eclipse[tra_ps].dep1_eff = eclipse[tra_ps].dep1*dil_ffi_eclip(dur1_min, float(ps_len), $
                                ffis=nps, randSeed=randSeed)
    eclipse[tra_ps].dur1_eff = (nps*ps_len)/(24.0*60.0)
    eclipse[tra_ps].dep2_eff = eclipse[tra_ps].dep2*dil_ffi_eclip(dur2_min, float(ps_len), $
                                ffis=nps, randSeed=randSeed+9001) ; diff seeds for each
    eclipse[tra_ps].dur2_eff = (nps*ps_len)/(24.0*60.0)
  endif

  if ~extMission then tra_ffi = where(star[eclipse.hostid].ffi gt 0)
  if extMission and ~ps_only then tra_ffi = where(star[eclipse.hostid].ffi eq 1 or $
                                     star[eclipse.hostid].nPntgs eq 0)
  if extMission and ps_only then tra_ffi = -1

  if (tra_ffi[0] ne -1) then begin
    assert, 0, 'warning: should not be called until extended mission is adapted to ffis'
    dur1_min = eclipse[tra_ffi].dur1*24.0*60.0
    dur2_min = eclipse[tra_ffi].dur2*24.0*60.0
    eclipse[tra_ffi].dep1_eff = eclipse[tra_ffi].dep1*dil_ffi_eclip(dur1_min, float(ffi_len), $
                                  ffis=ffis, randSeed=randSeed+9002)
    eclipse[tra_ffi].dur1_eff = (ffis*ffi_len)/(24.0*60.0)
    eclipse[tra_ffi].dep2_eff = eclipse[tra_ffi].dep2*dil_ffi_eclip(dur2_min, float(ffi_len), $
                                  ffis=ffis, randSeed=randSeed+9003)
    eclipse[tra_ffi].dur2_eff = (ffis*ffi_len)/(24.0*60.0)
  endif

  ; for each observed transiting eclipse, calculate preliminary snr
  obs = where((eclipse.neclip_obs1 + eclipse.neclip_obs2) gt NTRA_OBS_MIN)
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
    beb_ph = dblarr(total(mask1d), nobs)
    if (bebdil[0] ne -1) then begin
      nbeb = n_elements(bebdil)
      ;beb_ph = dblarr(total(mask1d), nbeb)
      dilute_beb, eclipse[obs[bebdil]], frac, ph_p, $
          dx[obs[bebdil]], dy[obs[bebdil]], bebvec, aspix=aspix, radmax=6.0, randSeed=randSeed+9008
      for jj=0,nbeb-1 do beb_ph[*,bebdil[jj]] = bebvec[*,jj]
    end   
    for jj=0,nobs-1 do begin
      star_ph[*,jj] = total(star_ph[sind[*,jj],jj], /cumulative)
      beb_ph[*,jj] = total(beb_ph[sind[*,jj],jj], /cumulative)
    end      
    ;dil_ph = dblarr(total(mask1d), nobs)

    zodi_flux, eclipse[obs].coord.elat, aspix, zodi_ph
    eclipse[obs].zodi_ph = zodi_ph

    print, 'Calculating noise'
    noises = dblarr(n_elements(obs), npix_max)
    dilution = dblarr(n_elements(obs), npix_max)
    ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
    exptime = dblarr(n_elements(obs)) + 3600.
    for ii=0,(npix_max-1) do begin ; see p13 S+15. This gets noises for all apertures; takes min below
      thiscr = cr[*,ii]
      calc_noise_eclip, star_ph, beb_ph, exptime, $
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
      ;shot_noises[*,ii] = shot_noise*1d6
      ;print, median(shot_noise*1d6)
      if (ii eq 0) then eclipse[obs].sat = (estar gt SATURATION)
    end
    noises = noises[*,(npix_min-1):(npix_max-1)]
    dilution = dilution[*,(npix_min-1):(npix_max-1)]
    minnoise = min(noises, ind, dimension=2)
    eclipse[obs].npix = ind / n_elements(obs) + npix_min
    eclipse[obs].dil = dilution[ind]
    ;stop
    for ss=0,nobs-1 do eclipse[obs[ss]].star_ph = star_ph[eclipse[obs[ss]].npix-1,ss]
    if (bebdil[0] ne -1) then begin
      beb_pixph = dblarr(nbeb)
      for tt=0,nbeb-1 do beb_pixph[tt] = beb_ph[eclipse[obs[bebdil[tt]]].npix-1,bebdil[tt]]
      beb_starph = eclipse[obs[bebdil]].star_ph
      minnoise[bebdil] = minnoise[bebdil]*beb_starph/beb_pixph
    end  
    ; Calculate SNR in phase-folded lightcurve
    et1_eclp = eclipse[obs].dur1_eff * 24.0
    et2_eclp = eclipse[obs].dur2_eff * 24.0
    eclipse[obs].snreclp1 = eclipse[obs].dep1_eff / sqrt(minnoise^2./et1_eclp + eclipse[obs].var^2.)
    eclipse[obs].snreclp2 = eclipse[obs].dep2_eff / sqrt(minnoise^2./et2_eclp + eclipse[obs].var^2.)
    eclipse[obs].snr1 = eclipse[obs].snreclp1 * sqrt(double(eclipse[obs].neclip_obs1))
    eclipse[obs].snr2 = eclipse[obs].snreclp2 * sqrt(double(eclipse[obs].neclip_obs2))
    eclipse[obs].snr  = sqrt(eclipse[obs].snr1^2. + eclipse[obs].snr2^2.)

    ; decide if it is 'detected'.
    det = WHERE(((eclipse.neclip_obs1 + eclipse.neclip_obs2) ge NTRA_OBS_MIN) and $
          (eclipse.snr ge SNR_MIN))

    ; Dilute if there are planets detected (preliminarily) with SNR > 5 (nb dilution makes SNR decrease)
    if (det[0] ne -1) then begin
      ;print, 'Calculating noise again'
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
      ; Binaries dilute planets, BEB targets. Not EBs or HEBs
      bindil = where(eclipse[det].class eq 1 or eclipse[det].class eq 3 or eclipse[det].class eq 5)
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
   
      ; BEBs, HEBs, and BTPs
      bebdil = where(eclipse[det].class eq 3 or eclipse[det].class eq 4 or eclipse[det].class eq 5)
      if (bebdil[0] ne -1) then begin
        nbeb = n_elements(bebdil)
        ;beb_ph = dblarr(total(mask1d), nbeb)
        dilute_beb, eclipse[det[bebdil]], frac, ph_p, $
            dx[det[bebdil]], dy[det[bebdil]], bebvec, aspix=aspix, radmax=6.0, randSeed=randSeed+9007
        for jj=0,nbeb-1 do beb_ph[*,bebdil[jj]] = bebvec[*,jj]
      end   
      
      ; Everything
      print, "Diluting with other target stars"
      dilute_eclipse_img, eclipse[det], star, frac, ph_p, $
          dx[det], dy[det], dilvec, aspix=aspix, sq_deg=13.4, radmax=6.0, randSeed=randSeed+9006
      for jj=0,ndet-1 do dil_ph[*,jj] += dilvec[*,jj] ; 160423: uh.. not +=?
   
      ; Calculate centroid shift and uncertainty (section 8 of S15)
      dep1 = eclipse[det].dep1_eff
      dep2 = eclipse[det].dep2_eff
      dur1 = eclipse[det].dur1_eff * eclipse[det].neclip_obs1 
      dur2 = eclipse[det].dur2_eff * eclipse[det].neclip_obs2 
 
      cennoises   = fltarr(ndet, npix_max-npix_min+1)
      xcennoise1s = fltarr(ndet, npix_max-npix_min+1)
      xcennoise2s = fltarr(ndet, npix_max-npix_min+1)
      ycennoise1s = fltarr(ndet, npix_max-npix_min+1)
      ycennoise2s = fltarr(ndet, npix_max-npix_min+1)
      xcenshift1s = fltarr(ndet, npix_max-npix_min+1)
      xcenshift2s = fltarr(ndet, npix_max-npix_min+1)
      ycenshift1s = fltarr(ndet, npix_max-npix_min+1)
      ycenshift2s = fltarr(ndet, npix_max-npix_min+1)
      for ii=0,(npix_max-npix_min) do begin
        thiscr = cr[*,ii]
        calc_noise_cen, star_ph, dil_ph, beb_ph, bebdil, $
            dur1, dur2, $
	          dep1, dep2, $
            xx, yy, sind, $
            readnoise, sys_limit, $
            xcen, ycen, $
            xcenshift1, xcenshift2, $
            ycenshift1, ycenshift2, $
            xcennoise1, xcennoise2, $
            ycennoise1, ycennoise2, $
            npix_aper=(ii+npix_min), $
            field_angle=eclipse[det].coord.field_angle, $
            cr_noise = star[detid].ffi*0., $
            subexptime=subexptime, $
            geom_area = effarea, $
            aspix=aspix, $
            zodi_ph=zodi_ph
        cennoises[*,ii] = sqrt(xcennoise1^2.+ycennoise1^2.)
        xcennoise1s[*,ii] = xcennoise1
        xcennoise2s[*,ii] = xcennoise2
        ycennoise1s[*,ii] = ycennoise1
        ycennoise2s[*,ii] = ycennoise2
        xcenshift1s[*,ii] = xcenshift1
        xcenshift2s[*,ii] = xcenshift2
        ycenshift1s[*,ii] = ycenshift1
        ycenshift2s[*,ii] = ycenshift2
      endfor
      mincennoise = min(cennoises, ind, dimension=2)
      ;ind = indgen(ndet)
      mincennoise = cennoises[ind]
      cenind = ind / ndet
      cenpix = cenind + npix_min
      censhift1 = sqrt(xcenshift1s[ind]^2. + ycenshift1s[ind]^2.)
      censhift2 = sqrt(xcenshift2s[ind]^2. + ycenshift2s[ind]^2.)
      eclipse[det].censhift1 = censhift1
      eclipse[det].censhift2 = censhift2
      eclipse[det].cenerr1 = sqrt((xcennoise1s[ind]*xcenshift1s[ind]/censhift1)^2. + $
          (ycennoise1s[ind]*ycenshift1s[ind]/censhift1)^2.)
      ; eclipse[det].cenerr1 = xcennoise1s[ind]*ycennoise1s[ind]/den1
      eclipse[det].cenerr2 = sqrt((xcennoise2s[ind]*xcenshift2s[ind]/censhift2)^2. + $
          (ycennoise2s[ind]*ycenshift2s[ind]/censhift2)^2.)
      ;eclipse[det].cenerr2 = xcennoise2s[ind]*ycennoise2s[ind]/den2
      
      ; Sort into the same pixel order as target star flux
      for jj=0,ndet-1 do begin
        star_ph[*,jj] = total(star_ph[sind[*,jj],jj], /cumulative)
        dil_ph[*,jj] = dil_ph[*,jj] + beb_ph[*,jj]; + tgt_ph[*,jj]
        dil_ph[*,jj] = total(dil_ph[sind[*,jj],jj], /cumulative)
        beb_ph[*,jj] = total(beb_ph[sind[*,jj],jj], /cumulative)
      end      

      noises = dblarr(n_elements(det), npix_max)
      dilution = dblarr(n_elements(det), npix_max)
      ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
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
      eclipse[det].snrhr = 1. / minnoise
      for ss=0,ndet-1 do eclipse[det[ss]].star_ph = star_ph[eclipse[det[ss]].npix-1,ss]
      
      ; Dilute the BEBs (and HEBs) now that we tricked aperture into excluding the bebs
      if (bebdil[0] ne -1) then begin
        beb_pixph = dblarr(nbeb)
        for tt=0,nbeb-1 do beb_pixph[tt] = beb_ph[eclipse[det[bebdil[tt]]].npix-1,bebdil[tt]]
        beb_starph = eclipse[det[bebdil]].star_ph
        minnoise[bebdil] = minnoise[bebdil]*beb_starph/beb_pixph
        eclipse[det[bebdil]].snrhr = 1. / minnoise[bebdil]
        eclipse[det[bebdil]].dil = eclipse[det[bebdil]].dil*beb_starph/beb_pixph
        ;eclipse[det[bebdil]].dep1_eff =  eclipse[det[bebdil]].dep1_eff*beb_pixph/beb_starph
        ;eclipse[det[bebdil]].dep2_eff =  eclipse[det[bebdil]].dep2_eff*beb_pixph/beb_starph
      end  
      ; Calculate SNR in phase-folded lightcurve
      et1_eclp = eclipse[det].dur1_eff * 24.0
      et2_eclp = eclipse[det].dur2_eff * 24.0 
      eclipse[det].snreclp1 = eclipse[det].dep1_eff / sqrt(minnoise^2./et1_eclp + eclipse[det].var^2.)
      eclipse[det].snreclp2 = eclipse[det].dep2_eff / sqrt(minnoise^2./et2_eclp + eclipse[det].var^2.)
      eclipse[det].snr1 = eclipse[det].snreclp1 * sqrt(double(eclipse[det].neclip_obs1))
      eclipse[det].snr2 = eclipse[det].snreclp2 * sqrt(double(eclipse[det].neclip_obs2))
      eclipse[det].snr  = sqrt(eclipse[det].snr1^2. + eclipse[det].snr2^2.)
   
      det1 = where((eclipse.neclip_obs1 ge NTRA_OBS_MIN) and $
             (eclipse.snr1 ge SNR_MIN))
      det2 = where((eclipse.neclip_obs2 ge NTRA_OBS_MIN) and $
             (eclipse.snr2 ge SNR_MIN))
      det = WHERE(((eclipse.neclip_obs1 + eclipse.neclip_obs2) ge NTRA_OBS_MIN) and $
            (eclipse.snr ge SNR_MIN))
      if (det1[0] ne -1) then eclipse[det1].det1 = 1
      if (det2[0] ne -1) then eclipse[det2].det2 = 1
      if (det[0] ne -1)  then begin
        eclipse[det].det = 1
        print, 'Detected ', n_elements(det), ' eclipses.'
      endif
    end ;det if
  end ; obs if
  assert, max(cr) eq 0, 'error: nonzero cosmic ray flux. means blowback to stack_prf'
end

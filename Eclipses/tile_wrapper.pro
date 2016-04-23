PRO tile_wrapper, fpath, fnums, outname, ps_only=ps_only, detmag=detmag, $
eclip=eclip, n_trial=n_trial, eclass=eclass, pla_err=pla_err, prf_file=prf_file, $
prototypeMode=prototypeMode,fCamCoordPri=fCamCoordPri,fCamCoordExt=fCamCoordExt, $
psPriFile=psPriFile,psExtFile=psExtFile, burtCatalog=burtCatalog

TIC ; Grab initial system time

; Assign prototype & extMission
numfil = N_ELEMENTS(fnums)   
randomNumbers = RANDOMU(seed, numfil)
randomIndices = SORT(ROUND(randomNumbers * (numfil-1)))
randomIndices = randomIndices[0:ROUND(numfil/10)]
assert, (prototypeMode eq 3) or (prototypeMode eq 2) or (prototypeMode eq 1) $
        or (prototypeMode eq 0), 'Invalid prototype mode assignment.'
case prototypeMode of
  0: fnums = fnums			; look at all tiles
  1: fnums = fnums[randomIndices[0]]	; look at a single random tile
  2: fnums = fnums[randomIndices[0:9]] ; look at 10 random tiles
  3: fnums = fnums[randomIndices] 	; look at ~300 random tiles
endcase

numfil = N_ELEMENTS(fnums)
if fCamCoordExt ne '' then extMission = 1 else extMission = 0

; Input files
ph_file = 'ph_T_filt.fits' ; photon fluxes for T=10 vs Teff
cr_file = 'crnoise.fits'   
tic_file = 'tic_teff.fits'
dart_file = 'dartmouth_grid.sav'
var_file = 'starvar.fits'
tband_file = 'tband.csv'

; User-adjustable settings (yes, that's you!)
assert, ~burtCatalog, 'b/c this was a one-timer'
nparam = 74 ; output table width
fov = 24. ; degrees
seg = 13  ; number of segments per hemisphere
effarea = 69.1 ;43.9 ;54.9 ;100. ;54.9 ;69.1 ; in cm^2. 
readnoise = 10. ;10. ;10.0 ; in e- per subexposure
subexptime = 2.0 ; sec in subexposure
thresh = 5.0 ; detection threshold in phase-folded lightcurve
tranmin = 1.0 ; minimum number of eclipses for detection
sys_limit=60. ;60. in ppm/hr
ffi_len=30. ; in minutes
ps_len=2. ; in minutes
duty_cycle=100.+fltarr(numfil) ; Time blanked around perigee (in minutes)
min_depth=1D-6 ; minimum transit depth to retain from eclipses
max_depth=1.0; maximum transit depth to retain from EBs
if (keyword_set(n_trial)) then n_trial=n_trial else n_trial = 10 ; number of trials in this run
if (keyword_set(prf_file)) then frac_file=prf_file $
  else frac_file = 'psfs/dfrac_t75_f3p31.fits' ; prf file 
if (keyword_set(ps_only)) then ps_only = 1 else ps_only = 0 ; Only run postage stamps?
if (keyword_set(detmag)) then detmag=detmag else detmag = 0 ; Only run postage stamps?
saturation=150000. ; e-
CCD_PIX = 4096. ; entire camera
GAP_PIX = 2.0/0.015 ; for 2 mm gap
orbit_period = 13.66d0 ; days per orbit
downlink = 16.0d0/24.0d0 ;3.66 for level1 ; downlink time in days
aspix = fov*3600/(CCD_PIX+GAP_PIX) ;20.43 ; arcseconds per pixel
if (keyword_set(eclass)) then eclass = eclass else $
eclass = [	1, $ ; Planets
            0, $ ; EBs
            0, $ ; BEBs
            0, $ ; HEBs
            0  ] ; BTPs

REARTH_IN_RSUN = 0.0091705248
AU_IN_RSUN = 215.093990942

numtargets = lonarr(numfil) ; numfil is n_elements(fnums), or ~3000 for standard all-sky.
numbkgnd = lonarr(numfil)
numdeeps = lonarr(numfil)
numPriPs = lonarr(numfil)
numExtPs = lonarr(numfil)
; Open the fits files
frac_fits = mrdfits(frac_file)
ph_fits = mrdfits(ph_file)
var_fits = mrdfits(var_file)
readcol, 'tband.csv', lam, t
tband = [[lam],[t]]
cr_fits = fltarr(100,64)
restore, dart_file
dartstruct = ss
tic_fits = mrdfits(tic_file)
restore, psPriFile ; Selected postage stamps: 0th column is stat, 1st healpix number, 2nd star ID
psPriStars = outPri
psHpNum = psPriStars[*,1]
psStarID = psPriStars[*,2]
delvarx, outPri
if extMission then begin
  restore, psExtFile
  psExtStars = outExt
  psExtHpNum = psExtStars[*,1]
  psExtStarID = psExtStars[*,2]
  delvarx, outExt
endif

SPAWN, 'cat main.pro' ; Print input file to log so you know what you did

totPriDet = 0L
totExtDet = 0L
totEclCounter = 0UL
star_out = DBLARR(5E4*n_trial,nparam) ; 35 million
ext_out = DBLARR(5e4*n_trial,nparam) ; *2

for ii=848, 858 do begin ;TODO switch back for testing
;for ii=0, numfil-1 do begin
  tileClock = TIC('tileNumber-' + STRTRIM(ii, 2) + '-' + STRING(fnums[ii]))
  fopenClock = TIC('fileOpen-' + STRTRIM(ii, 2))

  ; sav files are outStarStructs for bright catalog, starStructs for backgrounds and faints.
  print, 'Restoring files for tile ', fnums[ii]
  fname = fpath+'outStar-hp'+string(fnums[ii], format='(I04)')+'.sav'
  print, fname
  restore, fname
  targets = outStar ;"bright catalog", 1.58*10^8 stars. 
  numtargets[ii] = n_elements(targets)
  fname = fpath+'bk'+string(fnums[ii], format='(I04)')+'.sav'
  restore, fname
  bkgnds = star ;[where(star.mag.ksys gt 15)], at least roughly
  numbkgnd[ii] = n_elements(bkgnds)
  fname = fpath+'dp'+string(fnums[ii], format='(I04)')+'.sav'
  restore, fname
  deeps = star ;[where(star.mag.tsys gt 21)], at least roughly
  numdeeps[ii] = n_elements(deeps)
  print, 'File #', ii, 'tile #', fnums[ii], ' numtargets', numtargets[ii], ' numbkgnd', $ 
          numbkgnd[ii], ' numdeeps', numdeeps[ii]
  TOC, fopenClock
  ; previously had a first-cut for nPointings skip here. todo might be good for ffis.

  ; Selecting postage stamp / ffi stars
  psSelClock = TIC('psSel-' + STRTRIM(ii, 2))
  targets.ffi = 1 ; by definition all "target" stars will be in FFIs
  psInd = where(psHpNum eq fnums[ii], nPriPs)
  if nPriPs gt 0 then begin
    targets[psStarID[psInd]].ffi = 0 ; these are postage stamps in primary
    psSecInd = targets[psStarID[psInd]].companion.ind
    targets[psSecInd].ffi = 0
    assert, min(targets[psStarID[psInd]].starID eq psStarID[psInd]) eq 1, 'error: postage stamp matching'
    assert, min(targets.nPntgs eq make_array(n_elements(targets))) eq 1, 'error: input target pntgs'
    numPriPs[ii] += nPriPs
  endif
  if extMission then begin
    extPsInd = where(psExtHpNum eq fnums[ii], nExtPs)
    if nExtPs gt 0 then begin
      targets[psExtStarID[extPsInd]].nPntgs = 1 ; misnamed, but keep existing libs.
      extPsSecInd = targets[psExtStarID[extPsInd]].companion.ind
      targets[extPsSecInd].nPntgs = 1
      assert, min(targets[psExtStarID[extPsInd]].starID eq psExtStarID[extPsInd]) eq 1, 'extPs must match'
      numExtPs[ii] += nExtPs
    endif
  endif
  if ps_only eq 1 and extMission then $
    targets = targets[where(targets.ffi eq 0 or targets.nPntgs eq 1, targetCount)]
  if targetCount eq 0 then begin
    print, 'Skipping loop number', ii, 'hpNum ', fnums[ii], 'b/c gets no targets ever'
    continue
  endif
  assert, targetCount gt 0
  TOC, psSelClock

  ; Loop over each trial to generate eclipses. Each eclip (as long as not in multi system)
  ; gets unique coordinates for same trial (from host star). Not true for diff trials.
  ecliplen_tot = 0L
  for jj=0,n_trial-1 do begin
    targets.cosi = -1 + 2.0*randomu(seed, n_elements(targets)) ; randomize inclination

    makeEclipseClock = TIC('makeEclipse-' + STRTRIM(ii, 2) + '-' + STRTRIM(jj, 2))
    ecliplen =  make_eclipse(targets, bkgnds, eclip_trial, frac_fits, $
          ph_fits, dartstruct, tic_fits, eclass, tband, noEclComp, pla_err=pla_err, $
          min_depth=min_depth, max_depth=max_depth, ps_only=ps_only, extMission=extMission, $
          burtCatalog=burtCatalog)
    TOC, makeEclipseClock

    if (ecliplen gt 0) then begin
      ecIdx = ulindgen(ecliplen) + totEclCounter ; max at 2^32 ~= 4e9. Should be fine.
      eclip_trial.trial = jj + 1
      eclip_trial.uniqEclipID = ecIdx
      if (ecliplen_tot gt 0) then eclip = struct_append(eclip, eclip_trial) $
      else eclip = eclip_trial
      ecliplen_tot += ecliplen
      totEclCounter += ulong(ecliplen) ; saved over tiles
    endif
  endfor

  if (ecliplen_tot gt 0) then begin
    if (detmag eq 0) then begin
      assert, where(eclip.npointings ne 0) eq -1, 'error: eclipse pointings should be zero'
      for run=0,1 do begin ; interior loop for ext missions
        if run eq 0 then fCamCoord = fCamCoordPri
        if run eq 1 and ~extMission then continue
        if run eq 1 and extMission then fCamCoord = fCamCoordExt

        ; Survey: add this run's npointings and calc field angle (with dead ccd pixels)
        eclipSurveyClock = TIC('eclipSurvey-' + STRTRIM(ii, 2))
        eclip_survey, fov, eclip, fCamCoord ; eclip has eclipses from all trials (& runs)
        assert, max(eclip[*].coord.fov_ind) eq 0, 'fov_ind is zeroed for snr consistency'
        if run eq 0 then eclipCopy = eclip ; has eclip.npointings and CCD angles. No observed params.
        if run eq 1 then begin
          new_npointings = eclip.npointings
          new_field_angle = eclip.coord.field_angle
          new_fov_ind = eclip.coord.fov_ind
          new_fov_r = eclip.coord.fov_r
          eclip = eclipCopy ; revert all transit parameters back.
          eclip.npointings = new_npointings
          eclip.coord.field_angle = new_field_angle
          eclip.coord.fov_ind = new_fov_ind
          eclip.coord.fov_r = new_fov_r
        endif
        TOC, eclipSurveyClock
        if max(eclip.npointings) eq 0 then begin
          print, 'Skipping tileNum', fnums[ii], ' in run', run, ' because no stars pointings.'
          continue
        endif

        eclipObserveClock = TIC('eclipObserve-' + STRTRIM(ii, 2))
        assert, ps_only, 'todo: add observations for FFI allowed ext missions'
        if run eq 1 and ~extMission then continue
        if run eq 0 and ~extMission then assert, 0, 'todo: add this case'
          
        print, 'Entering eclip_observe with nelements eclip:', N_ELEMENTS(eclip)
        randSeed = ii ; seed randomized aspects of eclip_observe per-tile, not per-run
        eclip_observe, eclip, targets, bkgnds, deeps, $
                frac_fits, ph_fits, cr_fits, var_fits, $
                aspix=aspix, effarea=effarea, sys_limit=sys_limit, $
                readnoise=readnoise, thresh=thresh, tranmin=tranmin, ps_len=ps_len, $
                duty_cycle=duty_cycle[ii], ffi_len=ffi_len, saturation=saturation, $
                subexptime=subexptime, dwell_time=orbit_period, downlink=downlink, $
                extMission=extMission, randSeed=randSeed, ps_only=ps_only
        print, 'Exiting eclip_observing with nelements eclip:', N_ELEMENTS(eclip)
        TOC, eclipObserveClock

        ; only save ffiClass 2 and 4 for the primary mission, and only save ffiClass 3 and
        ; 4 for the extended. Needed for random number generation in eclip_observe
        if run eq 0 and extMission then begin
          primaryEclips = where(eclip.ffiClass eq 2 or eclip.ffiClass eq 4, nPrimaryEclips)
          if nPrimaryEclips gt 0 then eclipToObs = eclip[primaryEclips] else continue
        endif
        if run eq 1 and extMission then begin
          extEclips = where(eclip.ffiClass eq 3 or eclip.ffiClass eq 4, nExtEclips)
          if nExtEclips gt 0 then eclipToObs = eclip[extEclips] else continue
        endif

        det = where(eclipToObs.det1 or eclipToObs.det2 or eclipToObs.det)
        ASSERT, ecliplen_tot gt 0, 'ecliplen_tot should be gt 0.'
        endClock = TIC('endClock-' + STRTRIM(ii, 2))

        if (det[0] ne -1) then begin
          detid = eclipToObs[det].hostid
          ndet = n_elements(det)
          companionInd = MAKE_ARRAY(ndet, /integer, value=0)
          companionIndList = targets[detid].companion.ind
          for qq=0, ndet-1 do begin
            companionInd[qq] = WHERE(targets.starid eq companionIndList[qq])
          endfor

          bins = targets[detid].pri + 2*targets[detid].sec
          tmp_star = [[eclipToObs[det].trial], [targets[detid].mag.v], [targets[detid].mag.ic], $
                    [targets[detid].mag.t], [targets[detid].mag.j], [targets[detid].mag.h], $
                    [targets[detid].mag.k], [targets[detid].teff], [eclipToObs[det].coord.elon], $
                    [eclipToObs[det].coord.elat], [eclipToObs[det].coord.glon], $
                    [eclipToObs[det].coord.glat], $
                    [eclipToObs[det].coord.ra], [eclipToObs[det].coord.dec], [eclipToObs[det].p], $
                    [eclipToObs[det].a], [eclipToObs[det].s], [eclipToObs[det].cosi], $
                    [eclipToObs[det].teff2], [eclipToObs[det].m2], [eclipToObs[det].r2], $
                    [eclipToObs[det].dep1_eff], [eclipToObs[det].dur1], [eclipToObs[det].neclip_obs1], $
                    [eclipToObs[det].teff1], [eclipToObs[det].m1], [eclipToObs[det].r1], $
                    [eclipToObs[det].dep2_eff], [eclipToObs[det].dur2], [eclipToObs[det].neclip_obs2], $
                    [eclipToObs[det].snreclp1], [eclipToObs[det].gress1], [eclipToObs[det].snreclp2], $
                    [eclipToObs[det].gress2], [eclipToObs[det].k], [eclipToObs[det].snrhr], $
                    [eclipToObs[det].star_ph], [eclipToObs[det].bk_ph], [eclipToObs[det].zodi_ph], $
                    [eclipToObs[det].npix], [eclipToObs[det].dil], [targets[detid].ffi], $
                    [eclipToObs[det].npointings], [eclipToObs[det].sat], [eclipToObs[det].coord.fov_r], $
                    [eclipToObs[det].class], [eclipToObs[det].sep], [eclipToObs[det].icsys], $
                    [eclipToObs[det].tsys],  [eclipToObs[det].jsys], [eclipToObs[det].kpsys], $
                    [eclipToObs[det].censhift1], [eclipToObs[det].censhift2], [eclipToObs[det].cenerr1], $
                    [eclipToObs[det].cenerr2], [eclipToObs[det].var], [eclipToObs[det].coord.healpix_n], $
                    [eclipToObs[det].mult], [eclipToObs[det].tmult], [eclipToObs[det].pr], $
                    [bins], [targets[detid].companion.sep], [targets[companionInd].mag.t], $
                    [targets[detid].mag.dm], [targets[detid].age], [eclipToObs[det].det], $
                    [eclipToObs[det].det1], [eclipToObs[det].det2], [eclipToObs[det].hostid], $
                    [targets[detid].mag.micsys], [eclipToObs[det].ffiClass], $
                    [eclipToObs[det].uniqEclipID], [eclipToObs[det].dur1_eff], [eclipToObs[det].dur2_eff]]
          if run eq 0 then begin
            idx = lindgen(ndet) + totPriDet
            star_out[idx,*] = tmp_star
            totPriDet += ndet
          endif
          if run eq 1 then begin
            idx = lindgen(ndet) + totExtDet
            ext_out[idx,*] = tmp_star
            totExtDet += ndet
          endif
        endif
        TOC, endClock
        TOC, tileClock
      endfor ; ext mission loop (run per tile)
    endif
  endif
endfor ; end tile loop 
TOC 
print, 'Reached end at totPriDet = ', totPriDet
if extMission then print, 'Reached end at totExtDet = ', totExtDet
if totPriDet gt 0 and ~extMission then mwrfits, star_out[0:(totPriDet-1),*], outname ; write
if totPriDet gt 0 and totExtDet gt 0 and extMission then begin
  periodLoc = strpos(outname, '.')
  outNameNoExt = strmid(outname, 0, periodLoc)
  mwrfits, star_out[0:(totPriDet-1),*], outNameNoExt+'-pri.fits'
  mwrfits, ext_out[0:(totExtDet-1),*], outNameNoExt+'-ext.fits'
endif
print, total(numPriPs), ' primary postage stamps assigned'
assert, total(numPriPs) eq 200000
if extMission then print, total(numExtPs), ' extended postage stamps assigned'
assert, total(numExtPs) eq 200000
END

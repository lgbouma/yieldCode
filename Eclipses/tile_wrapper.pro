PRO tile_wrapper, fpath, fnums, outname, ps_only=ps_only, detmag=detmag, $
eclip=eclip, n_trial=n_trial, eclass=eclass, pla_err=pla_err, prf_file=prf_file, $
prototypeMode=prototypeMode,fCamCoord=fCamCoord,psFile=psFile,fTilesCounts=fTilesCounts, $
burtCatalog=burtCatalog

TIC ; Grab initial system time

; Assign prototype mode
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

; Assign extended mission mode
MT = 'I,I,F,F' ; get camera pointing info for observing specification 
READCOL, fCamCoord, F=FMT, pointingNumber, camNumber, camElat, camElong
nPointings = MAX(pointingNumber)+1 ; 26 pntgs 2yr, 52 4yr
nCams = 4
nCamPointings = nPointings * nCams
assert, nCamPointings eq N_ELEMENTS(pointingNumber), $
  'Need as many camPointings as elements passed in camPointingFile'
assert, (nPointings eq 26) or (nPointings eq 52), $
  'Invalid camera pointings specified. Need 2 years or 4 years.'
if nPointings eq 52 then extMission = 1 else extMission = 0
if extMission then missionCount = 2 else missionCount = 1 ; count used in for loop

; Input files
ph_file = 'ph_T_filt.fits' ; photon fluxes for T=10 vs Teff
cr_file = 'crnoise.fits'   
tic_file = 'tic_teff.fits'
dart_file = 'dartmouth_grid.sav'
var_file = 'starvar.fits'
tband_file = 'tband.csv'

; User-adjustable settings (yes, that's you!)
if burtCatalog eq 1 then nparam = 70 ; output table width
if burtCatalog eq 0 then nparam = 71
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

; Here we go!
numtargets = lonarr(numfil) ; numfil is n_elements(fnums), or ~3000 for standard all-sky.
numbkgnd = lonarr(numfil)
numdeeps = lonarr(numfil)
numps = lonarr(numfil)
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
; Get selected postage stamp data
restore, psFile
psStars = out
psStat = psStars[*,0]
psHpNum = psStars[*,1]
psStarID = psStars[*,2]
delvarx, out

SPAWN, 'cat main.pro' ; Print input file to log so you know what you did

; Pre-allocate for output: TODO possible to make preallocation less hacky?
assert, burtCatalog ne 1; TODO make less hacky
bigOutNumber = 1e6
totdet = 0L
star_out = DBLARR(bigOutNumber+1E4*n_trial,nparam)

for ii=0, numfil-1 do begin
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
  ;delvarx, star ; slow, rewrites anyway
  ;delvarx, outStar
  TOC, fopenClock

  print, 'Getting number of pointings each star receives'
  get_star_pntgs, fCamCoord, targets, 'placeholder.foo'
  print, 'Got number of pointings each star receives'

  if max(targets.npntgs) eq 0 then begin
    print, 'Skipping tileNum', fnums[ii], ' because no stars (ps or ffi) get any pointings.'
    continue
  endif

  ; Selected postage stamp / ffi stars
  psSelClock = TIC('psSel-' + STRTRIM(ii, 2))
  targets.ffi = 1 ; by definition all "target" stars will be in FFIs
  psInd = WHERE(psHpNum eq fnums[ii])
  targets[psStarID[psInd]].ffi = 0
  psSecInd = targets[psStarID[psInd]].companion.ind
  targets[psSecInd].ffi = 0
  assert, MIN(targets[psStarID[psInd]].starID eq psStarID[psInd]) eq 1, 'postage stamp matching'
  numps[ii] += N_ELEMENTS(WHERE(targets[psStarID[psInd]].ffi eq 0))
  TOC, psSelClock  

  if ps_only eq 1 then targets = targets[where(targets.ffi eq 0)]

  ecliplen_tot = 0L

	; At this point in the sim, "targets" is an array of outStarStruct objs. Smaller if ps_only.
  ; We loop over each trial to generate eclipses. Each eclip (as long as not in multi system) 
  ; gets unique coordinates for same trial (from host star). Not true for diff trials.

  for jj=0,n_trial-1 do begin
    targets.cosi = -1 + 2.0*randomu(seed, n_elements(targets)) ; randomize inclination
    ; Add eclipses. E.g., eclip_trial will come back with 2 eclipses, from 80 planets
    ; created about 7000 stars on the tile.
    makeEclipseClock = TIC('makeEclipse-' + STRTRIM(ii, 2) + '-' + STRTRIM(jj, 2))
    ecliplen =  make_eclipse(targets, bkgnds, eclip_trial, frac_fits, $
          ph_fits, dartstruct, tic_fits, eclass, tband, noEclComp, pla_err=pla_err, $
          min_depth=min_depth, max_depth=max_depth, ps_only=ps_only, $
          burtCatalog=burtCatalog)
    TOC, makeEclipseClock

    if (ecliplen gt 0) then begin
      eclip_trial.trial = jj + 1
      if (ecliplen_tot gt 0) then eclip = struct_append(eclip, eclip_trial) $
      else eclip = eclip_trial
      ecliplen_tot += ecliplen
    endif
  endfor

  if (ecliplen_tot gt 0) then begin
    if (detmag eq 0) then begin
      nPntgInit = eclip.npointings
      eclip.npointings = 0
      ; Survey: figure out npointings and field angle (nPntgs recomputed with dead ccd pixels)
      eclipSurveyClock = TIC('eclipSurvey-' + STRTRIM(ii, 2))
      eclip_survey, fov, eclip, fCamCoord
      TOC, eclipSurveyClock 
      
      ; Observe      
      eclipObserveClock = TIC('eclipObserve-' + STRTRIM(ii, 2))
      print, 'Entering eclip_observe with nelements eclip:', N_ELEMENTS(eclip)
      eclip_observe, eclip, targets, bkgnds, deeps, $
              frac_fits, ph_fits, cr_fits, var_fits, $
              aspix=aspix, effarea=effarea, sys_limit=sys_limit, $
              readnoise=readnoise, thresh=thresh, tranmin=tranmin, ps_len=ps_len, $
              duty_cycle=duty_cycle[ii], ffi_len=ffi_len, saturation=saturation, $
              subexptime=subexptime, dwell_time=orbit_period, downlink=downlink, $
              burtCatalog=burtCatalog
      print, 'Exiting eclip_observing with nelements eclip:', N_ELEMENTS(eclip)
      TOC, eclipObserveClock 

      det = where(eclip.det1 or eclip.det2 or eclip.det) ; only saves detected transits/eclipses
    endif 
	  ASSERT, ecliplen_tot gt 0, 'ecliplen_tot should be gt 0.'
    endClock = TIC('endClock-' + STRTRIM(ii, 2))

    ; Write to fits
    if (det[0] ne -1) then begin
      detid = eclip[det].hostid
      ndet = n_elements(det)
      companionInd = MAKE_ARRAY(ndet, /integer, value=0) ; TODO clean up (e.g., a write function)
      companionIndList = targets[detid].companion.ind
      for qq=0, ndet-1 do begin
        companionInd[qq] = WHERE(targets.starid eq companionIndList[qq])
      endfor

      bins = targets[detid].pri + 2*targets[detid].sec
      tmp_star = [[eclip[det].trial], [targets[detid].mag.v], [targets[detid].mag.ic], $
                [targets[detid].mag.t], [targets[detid].mag.j], [targets[detid].mag.h], $
                [targets[detid].mag.k], [targets[detid].teff], [eclip[det].coord.elon], $
                [eclip[det].coord.elat], [eclip[det].coord.glon], [eclip[det].coord.glat], $
                [eclip[det].coord.ra], [eclip[det].coord.dec], [eclip[det].p], $
                [eclip[det].a], [eclip[det].s], [eclip[det].cosi], $
                [eclip[det].teff2], [eclip[det].m2], [eclip[det].r2], $
                [eclip[det].dep1_eff], [eclip[det].dur1], [eclip[det].neclip_obs1], $
                [eclip[det].teff1], [eclip[det].m1], [eclip[det].r1], $ 
                [eclip[det].dep2_eff], [eclip[det].dur2], [eclip[det].neclip_obs2], $
                [eclip[det].snreclp1], [eclip[det].gress1], [eclip[det].snreclp2], $
                [eclip[det].gress2], [eclip[det].k], [eclip[det].snrhr], $
                [eclip[det].star_ph], [eclip[det].bk_ph], [eclip[det].zodi_ph], $
                [eclip[det].npix], [eclip[det].dil], [targets[detid].ffi], $
                [eclip[det].npointings], [eclip[det].sat], [eclip[det].coord.fov_r], $
                [eclip[det].class], [eclip[det].sep], [eclip[det].icsys], $
                [eclip[det].tsys],  [eclip[det].jsys], [eclip[det].kpsys], $ 
                [eclip[det].censhift1], [eclip[det].censhift2], [eclip[det].cenerr1], $
                [eclip[det].cenerr2], [eclip[det].var], [eclip[det].coord.healpix_n], $
                [eclip[det].mult], [eclip[det].tmult], [eclip[det].pr], $
                [bins], [targets[detid].companion.sep], [targets[companionInd].mag.t], $
                [targets[detid].mag.dm], [targets[detid].age], [eclip[det].det], $
                [eclip[det].det1], [eclip[det].det2], [eclip[det].hostid], $
                [targets[detid].mag.micsys], [targets[detid].mag.mvsys]]
      idx = lindgen(ndet) + totdet
      star_out[idx,*] = tmp_star
      totdet += ndet
    endif
    TOC, endClock 
    TOC, tileClock 
  endif
endfor ; end tile loop 
TOC 
print, 'Reached end at totdet = ', totdet
if (totdet gt 0) then mwrfits, star_out[0:(totdet-1),*], outname
print, total(numps), ' postage stamps assigned'
END

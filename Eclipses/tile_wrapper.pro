PRO tile_wrapper, fpath, fnums, outname, ps_only=ps_only, detmag=detmag, $
	eclip=eclip, n_trial=n_trial, eclass=eclass, pla_err=pla_err, prf_file=prf_file, $
	prototypeMode=prototypeMode,fCamCoord=fCamCoord,fTilesCounts=fTilesCounts, $
	radCutoff=radCutoff,burtCatalog=burtCatalog,pepperCatalog=pepperCatalog, $
	asteroseisCalc=asteroseisCalc

  TIC ; Grab initial system time

  ASSERT, asteroseisCalc eq 1, 'Why else would you be on this branch??'

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
  if burtCatalog eq 0 then nparam = 69
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
  duty_cycle=100.+fltarr(numfil) ; Time blanked around apogee (in minutes)
  ;min_depth=1D-6 ; minimum transit depth to retain from eclipses
  min_depth=1D-9 ; minimum transit depth to retain from eclipses TODO modd'd for asteroseis trials
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

  SPAWN, 'cat main.pro' ; Print input file to log so you know what you did

  ; Pre-allocate for output
  smallOutNumber = 1
  totdet = 0L
  star_out = dblarr(smallOutNumber+1E6*n_trial,nparam)
  nPepperEls = 9
  pepperDat = DBLARR(1, nPepperEls)

  inDat = read_csv_col('../asteroseis/asteroseisStars-2.csv')
  nSG = n_elements(indat.glon)
  PRINT, 'Total ', nSG, ' possible subgiant stars', STRING(10B)

  ; Read into starStruct:
  starDat = REPLICATE({starstruct}, nSG)
  starDat.dart = inDat.dart 
  starDat.r = inDat.r 
  starDat.m = inDat.mass 
  starDat.teff = inDat.teff 
  starDat.cosi = inDat.cosi
  starDat.age = inDat.age
  starDat.mini = inDat.mini
  starDat.feh = inDat.feh
  starDat.logg = inDat.logg
  starDat.pri = inDat.pri
  starDat.sec = inDat.sec
  starDat.spl = inDat.spl
  starDat.ffi = inDat.ffi
  starDat.gc = inDat.gc
  starDat.mag.av = inDat.av
  starDat.mag.dm = inDat.dm
  starDat.mag.mv = inDat.mv
  starDat.mag.mic = inDat.mic
  starDat.mag.mj = inDat.mj
  starDat.mag.v = inDat.v
  starDat.mag.g = inDat.g
  starDat.mag.r = inDat.r
  starDat.mag.i = inDat.i
  ; omitting ic magnitudes. not given, for some reason.
  starDat.mag.z = inDat.z
  starDat.mag.kp = inDat.kp
  starDat.mag.t = inDat.t
  starDat.mag.j = inDat.j
  starDat.mag.icsys = inDat.icsys
  starDat.mag.jsys = inDat.jsys
  starDat.mag.ksys = inDat.ksys
  starDat.mag.tsys = inDat.tsys
  starDat.mag.kpsys = inDat.kpsys
  starDat.mag.mvsys = inDat.mvsys
  starDat.mag.micsys = inDat.micsys
  starDat.mag.mjsys = inDat.mjsys
  starDat.mag.h = inDat.h
  starDat.mag.k = inDat.k
  ;companion data thrown out (so binary contamination will be a real issue)

  phi = inDat.glon*!dpi/180. ; convert to rad
  theta = inDat.glat*!dpi/180.
  ang2pix_ring, 16, theta, phi, sgTileNums
  ; now sgTileNums is list of tile numbers corresponding to all the subgiants

  RESTORE, fTilesCounts ; restores a "catalog" of tile #s and their npointings.
  counter = 0.
  for ii=0, numfil-1 do begin ; input data (as is) _misses_ many tiles. Skip them.
    sgInds = where(sgTileNums eq fnums[ii])
    if N_ELEMENTS(sgInds) eq 1 then begin
      if sgInds eq -1 then begin
        print, 'Skipping tile number', fnums[ii], 'because no subgiants here'
        CONTINUE
      endif
    endif
    targets = starDat[sgInds]
    print, targets
    numtargets[ii] = n_elements(targets)

    fname = fpath+'bk'+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    bkgnds = star ;[where(star.mag.ksys gt 15)]
    numbkgnd[ii] = n_elements(bkgnds)
    fname = fpath+'dp'+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    deeps = star ;[where(star.mag.tsys gt 21)]
    numdeeps[ii] = n_elements(deeps)
    ; Debugging
    print, 'File #', ii, 'tile #', fnums[ii], ' numtargets', numtargets[ii], ' numbkgnd', $ 
	    numbkgnd[ii], ' numdeeps', numdeeps[ii]
    delvarx, star

    ; Restore ~201k random coordinates for this tile from coordLib.
    RESTORE, '../../coordLib/coordHPnum'+STRTRIM(STRING(fnums[ii]),2)+'.sav'
    ipring = coordNum[0,*]
    glon = coordNum[1,*]
    glat = coordNum[2,*]
    elon = coordNum[3,*]
    elat = coordNum[4,*]
    ra = coordNum[5,*]
    dec = coordNum[6,*]
    DELVARX, coordNum

    ; if ps_only, then all the target stars become postage stamps. Else they're FFIs.
    ; no fancy primary/secondary star cuts. 
    ; 4047 of the subgiants are in single star systems
    ; 5495 are the _primary_ in their binary system
    ; 458 are the _secondary_ in their binary system 
    targets.ffi = 1
    if ps_only eq 1 then begin
      targets.ffi = 0 
      numps[ii] += N_ELEMENTS(targets)
    endif
    ecliplen_tot = 0L

    ; Now create planets for different trials.
    for jj=0,n_trial-1 do begin
      ; Re-radomize the inclination
      targets.cosi = -1 + 2.0*RANDOMU(seed, N_ELEMENTS(targets))
      ;targets.cosi = 0.0001

      ; Add eclipses. E.g., eclip_trial will come back with 2 eclipses, from 80 planets
      ; created about 7000 stars on the tile.
      makeEclipseClock = TIC('makeEclipse-' + STRTRIM(ii, 2) + '-' + STRTRIM(jj, 2))
      eclipOut =  make_eclipse(targets, bkgnds, eclip_trial, frac_fits, $
	  			  ph_fits, dartstruct, tic_fits, eclass, tband, noEclComp, pla_err=pla_err, $
            min_depth=min_depth, max_depth=max_depth, ps_only=ps_only, $
            burtCatalog=burtCatalog, asteroseisCalc=asteroseisCalc, startNum=counter)
      ecliplen = eclipOut[0]
      nPlanetsThisTrial = eclipOut[1]
      counter += nPlanetsThisTrial ; this is assuming we have <100000 total planets.
      TOC, makeEclipseClock 

      ; Add coordinates to the eclipses
      if (ecliplen gt 0) then begin
			eclip_trial.trial = jj + 1

      ; Get indices corresponding to random coordinates at this tile.
      thispix = where(ipring eq fnums[ii])
      ncoord = n_elements(thispix)
      coordind = lindgen(ecliplen)
      assert, (ecliplen lt ncoord), 'You need unique coords for each eclip.' ; LB 16/02/05

      ; Assign every eclipse object in eclip_trial unique coords
      eclip_trial.coord.elon = elon[thispix[coordind]]
      eclip_trial.coord.elat = elat[thispix[coordind]]
      eclip_trial.coord.ra = ra[thispix[coordind]]
      eclip_trial.coord.dec = dec[thispix[coordind]]
      eclip_trial.coord.glon = glon[thispix[coordind]]
      eclip_trial.coord.glat = glat[thispix[coordind]]
      eclip_trial.coord.healpix_n = fnums[ii]

      ; If you assigned different coords to eclipses w/ same host id, no good.
      ; loop over eclip_trial & set coord to whatever was first assigned to first eclip.
      if (ecliplen gt 0) then begin
        hostList = eclip_trial.hostid
        for kk=0, ecliplen-1 do begin
          thisEclipHostID = hostList[kk]
          sameStar = WHERE(thisEclipHostID eq hostList)
          nMultiEclipse = N_ELEMENTS(sameStar) ; # of eclipsing planets sameStar has
          for ll=0, nMultiEclipse-1 do begin
            eclip_trial[sameStar[ll]].coord = eclip_trial[sameStar[0]].coord
          endfor
        endfor
      endif

      ;assert, (ecliplen gt 0), 'Want eclips to be seen. Interesting for FFI multi questions'
      if (ecliplen_tot gt 0) then eclip = struct_append(eclip, eclip_trial) $
      else eclip = eclip_trial
      ecliplen_tot += ecliplen
      endif
    endfor

    ;Loop over each trial to generate eclipses. We guarantee that each eclip (as long as not
    ;in multi system) gets unique coordinates for same _trial_. Not true for diff trials.
    ;In fact, diff trials get the same coords (but diff inclinations -> diff transiters)
    if (ecliplen_tot gt 0) then begin
      if (detmag eq 0) then begin
        ; Survey: figure out npointings and field angles
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

        if burtCatalog eq 0 then begin
          det = where(eclip.det1 or eclip.det2 or eclip.det) ; only saves detected transits/eclipses
        endif 
        if burtCatalog eq 1 then begin
          det = WHERE(eclip.det or not(eclip.det)) ; any sysm w/ >=1 transiting planet
        endif
      endif 
      ASSERT, ecliplen_tot gt 0, 'ecliplen_tot should be gt 0.'
      if (ecliplen_tot eq 0) then begin
        det = where((targets[eclip.hostid].mag.ic lt detmag) or $ 
              (targets[eclip.hostid].mag.k lt detmag) or $
              (targets[eclip.hostid].mag.v lt detmag) or $
              (eclip.icsys lt detmag) or (eclip.kpsys lt detmag))
        print, "This never should be printed!"
      endif
        
      endClock = TIC('endClock-' + STRTRIM(ii, 2))
      if (det[0] ne -1) then begin
        detid = eclip[det].hostid
        ndet = n_elements(det)
        bins = targets[detid].pri + 2*targets[detid].sec
        if burtCatalog eq 0 then begin ; standard non-transiting detected planets
           tmp_star = [[eclip[det].trial], [targets[detid].mag.v], [targets[detid].mag.ic], $
                [targets[detid].mag.t], [targets[detid].mag.j], $
                [targets[detid].mag.h], [targets[detid].mag.k], [targets[detid].teff], $
                [eclip[det].coord.elon], [eclip[det].coord.elat], $
                [eclip[det].coord.glon], [eclip[det].coord.glat], $
                [eclip[det].coord.ra], [eclip[det].coord.dec], $
                [eclip[det].p], [eclip[det].a], [eclip[det].s], [eclip[det].cosi], $
                [eclip[det].teff2], [eclip[det].m2], [eclip[det].r2], $
                [eclip[det].dep1_eff], [eclip[det].dur1], [eclip[det].neclip_obs1], $
                [eclip[det].teff1], [eclip[det].m1], [eclip[det].r1], $ 
                [eclip[det].dep2_eff], [eclip[det].dur2], [eclip[det].neclip_obs2], $
                [eclip[det].snreclp1], [eclip[det].gress1], $
                [eclip[det].snreclp2], [eclip[det].gress2], $
                [eclip[det].k], [eclip[det].snrhr], $
                [eclip[det].star_ph], [eclip[det].bk_ph], [eclip[det].zodi_ph], $
                [eclip[det].npix], [eclip[det].dil], [targets[detid].ffi], [eclip[det].npointings] ,$
                [eclip[det].sat], [eclip[det].coord.fov_r], $
                [eclip[det].class], [eclip[det].sep], $
                [eclip[det].icsys],  [eclip[det].tsys],  [eclip[det].jsys], [eclip[det].kpsys], $ 
                [eclip[det].censhift1], [eclip[det].censhift2], $
                [eclip[det].cenerr1], [eclip[det].cenerr2], $
                [eclip[det].var], [eclip[det].coord.healpix_n], $
                [eclip[det].mult], [eclip[det].tmult], [eclip[det].pr], $
                [bins], [targets[detid].companion.sep], $ 
                [targets[targets[detid].companion.ind].mag.t], $
                [targets[detid].mag.dm], [targets[detid].age], [eclip[det].det], [eclip[det].det1], $
                [eclip[det].det2], eclip[det].hostid]
        endif
        idx = lindgen(ndet) + totdet
        star_out[idx,*] = tmp_star
        totdet += ndet
      endif
      TOC, endClock 
      TOC, tileClock 
    endif
    TOC 
  endfor ; end tile loop
  print, 'Reached end at totdet = ', totdet
  if (totdet gt 0) then mwrfits, star_out[0:(totdet-1),*], outname
  print, total(numps), ' postage stamps assigned'
END

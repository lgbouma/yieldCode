PRO tile_wrapper, fpath, fnums, outname, ps_only=ps_only, detmag=detmag, $
	eclip=eclip, n_trial=n_trial, eclass=eclass, pla_err=pla_err, prf_file=prf_file, $
	prototypeMode=prototypeMode,fCamCoord=fCamCoord,fTilesCounts=fTilesCounts

  TIC ; Grab initial system time
  SPAWN, 'cat main.pro' ; Print input file to log so you know what you did

  ; Assign prototype mode
  numfil = N_ELEMENTS(fnums)   
  randomNumbers = RANDOMU(seed, numfil)
  randomIndices = SORT(ROUND(randomNumbers * (numfil-1)))
  randomIndices = randomIndices[0:ROUND(numfil/100)]
  assert, (prototypeMode eq 2) or (prototypeMode eq 1) or (prototypeMode eq 0), $
	  	'Invalid prototype mode assignment.'
  case prototypeMode of
  	0: fnums = fnums			; look at all tiles
	1: fnums = fnums[randomIndices[0]]	; look at a single random tile
	2: fnums = fnums[randomIndices[0:9]] ; look at 10 random tiles
	3: fnums = fnums[randomIndices] 	; look at 1/100th of tiles (~290), randomly selected
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

  ; Don't phuck with physics, though
  REARTH_IN_RSUN = 0.0091705248
  AU_IN_RSUN = 215.093990942
  nparam = 65 ; output table width

  ; Here we go!
  numtargets = lonarr(numfil) ; numfil is n_elements(fnums), or ~3000 for standard all-sky.
  numbkgnd = lonarr(numfil)
  numdeeps = lonarr(numfil)
  numps = lonarr(numfil)
  ; Open the fits files
  frac_fits = mrdfits(frac_file)
;  rad_fits = mrdfits(rad_file)/60. ; put into pixels
  ph_fits = mrdfits(ph_file)
  var_fits = mrdfits(var_file)
  readcol, 'tband.csv', lam, t
  tband = [[lam],[t]]
  cr_fits = fltarr(100,64)
;  cr_fits = mrdfits(cr_file)
  restore, dart_file
  dartstruct = ss
  tic_fits = mrdfits(tic_file)

  ; Make random spherical coords
  tempBigStarNumber = 5e6	; WARNING: no good for FFIs. Size set by available local memory.
  print, 'Using ', tempBigStarNumber, ' as tempBigStarNumber to make coordinate list (WARNING).'
  u = randomu(seed, tempBigStarNumber)
  v = randomu(seed, tempBigStarNumber)
  phi = 2.*!dpi*u
  theta = acos(2.*v-1.)
  ang2pix_ring, 16, theta, phi, ipring

  totdet = 0L
  star_out = dblarr(tempBigStarNumber+1E6*n_trial,nparam)

  for mc=0, missionCount-1 do begin ; for nominal mission, then possible ext mission
  RESTORE, fTilesCounts ; restores a "catalog" of tile #s and their npointings.
  for ii=0, numfil-1 do begin ; for each healpix tile
    tileClock = TIC('tileNumber-' + STRTRIM(ii, 2) + '-' + STRING(fnums[ii]))
    fopenClock = TIC('fileOpen-' + STRTRIM(ii, 2))

	; If in postage stamp only mode and tile gets no pointings, skip tile.
	if ((cat[ii].npointings eq 0) and (ps_only eq 1)) then begin
		print, 'Skipping tileNum', fnums[ii], ' because it gets no pointings.'
		continue
	endif
	assert, cat[ii].npointings ne 0, 'Did not skip tile that is not observed.'

    ; Gather the .sav files. These are starstructs. 
    print, 'Restoring files for tile ', fnums[ii]
    fname = fpath+'hp'+string(fnums[ii], format='(I04)')+'.sav'
    print, fname
    restore, fname
    targets = star ;"bright catalog", 2.11e7 total stars. 
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
    TOC, fopenClock
    
    ; Choose which stars are postage stamps vs ffis
	; With 16/02/01 debug: postage stamps are only allowed on tiles that get >=1 pnting.
	; Their distributions are weighted according to the # of pointings their _tile_ gets.
    psSelClock = TIC('psSel-' + STRTRIM(ii, 2))
    targets.ffi = 1
    pri = where(targets.pri eq 1) ; targets.pri true if target is primary star in their system
    selpri = ps_sel(targets[pri].mag.t, targets[pri].teff, targets[pri].m, targets[pri].r, ph_fits, $
			rn_pix=15., npnt=cat[ii].npointings)
    if (selpri[0] ne -1) then begin 
      targets[pri[selpri]].ffi=0
      secffi = targets[pri[selpri]].companion.ind
      targets[secffi].ffi=0
      numps[ii] = numps[ii]+n_elements(selpri)
    endif
    sing = where((targets.pri eq 0) and (targets.sec eq 0)) ; targets.sec true are secondary of binary.
	; LB 16: unclear why this isn't targets.sec eq 1...
    selsing = ps_sel(targets[sing].mag.t, targets[sing].teff, targets[sing].m, $
		targets[sing].r, ph_fits, rn_pix=15., npnt=cat[ii].npointings)
    if (selsing[0] ne -1) then begin 
      targets[sing[selsing]].ffi=0
      numps[ii] = numps[ii]+n_elements(selsing)
    endif
    TOC, psSelClock  

    ecliplen_tot = 0L

	; At this point in the sim, "targets" is an array of starStruct objs. E.g., if 
	; PS only, it'll be of length about 7700 (->200k over all obs'd tiles).

    ;Loop over each trial to generate eclipses
    for jj=0,n_trial-1 do begin
      ; Re-radomize the inclination
      targets.cosi = -1 + 2.0*randomu(seed, n_elements(targets))

      ; Add eclipses
      makeEclipseClock = TIC('makeEclipse-' + STRTRIM(ii, 2) + '-' + STRTRIM(jj, 2))
      ecliplen =  make_eclipse(targets, bkgnds, eclip_trial, frac_fits, $
	  			  ph_fits, dartstruct, tic_fits, eclass, tband, pla_err=pla_err, $
	 	          min_depth=min_depth, max_depth=max_depth, ps_only=ps_only)
      TOC, makeEclipseClock 

      if (ecliplen gt 0) then begin
      	eclip_trial.trial = jj + 1
        ; Add coordinates to the eclipses
        thispix = where(ipring eq fnums[ii])
        ncoord = n_elements(thispix)
        coordind = lindgen(ecliplen) mod ncoord

		assert, (ecliplen lt ncoord), 'You need unique coords for each eclip.' ; LB 16/01/29
		print, 'nFile', ii, ' nTrial', jj, ' nTile', fnums[ii], $
			' Ecliplen', ecliplen, ' ncoords', ncoord

        glon = phi[thispix[coordind]]*180./!dpi
        glat = (theta[thispix[coordind]]-!dpi/2.)*180./!dpi
        ; Transform from galactic healpix to ecliptic observations
        euler, glon, glat, elon, elat, select=6
        euler, glon, glat, ra, dec, select=2
        eclip_trial.coord.elon = elon
        eclip_trial.coord.elat = elat
        eclip_trial.coord.ra = ra
        eclip_trial.coord.dec = dec
        eclip_trial.coord.glon = glon
        eclip_trial.coord.glat = glat
        eclip_trial.coord.healpix_n = fnums[ii]
      
        if (ecliplen_tot gt 0) then eclip = struct_append(eclip, eclip_trial) $
        else eclip = eclip_trial
        ecliplen_tot = ecliplen_tot + ecliplen
      endif
    endfor

    if (ecliplen_tot gt 0) then begin
      if (detmag eq 0) then begin
        ; Survey: figure out npointings and field angles
		eclipSurveyClock = TIC('eclipSurvey-' + STRTRIM(ii, 2))
        eclip_survey, fov, eclip, fCamCoord
		TOC, eclipSurveyClock 
      
        ; Observe      
		eclipObserveClock = TIC('eclipObserve-' + STRTRIM(ii, 2))
        eclip_observe, eclip, targets, bkgnds, deeps, $
          frac_fits, ph_fits, cr_fits, var_fits, $
          aspix=aspix, effarea=effarea, sys_limit=sys_limit, $ ;infil=sp_name,outfile=spo_name
          readnoise=readnoise, thresh=thresh, tranmin=tranmin, ps_len=ps_len, $
          duty_cycle=duty_cycle[ii], ffi_len=ffi_len, saturation=saturation, $
          subexptime=subexptime, dwell_time=orbit_period, downlink=downlink
		TOC, eclipObserveClock 

        det = where(eclip.det1 or eclip.det2 or eclip.det)
      endif else det = where((targets[eclip.hostid].mag.ic lt detmag) or $ 
	      (targets[eclip.hostid].mag.k lt detmag) or $
              (targets[eclip.hostid].mag.v lt detmag) or $
	      (eclip.icsys lt detmag) or (eclip.kpsys lt detmag))

      endClock = TIC('endClock-' + STRTRIM(ii, 2))
      if (det[0] ne -1) then begin
        detid = eclip[det].hostid
        ndet = n_elements(det)
        bins = targets[detid].pri + 2*targets[detid].sec
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
                [bins], [targets[detid].companion.sep], [targets[targets[detid].companion.ind].mag.t], $
 		[targets[detid].mag.dm], [targets[detid].age]]
        idx = lindgen(ndet) + totdet
        star_out[idx,*] = tmp_star
        totdet = totdet+ndet
      ;stard = star[det]
      ;if (keyword_set(sav)) then save, filen=spo_name, stard
      endif
      TOC, endClock 
      TOC, tileClock 
    endif
  endfor ; end tile loop
  endfor ; end mission count loop
  TOC 
  print, 'Reached end at totdet = ', totdet
  if (totdet gt 0) then mwrfits, star_out[0:(totdet-1),*], outname
  print, total(numps), ' postage stamps assigned'
END

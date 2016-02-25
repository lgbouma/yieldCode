function addSubGiantPlanets, star, pstruct, frac, ph_p, tband, noEclComp, err=err, $
	aspix=aspix, fov=fov, dressing=dressing, min_depth=min_depth, ps_only=ps_only, $
	burtCatalog=burtCatalog, startNum=startNum
;+
; NAME: addSubGiantPlanets
; PURPOSE: Over each trial (usually 1-10), over each tile (1-2908), populate
;	the stars on this trial with ecipsing objects.
; INPUTS: 
; 1. star: a starStruct object (called from targets in tile_wrapper), selected for postage 
; 	stamps and FFIs.
; 2. pstruct: passed a dummy name, will be populated with an eclipStruct of eclipsing planet
;	data.
; 3. frac: PRF data from Deb Woods.
; 4. ph_p: Photon fluxes for T=10 vs T_eff.
; 5. tband: TESS transmission bandpass vs wavelength
; 6. err: Option from main for "nominal" or "boosted" occurrence rates. Usually 0 (for nom).
; 7. aspix: arcsecond per pixel. Defaults to 21.1
; 8. fov: field of view in degrees. Defaults to 24.
; 9. dressing: called with Dressing = True. This means split out the bins as in Fig 8,
;		Sullivan+ 2015. Otherwise, no temperature dependent occurrences.
; 10. min_depth: minimum transit depth. Called with 0.
; 11. ps_only: option from main for only giving planets about PS stars, or all stars.
; 12. burtCatalog: do you want an output eclip_trial with some objects that do NOT
;		eclipse? (specifically: if at least one planet in system is transiting, as-written
;		we return all other planets in that system + that planet. Detection fine-tuning later)
; 13. noEclComp: companions of transiting planets that do not transit (burt catalog). Is
;		eclipCompStruct type (to keep it straight)
; 14. starting index to go thru period / radius draw (input: # of already existing stars with
;   planets that have gone thru). This is necessary to keep the draw random.
; OUTPUTS:
; 1. ntra: number of transiting planets
; 2. pstruct: an _eclip_Struct that gets passed back to `make_eclipse` with all added 
;		transit parameters.
;-
  ntra = 0
  planet_multi=0
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 0.6395 ; in m/s
 
  if (keyword_set(fov)) then fov=fov else fov=24.0 
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth=0.0
  ccd_pix = 4096.0
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(err)) then err=1 else err=0
  if (KEYWORD_SET(burtCatalog)) then burtCatalog=burtCatalog else burtCatalog=0

  nStars = n_elements(star)
  if (keyword_set(ps_only)) then ps = (star.ffi ne 1) else ps = intarr(nstars) + 1

  randList = RANDOMU(seed, nStars)
  randList2 = RANDOMU(seed, nStars) ; for dupes
  randList3 = RANDOMU(seed, nStars) ; for triples
  randList4 = RANDOMU(seed, nStars) ; for quads (the end)
  FMT='F'
  READCOL, '../asteroseis/periodVals.dat', F=FMT, periods ; log-normals in IDL requires a license
  READCOL, '../asteroseis/radiiVals.dat', F=FMT, radii ; so instead draw from 1e5 python-generated vals

  ; Count total number planets each star gets, allowing maximum of 3.
  nPlanets = 0
  starIndArr = []
  nPlanetsStarArr = []
  for i=0,nStars-1 do begin ; loop over all subgiants (at this tile)
    nPlanetsThisStar = 0
    if randList[i] lt 0.1 then begin ; first draw: do you get one planet?
      nPlanetsThisStar += 1
      nPlanets += 1
      if randList2[i] lt 0.1 then begin 
        nPlanetsThisStar += 1
        nPlanets +=1
        if randList3[i] lt 0.1 then begin 
          nPlanetsThisStar += 1
          nPlanets += 1
        endif
      endif
    endif
    starIndArr = [starIndArr, i]
    nPlanetsStarArr = [nPlanetsStarArr, nPlanetsThisStar]
  endfor
  
  print,'Making ', nPlanets, ' planets...'

  if (nplanets gt 0) then begin ; then make some planets
    planet_rad = DBLARR(nPlanets)
    planet_per = DBLARR(nPlanets)
    planet_hid = LONARR(nPlanets)

    ; get host IDs
    hIDlist = []
    singInd = where(nPlanetsStarArr eq 1)
    doubInd = where(nPlanetsStarArr eq 2)
    tripInd = where(nPlanetsStarArr eq 3)
    hIDlist = [hIDlist, singInd]
    hIDlist = [hIDlist, doubInd, doubInd]
    hIDlist = [hIDlist, tripInd, tripInd, tripInd]
    hIDlist = hIDlist[sort(hIDlist)]

    for i=0,nPlanets-1 do begin ; get radii, periods, and hostIDs for planets
      planet_rad[i] = radii[i+startNum]
      planet_per[i] = periods[i+startNum]
      planet_hid[i] = hIDlist[i]
    endfor

    ; if same host has planet periods with period ratio < 1.2, fix.
    ; loop over systems
    nSystems = n_elements(uniq(hIDlist))
    uniqSysList = hIDlist[uniq(hIDlist)]
    ;for i=0,nSystems-1 do begin
    ;  thisSysInd = where(hIDlist eq uniqSysList[i])
    ;  thisSysPeriods = planet_per[thisSysInd]
      
    ;endfor

    allid = planet_hid
    planet_a = (star[allid].m)^(1./3.) * (planet_per/365.25)^(2./3.); in AU [Kepler 3]
    planet_s = (star[allid].r)^2. * $
               (star[allid].teff/5777.0)^4. / (planet_a)^2. ; incident flux wrt sun-earth value
    ; Equilibrium temp
    planet_teq = (star[allid].teff)*sqrt(star[allid].r/(2.*planet_a*AU_IN_RSUN))
    ; Impact parameter
    planet_b = (planet_a*AU_IN_RSUN / star[allid].r) * star[allid].cosi; assumes circular orbit
    ; Stellar radius in AU
    min_a = 2.*star[allid].r/AU_IN_RSUN 

    ; Weiss & Marcy 2014 for planet mass assignment
    planet_m = dblarr(nplanets)
    lomass = where(planet_rad lt 1.5)
    if (lomass[0] ne -1) then planet_m[lomass] = 0.440*(planet_rad[lomass])^3. + $
                              0.614*(planet_rad[lomass])^4.
    himass = where(planet_rad ge 1.5)
    if (himass[0] ne -1) then planet_m[himass] = 2.69*(planet_rad[himass])^0.93
    ; RV amplitude
    planet_k = RV_AMP*planet_per^(-1./3.) * planet_m * $ 
               sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 

    ; Select transiting planets, and then work out transit properties
    dep1 = (planet_rad*REARTH_IN_RSUN / star[allid].r )^2.0
    tra = where((abs(planet_b) lt 1.0) and (planet_a gt min_a) and (dep1 gt min_depth))
    print, 'planet_a', planet_a, ' min_a', min_a

    if (tra[0] ne -1) then begin ; if you have at least one transiting planet
      traid = planet_hid[tra] ; hostIDs of transiting planets
      ntra = n_elements(tra)
      planet_multi = intarr(ntra)
      planet_tmulti = intarr(ntra)
      planet_pr = fltarr(ntra)
      for t=0,(ntra-1) do begin ; for each transiting system
        pos = where(planet_hid eq traid[t]) ; index of planet_per where other planets live
        tos = where(traid eq traid[t]) ; other trasiting planets live
        if (n_elements(pos) gt 1) then begin
          planet_multi[t] = n_elements(pos) - 1  ; Number of other planets in the system
          planet_tmulti[t] = n_elements(tos) - 1 ; Number of other transiting planets
          pp = [tra[t]] ; index of transit indices
          opos = setdifference(pos, pp) ; index of the periods other than this transiting one
          if (opos[0] ne -1) then begin
            planet_pr[t] = min(planet_per[opos]/planet_per[tra[t]] > $ 
                           planet_per[opos]^(-1)*planet_per[tra[t]])
          endif 
        endif
      endfor
    
      planet_eclip = replicate({eclipstruct}, ntra)
          
      r2 = planet_rad[tra]*REARTH_IN_RSUN
      dep2 = phot_ratio_planet(star[traid].teff, planet_teq[tra], star[traid].mag.t, $ 
                 star[traid].mag.dm, r2, ph_p, tband)
      planet_eclip.class=1
      planet_eclip.m1 = star[traid].m
      planet_eclip.m2 = planet_m[tra]/MSUN_IN_MEARTH
      planet_eclip.k = planet_k[tra]
      planet_eclip.r1 = star[traid].r
      planet_eclip.r2 = r2
      planet_eclip.teff1 = star[traid].teff
      planet_eclip.teff2 = planet_teq[tra]
      planet_eclip.a = planet_a[tra]
      planet_eclip.s = planet_s[tra]
      planet_eclip.p = planet_per[tra]
      planet_eclip.b = planet_b[tra]
      planet_eclip.tsys = star[traid].mag.tsys
      planet_eclip.kpsys = star[traid].mag.kpsys
      planet_eclip.icsys = star[traid].mag.icsys
      planet_eclip.jsys = star[traid].mag.jsys
      planet_eclip.hostid = planet_hid[tra]
      planet_eclip.cosi = star[traid].cosi ; LB 16/02/08
      planet_eclip.dep1 = (planet_eclip.r2 / planet_eclip.r1 )^2.0
      toodeep = where(planet_eclip.dep1 gt 1.0)
      if (toodeep[0] ne -1) then planet_eclip[toodeep].dep1 = 1.0
      planet_eclip.dep2 = dep2 
      planet_eclip.mult = planet_multi
      planet_eclip.tmult = planet_tmulti
      planet_eclip.pr   = planet_pr
      planet_eclip.dur1 = planet_eclip.r1 * planet_eclip.p * $ 
              sqrt(1.-(planet_eclip.b)^2.) / (!PI*planet_eclip.a*AU_IN_RSUN)
      planet_eclip.dur2 = planet_eclip.r1 * planet_eclip.p * $ 
              sqrt(1.-(planet_eclip.b)^2.) / (!PI*planet_eclip.a*AU_IN_RSUN)
      planet_eclip.gress1 = planet_eclip.r2 * planet_eclip.p / $
                (!PI*planet_eclip.a*AU_IN_RSUN*sqrt(1.-(planet_eclip.b)^2.))
      planet_eclip.gress2 = planet_eclip.r2 * planet_eclip.p / $
                (!PI*planet_eclip.a*AU_IN_RSUN*sqrt(1.-(planet_eclip.b)^2.))
 
      pstruct=planet_eclip
    endif
  endif
  print, 'Created ', ntra, ' transiting planets out of ', nplanets, ' around ', $
         nstars, ' stars.'
  return, [ntra, nPlanets]

  ; physics: want period radios >1.2
  ; really, binary info would be important here too. but they threw it out, so i will as well.
end

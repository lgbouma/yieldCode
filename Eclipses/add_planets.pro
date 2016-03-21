function add_planets, star, pstruct, frac, ph_p, tband, noEclComp, err=err, $
	aspix=aspix, fov=fov, dressing=dressing, min_depth=min_depth, ps_only=ps_only, $
	burtCatalog=burtCatalog
;+
; NAME: add_planets
; PURPOSE: Over each trial (usually 1-10), over each tile (1-2908), populate
;	the stars on this trial with ecipsing objects.
; INPUTS: 
; 1. star: an outStarStruct object (called from targets in tile_wrapper), selected for postage 
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
  G_CM_S2 = 98.1 ; cm/s^2
 
  if (keyword_set(fov)) then fov=fov else fov=24.0 
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth=0.0
  ccd_pix = 4096.0
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(err)) then err=1 else err=0
  if (KEYWORD_SET(burtCatalog)) then burtCatalog=burtCatalog else burtCatalog=0

  nstars = n_elements(star)
  if (keyword_set(ps_only)) then ps = (star.ffi ne 1) else ps = intarr(nstars) + 1

  if(keyword_set(dressing)) then begin
    hotstars = where(star.teff ge 4000. and ps)
    nhotstars = total(star.teff ge 4000. and ps)
    coolstars = where(star.teff lt 4000. and ps) 
    ncoolstars = total(star.teff lt 4000. and ps)
  endif else begin
    hotstars = where(ps)
    nhotstars = total(ps)
    coolstars = -1
    ncoolstars = 0
  endelse
  hot_companion_p = star[hotstars].companion.p
  cool_companion_p = star[coolstars].companion.p

  if (keyword_set(verbose)) then verbose=1 else verbose=0     
 
  fressin_period = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
  fressin_radius = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0]
  planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']

  rate_fressin = dblarr(11,5) ; period bin, radius bin
  rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
  rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
  rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
  rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
  rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
  rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
  rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
  rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
  rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
  rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
  rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
  
  sig_fressin = dblarr(11,5) ; period bin, radius bin
  sig_fressin[0,*] = [0.04, 0.03, 0.011, 0.003, 0.007]
  sig_fressin[1,*] = [0.15, 0.13, 0.03,  0.006, 0.018]
  sig_fressin[2,*] = [0.43, 0.23, 0.09,  0.03,  0.03]
  sig_fressin[3,*] = [0.60, 0.56, 0.19,  0.03,  0.04]
  sig_fressin[4,*] = [0.83, 0.73, 0.39,  0.07,  0.06]
  sig_fressin[5,*] = [1.05, 1.00, 0.64,  0.08,  0.06]
  sig_fressin[6,*] = [1.88, 1.48, 1.01,  0.12,  0.10]
  sig_fressin[7,*] = [2.81, 1.21, 1.05,  0.16,  0.17]
  sig_fressin[8,*] = [0.0,  2.20, 1.03,  0.17,  0.29]
  sig_fressin[9,*] = [0.0,  0.0,  0.90,  0.21,  0.28]
  sig_fressin[10,*] =[0.0,  0.0,  0.00,   0.15,  0.30]

  dressing_period = [0.5, 1.7, 5.5, 18.2, 60.3, 200.]
  dressing_radius = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
  rate_dressing = dblarr(5,7)
  rate_dressing[0,*] = [0.87, 1.70, 0.46, 0.10, 0.0, 0.0, 0.0]
  rate_dressing[1,*] = [7.25, 8.03, 3.16, 1.99, 1.13, 0.51, 0.31]
  rate_dressing[2,*] = [45.27, 18.17, 20.07, 15.89, 5.77, 2.45, 0.85]
  rate_dressing[3,*] = [0.0, 25.79, 22.90, 24.68, 12.16, 2.94, 0.60]
  rate_dressing[4,*] = [0.0, 34.46, 19.00, 16.08, 8.05, 1.91, 0.0]

  psig_dressing = dblarr(5,7)
  psig_dressing[0,*] = [0.78, 0.87, 0.50, 0.29, 0.15, 0.13, 0.13]
  psig_dressing[1,*] = [3.26, 2.50, 1.70, 1.47, 1.18, 0.87, 0.70]
  psig_dressing[2,*] = [15.87, 5.40, 5.48, 4.96, 3.36, 2.40, 1.66]
  psig_dressing[3,*] = [0.0, 10.14, 8.70, 8.91, 6.99, 4.30, 2.32]
  psig_dressing[4,*] = [0.0, 32.79, 16.43, 12.78, 9.32, 5.59, 2.18]
 
  nsig_dressing = dblarr(5,7)
  nsig_dressing[0,*] = [0.38, 0.56, 0.22, 0.05, 0.0, 0.0, 0.0]
  nsig_dressing[1,*] = [2.16, 1.85, 1.07, 0.80, 0.52, 0.25, 0.15]
  nsig_dressing[2,*] = [9.41, 3.93, 4.06, 3.59, 2.02, 1.10, 0.42]
  nsig_dressing[3,*] = [0.0, 6.51, 5.74, 5.94, 4.09, 1.43, 0.22]
  nsig_dressing[4,*] = [0.0, 11.46, 7.28, 6.12, 3.62, 0.83, 0.0]

  rate_fressin = rate_fressin/100.;
  sig_fressin = sig_fressin/100.;
  rate_dressing = rate_dressing/100.;
  nsig_dressing = nsig_dressing/100.;
  psig_dressing = psig_dressing/100.;

  nplanets = total(round((rate_fressin + sig_fressin) * nhotstars)) + $
             total(round((rate_dressing + psig_dressing) * ncoolstars))
  print, 'making ', nplanets, ' planets.'
  if (nplanets gt 0) then begin
    ; Pre-allocate for speed
    ;planet = replicate(template_planet, total(nplanets))
    planet_rad = dblarr(nplanets)
    planet_per = dblarr(nplanets)
    planet_hid = lonarr(nplanets) - 99 ; array of "-99"'s, length nPlanets
    idx0 = 0L
    for i=0,(n_elements(fressin_period)-2) do begin
      for j=(n_elements(fressin_radius)-2),0,-1 do begin
        thisrate = rate_fressin[i,j]
        if err then thisrate = thisrate + sig_fressin[i,j]*randomn(seed);
        binplanets = round(thisrate * nhotstars)
        if (binplanets gt 0) then begin 
           ; planet_per is either 0 or a smaller number since period is increasing.
           ; 1. Find the bad periods
           bd_ind = where((planet_per ge fressin_period[i]/1.2))
           ; 2. Find the remaining periods
           bd_hid = where((hot_companion_p lt fressin_period[i+1]*5.0) and $
           (hot_companion_p gt fressin_period[i]/5.0))
           if (bd_ind[0] ne -1) then gd_hid = setdifference(hotstars, planet_hid[bd_ind]) $
           else gd_hid = hotstars
           if (bd_hid[0] ne -1) then gd_hid = setdifference(gd_hid, bd_hid)
           if (gd_hid[0] ne -1) then begin
             ngdhot = n_elements(gd_hid)
             if (binplanets lt ngdhot) then begin
               rand_inds = sort(randomu(seed, ngdhot))
               use_inds = rand_inds[0:binplanets-1]
               hostids = gd_hid[use_inds]
             endif else begin 
               print, 'needed ', binplanets, ' but got ', ngdhot, ' hot ones' 
               hostids = gd_hid
               binplanets = ngdhot
             endelse
             if(j eq 0) then radpow = 0.0 else radpow = -1.7
             randomp, periods, -1.0, binplanets, $
               range_x = [fressin_period[i], fressin_period[i+1]], seed=seed
             randomp, radii, radpow, binplanets, $
               range_x = [fressin_radius[j], fressin_radius[j+1]], seed=seed
             ; Period check
             for k=0, (binplanets-1) do begin
               same_host = where(planet_hid eq hostids[k])
               if(total(planet_per[same_host] gt periods[k]/1.2 )) then begin
                 print, 'added ', periods[k], ' to ', planet_per[same_host]
                 stop
               endif
             endfor
             idx = lindgen(binplanets) + idx0
             planet_per[idx] = periods
             planet_rad[idx] = radii
             planet_hid[idx] = hostids
             
             idx0 = max(idx) + 1
           endif
        endif
      endfor
    endfor
    if (keyword_set(dressing)) then begin
      for i=0,(n_elements(dressing_period)-2) do begin
        for j=(n_elements(dressing_radius)-2),0,-1 do begin
          thisrate = rate_dressing[i,j]
          if err then begin
            thisrand = randomn(seed);
            if (thisrand gt 0) then thisrate = thisrate + psig_dressing[i,j]*thisrand $
            else thisrate = thisrate + nsig_dressing[i,j]*thisrand
          endif
          binplanets = round(thisrate * ncoolstars)
          ;print, 'binplanets in:', binplanets, ' per:', dressing_period[i], ' rad:', dressing_radius[j]
          if (binplanets gt 0) then begin
            bd_ind = where(planet_per ge dressing_period[i]/1.2)
            bd_hid = where((cool_companion_p lt dressing_period[i+1]*5.0) and $ 
                     (cool_companion_p gt dressing_period[i]/5.0))
            if (bd_ind[0] ne -1) then gd_hid = setdifference(coolstars, planet_hid[bd_ind]) $
            else gd_hid = coolstars
            if (bd_hid[0] ne -1) then gd_hid = setdifference(gd_hid, bd_hid)
            if (gd_hid[0] ne -1) then begin
              ngdcool = n_elements(gd_hid)
              if (binplanets lt ngdcool) then begin
                rand_inds = sort(randomu(seed, ngdcool))
                hostids = gd_hid[rand_inds[0:binplanets-1]]
              endif else begin
                print, 'needed ', binplanets, ' but got ', ngdcool, ' cool ones'
                hostids = gd_hid
                binplanets = ngdcool
              endelse
              ;endif else hostids = gd_hid[floor(double(ngdcool)*randomu(seed, binplanets))]
              if(j lt 2) then radpow = 0.0 else if (j eq 2) then radpow = -1 else radpow = -1.7
              randomp, periods, -1.0, binplanets, $
                range_x = [dressing_period[i], dressing_period[i+1]], seed=seed
              randomp, radii, radpow, binplanets, $
                range_x = [dressing_radius[j], dressing_radius[j+1]], seed=seed
              idx = lindgen(binplanets) + idx0
              planet_per[idx] = periods
              planet_rad[idx] = radii
              planet_hid[idx] = hostids
              idx0 = max(idx) + 1
            endif
          endif
        endfor
      endfor
    endif
    planet_per = planet_per[0:idx0-1]
    planet_rad = planet_rad[0:idx0-1]
    planet_hid = planet_hid[0:idx0-1]
    if idx0 gt 0 then nplanets = idx0 ;16/03/19 kind of surprising bug
    if idx0 eq 0 then begin
      allZeros = make_array(n_elements(planet_radi), value=0)
      assert, array_equal(planet_rad, allZeros)
    endif
    ; Tally multi-planet systems
    ; Work out orbital distance and impact parameter
    allid = planet_hid
    planet_a = (star[allid].m)^(1./3.) * (planet_per/365.25)^(2./3.); in AU
    planet_s = (star[allid].r)^2.0 * $
               (star[allid].teff/5777.0)^4. / (planet_a)^2. ; incident flux wrt sun-earth value
    ; Equilibrium Temp.
    planet_teq = (star[allid].teff)*sqrt(star[allid].r/(2.*planet_a*AU_IN_RSUN))
    ; Impact parameter
    planet_b = (planet_a*AU_IN_RSUN / star[allid].r) * star[allid].cosi; assumes circular orbit
    ; Stellar radius in AU
    min_a = 2.*star[allid].r/AU_IN_RSUN
  
    ; Weiss & Marcy 2014:
    planet_m = dblarr(nplanets)
    lomass = where(planet_rad lt 1.5)
    if (lomass[0] ne -1) then planet_m[lomass] = 0.440*(planet_rad[lomass])^3. + $
                              0.614*(planet_rad[lomass])^4.
    himass = where(planet_rad ge 1.5)
    if (himass[0] ne -1) then planet_m[himass] = 2.69*(planet_rad[himass])^0.93
    ; RV amplitude
    planet_k = RV_AMP*planet_per^(-1./3.) * planet_m * $ 
               sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 

    ; compPlaInd holds indices of nontransiting planets in systems w/ at least one transiting planet
    ; transPlaInd holds indices of transiting planets
    assert, burtCatalog ne 1, 'so i can clean up this hacky mess'

    ; Select transiting planets, and then work out transit properties
    dep1 = (planet_rad*REARTH_IN_RSUN / star[allid].r )^2.0
    tra = where((abs(planet_b) lt 1.0) and (planet_a gt min_a) and (dep1 gt min_depth))
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
      planet_eclip.coord = star[traid].coord ; outStarStruct benefit
      planet_eclip.npointings = star[traid].npntgs ; n.b. different from neclip_obs!
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
  return, ntra
end

function add_btps, star, bkgnd, estruct, frac, ph_p, tband, $ ;input
           mult, ps_only=ps_only, dressing=dressing, aspix=aspix, radmax=radmax, fov=fov, min_depth=min_depth

  keep = 0

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  if (keyword_set(fov)) then fov=fov else fov=24.0
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth=0.0
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(radmax)) then radmax=radmax else radmax=2.0
  if (keyword_set(sq_deg)) then sq_deg=sq_deg else sq_deg=13.54
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels
  radas = sqrt(sq_deg/!dpi)*3600.
  radpix = radas/aspix
  thresh = (double(radmax)/double(radpix))^2.

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 1731446.8 ; in m/s
  G_CM_S2 = 98.1 ; cm/s^2
 
;  if (keyword_set(fov)) then fov=fov else fov=24.0
;  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
;  ccd_pix = 4096.0
;  PIX_SCALE = fov*3600./ccd_pix

  nstars = n_elements(star)
  nbks = n_elements(bkgnd)
  nebtot = 0L

  if(keyword_set(dressing)) then begin
    hotstars = where(bkgnd.teff ge 4000.)
    nhotstars = total(bkgnd.teff ge 4000.)
    coolstars = where(bkgnd.teff lt 4000.)
    ncoolstars = total(bkgnd.teff lt 4000.)
  endif else begin
    hotstars = lindgen(nbks)
    nhotstars = nbks
    coolstars = -1
    ncoolstars = 0
  endelse
  
  companion_p = bkgnd.companion.p

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
  rate_fressin = rate_fressin/100.

  dressing_period = [0.68, 1.2, 2.0, 3.4, 5.9, 10., 17., 29., 50.]
  dressing_radius = [0.5, 0.7, 1.0, 1.4, 2.0, 2.8]
  rate_dressing = dblarr(8,5)
  rate_dressing[0,*] = [0.0026, 0.0020, 0.0039, 0.0, 0.0]
  rate_dressing[1,*] = [0.0, 0.0068, 0.011, 0.0019, 0.0023]
  rate_dressing[2,*] = [0.011, 0.013, 0.015, 0.011, 0.0]
  rate_dressing[3,*] = [0.0, 0.037, 0.028, 0.016, 0.0046]
  rate_dressing[4,*] = [0.0, 0.051, 0.050, 0.051, 0.031]
  rate_dressing[5,*] = [0.0, 0.038, 0.045, 0.030, 0.039]
  rate_dressing[6,*] = [0.0, 0.065, 0.048, 0.066, 0.11]
  rate_dressing[7,*] = [0.0, 0.0, 0.084, 0.027, 0.0]

  nplanets = total(round(rate_fressin * nhotstars)) + $
                total(round(rate_dressing * ncoolstars))
     ; Pre-allocate for speed
     ;planet = replicate(template_planet, total(nplanets))

;  pris = where(bkgnd.pri eq 1)
;  secs = where(bkgnd.sec eq 1)

   ; npri = n_elements(pris)
   ; r1 = bkgnd[pris].r
   ; r2 = bkgnd[bkgnd[pris].companion.ind].r
   ; m1 = bkgnd[pris].m
   ; m2 = bkgnd[bkgnd[pris].companion.ind].m
   ; teff1 = bkgnd[pris].teff
   ; teff2 = bkgnd[bkgnd[pris].companion.ind].teff
   ; tmag1 = bkgnd[pris].mag.t
   ; tsys = bkgnd[pris].mag.tsys
   ; kpsys = bkgnd[pris].mag.kpsys
   ; icsys = bkgnd[pris].mag.icsys
   ; jsys = bkgnd[pris].mag.jsys
   ; tmag2 = bkgnd[bkgnd[pris].companion.ind].mag.t
   ; ars = bkgnd[pris].companion.a*AU_IN_RSUN
   ; a   = bkgnd[pris].companion.a
   ; p = bkgnd[pris].companion.p
   ; ecc = bkgnd[pris].companion.ecc
    ;f = bkgnd[pris].companion.f
     
  for ii=0,mult-1 do begin
    ; Start over with new arrays
    planet_rad = dblarr(nplanets)
    planet_per = dblarr(nplanets)
    planet_hid = lonarr(nplanets)
    idx0 = 0L
    ; Re-randomize the inclination and w
    cosi = -1.0 + 2.0*randomu(seed, nbks)
    for i=0,(n_elements(fressin_period)-2) do begin
      for j=(n_elements(fressin_radius)-2),0,-1 do begin
        binplanets = round(rate_fressin[i,j] * nhotstars)
        if (binplanets gt 0) then begin
              ; planet_per is either 0 or a smaller number since period is increasing.
              ; 1. Find the bad periods
              bd_ind = where((planet_per ge fressin_period[i]/1.2) or $
		((companion_p lt fressin_period[i+1]*5.0) and (companion_p gt fressin_period[i]/5.0)))
              ; 2. Find the remaining periods
              gd_hid = hotstars
              if (bd_ind[0] ne -1) then gd_hid = setdifference(hotstars, planet_hid[bd_ind])
              if (gd_hid[0] ne -1) then begin
                ngdhot = n_elements(gd_hid)
                if (binplanets lt ngdhot) then begin
                  rand_inds = sort(randomu(seed, ngdhot))
                  hostids = gd_hid[rand_inds[0:binplanets-1]]
                endif else begin
                  hostids = gd_hid[rand_inds]
                  binplanets = ngdhot
                endelse
                ;endif else hostids = gd_hid[floor(double(ngdhot)*randomu(seed, binplanets))]
                if(j eq 0) then radpow = 0.0 else radpow = -1.7
                randomp, periods, -1.0, binplanets, $
		  range_x = [fressin_period[i], fressin_period[i+1]], seed=seed
                randomp, radii, radpow, binplanets, $
		  range_x = [fressin_radius[j], fressin_radius[j+1]], seed=seed
                ;tmp_planet.r = radii
                ;tmp_planet.p = periods
                idx = lindgen(binplanets) + idx0
                planet_per[idx] = periods
                planet_rad[idx] = radii
                planet_hid[idx] = hostids
                ;planet[idx] = tmp_planet
                ;delvar, tmp_planet
                ;if (binplanets gt 3) then stop  

                idx0 = max(idx) + 1
              endif ; gd_hid
        endif ; binplanets
      endfor ; j loop
    endfor ; i loop 
    if (keyword_set(dressing)) then begin
      for i=0,(n_elements(dressing_period)-2) do begin
        for j=(n_elements(dressing_radius)-2),0,-1 do begin
          binplanets = round(rate_dressing[i,j] * ncoolstars)
          if (binplanets gt 0) then begin
              bd_ind = where((planet_per ge dressing_period[i]/1.2) or $
	        ((companion_p lt dressing_period[i+1]*5.0) and (companion_p gt dressing_period[i]/5.0)))
              gd_hid = coolstars
              if (bd_ind[0] ne -1) then gd_hid = setdifference(coolstars, planet_hid[bd_ind])
              if (gd_hid[0] ne -1) then begin
                ngdcool = n_elements(gd_hid)
                if (binplanets lt ngdcool) then begin
                  rand_inds = sort(randomu(seed, ngdcool))
                  hostids = gd_hid[rand_inds[0:binplanets-1]]
                endif else begin
                  hostids = gd_hid[rand_inds]
                  binplanets = ngcool
                endelse
                ;endif else hostids = gd_hid[floor(double(ngdcool)*randomu(seed, binplanets))]
                if(j lt 2) then radpow = 0.0 else radpow = -1.7
                randomp, periods, -1.0, binplanets, $
	          range_x = [dressing_period[i], dressing_period[i+1]], seed=seed
                randomp, radii, radpow, binplanets, $
		  range_x = [dressing_radius[j], dressing_radius[j+1]], seed=seed
                ;tmp_planet.r = radii
                ;tmp_planet.p = periods
                idx = lindgen(binplanets) + idx0
                planet_per[idx] = periods
                planet_rad[idx] = radii
                planet_hid[idx] = hostids
                ;print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
                ;planet[idx] = tmp_planet
                ;delvar, tmp_planet
                idx0 = max(idx) + 1
              endif ; gd_hid
          endif ; binplanets
        endfor ; j loop
      endfor ; i loop
    endif ; dressing
    planet_per = planet_per[0:idx0-1]
    planet_rad = planet_rad[0:idx0-1]
    planet_hid = planet_hid[0:idx0-1]
    nplanets = idx0
    ; Tally multi-planet systems
    ; Work out orbital distance and impact parameter
    allid = planet_hid
    planet_a = (bkgnd[allid].m)^(1./3.) * (planet_per/365.25)^(2./3.); in AU
    planet_s = (bkgnd[allid].r)^2.0 * (bkgnd[allid].teff/5777.0)^4. / (planet_a)^2. ; incident flux wrt sun-earth value
    ; Equilibrium Temp.
    planet_teq = (bkgnd[allid].teff)*sqrt(bkgnd[allid].r/(2.*planet_a*AU_IN_RSUN))
    ; Impact parameter
    planet_b = (planet_a*AU_IN_RSUN / bkgnd[allid].r) * cosi[allid]; assumes circular orbit
    ; Stellar radius in AU
    min_a = 2.*star[allid].r/AU_IN_RSUN

    ; Weiss & Marcy 2014:
    planet_m = dblarr(nplanets)
    lomass = where(planet_rad lt 1.5)
    if (lomass[0] ne -1) then planet_m[lomass] = 0.440*(planet_rad[lomass])^3. + 0.614*(planet_rad[lomass])^4.
    himass = where(planet_rad ge 1.5)
    if (himass[0] ne -1) then planet_m[himass] = 2.69*(planet_rad[himass])^0.93
    ; RV amplitude
    planet_k = RV_AMP*planet_per^(-1./3.) * planet_m * $
        sqrt(1.0-cosi[allid]^2.) * (bkgnd[allid].m)^(-2./3.)

    ; Work out transit properties
    dep1 = (planet_rad*REARTH_IN_RSUN / bkgnd[allid].r )^2.0
    tra = where((abs(planet_b) lt 1.0) and (planet_a gt min_a) and (dep1 gt min_depth))
    ntra = 0
    if (tra[0] ne -1)  then begin
      traid = planet_hid[tra]
      ntra = n_elements(tra)
      planet_multi = intarr(ntra)
      planet_tmulti = intarr(ntra)
      planet_pr = fltarr(ntra)
      for t=0,(ntra-1) do begin
        pos = where(planet_hid eq traid[t])
        tos = where(traid eq traid[t])
        if (n_elements(pos) gt 1) then begin
          planet_multi[t] = n_elements(pos) - 1  ; Number of other planets in the system
          planet_tmulti[t] = n_elements(tos) - 1 ; Number of other transiting planets
          pp = [tra[t]]
          opos = setdifference(pos, pp)
          if (opos[0] ne -1) then begin
            planet_pr[t] = min(planet_per[opos]/planet_per[tra[t]] > planet_per[opos]^(-1)*planet_per[tra[t]])
          endif
        endif
      endfor
      planet_eclip = replicate({eclipstruct}, ntra)
      r2 = planet_rad[tra]*REARTH_IN_RSUN
      dep2 = phot_ratio_planet(bkgnd[traid].teff, planet_teq[tra], bkgnd[traid].mag.t, bkgnd[traid].mag.dm, r2, ph_p, tband)
      planet_eclip.class=5
      planet_eclip.m1 = bkgnd[traid].m
      planet_eclip.m2 = planet_m[tra]/MSUN_IN_MEARTH
      planet_eclip.k = planet_k[tra]
      planet_eclip.r1 = bkgnd[traid].r
      planet_eclip.r2 = r2
      planet_eclip.teff1 = bkgnd[traid].teff
      planet_eclip.teff2 = planet_teq[tra]
      planet_eclip.a = planet_a[tra]
      planet_eclip.s = planet_s[tra]
      planet_eclip.p = planet_per[tra]
      planet_eclip.b = planet_b[tra]
      planet_eclip.cosi = cosi[traid]
      planet_eclip.tsys = bkgnd[traid].mag.tsys
      planet_eclip.kpsys = bkgnd[traid].mag.kpsys
      planet_eclip.icsys = bkgnd[traid].mag.icsys
      planet_eclip.jsys = bkgnd[traid].mag.jsys
      planet_eclip.hostid = planet_hid[tra]
      planet_eclip.dep1 = (planet_eclip.r2 / planet_eclip.r1 )^2.0 * 10.^(-0.4*(bkgnd[traid].mag.t-bkgnd[traid].mag.tsys))
      toodeep = where(planet_eclip.dep1 gt 1.0)
      if (toodeep[0] ne -1) then planet_eclip[toodeep].dep1 = 1.0
      planet_eclip.dep2 = dep2 * 10.^(-0.4*(bkgnd[traid].mag.t-bkgnd[traid].mag.tsys))
      planet_eclip.mult = planet_multi
      planet_eclip.tmult = planet_tmulti
      planet_eclip.pr   = planet_pr
      planet_eclip.dur1 = planet_eclip.r1 * planet_eclip.p * sqrt(1.-(planet_eclip.b)^2.) / (!PI*planet_eclip.a*AU_IN_RSUN)
      planet_eclip.dur2 = planet_eclip.r1 * planet_eclip.p * sqrt(1.-(planet_eclip.b)^2.) / (!PI*planet_eclip.a*AU_IN_RSUN)
      planet_eclip.gress1 = planet_eclip.r2 * planet_eclip.p / $
        (!PI*planet_eclip.a*AU_IN_RSUN*sqrt(1.-(planet_eclip.b)^2.))
      planet_eclip.gress2 = planet_eclip.r2 * planet_eclip.p / $
        (!PI*planet_eclip.a*AU_IN_RSUN*sqrt(1.-(planet_eclip.b)^2.))



      if (ii gt 0) then estruct=struct_append(estruct, planet_eclip) $
      else estruct = planet_eclip
      nebtot = nebtot + ntra
    end ; if tra
  end ; mult loop

  keep = 0L
  if (nebtot gt 0) then begin
    randomp, sep, 1., nebtot, range_x=[0., radmax]
    print, 'Created ', nebtot, ' transiting planets in all.'
    keep = lonarr(nebtot)
    ; Blend with target stars
    for ii=0, nebtot-1 do begin
      ; draw star from random
      r = randomu(seed, nstars)
      ; How many of these are primaries and fall within radmax?
      if (keyword_set(ps_only)) then gd = where((r lt thresh) and (star.sec ne 1) and (star.ffi ne 1)) $
      else gd = where((r lt thresh) and (star.sec ne 1))
      ;print, 'found goods'
      if (gd[0] ne -1) then begin
        ; ID the brightest star
        brightt = min(star[gd].mag.t, ind)
        ;print, 'done with min'
        estruct[ii].hostid = gd[ind]
        estruct[ii].sep= sep[ii] ; in pixels
        keep[ii] = 1 ; set the keep flag
      end
    end

    if (total(keep) gt 0) then begin
      print, 'Keeping ', round(total(keep)), ' blended planets'
      keepers = where(keep)  
      estruct = estruct[keepers] 
    ;bkteff = estruct.teff1
    ;bkmagt = tsys[keepers]
    ;starteff = star[gd].teff
    ;starmagt = star[gd].mag.t
      
    ;phr1 = phot_ratio(bkteff, starteff, bkmagt, starmagt, ph_p) ; Flux ratios
    ;phr2 = 1.0-phr1
    ;phr = phr1/phr2
      
   ; bkrad  = r[gd]
    endif
  endif
return, total(keep)
end

function add_bebs, star, bkgnd, estruct, frac, ph_p, mult, $ ;input
  	aspix=aspix, radmax=radmax, ps_only=ps_only

  keep = 0

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
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

  pris = where(bkgnd.pri eq 1)
  secs = where(bkgnd.sec eq 1)

  for ii=0,mult-1 do begin
    npri = n_elements(pris)
    r1 = bkgnd[pris].r
    r2 = bkgnd[bkgnd[pris].companion.ind].r
    m1 = bkgnd[pris].m
    m2 = bkgnd[bkgnd[pris].companion.ind].m
    teff1 = bkgnd[pris].teff
    teff2 = bkgnd[bkgnd[pris].companion.ind].teff
    tmag1 = bkgnd[pris].mag.t
    tsys = bkgnd[pris].mag.tsys
    kpsys = bkgnd[pris].mag.kpsys
    icsys = bkgnd[pris].mag.icsys
    jsys = bkgnd[pris].mag.jsys
    tmag2 = bkgnd[bkgnd[pris].companion.ind].mag.t
    ars = bkgnd[pris].companion.a*AU_IN_RSUN
    a   = bkgnd[pris].companion.a
    p = bkgnd[pris].companion.p
    ecc = bkgnd[pris].companion.ecc
    f = bkgnd[pris].companion.f

    ; Re-randomize the inclination and w
    cosi = -1.0 + 2.0*randomu(seed, npri)
    w = 2.0*!dpi*randomu(seed, npri)
    b1 = ars*cosi/r1*(1.0-ecc^2.)/(1.0+ecc*sin(w))
    b2 = ars*cosi/r2*(1.0-ecc^2.)/(1.0-ecc*sin(w))

    roche1 = (3.*m1/m2)^(1./3.)*r2
    roche2 = (3.*m2/m1)^(1./3.)*r1
    ; Where are the (non-contact) eclipsing systems? 
;    bin_ecl = where((r1*abs(b1) lt (r1+r2)) and (ars gt roche))
    bin_ecl = where(((r1*abs(b1) lt (r1+r2)) or (r2*abs(b2) lt (r1+r2))) and $
                    (ars gt (r1 > r2)) and (ars gt (roche1 > roche2)))
 
    if (bin_ecl[0] ne -1) then begin
      neb = n_elements(bin_ecl)
      pdur14 = dblarr(neb)
      pdur23 = dblarr(neb)
      sdur14 = dblarr(neb) 
      sdur23 = dblarr(neb) 
      gress1 = dblarr(neb) 
      gress2 = dblarr(neb) 
      dep1 = dblarr(neb) 
      dep2 = dblarr(neb) 
      dur1 = dblarr(neb) 
      dur2 = dblarr(neb) 
      a1 = dblarr(neb)
      a2 = dblarr(neb)
      r1 = r1[bin_ecl]
      m1 = m1[bin_ecl]
      teff1 = teff1[bin_ecl]
      tmag1 = tmag1[bin_ecl]
      r2 = r2[bin_ecl]
      m2 = m2[bin_ecl]
      teff2 = teff2[bin_ecl]
      tmag2 = tmag2[bin_ecl]
      tsys = tsys[bin_ecl]
      kpsys = kpsys[bin_ecl]
      icsys = icsys[bin_ecl]
      jsys = jsys[bin_ecl]
      a = a[bin_ecl]
      ars = ars[bin_ecl]
      p = p[bin_ecl]
      b1 = b1[bin_ecl]
      b2 = b2[bin_ecl]
      cosi = cosi[bin_ecl]
      w = w[bin_ecl]
      ecc = ecc[bin_ecl]
      f = f[bin_ecl]

      ; Totals and grazes
      te1 = where((r1*abs(b1) lt abs(r1-r2)) and (ars gt (r1+r2)))
      te2 = where((r2*abs(b2) lt abs(r2-r1)) and (ars gt (r1+r2)))
      ge1 = where((r1*abs(b1) ge abs(r1-r2)) $
         and (r1*abs(b1) lt (r1+r2)) and (ars gt (r1+r2)))
      ge2 = where((r2*abs(b2) ge abs(r2-r1)) $
        and (r2*abs(b2) lt (r1+r2)) and (ars gt (r1+r2)))

      ; Same for all eclipse types
      sini = sqrt(1. - cosi^2.)
      pdur14 = (p/!dpi)*asin(r1/ars*sqrt((1.+r2/r1)^2.-b1^2.)/sini) ; Eqn 14 of Winn
      sdur14 = (p/!dpi)*asin(r2/ars*sqrt((1.+r1/r2)^2.-b2^2.)/sini)
      ; change for eccentricity
      pdur14 = pdur14*sqrt(1.0-ecc^2.)/(1.0 + ecc*sin(w)) ; Eqn 16 of Winn
      sdur14 = sdur14*sqrt(1.0-ecc^2.)/(1.0 - ecc*sin(w))
      ; flux ratios
      phr1 = phot_ratio(teff1, teff2, tmag1, tmag2, ph_p) ; Flux ratios
      phr2 = 1.0-phr1
      
      ; Only for total eclipses
      if (te1[0] ne -1) then begin
        pdur23[te1] = (p[te1]/!dpi)*asin(r1[te1]/ars[te1]*$ ; Eqn 15 of Winn
          sqrt((1.-r2[te1]/r1[te1])^2.-b1[te1]^2)/sini[te1]) ; 
        pdur23[te1] = pdur23[te1]*sqrt(1.0-ecc[te1]^2.)/(1.0 + ecc[te1]*sin(w[te1]))
        dur1[te1] = (pdur14[te1] + pdur23[te1])/2. ; Trapezoidal area
        a1[te1] = (r2[te1]/r1[te1])^2. < 1.0
      end
      if (te2[0] ne -1) then begin
        sdur23[te2] = (p[te2]/!dpi)*asin(r2[te2]/ars[te2]*$ ; Eqn 15 of Winn
          sqrt((1.-r1[te2]/r2[te2])^2.-b2[te2]^2)/sini[te2]) ; 
        sdur23[te2] = sdur23[te2]*sqrt(1.0-ecc[te2]^2.)/(1.0 - ecc[te2]*sin(w[te2]))
        dur2[te2] = (sdur14[te2] + sdur23[te2])/2.
        a2[te2] = (r1[te2]/r2[te2])^2. < 1.0
      end
      if (ge1[0] ne -1) then begin
        dur1[ge1] = pdur14[ge1]/2. ; FWHM of "V" shaped eclipse
        delt = r1[ge1]*abs(b1[ge1]) ; All defined on ph. 24-26 of Kopal (1979)
        ph1  = acos((delt^2. + r1[ge1]^2. - r2[ge1]^2.)/(2.*r1[ge1]*delt))
        ph2  = acos((delt^2. - r1[ge1]^2. + r2[ge1]^2.)/(2.*r2[ge1]*delt))
        da1  = r1[ge1]^2.*(ph1-0.5*sin(2*ph1))
        da2  = r2[ge1]^2.*(ph2-0.5*sin(2*ph2))
        a1[ge1] = (da1+da2)/(!dpi*r1[ge1]^2.) < 1.0
      end
      if (ge2[0] ne -1) then begin
        dur2[ge2] = sdur14[ge2]/2.
        delt = r2[ge2]*abs(b2[ge2]) ; All defined on ph. 24-26 of Kopal (1979)
        ph1  = acos((delt^2. + r1[ge2]^2. - r2[ge2]^2.)/(2.*r1[ge2]*delt))
        ph2  = acos((delt^2. - r1[ge2]^2. + r2[ge2]^2.)/(2.*r2[ge2]*delt))
        da1  = r1[ge2]^2.*(ph1-0.5*sin(2*ph1))
        da2  = r2[ge2]^2.*(ph2-0.5*sin(2*ph2))
        a2[ge2] = (da1+da2)/(!dpi*r2[ge2]^2.) < 1.0
      end


      dep1 = phr1*a1
     ; toodeep = where(a1 gt 1.0)
      ;if (toodeep[0] ne -1) then dep1[toodeep] = phr1[toodeep]
      dep2 = phr2*a2
      ;toodeep = where(a2 gt 1.0)
      ;if (toodeep[0] ne -1) then dep2[toodeep] = phr2[toodeep]
      eclipse_hid = pris[bin_ecl]
      
      ; RV amplitude
      ;eclipse_k = RV_AMP*eclipse_p^(-1./3.) * eclipse_m * $ 
      ;	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 
      gress1 = (pdur14-pdur23)/2.0
      gress2 = (sdur14-sdur23)/2.0
      ; Work out transit properties
      eclip = replicate({eclipstruct}, neb)
      eclip.class=3
      eclip.m1 = m1
      eclip.m2 = m2
      eclip.k = RV_AMP*2.0*!dpi*a*m2 * $ 
	sqrt(1.0-cosi^2.)/(p*m1)
      eclip.r1 = r1
      eclip.r2 = r2
      eclip.teff1 = teff1
      eclip.teff2 = teff2
      eclip.a = a
      ;planet_eclip.s = s
      eclip.p = p
      eclip.b = b1
      eclip.cosi = cosi
      eclip.w = w
      eclip.ecc = ecc
      eclip.f = f
      eclip.dep1 = dep1
      eclip.dep2 = dep2
      eclip.dur1 = dur1
      eclip.dur2 = dur2
      eclip.gress1 = gress1
      eclip.gress2 = gress2
      eclip.tsys = tsys
      eclip.kpsys = kpsys
      eclip.icsys = icsys
      eclip.jsys = jsys
      ;print, 'Created ', neb, ' eclipsing binaries out of ', n_elements(pris), ' primaries.'
      if (ii gt 0) then estruct=struct_append(estruct, eclip) $
	else estruct = eclip
      nebtot = nebtot + neb
    end ; if ebs
  end ; mult loop
  keep = 0L
  if (nebtot gt 0) then begin
    randomp, sep, 1., nebtot, range_x=[0., radmax]
    print, 'Created ', nebtot, ' eclipsing binaries in all.'
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
      print, 'Keeping ', round(total(keep)), ' blended binaries'
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

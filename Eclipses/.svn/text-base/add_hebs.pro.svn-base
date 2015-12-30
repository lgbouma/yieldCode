function add_hebs, star, eclip, $
        frac, ph_p, dartstruct, tefftic, $ ;input
  	aspix=aspix, radmax=radmax, ps_only=ps_only
  
  neb = 0.
  nspl = 0.

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

  if keyword_set(ps_only) then spl = where(star.spl eq 1 and star.ffi ne 1) $
  else spl = where(star.spl eq 1)
  if (spl[0] ne -1) then begin
    hostid = star[spl].companion.ind
    hostp = star[spl].companion.p
    hostm = star[hostid].m
    hostsep = star[spl].companion.sep/aspix
    nspl = n_elements(spl)
    mtot = star[spl].m
    age = star[spl].age
    feh = star[spl].feh
    ;q = randomu(seed, nspl)*0.9 + 0.1
    q = randomn(seed, nspl)*0.2 + 0.5
    q[where((randomu(seed, nspl) lt 0.20) or (q lt 0.1) or (q gt 1.0))] = 1.0
    m1 = mtot/(1.0+q)
    m2 = mtot*q/(1.0+q)
    pu = randomu(seed, nspl)
    p = 0.2*hostp*10.^(-2.*pu)
    a = (m1+m2)^(1./3.)*(p/365.25)^(2./3.)
    ; Fill in stellar properties
    dm = star[spl].mag.dm
    av = star[spl].mag.av
    dartmouth_interpolate, dartstruct, m1, age, feh, $
        rad=r1, ic=ic1, teff=teff1, v=v1, rc=rc1, j=j1, h=h1, ks=ks1
    dartmouth_interpolate, dartstruct, m2, age, feh, $
        rad=r2, ic=ic2, teff=teff2, v=v2, rc=rc2, j=j2, h=h2, ks=ks2
    tmag1 = interpol(tefftic[*,1], tefftic[*,0], teff1) + ic1 + dm + 0.479*av
    tmag2 = interpol(tefftic[*,1], tefftic[*,0], teff2) + ic1 + dm + 0.479*av
    ;tmag2 = tmag1
    ;r2 = r1
    ;teff2 = teff1
    vr1 = v1-rc1
    vr2 = v2-rc2
    ; Jordi (2008)
    newr1 = rc1
    newr2 = rc2
    ; Re-create Sloan r
    vr11 = where(vr1 le 0.93)
    vr21 = where(vr2 le 0.93)
    if (vr11[0] ne -1) then newr1[vr11] = rc1[vr11] + 0.275*vr1[vr11] + 0.086
    if (vr21[0] ne -1) then newr2[vr21] = rc2[vr21] + 0.275*vr2[vr21] + 0.086
    vr12 = where(vr1 gt 0.93)
    vr22 = where(vr2 gt 0.93)
    if (vr12[0] ne -1) then newr1[vr12] = rc1[vr12] +  0.71*vr1[vr12] - 0.31
    if (vr22[0] ne -1) then newr2[vr22] = rc2[vr22] +  0.71*vr2[vr22] - 0.31
    ; Re-create Sloan i
    newi1 = ic1 + 0.251*(rc1-ic1) + 0.325
    newi2 = ic2 + 0.251*(rc2-ic2) + 0.325
    ; Re-create Kp
    kp1 = newi1
    kp2 = newi2
    ri1 = newr1-newi1
    ri2 = newr2-newi2
    ri11 = where(ri1 le 0.673)
    ri21 = where(ri2 le 0.673)
    if (ri11[0] ne -1) then kp1[ri11] = 0.65*newr1[ri11] + 0.36*newi1[ri11]
    if (ri21[0] ne -1) then kp2[ri21] = 0.65*newr1[ri21] + 0.36*newi1[ri21]
    ri12 = where(ri1 gt 0.673)
    ri22 = where(ri2 gt 0.673)
    if (ri12[0] ne -1) then kp1[ri12] = 1.2*newr1[ri12] + 0.2*newi1[ri12]
    if (ri22[0] ne -1) then kp2[ri22] = 1.2*newr2[ri22] + 0.2*newi2[ri22]

    tsys = -2.5*alog10(10.^(-0.4*tmag1) + 10.^(-0.4*tmag2)) 
    icsys = -2.5*alog10(10.^(-0.4*ic1) + 10.^(-0.4*ic2)) + dm + 0.479*av
    kpsys = -2.5*alog10(10.^(-0.4*kp1) + 10.^(-0.4*kp2)) + dm + 0.615*av
    jsys = -2.5*alog10(10.^(-0.4*j1) + 10.^(-0.4*j2)) + dm + 0.282*av

    ; over-write sys mags
    star[spl].mag.tsys = tsys
    star[spl].mag.jsys = jsys
    star[spl].mag.icsys = icsys
    star[spl].mag.kpsys = kpsys

    ; Separation in Rsun
    ars = a*AU_IN_RSUN

    ; eccentricity
    ecc = randomu(seed, nspl)*(atan((alog10(p)-1.5)*2.)+!dpi/2.)/!dpi
    ecc = ecc > 0.0
    ecc = ecc < 1.0
    ; Argument of periastron: 
    w = 2.*!dpi*randomu(seed, nspl)

    ; Mean anomaly is linear in time
    ma = 2.*!dpi*randomu(seed, nspl)

    ; Find eccentric anomaly
    e = fltarr(nspl)
    for jj=0, nspl-1 do begin
        e[jj] = keplereq(ma[jj],ecc[jj])
    end
      
    ; Find true anomaly
    f = 2.0*atan(sqrt(1.0+ecc)*sin(e/2.0), sqrt(1.0-ecc)*cos(e/2.0))

    ; Re-randomize the inclination. These are circular orbits
    cosi = -1.0 + 2.0*randomu(seed, nspl)
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

      hostid = hostid[bin_ecl]
      hostsep = hostsep[bin_ecl]

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
      phr1 = phot_ratio(teff1, teff2, tmag1, tmag2, ph_p)
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
      ;toodeep = where(a1 gt 1.0)
      ;if (toodeep[0] ne -1) then dep1[toodeep] = phr1[toodeep]
      dep2 = phr2*a2
      ;toodeep = where(a2 gt 1.0)
      ;if (toodeep[0] ne -1) then dep2[toodeep] = phr2[toodeep]
      
      ; RV amplitude
      ;eclipse_k = RV_AMP*eclipse_p^(-1./3.) * eclipse_m * $ 
      ;	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 
      gress1 = (pdur14-pdur23)/2.0
      gress2 = (sdur14-sdur23)/2.0
      ; Work out transit properties
      eclip = replicate({eclipstruct}, neb)
      eclip.class=4
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
      eclip.hostid = hostid
      eclip.sep = hostsep ; in pixels
      print, 'Created ', neb, ' hierarchical eclipsing binaries out of ', nspl
    end ; if ebs
  end ; spl loop
  return, neb
end

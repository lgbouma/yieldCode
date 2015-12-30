PRO fits2sav, fname, dartmouth, tefftic, jlfr=jlfr, nstar=nstar, $
	icmax=icmax, dmax=dmax, homies=homies, dbl=dbl, kmin=kmin, tmin=tmin, $
        dartcor=dartcor, suffix=suffix
  dat = mrdfits(fname, 0, h, /SILENT)
  print, fname

  if keyword_set(jlfr) then begin
    sz = size(jlfr)
    nmj = sz[1]
    lf_mj = jlfr[*,0]
    lf_dm = lf_mj[1]-lf_mj[0]
    lf_r  = jlfr[*,1]
    ;print, "Applying LF from J=",min(lf_mj)-lf_dm/2.," to J=",max(lf_mj)+lf_dm/2.," with a median change of ",median(lf_r)
  endif
       
  dartmin = 0.14
  dartmax = 0.78
  ; Abs. mag ranges for binary properties
  mp_min = [0.0, 0.1, 0.6, 0.8, 1.0, 1.4]  
  mp_max = [     0.1, 0.6, 0.8, 1.0, 1.4, 99.]  
  nmp = n_elements(mp_min)
  mf = [0.22, 0.26, 0.34, 0.41, 0.50, 0.75]
  vtf = [100., 4500., 5000., 5500., 6000., 1D6]
  var = [0.03, 0.02, 0.01, 0.01, 0.01] ; Basri 2013
  abar = [4.5, 5.3, 20.,  45.,  45., 350]
  psig = [0.5, 1.3, 2.0,  2.3,  2.3, 3.0]
  qgam = [4.0, 0.4, 0.35, 0.3, 0.3, -0.5]
  homf = [0.0, 3.9, 3.8, 3.7, 3.7, 3.7]
  qnorm = 0.9*(qgam+1.0)*mf/(1.0-0.1^(qgam+1.0))
;  print, qnorm
;  readcol, fname, gc, logA, z, mini, logL, logT, logG, dm, av, $
;	comp, bol, t, j, h, ks, kp, g, r, i, z, dd, mnow 
  comp = dat[*,9]
  t    = dat[*,11] 
  ks   = dat[*,14] 
  if (keyword_set(kmin)) then gd = where(ks gt kmin and comp eq 1) $
  else if (keyword_set(tmin)) then gd = where(t gt tmin and comp eq 1) $
  else gd = where(comp eq 1)
  gd = [gd, gd+1] ; get the secondaries

  gc   = dat[gd,0]
  age  = 10.^(dat[gd,1]-9.0)
  feh  = alog10(dat[gd,2]/19.0)
  mini = dat[gd,3] 
  logL = dat[gd,4] 
  teff = 10.^(dat[gd,5])
  logG = dat[gd,6] 
  dm   = dat[gd,7] 
  av   = dat[gd,8]
  comp = dat[gd,9]
  bol  = dat[gd,10]
  t    = dat[gd,11] 
  j    = dat[gd,12] 
  h    = dat[gd,13] 
  ks   = dat[gd,14] 
  kp   = dat[gd,15] 
  g    = dat[gd,16] 
  r    = dat[gd,17] 
  i    = dat[gd,18] 
  z    = dat[gd,19] 
  dd   = dat[gd,20] 
  mnow = dat[gd,21] 

  ndat = n_elements(mini)
  
 ; Derived quantities
  m = mnow
  rad = sqrt(m)/sqrt(10.^logG/27542.3)
  gr = g-r
  gr1 = where(gr le 1.3811)
  gr2 = where(gr gt 1.3811)
  v = g
  if (gr1[0] ne -1) then v[gr1] = 0.4534*g[gr1] + 0.5466*r[gr1] - 0.03794
  if (gr2[0] ne -1) then v[gr2] = 0.1814*g[gr2] + 0.8186*r[gr2] - 0.32922
  ;v = g - 0.5784*(g - r) - 0.0038 ;Lupton 2005
  ic = i - 0.16*(r - z) - 0.3842 
  ;ic = i - 0.3780*(i - z) - 0.3974 ; Lupton 2005
  ;icsys = ic
  ;jsys = j
  ;ksys = ks
  ;tsys = t
  mv = v - av - dm
  mic = ic - 0.479*av - dm
  mj = j - 0.282*av - dm
  
  dart = intarr(ndat)
  darts = where(m gt dartmin and m lt dartmax)
  if ((darts[0] ne -1) and keyword_set(dartcor)) then begin
    newdm = dm[darts]
    newav = av[darts]
    dartmouth_interpolate, dartmouth, m[darts], age[darts], feh[darts], $
	rad=newrad, ic=newic, teff=newteff, v=newv, rc=newrc, j=newj, h=newh, ks=newks
    newt = interpol(tefftic[*,1], tefftic[*,0], newteff) + newic
    dart[darts] = 1
    rad[darts] = newrad
    mv[darts] = newv
    mic[darts] = newic
    mj[darts] = newj
    ic[darts] = newic + newdm + 0.479*newav
    t[darts] = newt + newdm + 0.479*newav
    teff[darts] = newteff
    v[darts] = newv + newdm + newav
    vr = newv-newrc
    ; Jordi (2008)
    newr = newrc
    ; Re-create Sloan r
    vr1 = where(vr le 0.93)
    if (vr1[0] ne -1) then newr[vr1] = newrc[vr1] + 0.275*vr[vr1] + 0.086
    vr2 = where(vr gt 0.93)
    if (vr2[0] ne -1) then newr[vr2] = newrc[vr2] +  0.71*vr[vr2] - 0.31
    ; Re-create Sloan i
    newi = newic + 0.251*(newrc-newic) + 0.325
    ; Re-create Kp
    newkp = newi
    ri = newr-newi
    ri1 = where(ri le 0.673)
    if (ri1[0] ne -1) then newkp[ri1] = 0.65*newr[ri1] + 0.36*newi[ri1]
    ri2 = where(ri gt 0.673)
    if (ri2[0] ne -1) then newkp[ri2] = 1.2*newr[ri2] + 0.2*newi[ri2]
    r[darts] = newr + newdm + 0.751*newav
    i[darts] = newi + newdm + 0.479*newav
    kp[darts] = newkp + newdm + 0.615*newav
    j[darts] = newj + newdm + 0.282*newav
    h[darts] = newh + newdm + 0.190*newav
    ks[darts] = newks + newdm + 0.114*newav
    ;icsys[darts] = ic[darts]
    ;jsys[darts] = j[darts]
    ;ksys[darts] = ks[darts]
    ;tsys[darts] = t[darts]
  endif

  pris = where(comp eq 1)
  secs = where(comp eq 2)
  
  if (keyword_set(dbl)) then begin
     pris = [pris,pris]
     secs = [secs,secs]
  end
 
  ; Enforce luminosity function on primaries + singles
  if (keyword_set(jlfr)) then begin
    rm = []
    ad = []
    for ii=0,nmj-1 do begin
      thismj = where((mj[pris] gt lf_mj[ii]-lf_dm/2.) and (mj[pris] le lf_mj[ii]+lf_dm/2.) and (gc[pris] le 2))
      ; If there are no stars in this mag bin, then do nothing
      if (thismj[0] ne -1) then begin
        ; Select stars for removal if over-abundant
        if (lf_r[ii] lt 1.0) then begin
          rmp = where(randomu(seed,n_elements(thismj)) gt (lf_r[ii]))
          if (rmp[0] ne -1) then rm = [rm,thismj[rmp]]
        ; select stars for duplication if under-abundant
	endif else if (lf_r[ii] gt 1.0) then begin
          if (lf_r[ii] ge 2.0) then begin
            addp = lindgen(n_elements(thismj)*floor(lf_r[ii])) mod n_elements(thismj)
            ad = [ad, thismj[addp]]
          end
          addp = where(randomu(seed, n_elements(thismj)) lt (lf_r[ii] mod 1)) 
          if (addp[0] ne -1) then ad = [ad,thismj[addp]]
        end
      endif
    endfor ; M_J loop
    ;if (n_elements(rm) gt 0) then begin
    print, 'Removing ', n_elements(rm), ' stars and duplicating ', n_elements(ad), ' stars'
    pris = [pris, pris[ad]]
    secs = [secs, secs[ad]]
    remove, rm, pris ; also remove the indices
    remove, rm, secs ; also remove the indices
  endif ; LF enforcement

  npri = n_elements(pris)
  star = replicate({starstruct}, npri)
  
  m1 = mini[pris]
  m2 = mini[secs]
  q = m2/m1

  ; error checking on feh
  ;mz = where(z lt 0)
  ;if (mz[0] ne -1) then z[mz] = 19.0
  
  star.age = age[pris]
  star.feh = feh[pris]
  star.mini = mini[pris]
  star.teff = teff[pris]
  star.logg = logG[pris]
  star.mag.dm = dm[pris]
  star.mag.av = av[pris]
  star.mag.g = g[pris]
  star.mag.v = v[pris] 
  star.mag.r = r[pris]
  star.mag.i = i[pris]
  star.mag.ic = ic[pris]  
  star.mag.z = z[pris]
  star.mag.kp = kp[pris]
  star.mag.j = j[pris]
  star.mag.h = h[pris]
  star.mag.k = ks[pris]
  star.mag.t = t[pris]
  star.mag.mv = mv[pris] ;star.mag.v - star.mag.av - star.mag.dm
  star.mag.mic = mic[pris] ;star.mag.ic - 0.479*star.mag.av - star.mag.dm
  star.mag.mj = mj[pris] ;star.mag.j - 0.282*star.mag.av - star.mag.dm
  star.mag.icsys = ic[pris] ; star.mag.ic
  star.mag.kpsys = kp[pris] ; star.mag.ic
  star.mag.jsys = j[pris] ;star.mag.j
  star.mag.ksys = ks[pris] ;star.mag.j
  star.mag.tsys = t[pris] ;star.mag.t
  star.mag.mvsys = mv[pris] ;bin_star.mag.v - bin_star.mag.av - bin_star.mag.dm
  star.mag.micsys = mic[pris] ; bin_star.mag.ic - 0.479*bin_star.mag.av - bin_star.mag.dm
  star.mag.mjsys = mj[pris] ; bin_star.mag.j - 0.282*bin_star.mag.av - bin_star.mag.dm
  star.m = m[pris]
  star.r = rad[pris]
  star.gc = gc[pris]

 ;  darts = where(star.m gt dartmin and star.m lt dartmax)
 ;  if (darts[0] ne -1) then begin
 ;    newdm = star[darts].mag.dm
 ;    newav = star[darts].mag.av
 ;   dartmouth_interpolate, dartmouth, star[darts].m, star[darts].age, star[darts].feh, $
 ;	rad=newrad, ic=newic, teff=newteff, v=newv, rc=newr, j=newj, h=newh, ks=newks
 ;    newt = interpol(tefftic[*,1], tefftic[*,0], newteff) + newic
 ;    star[darts].dart = 1
 ;   star[darts].r = newrad
 ;   star[darts].mag.mv = newv
 ;   star[darts].mag.mic = newic
 ;   star[darts].mag.mj = newj
 ;   star[darts].mag.ic = newic + newdm + 0.479*newav
 ;   star[darts].mag.t = newt + newdm + 0.479*newav
 ;   star[darts].teff = newteff
 ;   star[darts].mag.v = newv + newdm + newav
 ;   star[darts].mag.r = newr + newdm + 0.751*newav
 ;   star[darts].mag.j = newj + newdm + 0.282*newav
 ;   star[darts].mag.h = newh + newdm + 0.190*newav
 ;   star[darts].mag.k = newks + newdm + 0.114*newav
 ;   star[darts].mag.icsys = star[darts].mag.ic
 ;   star[darts].mag.jsys = star[darts].mag.j
 ;   star[darts].mag.tsys = star[darts].mag.t
 ;   ; don't need to over-write mass, age, feh. Maybe add back extinction? 
 ; end
 ;

  idx0 = long(npri)

  for ii=0,nmp-1 do begin
    sind = where((randomu(seed, npri) lt qnorm[ii]*q^qgam[ii]) and (star.m gt mp_min[ii]) and (star.m le mp_max[ii]))   
    if (sind[0] ne 0) then begin 
      nsec = n_elements(sind)
      bin_star = replicate({starstruct}, nsec)
      bin_star.age = age[secs[sind]] ;10.^(logA[secs[sind]]-9.0)
      bin_star.feh = feh[secs[sind]] ; alog10(z[secs[sind]]/19.0)
      bin_star.mini = mini[secs[sind]]
      bin_star.teff = teff[secs[sind]] ;10.^(logT[secs[sind]])
      bin_star.logg = logG[secs[sind]]
      bin_star.mag.dm = dm[secs[sind]]
      bin_star.mag.av = av[secs[sind]]
      bin_star.mag.g = g[secs[sind]]
      bin_star.mag.v = v[secs[sind]] ;g[secs[sind]] - 0.5784*(g[secs[sind]]-r[secs[sind]]) - 0.0038 ;Lupton 2005
      bin_star.mag.r = r[secs[sind]]
      bin_star.mag.i = i[secs[sind]]
      bin_star.mag.ic = ic[secs[sind]] ;i[secs[sind]] - 0.3780*(i[secs[sind]]-z[secs[sind]]) - 0.3974 ; Lupton 2005
      bin_star.mag.z = z[secs[sind]]
      bin_star.mag.kp = kp[secs[sind]]
      bin_star.mag.j = j[secs[sind]]
      bin_star.mag.h = h[secs[sind]]
      bin_star.mag.k = ks[secs[sind]]
      bin_star.mag.t = t[secs[sind]]
      bin_star.mag.mv = mv[secs[sind]] ;bin_star.mag.v - bin_star.mag.av - bin_star.mag.dm
      bin_star.mag.mic = mic[secs[sind]] ; bin_star.mag.ic - 0.479*bin_star.mag.av - bin_star.mag.dm
      bin_star.mag.mj = mj[secs[sind]] ; bin_star.mag.j - 0.282*bin_star.mag.av - bin_star.mag.dm
      bin_star.m = m[secs[sind]]
      bin_star.r = rad[secs[sind]] ;sqrt(bin_star.m)/sqrt(10.^bin_star.logg/27542.3)
      bin_star.gc = gc[secs[sind]] ;sqrt(bin_star.m)/sqrt(10.^bin_star.logg/27542.3)
   
;      darts = where(bin_star.m gt dartmin and bin_star.m lt dartmax)
;      if (darts[0] ne -1) then begin
;        newdm = bin_star[darts].mag.dm
;        newav = bin_star[darts].mag.av
;        dartmouth_interpolate, dartmouth, bin_star[darts].m, bin_star[darts].age, bin_star[darts].feh, $
;	  rad=newrad, ic=newic, teff=newteff, v=newv, rc=newr, j=newj, h=newh, ks=newks
 ;       newt = interpol(tefftic[*,1], tefftic[*,0], newteff) + newic
 ;       bin_star[darts].dart = 1
 ;       bin_star[darts].r = newrad
 ;       bin_star[darts].mag.mv = newv
 ;       bin_star[darts].mag.mic = newic
 ;       bin_star[darts].mag.mj = newj
 ;       bin_star[darts].mag.ic = newic + newdm + 0.479*newav
 ;       bin_star[darts].mag.t = newt + newdm + 0.479*newav
 ;       bin_star[darts].teff = newteff
 ;       bin_star[darts].mag.v = newv + newdm + newav
 ;       bin_star[darts].mag.r = newr + newdm + 0.751*newav
 ;       bin_star[darts].mag.j = newj + newdm + 0.282*newav
 ;       bin_star[darts].mag.h = newh + newdm + 0.190*newav
 ;       bin_star[darts].mag.k = newks + newdm + 0.114*newav
 ;       bin_star[darts].mag.icsys = bin_star[darts].mag.ic
 ;       bin_star[darts].mag.jsys = bin_star[darts].mag.j
 ;       bin_star[darts].mag.tsys = bin_star[darts].mag.t
 ;       ; don't need to over-write mass, age, feh. Maybe add back extinction? 
 ;     end
    
      bin_inds = lindgen(nsec)+idx0
      ;print, min(bin_inds), max(bin_inds)
      idx0 = idx0 + nsec
      ; Set the companionship flags
      star[sind].pri = 1
      bin_star.sec = 1
      ; Cross-reference the indices
      star[sind].companion.ind = bin_inds 
      bin_star.companion.ind = sind
      ; Cross-reference the indices
      star[sind].companion.m = bin_star.m 
      bin_star.companion.m = star[sind].m
     
      ; Add up the fluxes
      star[sind].mag.mvsys = -2.5*alog10(10.^(-0.4*bin_star.mag.mv) + 10.^(-0.4*star[sind].mag.mv))
      bin_star.mag.mvsys = star[sind].mag.mvsys
      star[sind].mag.micsys = -2.5*alog10(10.^(-0.4*bin_star.mag.mic) + 10.^(-0.4*star[sind].mag.mic))
      bin_star.mag.micsys = star[sind].mag.micsys
      star[sind].mag.mjsys =  -2.5*alog10(10.^(-0.4*bin_star.mag.mj) + 10.^(-0.4*star[sind].mag.mj))
      bin_star.mag.mjsys = star[sind].mag.mjsys
      star[sind].mag.icsys = -2.5*alog10(10.^(-0.4*bin_star.mag.ic) + 10.^(-0.4*star[sind].mag.ic))
      bin_star.mag.icsys = star[sind].mag.icsys
      star[sind].mag.kpsys = -2.5*alog10(10.^(-0.4*bin_star.mag.kp) + 10.^(-0.4*star[sind].mag.kp))
      bin_star.mag.kpsys = star[sind].mag.kpsys
      star[sind].mag.jsys =  -2.5*alog10(10.^(-0.4*bin_star.mag.j) + 10.^(-0.4*star[sind].mag.j))
      bin_star.mag.jsys = star[sind].mag.jsys
      star[sind].mag.ksys =  -2.5*alog10(10.^(-0.4*bin_star.mag.k) + 10.^(-0.4*star[sind].mag.k))
      bin_star.mag.ksys = star[sind].mag.ksys
      star[sind].mag.tsys =  -2.5*alog10(10.^(-0.4*bin_star.mag.t) + 10.^(-0.4*star[sind].mag.t))
      bin_star.mag.tsys = star[sind].mag.tsys
      
      ; Convert mean separation into mean period
      mtot = star[sind].m + star[sind].companion.m
      logpbar = alog10(365.25*abar[ii]^(1.5)*mtot^(-0.5))
      ;print, median(logpbar)
      
      ; Randomize the period
      logp = randomn(seed, nsec)*psig[ii] + logpbar
      bin_star.companion.p = 10.^logp
      
      ; Cross-reference p and a
      star[sind].companion.p = bin_star.companion.p
      bin_star.companion.a = mtot^(1./3.)*(star[sind].companion.p/365.25)^(2./3.)
      star[sind].companion.a = bin_star.companion.a
      
      ;  Inclination and phase of binary
      cosi = -1.0 + 2.0*randomu(seed, nsec)
      star[sind].companion.cosi = cosi
      bin_star.companion.cosi = cosi
      
      ; Randomize the eccentricity. Total chi-by-eye from Raghavan et al. 2010
      ecc = randomu(seed, nsec)*(atan((logp-1.5)*2.)+!dpi/2.)/!dpi
      ecc = ecc > 0.0
      ecc = ecc < 1.0
      star[sind].companion.ecc = ecc
      bin_star.companion.ecc = ecc
      
      ; Argument of periastron: need to know for angsep, but gets re-randomized for eclipses
      w = 2.*!dpi*randomu(seed, nsec) 
      star[sind].companion.w = w
      bin_star.companion.w = w
      
      ; Mean anomaly is linear in time
      ma = 2.*!dpi*randomu(seed, nsec) 

      ; Find eccentric anomaly
      e = fltarr(nsec)
      ;print, 'Entering keplereq loop'
      for jj=0, nsec-1 do begin
        e[jj] = keplereq(ma[jj],ecc[jj])
      end 
      ;print, 'Done with keplereq loop'

      ; Find true anomaly
      f = 2.0*atan(sqrt(1.0+ecc)*sin(e/2.0), sqrt(1.0-ecc)*cos(e/2.0))
      star[sind].companion.f = f
      bin_star.companion.f = f

      ; Angular separation on the sky
      aang = star[sind].companion.a/(10.*10.^(star[sind].mag.dm/5.)) ; AU/pc -> arcseconds!
      swfs = sin(w+f)^2.0
      sis  = 1.0-cosi^2. 
      angseps = aang*(1.0-ecc^2.)/(1.0+ecc*cos(f))*sqrt(1.0 - swfs*sis)
      star[sind].companion.sep = angseps
      bin_star.companion.sep = angseps
      
      ; Mark stars for triples and quadruples
      if (homf[ii] gt 0) then begin
        tqind = where((randomu(seed, nsec) lt homf[ii]^(-1.)))
        if (tqind[0] ne -1) then begin
          ntq = n_elements(tqind)
          ; For quadruples, BOTH pri and sec get split
          qind = where((randomu(seed, ntq) lt homf[ii]^(-1.)), complement=tind)
          if (qind[0] ne -1) then begin
            star[sind[tqind[qind]]].spl = 1
            bin_star[tqind[qind]].spl = 1
          endif
          ; For triples, pri OR sec gets split
          if (tind[0] ne -1) then begin
            nt = n_elements(tind)
            tpind = where((randomu(seed, nt) lt 0.5), complement=tsind)
            if (tpind[0] ne -1) then star[sind[tqind[tind[tpind]]]].spl = 1
            if (tsind[0] ne -1) then bin_star[tqind[tind[tsind]]].spl = 1
          endif
        endif
      endif 

      star = struct_append(star, bin_star)
      delvar, bin_star
    endif
  end

;  for ii=0, n_elements(vtf)-2 do begin
;    thisteff = where(star.teff ge vtf[ii] and star.teff lt vtf[ii+1])
;    if (thisteff[0] ne -1) then begin
;      u = randomu(seed, n_elements(thisteff))
;      star[thisteff].var = var[ii]*10.^u
;    end
;  end
  
  ; Output for making distance-limited or mag-limited catalogs
  if (keyword_set(dmax)) then begin
	near = where(star.mag.dm le dmax)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  end else if (keyword_set(icmax)) then begin
	near = where(star.mag.ic le icmax)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  end else if (keyword_set(homies)) then begin
	near = where(star.spl)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  end else nstar = n_elements(star)
  if (keyword_set(suffix)) then newfname = repstr(fname, '.fits', '_'+suffix+'.sav') else $
  	newfname = repstr(fname, 'fits', 'sav')
  save, star, filen=newfname

END

pro bin_stars, outfile=outfile, verbose=verbose, dst=dst, csr=csr, det=det, per=per, mag=mag, bin=bin

  if (keyword_set(verbose)) then v=1 else v=0
  if (keyword_set(outfile)) then outfile=outfile else outfile='s.sav'
  if (keyword_set(csr)) then begin
   readcol, '../PhotonFluxes/stars_table.txt', $
           comment='#', f='D,D,D,D,D,D,D,D,D', $
           vmag, m, r, teff, mbol, imag, kmag, ph2, number_density
           dr = dblarr(n_elements(r))
           dr[1:n_elements(r)-2] = (r[0:n_elements(r)-2]-r[2:n_elements(r)-1])/2.0
           dr[0] = r[0]-r[1]
           dr[n_elements(r)-1] = r[n_elements(r)-2]-r[n_elements(r)-1]
	   phi = number_density 
           dndr = (4./3.)*!dpi*10.^3*0.1*number_density/dr ; at 10 pc
  endif else begin
    restore, filen='../Stars/star_properties_vijk.sav'
    vol10 = (4./3.)*!PI*(10.0^3.)
    dr = r[1]-r[0] + dblarr(n_elements(r))
    dt = deriv(teff)
    phi = dndr*dr/0.1/vol10       ; now phi is in stars per cubic parsec per radius bin
  endelse
  
; what kind of star catalog?
  if (keyword_set(dst)) then begin
    dmax = float(dst)+0.0*r ;72.0 + 128.*r;
  endif else if (keyword_set(det)) then begin
    AU_IN_RSUN = 215.093990942
    REARTH_IN_RSUN = 0.0091705248 
    geom_area = 70 ;74.6
    a = m^(1./3.) * (float(per)/365.25)^(2./3.)   ; in AU
    dur = r * float(per) / (!DPI*a*AU_IN_RSUN)
    exptime = 2.*dur*24.*3600
    dep = (REARTH_IN_RSUN * float(det) / r)^2.0
    sig = dep/7.0
    rn = 10.*sqrt(4.0*exptime/2.0)
    minphot = (1.+sqrt(1.+4.*sig^2.*rn^2.))/(2.*sig^2.)
    megaph_s_cm2_0mag = 1.67
    ilim = -2.5*alog10(minphot/(megaph_s_cm2_0mag * 1D6 * geom_area * exptime))
    if keyword_set(mag) then begin
 	ilim[where((r le 0.5) and (ilim lt 12))] = 13.0
        ilim[where((r gt 0.5) and (ilim lt 13))] = 12.0
;      ilim[where(ilim lt 12.5)] = 12.5      
    endif
    dmax = 10.*10.0^(0.2*(ilim-imag))
  endif else if (keyword_set(mag)) then begin
    ilim = 12.0 + 0.*r
    q = where(r le 0.5)
    ilim[q] = 13.0
    dmax = 10.*10.0^(0.2*(ilim-imag))
  endif
  print, dmax
  h = 300.0 + 0.0*r


  nstars = number_of_stars_expz(dmax, phi, h)
  nstars = round(nstars,/L64)
  totstars = total(nstars) + $
		total(round(nstars[where(m le 0.5)]*0.26)) + $
		total(round(nstars[where(m gt 0.5)]*0.44))
  print, 'R_star: ', min(r), max(r)
 ; print, 'I: ', min(ilim), max(ilim)
  print, 'd: ', min(dmax), max(dmax)
  print, 'Nstars: ', min(nstars), max(nstars), total(nstars), total(totstars)

  template_coord = {$
                   ra: 0.0, $   ; RA (deg)
		   dec: 0.0, $  ; DEC (deg)
                   elon: 0.0, $ ; Ecliptic Long. (deg)
                   elat: 0.0, $ ; Ecliptic Lat. (deg)
                   glon: 0.0, $ ; Galactic Long. (deg)
                   glat: 0.0, $ ; Galactic Lat. (deg)
                   d: 0.0, $    ; Distance (pc)
                   fov_r: 0.0 $ ; Field angle (pixels from axis)
              }
  
  template_mag = {$
                 mv: 0.0, $ ; absolute
                 mi: 0.0, $ ; absolute
                 mj: 0.0, $ ; absolute
                 mk: 0.0, $ ; absolute
                 v: 0.0, $ ; apparent
                 i: 0.0, $ ; apparent
                 j: 0.0, $  ; apparent
                 k: 0.0 $  ; apparent
                 }
 
  template_companion = { $	; Mutual properties
  		 ind: 0L, $ 	; Cross-reference to companion
		 sep: 0.0, $	; Angular separation on sky (arcsec)
		 a: 0.0, $	; Semimajor axis (AU)
		 p: 0.0, $	; Period (days)
  		 m: 0.0 $	; Mass (Solar)
		 }

  template_star = {$
                  r: 0.0, $     ; Radius (solar)
                  m: 0.0, $     ; Mass (solar)
                  teff: 0.0, $  ; Effective T (Kelvin)
 		  cosi: 0.0, $  ; cos inclination of planets
                  npointings: 0, $ ; Number of TESS pointings this star gets
                  coord: template_coord, $ 
                  mag: template_mag, $
                  companion: template_companion, $
		  npix: 0, $       ; Optimal number of pix in aperture
		  snr: 0.0, $      ; SNR per hour
 		  sat: 0, $        ; saturation flag
		  dil: 0.0, $      ; dilution factor (0+)
                  pri: 0, $	   ; Primary of binary?
		  sec: 0, $	   ; Secondary of binary? 
  		  ffi: 0 $	   ; Full-frame only (default is postage stamp)
		}

  star = replicate(template_star, total(nstars))

  for i=0,n_elements(r)-1 do begin
     
     if (nstars[i] lt 1) then continue
     tmp_star = replicate(template_star, nstars[i])
     tmp_star.r = r[i] + dr[i]*(randomu(seed, nstars[i]) - 0.5)
     tmp_star.m = m[i]
     tmp_star.teff = teff[i] + dt[i]*(randomu(seed, nstars[i])-0.5)
     tmp_star.mag.mv = vmag[i]
     tmp_star.mag.mi = imag[i]
     tmp_star.mag.mk = kmag[i]
     
     if (keyword_set(csr)) then tmp_star.mag.mj = kmag[i] else tmp_star.mag.mj = jmag[i]

     assign_xyz, nstars[i], dmax[i], h[i], x, y, z, d

     tmp_star.coord.d = d
     dm = 5.0*alog10(tmp_star.coord.d/10.0D0)
     tmp_star.mag.v = tmp_star.mag.mv + dm
     tmp_star.mag.i = tmp_star.mag.mi + dm
     tmp_star.mag.j = tmp_star.mag.mj + dm
     tmp_star.mag.k = tmp_star.mag.mk + dm

     tmp_star.coord.glat = 180./!PI*atan(z,sqrt(x^2. + y^2.))
     tmp_star.coord.glon = 180./!PI*atan(y,x)

     euler, tmp_star.coord.glon, tmp_star.coord.glat, elon, elat, 6
     q = where(elon gt 180.)
     if (q[0] ne -1) then elon[q] = elon[q]-360.0
     euler, tmp_star.coord.glon, tmp_star.coord.glat, ra, dec, 2
     tmp_star.coord.elon = elon & tmp_star.coord.elat = elat
     tmp_star.coord.ra = ra & tmp_star.coord.dec = dec

     idx = lindgen(nstars[i])
     if (i gt 0) then idx = idx + total(nstars[0:i-1])
     print, i, r[i], n_elements(tmp_star), n_elements(idx), max(idx)/float(total(nstars))
     star[idx] = tmp_star
     delvar, tmp_star

  endfor

  ; Binarity
  if (keyword_set(bin)) then begin
  ; Indices of the hosts
  bin1 = where((randomu(seed, total(nstars)) lt 0.26) and (star.m le 0.6))
  ; Generate mass of companion
  randomp, q, 0.4, n_elements(bin1), range_x=[0.0, 1.0], seed=seed
  newm = star[bin1].m*q
  gdbin = where(newm ge min(m)) ; Enforce brown dwarf desert
  bin1 = bin1[gdbin]
  ; Struct to hold the companions
  bin_star = replicate(template_star, n_elements(bin1))
  bin_star.m = newm[gdbin]
  ; Indices of the binaries
  bin_inds = lindgen(n_elements(bin1))+max(idx)
  idx = max(idx) + n_elements(bin1)
  print, 'bin1',  n_elements(bin1), min(bin_inds), max(bin_inds)
  ; Set the companionship flags
  star[bin1].pri = 1
  bin_star.sec = 1
  ; Cross-reference the indices
  star[bin1].companion.ind = bin_inds
  bin_star.companion.ind = bin1
  ; Generate period of companion
  logpbar = alog10(365.25*5.3^(1.5)*(star[bin1].m*(1.0+q))^(-0.5))
  print, median(logpbar)
  logpsig = 1.3
  logp = randomn(seed, n_elements(bin1))*logpsig + logpbar
  bin_star.companion.p = 10.^logp
  star[bin1].companion.p = bin_star.companion.p
  bin_star.companion.a = (star[bin1].m*(1.0+q))^(1./3.)*(star[bin1].companion.p/365.25)^(2./3.)
  star[bin1].companion.a = bin_star.companion.a
  ; Now fill in the other stellar parameters (interpolate based on mass)
  bin_star.r = interpol(r, m, bin_star.m)
  ; Enforce minimum radius
  bin_star[where(bin_star.r lt 0.1)].r = 0.1
  bin_star.teff = interpol(teff, m, bin_star.m)
  bin_star.coord = star[bin1].coord
  bin_star.mag.mv = interpol(vmag, m, bin_star.m)
  bin_star.mag.mi = interpol(imag, m, bin_star.m)
  bin_star.mag.mj = interpol(jmag, m, bin_star.m)
  bin_star.mag.mk = interpol(kmag, m, bin_star.m)
  bin_star.mag.v = bin_star.mag.mv + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.i = bin_star.mag.mi + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.j = bin_star.mag.mj + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.k = bin_star.mag.mk + 5.0*alog10(bin_star.coord.d/10.)
  angseps = star[bin1].companion.a/star[bin1].coord.d 
  cosi = -1.0 + 2.0*randomu(seed, n_elements(bin1))
  phi = !DPI*2.0*randomu(seed, n_elements(bin1))
  angseps = angseps*sqrt(sin(phi)^2. + cosi^2.*cos(phi)^2.)
  star[bin1].companion.sep = angseps
  bin_star.companion.sep = angseps
  
  star = struct_append(star, bin_star)
  delvar, bin_star

  bin2 = where((randomu(seed, total(nstars)) lt 0.44) and (star.m gt 0.6))
  ; Generate mass of companion
  randomp, q, 0.3, n_elements(bin2), range_x=[0.0, 1.0], seed=seed
  newm = star[bin2].m*q
  gdbin = where(newm ge min(m)) ; Enforce brown dwarf desert
  bin2 = bin2[gdbin] 
  ; Struct to hold the companions
  bin_star = replicate(template_star, n_elements(bin2))
  bin_star.m = newm[gdbin]
  ; Indices of the binaries
  bin_inds = lindgen(n_elements(bin2))+idx
  print, 'bin2',  n_elements(bin2), min(bin_inds), max(bin_inds)
  ; Set the companionship flags
  star[bin2].pri = 1
  bin_star.sec = 1
  ; Cross-reference the indices
  star[bin2].companion.ind = bin_inds
  bin_star.companion.ind = bin2
  logpbar = alog10(365.25*45.^(1.5)*(star[bin2].m*(1.0+q))^(-0.5))
  print, median(logpbar)
  logpsig = 2.3
  logp = randomn(seed, n_elements(bin2))*logpsig + logpbar
  bin_star.companion.p = 10.^logp
  star[bin2].companion.p = bin_star.companion.p
  bin_star.companion.a = (star[bin2].m*(1.0+q))^(1./3.)*(star[bin2].companion.p/365.25)^(2./3.)
  star[bin2].companion.a = bin_star.companion.a
  ; Now fill in the other stellar parameters (interpolate based on mass)
  bin_star.r = interpol(r, m, bin_star.m)
  bin_star[where(bin_star.r lt 0.1)].r = 0.1
  bin_star.teff = interpol(teff, m, bin_star.m)
  bin_star.coord = star[bin2].coord
  bin_star.mag.mv = interpol(vmag, m, bin_star.m)
  bin_star.mag.mi = interpol(imag, m, bin_star.m)
  bin_star.mag.mj = interpol(jmag, m, bin_star.m)
  bin_star.mag.mk = interpol(kmag, m, bin_star.m)
  bin_star.mag.v = bin_star.mag.mv + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.i = bin_star.mag.mi + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.j = bin_star.mag.mj + 5.0*alog10(bin_star.coord.d/10.)
  bin_star.mag.k = bin_star.mag.mk + 5.0*alog10(bin_star.coord.d/10.)
  angseps = 2.*star[bin2].companion.a/star[bin2].coord.d
  cosi = -1.0 + 2.0*randomu(seed, n_elements(bin2))
  phi = !DPI*2.0*randomu(seed, n_elements(bin2))
  angseps = angseps*(1.0-cosi*sin(phi)^2.0)
  star[bin2].companion.sep = angseps
  bin_star.companion.sep = angseps
  
  star = struct_append(star, bin_star)
  delvar, bin_star
  endif
  save, filen=outfile, star

  print, 'Total number of stars = ', n_elements(star)
  if (v) then begin


     !p.multi=[0,3,3]
     !p.charsize=3
     plotsym,0,/fill

;  don't plot more than 10000 points
     if (n_elements(star) gt 10000) then begin
        q = round(double(n_elements(star)-1)*randomu(seed,10000))
     endif else begin
        q = indgen(n_elements(star))
     endelse

     aitoff, star[q].coord.glon, star[q].coord.glat, x, y
     plot, x, y, psym=3, /isotropic, xtit='glon', ytit='glat'

     aitoff, star[q].coord.elon, star[q].coord.elat, x, y
     plot, x, y, psym=3, /isotropic, xtit='elon', ytit='elat'

     aitoff, star[q].coord.ra, star[q].coord.dec, x, y
     plot, x, y, psym=3, /isotropic, xtit='ra', ytit='dec'

     q = where(star.coord.d lt dmax[0])
     vol = (4./3.)*!PI*(dmax[0])^3.
     plothist, star[q].mag.mv, x, y, bin=1.0, xra=[3,15], xsty=1, /noplot
     plot, x, y/vol, psym=8, $
           xtit='abs V mag', xra=[3,15], xsty=1, $
           ytit='number mag!E-1!N pc!E-3!N'

     readcol, '../Stars/recons_vlf.txt', f='D,D', comment='#', mv, n_recons, /silent
     oplot, mv, n_recons/vol10, psym=10, color=fsc_color('Red')

     plothist, star[q].mag.mj, x, y, bin=0.5, xra=[3,15], xsty=1, /noplot
     plot, x, y/vol, psym=8, $
           xtit='abs J mag', xra=[3,12], xsty=1, $
           ytit='number (0.5 mag)!E-1!N pc!E-3!N'

     readcol, '../Stars/jlfs.txt', f='D,D', comment='#', mj, nj, /silent
     oplot, mj, nj, psym=10, color=fsc_color('Red')

     plothist, star[q].r, x, y, bin=0.05, xra=[0.1,1], xsty=1, /noplot
     plot, x, y/vol, psym=10, $
           xtit='radius [R!Dsun!N]', xra=[0.1,1], xsty=1, $
           ytit='number pc!E-3!N per radius bin'

     q = where(abs(star.coord.glon-60.) lt 30.0)
     plothist, star[q].coord.glat, bin=2, x, y, /noplot
     y = y/cos(x*!PI/180.)
     plot, x, y, psym=10, xtit='galactic latitude', ytit='rel number of stars sr!E-1!N'

  endif


end

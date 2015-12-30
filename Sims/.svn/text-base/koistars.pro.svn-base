pro koistars, verbose=verbose, outfile=outfile

  if (keyword_set(verbose)) then v=1 else v=0
  if (keyword_set(outfile)) then fname=outfile else fname='spp_koi.sav'

;  restore, filen='../Stars/star_properties.sav'
  koidat = mrdfits('../Kepler/koi.fits')
  
  ra   =  reform(koidat[1,*])
  dec   = reform(koidat[2,*])
  srad  = reform(koidat[20,*])
  smass = reform(koidat[21,*])
  teff  = reform(koidat[17,*])
  imag  = reform(koidat[25,*])
  jmag  = reform(koidat[26,*])
  a     = reform(koidat[10,*])
  b     = reform(koidat[8,*])
  p     = reform(koidat[4,*])
  inc   = reform(koidat[9,*])
  dur   = reform(koidat[6,*]/24.0) ; hours to days
  dep   = reform(koidat[5,*]*1.0e-6) ; ppm to unity
  
  nkoi = n_elements(ra)
  print, 'Creating ', nkoi, ' KOIs'
; what kind of star catalog?

;  ilim = 12.0 + 0.*r
;  q = where(r le 0.5)
;  ilim[q] = 13.0
;  dmax = 10.*10.0^(0.2*(ilim-imag))

;  dmax = 72.0 + 128.*r

; Fressin bins
;  n_per = 11
;  n_rad = 5

;  h = 300.0 + 0.0*r

;  vol10 = (4./3.)*!PI*(10.0^3.)
;  dr = r[1]-r[0]
;  phi = dndr*dr/0.1/vol10       ; now phi is in stars per cubic parsec per radius bin

;  nstars = number_of_stars_expz(dmax, phi, h)
;  nstars = round(nstars)

  template_planet = {$
                    n: 0, $
                    r: 0.0, $ ; can eventually replace by dblarr(5), $
                    p: 0.0, $ 
                    a: 0.0, $
                    s: 0.0, $ ; incident flux in units of Sun-->Earth flux
                    b: 0.0, $
                    tra: 0, $
                    dep: 0.0, $
                    dur: 0.0, $
                    ntra_obs: 0, $
                    det: 0, $
                    snr: 0.0, $
 		    snrtran:0.0 $
                    }

 ; template_planets = { 	r: dblarr(n_rad), $
;			p: dblarr(n_per), $
;			a: dblarr(n_per), $
;			s: dblarr(n_per), $
;			b: dblarr(n_per), $
;			tra: intarr(n_per), $
;			dur: dblarr(n_per), $
;			ntra_obs: intarr(n_per), $
;			snr: dblarr(n_per), $
;			dep: dblarr(n_rad), $
;			det: intarr(n_per, n_rad)}

  template_coord = {$
                   ra: 0.0, $
                   dec: 0.0, $
                   elon: 0.0, $
                   elat: 0.0, $
                   glon: 0.0, $
                   glat: 0.0, $
                   d: 0.0, $
 		   fov_r: 0.0 $
              }
  
  template_mag = {$
                 mv: 0.0, $ ; absolute
                 mi: 0.0, $ ; absolute
                 mj: 0.0, $ ; absolute
                 v: 0.0, $ ; apparent
                 i: 0.0, $ ; apparent
                 j: 0.0 $ 
                 }

  template_star = {$
                  r: 0.0, $		; Radius in R_sun
                  m: 0.0, $		; Mass in M_sun
                  teff: 0.0, $		; Temperature in K
 		  cosi: 0.0, $
		  dx: 0, $		; Pixel x offset (0-1)
		  dy: 0, $		; Pixel y offset (0-1)
 		  npix: 0, $
 		  snr: 0.0, $
		  sat: 0, $
                  dil: 0.0, $
                  npointings: 0, $	; Number of TESS observation blocks
                  coord: template_coord, $
                  mag: template_mag, $
                  planet: template_planet}
		  ;planet_hz: template_planet}	; Single planet assignment
;                 planets: template_planets $   ; Planet detectability 


;  for i=0,n_elements(r)-1 do begin
     
;     if (nstars[i] lt 1) then continue

     star = replicate(template_star, nkoi)
     star.r = srad
     star.m = smass
     star.teff = teff
;    tmp_star.mag.mv = vmag[i]
     star.mag.i = imag
     star.mag.j = jmag

;    assign_xyz, nstars[i], dmax[i], h[i], x, y, z, d

;     tmp_star.coord.d = d
;     dm = 5.0*alog10(tmp_star.coord.d/10.0D0)
;     tmp_star.mag.v = tmp_star.mag.mv + dm
;     tmp_star.mag.i = tmp_star.mag.mi + dm
;     tmp_star.mag.j = tmp_star.mag.mj + dm
     
     glactc, ra, dec, 2000, gl, gb, 1, /DEGREE
     star.coord.glat = gb ;180./!PI*atan(z,sqrt(x^2. + y^2.))
     star.coord.glon = gl ;180./!PI*atan(y,x)

     euler, star.coord.glon, star.coord.glat, elon, elat, 6
     ;q = where(elon gt 180.)
     ;if (q[0] ne -1) then elon[q] = elon[q]-360.0
     ;euler, tmp_star.coord.glon, tmp_star.coord.glat, ra, dec, 2
     star.coord.elon = elon & star.coord.elat = elat
     star.coord.ra = ra & star.coord.dec = dec

;     if (n_elements(star) eq 0) then begin
;        star = tmp_star
;     endif else begin
;        star = struct_append(star, tmp_star)
;     endelse

;  endfor

  ; Random orientation of all stars
  star.cosi = cos((!DPI/2.0)*inc) ;cos((!DPI/2.0)*randomu(seed, n_elements(star)))
  ; Random pixel displacement of all stars
  star.dx   = floor(10.*randomu(seed, nkoi))
  star.dy   = floor(10.*randomu(seed, nkoi))
 ;  psf_x = mrdfits('../xPSF.fits')
 ;  psf_y = mrdfits('../yPSF.fits')
 ;  psf_f = mrdfits('../OnAxisPSF.fits')
  ; for ii = 0, n_elements(star)-1 do begin
  ;  calc_opt_npix_offset, star[ii].mag.i, 3600., npix, snr, $
  ;      psf_f = psf_f, psf_x = psf_x, psf_y = psf_y, $
  ;	dxdy=[star[ii].dx,star[ii].dy], $
  ;      elon=star[ii].coord.elon, $
  ;      elat=star[ii].coord.elat, $
  ;      teff=star[ii].teff, fluxfrac=frac
  ;  star[ii].npix = npix
  ;  star[ii].frac = frac
  ;  if(v) then print, 'Star ', ii, ' NPIX = ', npix
  ; end

; Planet params
  star.planet.n=1
  star.planet.p=p
  star.planet.a = a
  star.planet.s = (star.r)^2.0 * (star.teff/5777.0)^4. / (star.planet.a)^2. ;
  star.planet.b = b
  star.planet.tra = 1
  star.planet.dur = dur
  star.planet.dep = dep

 ; for ii=0, nkoi-1 do begin
 ;   if(total(star[ii].planet.tra) gt 0) then begin
 ;     calc_opt_npix_offset, star[ii].mag.i, 3600., npix, snr, $
 ;       psf_f = psf_f, psf_x = psf_x, psf_y = psf_y, $
 ;       dxdy=[star[ii].dx,star[ii].dy], $
 ;       ;saturation=SATURATION, $
 ;       verbose=0, $
 ;       elon=star[ii].coord.elon, $
 ;       elat=star[ii].coord.elat, $
 ;       teff=star[ii].teff, fluxfrac=frac
 ;     star[ii].npix = npix
 ;     star[ii].frac = frac
 ;     star[ii].snr  = snr
 ;     if(v) then print, 'Star ', ii, ' NPIX = ', npix, ' SNR = ', snr
 ;     if(v) then print, 'Ecliptic lon: ', star[ii].coord.elon, ' lat: ', star[ii].coord.elat
;    endif
;  endfor
  save, filen=fname, star


  if (v) then begin

     print, 'Total number of stars = ', n_elements(star)

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

pro observe, filen=filen

;;;;;; basic parameters here

  DWELL_TIME = 27.4d0                     ; days per field
  DWELL_TIME = DWELL_TIME - 16.0d0/24.0d0 ; correct for downlink time
;  DWELL_TIME = 27.4d0/2.0                ; days per field
;  DWELL_TIME = DWELL_TIME - 8.0d0/24.0d0  ; correct for downlink time
  SNR_MIN = 7.0
  NTRA_OBS_MIN = 2
  GEOM_AREA = 73.0; cm^2
  SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0
  PIX_SCALE = 20.0; arcsec per pixel side
  SURVEY_FILE = '../Survey/survey.sav'
  RMAX_HAB = 2.0 & SHZ_IN = 0.5 & SHZ_OUT = 1.5  ; radii, incident fluxes corresponding to inner/outer edge of HZ

;;;;;;; set up some plotting stuff

  !p.charsize=2
  syms_multiplier=1

  if (keyword_set(filen)) then begin

     set_plot,'ps'
     device,filen='sim_'+filen+'.eps',/encapsulated,xsize=6.5,ysize=9.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.5
     syms_multiplier=0.5

  endif

  !P.multi=[0,2,3]
  plotsym,8,/fill
  c=fsc_color(['Red','Orange Red','Orange','Yellow','Green','Blue','Violet'])

;;;;;;;; let's get to work

  restore, 'spp_koi.sav'
  nstars = n_elements(star)
  
; how long we have observed each star?

  restore, SURVEY_FILE
  si = size(Number_of_Pointings)
  n_pix_lon = si[1] & n_pix_lat = si[2]

  for i=0L,nstars-1 do begin
     i_mask = floor(n_pix_lon * star[i].coord.elon/360.0)
     j_mask = floor(n_pix_lat * (star[i].coord.elat+90.0)/180.0)
     star[i].npointings = Number_of_Pointings[i_mask,j_mask]
  endfor

  nplot = 10000
  plotme = round(double(nstars-1)*randomu(seed,nplot))
  color_index = star[plotme].npointings
  q = where(color_index gt 6)
  if (q[0] ne -1) then color_index[q] = 6
  
  aitoff, star[plotme].coord.elon, star[plotme].coord.elat, x, y
  plot, x, y, psym=3, /isotropic, /nodata, $
        tit='number of pointings', $
        xtit='ecliptic longitude', $
        ytit='ecliptic latitude'

  for i=0,nplot-1 do begin
     oplot, [x[i]], [y[i]], psym=3, color=c[color_index[i]]
  endfor

; to plot a band representing the galactic plane
;
;  m = 100
;  glon = -180. + 360.0*dindgen(m)/double(m-1)
;  glat = 0.0*glon
;  euler, glon, glat, elon, elat, select=6
;  aitoff, elon, elat, x_gal, y_gal
;  oplot, x_gal, y_gal, psym=8, syms=0.5*syms_multiplier, color=fsc_color('Black') 

; for each transiting planet, how many transits were observed?

  tra = where(star.planet.tra gt 0)
  star[tra].planet.ntra_obs = n_eclipses(star[tra].planet.p, DWELL_TIME*double(star[tra].npointings))

; for each observed transiting planet, calculate snr

  obs = where(star.planet.tra gt 0 and star.planet.ntra_obs gt 0)

;  for i=0,n_elements(obs)-1 do begin

     exptime =  star[obs].planet.dur * 24.0 * 3600.
     calc_noise, star[obs].mag.i, exptime, noise, $
;     calc_noise_hack, star[j].mag.i, exptime, noise, $
                 teff=star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat
     star[obs].planet.snrtran = star[obs].planet.dep/noise
;     j = obs[i]

     exptime = double(star[obs].planet.ntra_obs) * star[obs].planet.dur * 24.0 * 3600.
     calc_noise, star[obs].mag.i, exptime, noise, $
;     calc_noise_hack, star[j].mag.i, exptime, noise, $
                 teff=star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat
     star[obs].planet.snr = star[obs].planet.dep/noise

;  endfor

;   decide if it is 'detected'.

  detected = where(star.planet.tra gt 0 and star.planet.ntra_obs ge NTRA_OBS_MIN and star.planet.snr ge SNR_MIN)
  star[detected].planet.det = 1

; Report basic parameters

  if (keyword_set(filen)) then begin
     openw, lun, 'sim_' + filen + '.dat', /get_lun
  endif else begin
     lun=-1
  endelse

  fmt = '(A25, F10.2)'
  printf,lun,f=fmt,'dwell_time = ', DWELL_TIME
  printf,lun,f=fmt,'geom_area = ', geom_area
  printf,lun,f=fmt,'pix_scale = ', pix_scale
  printf,lun,f=fmt,'e_pix_ro = ', e_pix_ro
  printf,lun,f=fmt,'sys_limit = ', sys_limit
  printf,lun,f=fmt,'snr_min = ', SNR_MIN
  printf,lun,f=fmt,'ntra_obs_min = ', ntra_obs_min
  fmt = '(A25, A40)'
  printf,lun,f=fmt,'survey file = ', survey_file

; Report breakdown of different types of DETECTED TRANSITING planets

  readcol, 'planet_definitions.txt', /silent, comment='#', f='A,D,D,D,D', $
           planet_type, p_min, p_max, r_min, r_max

  fmt = '(A20, A9, I3, A3, I3, A9, F5.2, A3, F5.2, A4, I6)'
  for i=0,n_elements(planet_type)-1 do begin
     q = where(star.planet.det gt 0 and $
               star.planet.p ge p_min[i] and star.planet.p lt p_max[i] and $
               star.planet.r ge r_min[i] and star.planet.r lt r_max[i])
     if (q[0] ne -1) then npla=n_elements(q) else npla=0
     printf, lun, f=fmt, planet_type[i], ' (Period ', p_min[i], ' - ', p_max[i], ', Radius ', r_min[i], ' - ', r_max[i], ') = ', npla
  endfor

  hab = where(star.planet.det gt 0 and $
              star.planet.s ge SHZ_IN and star.planet.s le SHZ_OUT and $
              star.planet.r le RMAX_HAB)
  if (hab[0] eq -1) then nhab = 0 else nhab = n_elements(hab)
  printf, lun, f='(A20,I3)', 'R<2 HZ planets = ', nhab

  if (keyword_set(filen)) then begin
     close,lun
     free_lun, lun
  endif

; next make a bunch of plots of all the transiting planets.

  tra = where(star.planet.tra gt 0)
  aitoff, star[tra].coord.elon, star[tra].coord.elat, x, y
  plot, x, y, psym=3, /isotropic, $
        tit='transit observations', $
        xtit='ecliptic longitude', $
        ytit='ecliptic latitude'

  det = where(star.planet.det gt 0)
  aitoff, star[det].coord.elon, star[det].coord.elat, x, y

  plotsym,8,/fill
  for i=0,n_elements(det)-1 do begin
     j = det[i]
     syms = 0.3*alog10(star[j].planet.snr)
     color_index = star[j].planet.ntra_obs
     if (color_index gt 6) then color_index=6
     oplot, [x[i]], [y[i]], psym=8, syms=syms*syms_multiplier, color=c[color_index]
  endfor

  xvar = [ [star.mag.i],    [star.planet.p],   [star.r],        [star.planet.p] ]
  yvar = [ [star.planet.r], [star.planet.dep], [star.planet.r], [star.planet.r] ]

  xtit = [ 'apparent I mag',$
           'period [d]',$
           'star radius [R!Dsun!N]',$
           'period [d]' $
         ]

  ytit = [ 'planet radius [R!DE!N]',$
           'transit depth',$
           'planet radius [R!DE!N]',$
           'planet radius [R!DE!N]' $
         ]

  xlog = [0,1,0,1]
  ylog = [0,1,0,0]

  xra1 = [6, 0.3, 0.1, 0.3]
  xra2 = [14, 130, 1.0, 130]
  yra1 = [0.5, 5d-5, 0.5, 0.5]
  yra2 = [4.0, 5d-2, 4.0, 4.0]

  n = size(xvar) & n = n[2]
  
  for i=0,n-1 do begin

     plot, xvar[tra,i], yvar[tra,i], psym=3, $
           xtit=xtit[i], xra=[xra1[i],xra2[i]], xlog=xlog[i], xsty=1, $
           ytit=ytit[i], yra=[yra1[i],yra2[i]], ylog=ylog[i], ysty=1

     for j=0,n_elements(det)-1 do begin
        k = det[j]
        plotsym,8,/fill
        syms = 0.6*alog10(star[k].planet.snr)
        color_index = star[k].planet.ntra_obs
        if (color_index gt 6) then color_index=6
        oplot, [xvar[k,i]], [yvar[k,i]], psym=8, syms=syms*syms_multiplier, color=c[color_index]
     endfor

  endfor

  save, filen='spo_'+filen+'.sav', star

  if (keyword_set(filen)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

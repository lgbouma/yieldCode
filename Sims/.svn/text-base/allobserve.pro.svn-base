pro allobserve, filen=filen, verbose=verbose
  if (keyword_set(verbose)) then v=1 else v=0
;;;;;; basic parameters here

  DWELL_TIME = 13.7d0                     ; days per field
  DOWNLINK_TIME = 16.0d0/24.0d0 ; correct for downlink time
  DTY_CYC = 1.0
;  DWELL_TIME = 27.4d0/2.0                ; days per field
;  DWELL_TIME = DWELL_TIME - 8.0d0/24.0d0  ; correct for downlink time
  SNR_MIN = 7.0
  NTRA_OBS_MIN = 2
  GEOM_AREA = 73.0; cm^2
  SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0
  SATURATION = 150000.0 ; Saturation level in e-
  SUB_EXP_TIME = 2.0  ; sec
  PIX_SCALE = 20.0; arcsec per pixel side
  SURVEY_FILE = '../Survey/survey.sav'
  PSF_FILE = 'HalfEdgePSF.fits'
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

  restore, 'spp_halfedge.sav'
  n_stars = n_elements(star)
  n_rad = n_elements(star[0].planets.r)
  n_per = n_elements(star[0].planets.p)  

; how long we have observed each star?

  restore, SURVEY_FILE
  si = size(Number_of_Pointings)
  n_pix_lon = si[1] & n_pix_lat = si[2]

  for i=0L,n_stars-1 do begin
     i_mask = floor(n_pix_lon * star[i].coord.elon/360.0)
     j_mask = floor(n_pix_lat * (star[i].coord.elat+90.0)/180.0)
     star[i].npointings = Number_of_Pointings[i_mask,j_mask]
  endfor

  nplot = 10000
  plotme = round(double(n_stars-1)*randomu(seed,nplot))
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
; (this randomizes the planets' orbital phase 
 
; for each observed transiting planet, calculate snr
  obs = where(	(total(star.planets.tra,1) gt 0) and $
		(star.npix gt 0))

;  for ii=0,n_stars-1 do begin
    ;if (total(star[ii].planet.ntra_obs) gt 0) then begin
  for ii=0, n_per-1 do begin
    star[obs].planets.ntra_obs[ii] = $
			n_eclip_blank(star[obs].planets.p[ii], $ 
			DWELL_TIME, 2.0*double(star[obs].npointings), $
			periblank=DOWNLINK_TIME, $
			apoblank=(1.0-DTY_CYC)*(DWELL_TIME-DOWNLINK_TIME))
    exptime = double(star[obs].planets.ntra_obs[ii]) * star[obs].planets.dur[ii] * 24.0 * 3600.
    ;if (exptime gt 0) then begin
    calc_noise, star[obs].mag.i, exptime, noise, $
	 npix=star[obs].npix, $
      	 frac_aper=star[obs].frac, $ 
	 subexptime=SUB_EXP_TIME, $
         teff=star[obs].teff, $
         e_pix_ro = E_PIX_RO,$
         geom_area = GEOM_AREA, $
         sys_lim = SYS_LIMIT, $
         pix_scale = PIX_SCALE, $
         elon=star[obs].coord.elon, $
         elat=star[obs].coord.elat
    star[obs].planets.snr[ii] = 1.0/noise
    ; nonzero exposure time
    for jj=0, n_rad-1 do begin
      star[obs].planets.det[ii,jj] = ((star[obs].planets.dep[jj]/noise gt SNR_MIN) and $
				      (star[obs].planets.ntra_obs[ii] ge NTRA_OBS_MIN))
    endfor 	; radius loop
     
  endfor		; period loop

;   decide if it is 'detected'.

;  detected = where(star.planet.tra gt 0 and star.planet.ntra_obs ge NTRA_OBS_MIN and star.planet.snr ge SNR_MIN)
;  star[detected].planet.det = 1

; Report basic parameters
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

  nstars_detectable = total(star.planets.det,3)
  nplanets_detected = round(nstars_detectable*rate_fressin)

  print, total(nstars_detectable,1)
  print, total(nplanets_detected,1)
  
  hz_star_tally=intarr(n_rad)
  hz_planet_tally=intarr(n_rad)
  for ii=0,n_per-1 do begin
    s = star.planets.s[ii]
    hz = where((s ge SHZ_IN) and (s le SHZ_OUT))
    ;print, n_elements(hz)
    if (hz[0] ne -1) then begin
	hz_star_tally = hz_star_tally+ $ 
		round(total(star[hz].planets.det[ii,*],3))
        hz_planet_tally = hz_planet_tally+round(rate_fressin[ii,*] * $
		total(star[hz].planets.det[ii,*],3))
    endif
  endfor
  print, hz_star_tally
  print, hz_planet_tally
end

pro calc_opt_npix_offset, $
  ; Mandatory inputs:
  imag, $ 			; Cousins I mag
  exptime, $			; Exposure time in seconds
  ; Mandatory outputs:
  npix_opt, $ 			; Optimal (integer) number of pixels
  snr_opt, $			; Optimized SNR, 0 if saturated
  ; Other inputs:
  dxdy=dxdy, $ 			; Centroid offset, generated by randomu(seed,1)
  gridsize=gridsize, $		; Grid size parameter
  saturation=saturation, $ 	; Full well capacity (electrons)
  teff=teff, $                  ; effective temperature in Kelvins
  elon=elon, elat=elat, $       ; ecliptic coordinates in degrees
  subexptime=subexptime, $      ; subexposure time (n_exp = exptime/subexptime)
  npix_aper = npix_aper, $      ; number of pixels in photometric aperture
  frac_aper = frac_aper, $      ; fraction of flux enclosed in photometric aperture
  e_pix_ro = e_pix_ro, $        ; rms in no. photons/pixel from readout noise
  geom_area = geom_area, $      ; geometric collecting area
  pix_scale = pix_scale, $      ; arcsec per pixel
  sys_limit = sys_limit, $      ; minimum uncertainty in 1 hr of data, in ppm
  psf_f = psf_f, $		; e.g. mrdfits('OnAxisPSF.fits')
  psf_x = psf_x, psf_y = psf_y, $
  fluxfrac = fluxfrac, $	; Output the enclosed fraction
  verbose=verbose	        ; request verbose output

  !p.charsize=2

  if (keyword_set(dxdy)) then begin
     dx = dxdy[0]
     dy = dxdy[1]
  endif else begin
     dx = 0.5
     dy = 0.5
  endelse

  if (keyword_set(teff)) then teff=teff else teff=5000.
  if (keyword_set(elon)) then elon=elon else elon=0.0
  if (keyword_set(elat)) then elat=elat else elat=30.0
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.0
  if (keyword_set(e_pix_ro)) then e_pix_ro=e_pix_ro else e_pix_ro=10.
  if (keyword_set(geom_area)) then geom_area=geom_area else geom_area=73.0
  if (keyword_set(pix_scale)) then pix_scale=pix_scale else pix_scale=20.0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(gridsize)) then gridsize = gridsize else gridsize = 8
  if (keyword_set(verbose))  then v=1 else v=0
  if (keyword_set(psf_f)) then psf_f = psf_f else psf_f=mrdfits('../OnAxisPSF.fits')
  if (keyword_set(psf_x)) then psf_x = psf_x else psf_x=mrdfits('../xPSF.fits')
  if (keyword_set(psf_y)) then psf_y = psf_y else psf_y=mrdfits('../yPSF.fits')
 
  check_sat = keyword_set(saturation)

  ; Define the pixel boundaries
  x = (fltarr(gridsize)+1.0)##(findgen(gridsize)-ceil(gridsize/2.0)+1.0)
  y = (findgen(gridsize)-ceil(gridsize/2.0)+1.0)##(fltarr(gridsize)+1.0)
  
  ;psf_x = mrdfits('../xPSF.fits')
  ;psf_y = mrdfits('../yPSF.fits')
  ;psf   = mrdfits(psf_file)
  ;if (psf_file eq '../HalfEdgePSF.fits') then dy = dy-0.133333
;  tot = int_tabulated_2D(psf_x, psf_y, psf)
;  psf = psf/tot
;  if(v) then print, x
;  if(v) then print, y
  ; Pixel centers are [0.5, 0.5] away from boundaries. Calculate offset from these.
  r2 = (x+0.5-dx)^2 + (y+0.5-dy)^2

  sort_ind = sort(r2)
  r_sort  =  sqrt(r2[sort_ind])
  x_sort  =  x[sort_ind]
  y_sort  =  y[sort_ind]

;  if(v) then print, r_sort

  ii=0
  snrprev = 0.0
  snr     = 0.0
  accumfrac = 0.0
  accumprev = 0.0
  while (((snrprev lt snr) and (ii lt gridsize^2-2)) or (ii lt 3)) do begin 
    ;openw, lun, 'halfnhalf.pro', /get_lun
    ;printf, lun, 'function halfnhalf, y'
    ;printf, lun, '    return,  [', y_sort[ii]-dy, ', ', y_sort[ii]+1.0-dy, ']'
    ;printf, lun, 'end'
    ;free_lun, lun
    ; resolve_routine, 'halfnhalf', /IS_FUNCTION
    ; Do 2d integral over PSF for each additional pixel   
    ;accumfrac = accumfrac + int_2d('lorentz_psf', ([x_sort[ii]-dx, x_sort[ii]+1.0-dx]), 'halfnhalf', 96)
    accumprev = accumfrac
    thispix = where(	(psf_x ge x_sort[ii]-dx) and $
			(psf_x le x_sort[ii]+1.0-dx) and $
			(psf_y ge y_sort[ii]-dy) and $
			(psf_y le y_sort[ii]+1.0-dy))
    accumfrac = accumfrac + int_tabulated_2D(psf_x[thispix], psf_y[thispix], psf_f[thispix])
    calc_noise, imag, exptime, noise, npix_aper=(ii+1), frac_aper=(accumfrac), e_star_sub=estar, $
	teff=teff, elon=elon, elat=elat, subexptime=subexptime, e_pix_ro=e_pix_ro, geom_area=geom_area, $
	pix_scale=pix_scale, sys_limit=sys_limit
    if ((ii EQ 0) and (check_sat)) then begin
	if (estar GT saturation) then break ;outputs npix=-1, snr=0	
    endif
    snrprev = snr
    if (finite(noise) and noise ne 0.0) then $
	snr = 1.0/noise else snr = 0.0
    if (v) then begin
	 print, 'NPIX = ', ii, '  NOISE = ', noise, '  FRAC = ', accumfrac
	print, 'X = ', x_sort[ii], '  Y = ', y_sort[ii], '  R = ', r_sort[ii]
    end
    ii=ii+1  
endwhile
  
  npix_opt = ii-1	 ;Previous value of npix BEFORE the snr turned over
  snr_opt  = snrprev
  fluxfrac = accumprev
end

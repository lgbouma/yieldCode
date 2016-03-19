pro calc_noise_eclip, $
;
; mandatory inputs
;
   ph_star, $                       ; ph/s/cm^2 from (npixels x nstars)
   ph_dil, $                        ; ph/s/cm^2 per pixel
   exptime, $                       ; total exposure time in seconds
   e_pix_ro, $                      ; read noise per subexposure in e-
   sys_limit, $ 	                  ; noise floor (ppm)
;
; mandatory outputs
;
   noise, $                         ; fractional rms (= 1/SNR)
;
; optional inputs
;
   subexptime=subexptime, $         ; subexposure time (n_exp = exptime/subexptime)
   npix_aper = npix_aper, $         ; number of pixels in photometric aperture
   ;frac_aper = frac_aper, $         ; fraction of flux enclosed in photometric aperture
   geom_area = geom_area, $         ; geometric collecting area
   aspix = aspix, $                 ; arcsec per pixel
   zodi_ph = zodi_ph, $             ; zodiacal photons/s/cm^2
   verbose=verbose, $               ; request verbose output
   field_angle = field_angle, $	    ; Field angle for effective area
   bin_sys=bin_sys, $		    ; Is this a binary?
   bin_sep=bin_sep, $		    ; Separation to binary
   bin_ph=bin_ph, $		    ; imag of binary
   cr_noise=cr_noise, $             ; noise from cosmic rays

; optional outputs
;
   noise_star=noise_star, $         ; noise from star counts alone
   noise_dil=noise_dil, $           ; noise from sky counts only
   noise_zodi=noise_zodi, $           ; noise from sky counts only
   noise_ro=noise_ro,$              ; noise from readout only
   noise_cr=noise_cr,$              ; noise from readout only
   noise_sys=noise_sys, $           ; noise from systematic limit only
   e_tot_sub=e_tot_sub, $ 	    ; subexposure electron count (for saturation check)
   
   dilution=dilution		    ; background dilution factor
;
;
  if (keyword_set(zodi_ph)) then zodi_ph=zodi_ph else zodi_ph = 0.
  if (keyword_set(cr_noise)) then cr_noise=cr_noise else cr_noise = 0.
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.0
  if (keyword_set(geom_area)) then geom_area=geom_area else geom_area=69.1
  if (keyword_set(pix_scale)) then pix_scale=pix_scale else pix_scale=21.1
  if (keyword_set(verbose)) then v=1 else v=0
  if (keyword_set(field_angle)) then field_angle=field_angle else field_angle=0.0
  
  sz_ph = size(ph_star)
  npix_max = sz_ph[1]
  if (v) then print, 'npix_max = ', npix_max

  if (keyword_set(npix_aper)) then begin
     npix_aper=npix_aper
  endif else begin
     npix_aper=3+intarr(n_elements(ph_star))
  endelse
 
  if (v) then print, 'npix_aper = ', npix_aper
        
  omega_pix = aspix^2.
  n_exposures = exptime/subexptime

  ; electrons from the star
  e_star =   reform(ph_star[npix_aper-1,*]) * geom_area * cos(!DPI * field_angle/180.)* exptime
  e_pix_dil = reform(ph_dil[npix_aper-1,*]) * geom_area * cos(!DPI * field_angle/180.)* exptime
  
  if (keyword_set(noise_cr)) then noise_cr = reform(noise_cr[*,npix_aper-1]) else noise_cr = 0.0 
  
  if (v) then print, 'e_star = ', median(e_star)

  ; e- from zodi
  e_pix_zodi = zodi_ph * geom_area * cos(!DPI * field_angle/180.) *  exptime * npix_aper

  if (v) then print, 'vmag_zodi = ', median(vmag_zodi)
  if (v) then print, 'e_pix_zodi = ', median(e_pix_zodi)

; OLD  
; e/pix from background stars
;  euler, elon, elat, glon, glat, select=5
;  if (v) then print, 'glon, glat = ', glon, glat
;  dlon = glon
;  q = where(dlon gt 180.) & if(q[0] ne -1) then dlon[q] = 360.-dlon[q]
;  dlon = abs(dlon)/180.0D0
;  p = bk_p[fov_ind,*]
;  dlat = abs(glat)/90.0D0
;  imag_bgstars = p[*,0] + $
;	p[*,1]*(1.0-exp(-dlon/p[*,2])) + $
;	p[*,3]*(1.0-exp(-dlat/p[*,4])) + $
;	p[*,5]*sqrt(dlon*dlat)
;  e_pix_bgstars = 10.0^(-0.4*imag_bgstars) * 1.7D6 * $
;	geom_area * cos(!DPI * field_angle/180.) * omega_pix * exptime
;  if (v) then print, 'imag_bgstars = ', median(imag_bgstars)
;  if (v) then print, 'e_pix_bgstars = ', median(e_pix_bgstars)

  ; compute noise sources
  ;e_tot_sky = e_pix_zodi + e_pix_dil

;  e_bin = dblarr(n_elements(e_tot_sky))
;  if(keyword_set(bin_sys)) then nbin=total(bin_sys) else nbin=0
; e/pix from companion
;  if (nbin gt 0) then begin
;    if (v) then print, 'Dealing with ', nbin, ' binaries'
;    bins = where(bin_sys)
;    pix_sep = bin_sep[bins]/pix_scale  	; from arcsec to pixels
;    r = frac_off[fov_ind[bins],*]	; distance (in pixels)
;    di = frac_off[fov_ind[bins]+4,*]    ; imag attenuation
;    dibin = interpol(di, r, pix_sep)    ; interpolate over spline fit
;    e_bin[bins] = bin_ph[bins] * $
;        geom_area * cos(!DPI * field_angle/180.)* exptime * $
;        reform(ph_star[npix_aper-1,*])/reform(ph_star[npix_max-1,*]) ; frac
;    if (v) then print, 'e_pix_bin = ', median(e_bin[bins])
;  endif
  
  e_tot = e_pix_zodi + e_pix_dil + e_star ; noise
  e_phot =  e_pix_dil + e_star ; in photometric aperture
  ; electrons in cadence
  e_tot_sub = e_tot*subexptime/exptime
  
  noise_star = sqrt(e_star) / e_phot
  noise_dil = sqrt(e_pix_dil) / e_phot
  noise_zodi  = sqrt(e_pix_zodi) / e_phot
  noise_ro   = sqrt(npix_aper*n_exposures)*e_pix_ro / e_phot
  noise_cr   = cr_noise / e_phot
  dilution = e_phot / e_star
  noise_sys  = 0.0*noise_star + sys_limit/1d6/sqrt(exptime/3600.)

  noise = sqrt( noise_star^2. + noise_zodi^2. + noise_dil^2. + noise_ro^2. + noise_sys^2. + noise_cr^2.)

end

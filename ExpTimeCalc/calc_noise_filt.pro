pro calc_noise_filt, $
;
; mandatory inputs
;
   ph_star, $                       ; ph/s/cm^2 from (npixels x nstars)
   exptime, $                       ; total exposure time in seconds
;
; mandatory outputs
;
   noise, $                         ; fractional rms (= 1/SNR)
;
; optional inputs
;
   elon=elon, elat=elat, $          ; ecliptic coordinates in degrees
   subexptime=subexptime, $         ; subexposure time (n_exp = exptime/subexptime)
   npix_aper = npix_aper, $         ; number of pixels in photometric aperture
   ;frac_aper = frac_aper, $         ; fraction of flux enclosed in photometric aperture
   e_pix_ro = e_pix_ro, $           ; rms in no. photons/pixel from readout noise
   geom_area = geom_area, $         ; geometric collecting area
   pix_scale = pix_scale, $         ; arcsec per pixel
   sys_limit = sys_limit, $         ; minimum uncertainty in 1 hr of data, in ppm
   verbose=verbose, $               ; request verbose output
;   red=red, $		     	    ; consider a redder bandpass
;   al_bk = al_bk, $		    ; Use Al's sky background
;   al_phot = al_phot, $             ; Use Al's photomery
   fov_ind = fov_ind, $             ; 0-3
   field_angle = field_angle, $	    ; Field angle for effective area
   bk_p=bk_p, $		    	    ; Background polynomial fit
   ph_p=ph_p, $			    ; Photon flux fit
;   frac_off=frac_off, $	 	    ; Fractional offset file (mag vs. pix) 
   bin_sys=bin_sys, $		    ; Is this a binary?
   bin_sep=bin_sep, $		    ; Separation to binary
   bin_ph=bin_ph, $		    ; imag of binary

; optional outputs
;
   noise_star=noise_star, $         ; noise from star counts alone
   noise_sky=noise_sky, $           ; noise from sky counts only
   noise_ro=noise_ro,$              ; noise from readout only
   noise_sys=noise_sys, $           ; noise from systematic limit only
   noise_bin=noise_bin, $	    ; noise from binary companion
   e_star_sub=e_star_sub, $ 	    ; subexposure electron count (for saturation check)
   
   dilution=dilution		    ; background dilution factor
;
;
  if (keyword_set(teff)) then teff=teff else teff=5000.
  if (keyword_set(elon)) then elon=elon else elon=0.0
  if (keyword_set(elat)) then elat=elat else elat=30.0
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.0
  if (keyword_set(e_pix_ro)) then e_pix_ro=e_pix_ro else e_pix_ro=10.
  if (keyword_set(geom_area)) then geom_area=geom_area else geom_area=69.1
  if (keyword_set(pix_scale)) then pix_scale=pix_scale else pix_scale=21.0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
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
  ;if (keyword_set(frac_aper)) then begin
  ;   frac_aper=frac_aper[fov_ind]
  ;endif else begin
     ;frac_aper = 0.76
  ;   fwhm = 1.2
  ;   r_aper = sqrt(npix_aper/!PI)
  ;   frac_aper = 1.0 - (fwhm/(2.*!PI)) / ((fwhm/2.)^2 + r_aper^2.) ; from Roland 9 July 2013
  ;endelse
        
  omega_pix = pix_scale^2.
  n_exposures = exptime/subexptime

  ; electrons from the star
  e_star = reform(ph_star[npix_aper-1,*]) * geom_area * cos(!DPI * field_angle/180.)* exptime
  
  ; electrons in cadence
  e_star_sub = e_star*subexptime/exptime

  if (v) then print, 'e_star = ', median(e_star)

  ; e/pix from zodi
  dlat = (abs(elat)-90.)/90.
  vmag_zodi = 23.345 - 1.148*dlat^2.
  e_pix_zodi = 10.0^(-0.4*(vmag_zodi-22.8)) * 2.56D-3 * $
	geom_area * cos(!DPI * field_angle/180.) * omega_pix * exptime

  if (v) then print, 'vmag_zodi = ', median(vmag_zodi)
  if (v) then print, 'e_pix_zodi = ', median(e_pix_zodi)

  ; e/pix from background stars

  euler, elon, elat, glon, glat, select=5
;  if (v) then print, 'glon, glat = ', glon, glat

  dlon = glon
  q = where(dlon gt 180.) & if(q[0] ne -1) then dlon[q] = 360.-dlon[q]
  dlon = abs(dlon)/180.0D0
  p = bk_p[fov_ind,*]
  dlat = abs(glat)/90.0D0
  imag_bgstars = p[*,0] + $
	p[*,1]*(1.0-exp(-dlon/p[*,2])) + $
	p[*,3]*(1.0-exp(-dlat/p[*,4])) + $
	p[*,5]*sqrt(dlon*dlat)

  e_pix_bgstars = 10.0^(-0.4*imag_bgstars) * 1.7D6 * $
	geom_area * cos(!DPI * field_angle/180.) * omega_pix * exptime

  if (v) then print, 'imag_bgstars = ', median(imag_bgstars)
  if (v) then print, 'e_pix_bgstars = ', median(e_pix_bgstars)

  ; compute noise sources
  e_tot_sky = npix_aper*(e_pix_zodi + e_pix_bgstars) 

  e_bin = dblarr(n_elements(e_tot_sky))
  if(keyword_set(bin_sys)) then nbin=total(bin_sys) else nbin=0
; e/pix from companion
  if (nbin gt 0) then begin
    if (v) then print, 'Dealing with ', nbin, ' binaries'
    bins = where(bin_sys)
;    pix_sep = bin_sep[bins]/pix_scale  	; from arcsec to pixels
;    r = frac_off[fov_ind[bins],*]	; distance (in pixels)
;    di = frac_off[fov_ind[bins]+4,*]    ; imag attenuation
;    dibin = interpol(di, r, pix_sep)    ; interpolate over spline fit
    e_bin[bins] = bin_ph[bins] * $
        geom_area * cos(!DPI * field_angle/180.)* exptime * $
        reform(ph_star[npix_aper-1,*])/reform(ph_star[npix_max-1,*]) ; frac
    if (v) then print, 'e_pix_bin = ', median(e_bin[bins])
  endif
  
  e_tot = e_tot_sky + e_star + e_bin
  
  noise_star = sqrt(e_star) / e_tot
  noise_sky  = sqrt(e_tot_sky) / e_tot
  noise_ro   = sqrt(npix_aper*n_exposures)*e_pix_ro / e_tot
  noise_bin  = sqrt(e_bin) / e_tot
  noise_sys  = 0.0*noise_star + sys_limit/1d6/sqrt(exptime/3600.)

  dilution = e_tot / e_star
  ;if (al) then noise = sqrt(1.0/(e_star + e_tot_sky) + noise_sys^2.) else $
  noise = sqrt( noise_star^2. + noise_sky^2. + noise_ro^2. + noise_sys^2. + noise_bin^2. )
;  if (v) then begin
;     print, 'noise_star [ppm] = ', noise_star*1d6
;     print, 'noise_sky  [ppm] = ', noise_sky*1d6
;     print, 'noise_ro   [ppm] = ', noise_ro*1d6
;     print, 'noise_sys  [ppm] = ', noise_sys*1d6
;     print, 'noise      [ppm] = ', noise*1d6
;  endif

end

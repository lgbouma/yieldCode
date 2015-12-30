pro calc_noise_hack, $
;
; mandatory inputs
;
   imag, $                          ; apparent mag in Cousins I band
   exptime, $                       ; total exposure time in seconds
;
; mandatory outputs
;
   noise, $                         ; fractional rms (= 1/SNR)
;
; optional inputs
;
   teff=teff, $                     ; effective temperature in Kelvins
   elon=elon, elat=elat, $          ; ecliptic coordinates in degrees
   subexptime=subexptime, $         ; subexposure time (n_exp = exptime/subexptime)
   npix_aper = npix_aper, $         ; number of pixels in photometric aperture
   frac_aper = frac_aper, $         ; fraction of flux enclosed in photometric aperture
   e_pix_ro = e_pix_ro, $           ; variance in no. photons/pixel from readout noise
   geom_area = geom_area, $         ; geometric collecting area
   pix_scale = pix_scale, $         ; arcsec per pixel
   sys_limit = sys_limit, $         ; minimum uncertainty in 1 hr of data, in ppm
   verbose=verbose, $               ; request verbose output
;
; optional outputs
;
   noise_star=noise_star, $     ; noise from star counts alone
   noise_sky=noise_sky, $       ; noise from sky counts only
   noise_ro=noise_ro,$          ; noise from readout only
   noise_sys=noise_sys          ; noise from systematic limit only
;
;
  if (keyword_set(teff)) then teff=teff else teff=5000.
  if (keyword_set(elon)) then elon=elon else elon=0.0
  if (keyword_set(elat)) then elat=elat else elat=30.0
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.0
  if (keyword_set(e_pix_ro)) then e_pix_ro=e_pix_ro else e_pix_ro=10.
  if (keyword_set(geom_area)) then geom_area=geom_area else geom_area=73.0
  if (keyword_set(pix_scale)) then pix_scale=pix_scale else pix_scale=20.0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=0.0
  if (keyword_set(verbose)) then v=1 else v=0

  if (keyword_set(npix_aper)) then begin
     npix_aper=npix_aper
  endif else begin
     npix_aper=optimal_npix_hack(imag)
  endelse

  if (keyword_set(frac_aper)) then begin
     frac_aper=frac_aper
  endif else begin
     fwhm = 0.6
     r_aper = sqrt(npix_aper/!PI)
     frac_aper = 1.0 - (fwhm/(2.*!PI)) / ((fwhm/2.)^2 + r_aper^2.) ; from Roland 9 July 2013
  endelse
        
  omega_pix = pix_scale^2.
  n_exposures = exptime/subexptime

  if(v) then begin
     print, 'imag = ', imag
     print, 'exptime = ', exptime
     print, 'teff = ', teff
     print, 'elon = ', elon
     print, 'elat = ', elat
     print, 'npix_aper = ', npix_aper
     print, 'frac_aper = ', frac_aper
     print, 'subexptime = ', subexptime
     print, 'n_exposures = ', n_exposures
     print, 'e_pix_ro = ', e_pix_ro
     print, 'geom_area = ', geom_area
     print, 'pix_scale = ', pix_scale
     print, 'omega_pix = ', omega_pix
     print, 'sys_limit = ', sys_limit
  endif

; electrons from the star

  megaph_s_cm2_0mag = 1.6301336 + 0.14733937*(teff-5000.)/5000.
  e_star = 10.0^(-0.4*imag) * 1D6 * megaph_s_cm2_0mag * geom_area * exptime * frac_aper

  if (v) then print, 'e_star = ', e_star

; e/pix from zodi

  dlat = (abs(elat)-90.)/90.
  vmag_zodi = 23.345 - 1.148*dlat^2.
  e_pix_zodi = 10.0^(-0.4*(vmag_zodi-22.8)) * 2.39D-3 * geom_area * omega_pix * exptime

  if (v) then print, 'vmag_zodi = ', vmag_zodi
  if (v) then print, 'e_pix_zodi = ', e_pix_zodi

; e/pix from background stars

  euler, elon, elat, glon, glat, select=5
  if (v) then print, 'glon, glat = ', glon, glat

  dlat = abs(glat)/40.0D0
  dlon = glon
  q = where(dlon gt 180.) & if(q[0] ne -1) then dlon[q] = 360.-dlon[q]
  dlon = abs(dlon)/180.0D0
  p = [18.9733D0, 8.833D0, 4.007D0, 0.805D0]

  imag_bgstars = p[0] + p[1]*dlat + p[2]*dlon^(p[3])
  e_pix_bgstars = 10.0^(-0.4*imag_bgstars) * 1.7D6 * geom_area * omega_pix * exptime

  if (v) then print, 'imag_bgstars = ', imag_bgstars
  if (v) then print, 'e_pix_bgstars = ', e_pix_bgstars

; compute noise sources

  noise_star = sqrt(e_star) / e_star
  noise_sky  = sqrt(npix_aper*(e_pix_zodi + e_pix_bgstars)) / e_star
  noise_ro   = sqrt(npix_aper*n_exposures)*e_pix_ro / e_star
  noise_sys  = 0.0*noise_star + sys_limit/1d6/sqrt(exptime/3600.)

  noise = sqrt( noise_star^2. + noise_sky^2. + noise_ro^2. + noise_sys^2. )

  if (v) then begin
     print, 'noise_star [ppm] = ', noise_star*1d6
     print, 'noise_sky  [ppm] = ', noise_sky*1d6
     print, 'noise_ro   [ppm] = ', noise_ro*1d6
     print, 'noise_sys  [ppm] = ', noise_sys*1d6
     print, 'noise      [ppm] = ', noise*1d6
  endif

end

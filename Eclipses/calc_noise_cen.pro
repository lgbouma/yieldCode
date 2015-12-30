pro calc_noise_cen, $
;
; mandatory inputs
;
   ph_star, $                       ; ph/s/cm^2 from (npixels x nstars)
   ph_dil, $                        ; ph/s/cm^2 per pixel
   ph_beb, $                        ; ph/s/cm^2 per pixel from beb
                           ; ph/s/cm^2 from other target stars (subtracted)
   bebind, $
   dur1, $                      ; BEB/HEB dur 
   dur2, $                      ; BEB/HEB dur
   dep1, $                      ; BEB/HEB depth
   dep2, $                      ; BEB/HEB depth
   xx, $
   yy, $
   sind, $
   e_pix_ro, $                      ; read noise per subexposure in e-
   sys_limit, $ 	            ; noise floor (ppm)
;
; mandatory outputs
;
   xcen, $
   ycen, $
   xcenshift1, $       ; output centroid shifts for BEB/HEBs
   xcenshift2, $       ; 
   ycenshift1, $       ;
   ycenshift2, $       ;
   xcennoise1, $
   xcennoise2, $
   ycennoise1, $
   ycennoise2, $

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
   noise_sky=noise_sky, $           ; noise from sky counts only
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
  
  nstar = n_elements(dur1)
  if (bebind[0] ne -1) then nbeb = n_elements(beb_ind) $
  else nbeb = 0

  if (keyword_set(npix_aper)) then begin
     npix_aper=npix_aper
  endif else begin
     npix_aper=3+intarr(n_elements(ph_star))
  endelse
 
  if (v) then print, 'npix_aper = ', npix_aper
        
  exptime1 = dur1*24.0*3600.
  exptime2 = dur2*24.0*3600.
  exptime = 3600.
  
  omega_pix = aspix^2.
  n_exposures = exptime/subexptime
  rn_pix   = sqrt(n_exposures)*e_pix_ro

  ; truncate indices
  npix_sind = sind[0:npix_aper-1,*]

  if (keyword_set(noise_cr)) then noise_cr = reform(noise_cr[*,npix_aper-1]) else noise_cr = 0.0 
  
  if (v) then print, 'e_star = ', median(e_star)

  ; e/pix from zodi
  e_pix_zodi = zodi_ph * geom_area * cos(!DPI * field_angle/180.) * exptime

  if (v) then print, 'vmag_zodi = ', median(vmag_zodi)
  if (v) then print, 'e_pix_zodi = ', median(e_pix_zodi)


  xcenshift1 = fltarr(nstar) ; output vectors
  xcenshift2 = fltarr(nstar)
  ycenshift1 = fltarr(nstar)
  ycenshift2 = fltarr(nstar)
  xcennoise = fltarr(nstar)
  ycennoise = fltarr(nstar)
  xcen = fltarr(nstar)
  ycen = fltarr(nstar)
  for ii=0,nstar-1 do begin
  ; electrons from the star
    this_estar0 = ph_star[*,ii] * geom_area * cos(!DPI * field_angle[ii]/180.) * exptime
    this_estar1 = this_estar0*(1.0-dep1[ii])
    this_estar2 = this_estar0*(1.0-dep2[ii])
    this_edil  = ph_dil[*,ii]  * geom_area * cos(!DPI * field_angle[ii]/180.) * exptime
    this_ebeb0 = ph_beb[*,ii]  * geom_area * cos(!DPI * field_angle[ii]/180.) * exptime
    this_ebeb1 = this_ebeb0*(1.0-dep1[ii])
    this_ebeb2 = this_ebeb0*(1.0-dep2[ii])
    ;this_etgt  = ph_tgt[*,ii] * geom_area * cos(!DPI * field_angle[ii]/180.) * exptime
    this_sind = npix_sind[*,ii]                ; sorting indices
    this_epix_all = this_estar0 + this_edil + this_ebeb0 + e_pix_zodi[ii]  ; electrons per pixel, unsorted
    if total(ii eq bebind) then begin
      this_epix0 = this_estar0 + this_edil + this_ebeb0    ; electrons per pixel, unsorted
      this_epix1 = this_estar0 + this_edil + this_ebeb1    ; electrons per pixel, unsorted
      this_epix2 = this_estar0 + this_edil + this_ebeb2    ; electrons per pixel, unsorted
    endif else begin
      this_epix0 = this_estar0 + this_edil + this_ebeb0    ; electrons per pixel, unsorted
      this_epix1 = this_estar1 + this_edil + this_ebeb0    ; electrons per pixel, unsorted
      this_epix2 = this_estar2 + this_edil + this_ebeb0    ; electrons per pixel, unsorted
    endelse
    this_epix_all_sind = this_epix_all[this_sind]      ; sorted electrons per pixel
    this_epix0_sind = this_epix0[this_sind]      ; sorted electrons per pixel
    this_estar0_sind = this_estar0[this_sind]      ; sorted electrons per pixel
    this_epix1_sind = this_epix1[this_sind]      ; sorted electrons per pixel
    this_epix2_sind = this_epix2[this_sind]      ; sorted electrons per pixel
    
    this_etot_all = total(this_epix_all_sind)          ; total electrons
    this_etot0 = total(this_epix0_sind)          ; total electrons
    this_estartot = total(this_estar0_sind)          ; total electrons
    this_etot1 = total(this_epix1_sind)          ; total electrons
    this_etot2 = total(this_epix2_sind)          ; total electrons
    ;this_ntot = sqrt(this_etot_all + npix_aper*rn_pix^2. + cr_noise[ii]^2.)/this_estartot      ; total noise 
    ;this_ntot1 = sqrt(this_etot1 + npix_aper*rn_pix^2.)/this_etot1      ; total noise 
    ;this_ntot2 = sqrt(this_etot2 + npix_aper*rn_pix^2.)/this_etot2      ; total noise 
    this_epix_noise = sqrt(this_epix_all_sind + rn_pix^2. + cr_noise[ii]^2.)/this_etot0 ; noise per pixel

    ; calculate centroid
    xc0 = this_epix0_sind*xx[this_sind] ; counts*x
    yc0 = this_epix0_sind*yy[this_sind] ; counts*y
    xcen[ii] = total(xc0)/this_etot0    ; x centroid out-of-eclipse
    ycen[ii] = total(yc0)/this_etot0    ; y centroid out-of-eclipse
    xc1 = this_epix1_sind*xx[this_sind] ; in eclipse 1
    yc1 = this_epix1_sind*yy[this_sind]
    xcenshift1[ii] = total(xc1)/this_etot1 - xcen[ii]
    ycenshift1[ii] = total(yc1)/this_etot1 - ycen[ii]
    xc2 = this_epix2_sind*xx[this_sind] ; in eclipse 2
    yc2 = this_epix2_sind*yy[this_sind]
    xcenshift2[ii] = total(xc2)/this_etot2 - xcen[ii]
    ycenshift2[ii] = total(yc2)/this_etot2 - ycen[ii]
    ; calculate centroid noise
    xcn = this_epix_noise*(xx[this_sind] - xcen[ii])
    ycn = this_epix_noise*(yy[this_sind] - ycen[ii])
    xcennoise[ii] = sqrt(total(xcn^2.))
    ycennoise[ii] = sqrt(total(ycn^2.))
   ; if (npix_aper gt 50) then stop
  end
  xcennoise1 = xcennoise*sqrt(exptime/exptime1)
  xcennoise2 = xcennoise*sqrt(exptime/exptime2)
  ycennoise1 = ycennoise*sqrt(exptime/exptime1)
  ycennoise2 = ycennoise*sqrt(exptime/exptime2)

end

pro eclip_survey, n_segs, fov, star, offset=offset
;+
; NAME: eclip_survey
; PURPOSE: figure out the number of pointings and field angles each star (eclipse?)
; 	on a given tile gets.
; INPUTS: 
;	n_segs: number of observing segments (13) per hemisphere
;	fov: field of view for each camera (23deg)
;	star: ECLIPSE object with planets and eclipses assigned. 
;	In the standard run, this sub-routine is passed *eclip*
;	offset: off the ecliptic (nominal: passed "skirt" of 6deg)
; OUTPUTS: 
;	star object, with npointings and field angles added (needed
;	for subsequent "observing")
;-
  if (keyword_set(offset)) then offset=offset else offset=0.0
  ; LB 16/01/29: originally, this was labelled stars. It's not a "star". You're passing an
  ; eclip object. (??)
  print, 'Surveying ', n_elements(star), ' eclipses.'
  n_cams = 4
  ccd_pix = 4096.
  gap_pix = 2.0/0.015
  pix_scale = fov*3600./float(ccd_pix+gap_pix)
  pix_scale_deg = fov/float(ccd_pix+gap_pix)
  ccd_ctr = [1.0,1.0]*float(ccd_pix+gap_pix)/2.0
  delt = [1.0,1.0]*pix_scale_deg

  elat_cams = (indgen(n_cams)+0.5)*(fov) + offset
  elon_cams = intarr(n_cams)
  elon_segs = 360.0*(findgen(n_segs))/float(n_segs)

  old_onchip = intarr(n_elements(star))
 
  for hemi=-1,1,2 do begin
    for cam=0,n_cams-1 do begin
      for seg=0,n_segs-1 do begin
         make_astr, astr, crpix=ccd_ctr, $
                        crval=[elon_segs[seg],float(hemi)*elat_cams[cam]], $
                        delta=delt
         ad2xy, star.coord.elon, star.coord.elat, astr, x, y
         onchip =   ( $ ; Upper left
                       ((x gt 0.0) and $
                        (x lt (ccd_pix-gap_pix)/2.-1.) and $
                        (y gt 0.0) and $
                        (y lt (ccd_pix-gap_pix)/2.-1.)) $
                       or $ ; Upper Right
                       ((x gt (ccd_pix+gap_pix)/2.-1.) and $
                        (x lt (ccd_pix+gap_pix)-1.) and $
                        (y gt 0.0) and $
                        (y lt (ccd_pix-gap_pix)/2.-1.)) $
                       or $ ; Lower Right
                       ((x gt (ccd_pix+gap_pix)/2.-1.) and $
                        (x lt (ccd_pix+gap_pix)-1.) and $
                        (y gt (ccd_pix+gap_pix)/2.-1.) and $
                        (y lt (ccd_pix+gap_pix)-1.)) $
                       or $ ; Lower Left
                       ((x gt 0.0) and $
                        (x lt (ccd_pix-gap_pix)/2.-1.) and $
                        (y gt (ccd_pix+gap_pix)/2.-1.) and $
                        (y lt (ccd_pix+gap_pix)-1.)))
         if (total(onchip gt 0)) then begin
           prev_npointings = star.npointings
           ; only count consecutive observations:
           ; so it's a new observation OR the previous sector got it
           incr = onchip and ((not prev_npointings) or old_onchip) 
           new_npointings = prev_npointings + incr
           chipind = where(incr gt 0)
           r = sqrt((x[chipind]-ccd_ctr[0])^2. + (y[chipind]-ccd_ctr[1])^2.) 
           star[chipind].npointings = new_npointings[chipind]
           star[chipind].coord.fov_r = r/float(new_npointings[chipind]) + $
		float(prev_npointings[chipind])*star[chipind].coord.fov_r/float(new_npointings[chipind])
           old_onchip = onchip
           print, 'Hemi=',hemi,' Seg=',seg,' Cam=',cam, $
                ' ELon=',elon_segs[seg],' ELat=',float(hemi)*elat_cams[cam], $
                ' Nstars=',n_elements(onchip)
        end
      endfor
    endfor
  endfor
  ; Determine FOV index
  r = star.coord.fov_r
  field_angle = r * fov / (ccd_pix+gap_pix)
  fov_ind = intarr(n_elements(field_angle))
  ; 10/25/2015 new field angles: [0,7.14,14.26,20]
  ; old field angles were [0,6,12,17] deg
  fov_ind[where((field_angle ge 0.14875*fov) and $ ; was 0.125 
                (field_angle lt 0.44583*fov))] = 1 ; was 0.375
  fov_ind[where((field_angle ge 0.44583*fov) and $
                (field_angle lt 0.71375*fov))] = 2 ; was 0.604
  fov_ind[where( field_angle ge 0.71375*fov)] = 3
  star.coord.field_angle = field_angle
  star.coord.fov_ind=fov_ind

end

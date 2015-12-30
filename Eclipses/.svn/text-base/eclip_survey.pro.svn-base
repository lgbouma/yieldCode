PRO eclip_survey, n_segs, fov, star, offset=offset

  if (keyword_set(offset)) then offset=offset else offset=0.0
  print, 'Surveying ', n_elements(star), ' stars.'
  ;fov = 23.0
  ;n_segs = 13
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
           ;print, 'Hemi=',hemi,' Seg=',seg,' Cam=',cam, $
                ;' ELon=',elon_segs[seg],' ELat=',float(hemi)*elat_cams[cam], $
                ;' Nstars=',n_elements(onchip)
        end
      endfor
    endfor
  endfor
  ; Determine FOV index
  r = star.coord.fov_r
  field_angle = r * fov / (ccd_pix+gap_pix)
  fov_ind = intarr(n_elements(field_angle))
  fov_ind[where((field_angle ge 0.125*fov) and $ ; was 104
                (field_angle lt 0.375*fov))] = 1   ; was 365 
  fov_ind[where((field_angle ge 0.375*fov) and $
                (field_angle lt 0.604*fov))] = 2 ; was  591
  fov_ind[where( field_angle ge 0.604*fov)] = 3
  star.coord.field_angle = field_angle
  star.coord.fov_ind=fov_ind

end

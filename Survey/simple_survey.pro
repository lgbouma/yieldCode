PRO simple_survey, n_segs, fov, infile, outfile, offset=offset

  restore, infile
  if (keyword_set(offset)) then offset=offset else offset=0.0
  print, 'Surveying ', n_elements(star), ' stars.'
  ;fov = 23.0
  ;n_segs = 13
  n_cams = 4
  ccd_pix = 4096
  PIX_SCALE = fov*3600./float(ccd_pix)
  pix_scale_deg = fov/float(ccd_pix)
  ccd_ctr = [1.0,1.0]*float(ccd_pix)/2.0
  delt = [1.0,1.0]*pix_scale_deg

  elat_cams = (indgen(n_cams)+0.5)*(fov) + offset
  elon_cams = intarr(n_cams)
  elon_segs = 360.0*(findgen(n_segs))/float(n_segs)
 
 for hemi=-1,1,2 do begin
    for seg=0,n_segs-1 do begin
      for cam=0,n_cams-1 do begin
         make_astr, astr, crpix=ccd_ctr, $
                        crval=[elon_segs[seg],float(hemi)*elat_cams[cam]], $
                        delta=delt
         ad2xy, star.coord.elon, star.coord.elat, astr, x, y
         onchip = where((x gt 0.0) and $
                        (x lt ccd_pix-1) and $
                        (y gt 0.0) and $
                        (y lt ccd_pix-1))
         r = sqrt((x[onchip]-ccd_ctr[0])^2 + (y[onchip]-ccd_ctr[1])^2) 
         prev_npointings = star[onchip].npointings
         new_npointings = prev_npointings+1
         star[onchip].npointings = new_npointings
         star[onchip].coord.fov_r = r/float(new_npointings) + $
		float(prev_npointings)*star[onchip].coord.fov_r/float(new_npointings)
         print, 'Hemi=',hemi,' Seg=',seg,' Cam=',cam, $
                ' ELon=',elon_segs[seg],' ELat=',float(hemi)*elat_cams[cam], $
                ' Nstars=',n_elements(onchip)
      endfor
    endfor
  endfor

  save, star, filen=outfile
  return
end

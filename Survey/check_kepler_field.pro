pro check_kepler_field

  ra = 19.0 + 22.0/60.0 + 40./3600.
  ra = ra*15.
  dec = 44.0 + 30./60.
  
  euler, ra, dec, elon, elat, 3
  
  restore, 'survey_fov24_15segments.sav'
  si = size(Number_of_Pointings)
  n_pix_lon = si[1] & n_pix_lat = si[2]

  i_mask = floor(n_pix_lon * elon/360.0)
  j_mask = floor(n_pix_lat * (elat+90.0)/180.0)
  npointings = Number_of_Pointings[i_mask,j_mask]
  print, npointings

end

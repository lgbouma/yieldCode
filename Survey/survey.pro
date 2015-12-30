pro survey

; specify the FOV in degrees in a side (assumes square)
  fov = 23.0
;  fov = 25.7

; specify the centers of each camera field relative to the middle of
; the bottom edge of camera 0
;  delta_lon_camera = [0,0,0,0]           ; put all the cameras at a common ecliptic longitude
;  delta_lat_camera = ([0,1,2,3]+0.5)*fov ; stack them vertically in ecliptic latitude
  delta_lon_camera = [0,0,0,0]           ; put all the cameras at a common ecliptic longitude
  delta_lat_camera = ([0,1,2,3]+0.5)*(fov) ; stack them vertically in ecliptic latitude

; the "0" camera is the reference camera, and the coordinates
; are the coordinates relative to the middle of the bottom edge of camera 0
  lon_pointing = 360.0*dindgen(13)/13. ; evenly spaced in ecliptic longitude
;  lon_pointing = 360.0*dindgen(15)/15. ; evenly spaced in ecliptic longitude

  lat_pointing = 0*lon_pointing        ; bottom of bottom camera is at ecliptic lat=0

; now set the resolution of the maps.
; the 'all sky' maps will be 2NxN degrees, with each pixel
; representing a long extent of resol_allsky degrees and a lat extent
; of resol_allsky deg, centered on that point.
; The solid angle is obviously variable (less near the poles)

  resol_allsky = 0.25             ; degrees per pixel of the all sky map
  resol_camera = resol_allsky/10. ; degrees per pixel of the camera images

;;;; now get busy

  n_cameras = n_elements(delta_lat_camera)
  n_pointings = n_elements(lat_pointing)
  for i=0,n_cameras-1 do print, 'Camera number, dlon, dlat = ', i, delta_lon_camera[i], delta_lat_camera[i]
  for i=0,n_pointings-1 do print, 'Pointing number, lon, lat = ', i, lon_pointing[i], lat_pointing[i]

  n_pix_lon = round(360.0/resol_allsky); e.g. for 1 degree maps, use 360 pixels representing 0 to 359 deg
  n_pix_lat = round(180.0/resol_allsky)
  
  Number_of_Pointings = intarr(n_pix_lon, n_pix_lat)

  n_pix_x = round(fov/resol_camera)
  n_pix_y = n_pix_x

  gridmaker,n_pix_x,n_pix_y,1,1,ix,iy, type='float'
                                ; ix and iy camera-image arrays giving
                                ; the x and y indices of each pixel in
                                ; a camera image
   
  for hemisphere=-1,1,2 do begin

     delta_lat_camera = -delta_lat_camera ; this flipping has the effect of doing the south first, then the north

     for i_pointing=0,n_pointings-1 do begin

        This_Pointing_Bitmask = 0 * Number_of_Pointings
        
        for i_camera=0,n_cameras-1 do begin
           
           print, 'Hemi = ', hemisphere, ' Camera = ', i_camera, ' Pointing = ', i_pointing

           This_Camera_Bitmask = 0 * Number_of_Pointings

           lon_center = lon_pointing[i_pointing] + delta_lon_camera[i_camera]
           lat_center = lat_pointing[i_pointing] + delta_lat_camera[i_camera]

                                ; now convert the pixel coordinates of the camera image into lat and
                                ; lon on the sky using some IDLASTRO tools

           make_astr, astr, $
                      crpix=[n_pix_x/2,n_pix_y/2],$       ; indices of the reference pixel (center of image)
                      crval=[lon_center,lat_center], $    ; lon and lat of the reference pixel (center of camera pointing)
                      delta=[resol_camera,resol_camera],$ ; degrees per pixel at reference pixel location
                      ctype=['RA---TAN','DEC--TAN']

           mkhdr, h, 0*This_Pointing_Bitmask ; makes an empty FITS header
           putast, h, astr                   ; puts the astrometry info in that header
           xyad, h, ix, iy, lon, lat         ; now makes camera-image arrays with the lat and lon of each pixel

           oops = where(lon ge 360.,nk) & if nk gt 0 then lon[oops] = lon[oops] - 360.
           oops = where(lon lt 0.,nk) & if nk gt 0 then lon[oops] = lon[oops] + 360.

                                ; for each pixel in the camera image,
                                ; find the corresponding pixel in the
                                ; Bitmask, and set it to 1

           for i=0L,n_elements(lon)-1 do begin
              i_mask = round((n_pix_lon-1) * lon[i]/360.0)
              j_mask = round((n_pix_lat-1) * (lat[i]+90.0)/180.0)
              This_Camera_Bitmask[i_mask,j_mask] = (This_Camera_Bitmask[i_mask,j_mask] > 1)
           endfor

           This_Pointing_Bitmask = This_Pointing_Bitmask + This_Camera_Bitmask ; assumes FOVs are non overlapping at a given time

        endfor; done with this pointing

        Number_of_Pointings = Number_of_Pointings + This_Pointing_Bitmask

     endfor

  endfor

  save, Number_of_Pointings, filen='survey.sav'

  plot_survey

end

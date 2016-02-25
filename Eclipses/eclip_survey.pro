pro eclip_survey, fov, eclip, fCamCoord
;+
; NAME: eclip_survey
; PURPOSE: figure out the number of pointings and field angles each ECLIPSE (not star)
; 	on a given tile gets.
;	Quite similar to `starSurvey.pro` in cameraPointings, but here CCDs are more complex.
; INPUTS: 
;	n_segs: integer number of observing segments (13) per hemisphere.
;	fov: field of view for each camera (24deg).
;	eclip: ECLIPSE object with planets and eclipses assigned. 
;	fCamCoord: string pointing to file name of data file with camera
;		pointing info (pointing #, cam #, elat, elong for cameras)
; OUTPUTS: 
;	eclip object, with npointings and field angles added (needed
;	for subsequent "observing")
; COMMENTS:
;   Nomenclature: "pointing" == one grouping of 4 cameras. There are 26 
;   pointings. Let "camPointing" == one camera pointing, there are 26*4=104 (2yr).
;-

	nEcl = n_elements(eclip)
	print, 'Surveying', nEcl, ' possible eclipses over all trials (getting #pntgs and field angles).'

	MT = 'I,I,F,F' ; get camera pointing info for observing specification 
	READCOL, fCamCoord, F=FMT, pointingNumber, camNumber, camElat, camElong 
	nPointings = MAX(pointingNumber)+1 ; 26 pntgs 2yr, 52 4yr
	nCams = 4
    nCamPointings = nPointings * nCams
    assert, nCamPointings eq N_ELEMENTS(pointingNumber), $
	    'Need as many camPointings as elements passed in camPointingFile'

	;CCD info
	ccdPix = 4096.
	gapPix = 2.0/0.015	; \approx 133
	pixScale = fov*3600./float(ccdPix+gapPix)  ; pixel scale in arcseconds (about 20 arcsec / pixel)
	pixScaleDeg = fov/float(ccdPix+gapPix)    ; pixel scale in degrees (arcsec*3600)
	ccdCtr = [1.0,1.0]*float(ccdPix+gapPix)/2.0 ; center of CCD image in arcseconds 
	delt = [1.0,1.0]*pixScaleDeg		

	old_onChip = intarr(n_elements(eclip))

	; Below this, we determine whether eclipses/transits are observed, and npointings
	; they receive if so. This is back to geometry, in the style of /camPointings/*, but
	; fancier b/c we want CCD gaps and field angles.

	for i=0, nCamPointings-1 do begin ; loop over camera pointings
        make_astr, astr, $       ; build astrometry structure from input params
        crpix=ccdCtr, $      ; indices of reference pixel (center of image)
        crval=[camElong[i], camElat[i]], $ ; lon & lat of reference pixel
        delta=delt, $        ; degrees per pixel at reference pixel location
        ctype=['RA---TAN', 'DEC--TAN']  ; coordt type, RA & DEC is fine.

		ad2xy, eclip.coord.elon, eclip.coord.elat, astr, x, y

		; 16/02/02 LB: this logic seems fine. "-1." probably shouldn't be there, but
		; negligible effect for 4k*4k pixels.
		onChip =   ( $ ; Lower left
					   ((x gt 0.0) and $
						(x lt (ccdPix-gapPix)/2.-1.) and $
						(y gt 0.0) and $
						(y lt (ccdPix-gapPix)/2.-1.)) $
					   or $ ; Lower Right
					   ((x gt (ccdPix+gapPix)/2.-1.) and $
						(x lt (ccdPix+gapPix)-1.) and $
						(y gt 0.0) and $
						(y lt (ccdPix-gapPix)/2.-1.)) $
					   or $ ; Upper Right
					   ((x gt (ccdPix+gapPix)/2.-1.) and $
						(x lt (ccdPix+gapPix)-1.) and $
						(y gt (ccdPix+gapPix)/2.-1.) and $
						(y lt (ccdPix+gapPix)-1.)) $
					   or $ ; Upper Left
					   ((x gt 0.0) and $
						(x lt (ccdPix-gapPix)/2.-1.) and $
						(y gt (ccdPix+gapPix)/2.-1.) and $
						(y lt (ccdPix+gapPix)-1.)))
		
		;Remember you're calling this subroutine for every tile. Thus there will be
		;many pointings which do not see that tile. Logic below this is Peter's.	
		if TOTAL(onChip) gt 0 then begin
	    prev_npointings = eclip.npointings
	    ; only count consecutive observations:
	    ; so it's a new observation OR the previous sector got it
	    incr = onChip and ((not prev_npointings) or old_onChip) 
	    new_npointings = prev_npointings + incr
	    chipind = where(incr gt 0)
	    r = sqrt((x[chipind]-ccdCtr[0])^2. + (y[chipind]-ccdCtr[1])^2.) 
	    eclip[chipind].npointings = new_npointings[chipind]
	    eclip[chipind].coord.fov_r = r/float(new_npointings[chipind]) + $
	    float(prev_npointings[chipind])*eclip[chipind].coord.fov_r/float(new_npointings[chipind])
	    old_onChip = onChip
	   	endif
	endfor

	; Determine FOV index
	r = eclip.coord.fov_r
	field_angle = r * fov / (ccdPix+gapPix)
	fov_ind = intarr(n_elements(field_angle))
	; 10/25/2015 PS: new field angles: [0,7.14,14.26,20]
	; old field angles were [0,6,12,17] deg
	fov_ind[where((field_angle ge 0.14875*fov) and $ ; was 0.125 
				(field_angle lt 0.44583*fov))] = 1 ; was 0.375
	fov_ind[where((field_angle ge 0.44583*fov) and $
				(field_angle lt 0.71375*fov))] = 2 ; was 0.604
	fov_ind[where( field_angle ge 0.71375*fov)] = 3
	eclip.coord.field_angle = field_angle
	eclip.coord.fov_ind=fov_ind

end

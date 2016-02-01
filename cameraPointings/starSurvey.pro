pro starSurvey, camPointingFile, catFile
;+
; NAME: starSurvey
; PURPOSE: given an input list of camera pointings and a catalog of 
; 	stars (e.g., postage stamps), calculate how many pointings
; 	each star gets.
; INPUTS: 
;	camPointingFile: text file (e.g., "nominalCamPointings.dat")
;	formatted in columns as:
;	- pointing number (nominal run goes from 0 to 23)
;	- camera number (0, 1, 2, 3)
;	- ecliptic latitude (-90 to 90 deg)
; 	- ecliptic longitude (0 to 360)
;	catFile: .sav file, (e.g., "psStruct.sav") containing object of 
;	type pointingStruct, output from genPSstarStruct.pro. It's a 
;	catalog of PS stars, with tileNum, PSnum, elong, and elat for the star. 
;	CAUTION: it does NOT have any non-zero pointing numbers yet.
; OUTPUTS: 
;	Catalog, with star.npointings calculated. 
; COMMENTS:
;	Originally written to answer the question: "how many stars that are 
;	assigned as postage stamps are actually observed?" (since these are
;	different things in the Sullivan+ 2015 simulations). Modelled it from
; 	eclip_survey.pro in the Eclipses/* dir.
;	
; 	Nomenclature: "pointing" == one grouping of 4 cameras. There are 26 
; 	pointings. Let "camPointing" == one camera pointing, there are 26*4=104.
;-

	;camPointingFile = 'nominalCamPointings.dat'
	;catFile = 'psStruct.sav'

	FMT= 'I,I,F,F'
	READCOL, camPointingFile, F=FMT, pointingNum, camNum, elatCams, elongCams

	RESTORE, catFile 	
	nPS = N_ELEMENTS(cat)
	PRINT, 'Total ', nPS, ' possible postage stamps', STRING(10B)

	fov = 24.	; degrees
	nSeg = 13	; integer number of observing segments per hemisphere (presumably 13 is fixed)
	nPointings = MAX(pointingNum)+1 ; 26 pointings (but generated in python counting)
	nCams = 4
	nCamPointings = nPointings * nCams
	assert, nCamPointings eq N_ELEMENTS(pointingNum), $ 
		'Need as many camPointings as elements passed in camPointingFile'

	ccdPix = 4096.
	gapPix = 2.0/0.015	; which is \approx 133
	pixScale = fov*3600./FLOAT(ccdPix+gapPix)  ; pixel scale in arcseconds (about 20 arcsec / pixel)
	pixScaleDeg = fov/FLOAT(ccdPix+gapPix)     ; pixel scale in degrees (arcsec*3600)
	ccdCtr = [1.0,1.0]*FLOAT(ccdPix+gapPix)/2.0 ; center of CCD image in arcseconds 
	delt = [1.0,1.0]*pixScaleDeg		
	
	for i=0, nCamPointings-1 do begin ; loop over camera pointings
		make_astr, astr, $ 		 ; build astrometry structure from input params
			crpix=ccdCtr, $		 ; indices of reference pixel (center of image)
			crval=[elongCams[i], elatCams[i]], $ ; lon & lat of reference pixel
			delta=delt, $		 ; degrees per pixel at reference pixel location
			ctype=['RA---TAN', 'DEC--TAN']	; coordt type, RA & DEC is fine.

		ad2xy, cat.coord.elon, cat.coord.elat, astr, x, y	
		; x and y (units of arcsec) of stars in CCD image for this camera orientation

		; Is the star (known x,y coords) on CCD chip?
		onChip =   ($   ; Upper left
				   ((x gt 0.0) and $
					(x lt (ccdPix-gapPix)/2.-1.) and $
					(y gt 0.0) and $
					(y lt (ccdPix-gapPix)/2.-1.)) $
				   or $ ; Upper Right
				   ((x gt (ccdPix+gapPix)/2.-1.) and $
					(x lt (ccdPix+gapPix)-1.) and $
					(y gt 0.0) and $
					(y lt (ccdPix-gapPix)/2.-1.)) $
				   or $ ; Lower Right
				   ((x gt (ccdPix+gapPix)/2.-1.) and $
					(x lt (ccdPix+gapPix)-1.) and $
					(y gt (ccdPix+gapPix)/2.-1.) and $
					(y lt (ccdPix+gapPix)-1.)) $
				   or $ ; Lower Left
				   ((x gt 0.0) and $
					(x lt (ccdPix-gapPix)/2.-1.) and $
					(y gt (ccdPix+gapPix)/2.-1.) and $
					(y lt (ccdPix+gapPix)-1.)))
				
		assert, TOTAL(onChip) gt 0, 'Your pointings are off b/c some do not see *anything*'

		cat.nPointings += onChip

	    PRINT, 'nCamPointing', i, ' pointingNum', pointingNum[i], $
			  ' camNum', camNum[i], ' eLongCam', eLongCams[i], ' eLatCam', eLatCams[i], $
			  ' nStarsOnChip=', N_ELEMENTS(onchip)

	endfor
	
	;SAVE, cat, FILENAME='psWithPointings.sav'
	SAVE, cat, FILENAME='tilesMedWithNominalPointings.sav'
	tilePrint = 1
	tStr = ' tiles'
	psStr = ' postage stamps'
	if tilePrint then str=tStr else str=psStr
	
	PRINT, STRING(10B), '==========', STRING(10B), $
		'Total # assigned' + str, N_ELEMENTS(cat.nPointings), $
		'Total #' + str + ' with 0 pointings', N_ELEMENTS(WHERE(cat.nPointings eq 0)), $
		'Total #' + str + ' with >0 pointings', N_ELEMENTS(WHERE(cat.nPointings gt 0)), $
	    STRING(10B), '=========='
end

pro starSurvey, camPointingFile, catFile
;+
; NAME: starSurvey
; PURPOSE: given an input list of camera pointings and a catalog of 
; 	stars (e.g., postage stamps), calculate how many pointings
; 	each star gets.
; INPUTS: 
;	camPointingFile: text file formatted in columns as
;	- pointing number (nominal run goes from 0 to 23)
;	- camera number (0, 1, 2, 3)
;	- ecliptic latitude (-90 to 90 deg)
; 	- ecliptic longitude (0 to 360)
;	catFile: .sav file, containing object of type pointingStruct, 
;	output from genPSstarStruct.pro. 
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

	oldOnChip = INTARR(nPS)	; zeros array with size = # postage stamps
	
	for i=0, nCamPointings-1 do begin ; loop over camera pointings
		make_astr, astr, $ 		 ; build astrometry structure from input params
			crpix=ccdCtr, $		 ; indices of reference pixel (center of image)
			crval=[elongCams[i], elatCams[i]], $ ; lon & lat of reference pixel
			delta=delt, $		 ; degrees per pixel at reference pixel location
			ctype=['RA---TAN', 'DEC--TAN']	; coordt type, RA & DEC is fine.

		ad2xy, cat.coord.elon, cat.coord.elat, astr, x, y	
		; these x and y coordts are x and y (units of arcsec) of CCD image

		; Is the star on CCD chip?
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
				
		if (TOTAL(onChip gt 0)) then begin
		   nPointingsPrev = cat.nPointings
		   ; only count consecutive observations:
		   ; so it's a new observation OR the previous sector got it
		   incr = onChip and ((not nPointingsPrev) or oldOnChip) 
		   nPointingsNew = nPointingsPrev + incr
		   chipInd = WHERE(incr gt 0)
		   r = SQRT((x[chipInd]-ccdCtr[0])^2. + (y[chipInd]-ccdCtr[1])^2.) 
		   cat[chipInd].nPointings = nPointingsNew[chipInd]
		   oldOnChip = onChip
		   PRINT, 'nCamPointing', i, ' pointingNum', pointingNum[i], $
			      ' camNum', camNum[i], ' eLongCam', eLongCams[i], ' eLatCam', eLatCams[i], $
				  ' nStarsOnChip=', N_ELEMENTS(onchip)
		end
	endfor
	SAVE, cat, FILENAME='psWithPointings.sav'
	PRINT, STRING(10B), '==========', STRING(10B), $
		'Total # assigned postage stamps', N_ELEMENTS(cat.nPointings), $
		'Total # postage stamps with 0 pointings', N_ELEMENTS(WHERE(cat.nPointings eq 0)), $
		'Total # postage postage with >0 pointings', N_ELEMENTS(WHERE(cat.nPointings gt 0)), $
	    STRING(10B), '=========='
end

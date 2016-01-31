PRO genTileCoords
; NAME: genTileCoords
; PURPOSE: generate a list of coordinates, indexed by healpix tile, that give
; 	mean ecliptic latitude and longitude of that tile. This will then be used
;	with genPSstarStruct.pro and starSurvey.pro to calculate how many pointings
;	the center of each tile gets. This will be used to skip tiles that (at 1st
;	order) receive 0 pointings.
; The check on how accurate this is will be whether we get ~170k postage stamps
; by imposing the skips or not.

	fpath = '/home/luke/Dropbox/tessSimulation/trilegal/'
	fnumsPath = '/home/luke/Dropbox/tessSimulation/tesscode/Eclipses/'
	fnums = mrdfits(fnumsPath + 'fnums.fits')
	nTiles = N_ELEMENTS(fnums)

	bigNumber = 5e6
	u = RANDOMU(seed, bigNumber)
	v = RANDOMU(seed, bigNumber)
	phi = 2.*!dpi*u
	theta = ACOS(2.*v-1.)
	ang2pix_ring, 16, theta, phi, ipring

	for i=0, nTiles-1 do begin
		allowedInd = WHERE(ipring eq fnums[i])
		
		glon = phi[allowedInd]*180./!dpi
		glat = (theta[allowedInd]-!dpi/2.)*180./!dpi

		euler, glon, glat, elong, elat, select=6

		avgLat = MEAN(elat)
		avgLong = MEAN(elong)
		
		PRINT, 'tileNum:', fnums[i], ' psNum: 1', ' eclLong:', avgLong, ' eclLat:', avgLat
	endfor
		
END

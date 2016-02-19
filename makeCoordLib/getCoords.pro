PRO getCoords
; See README in this coordLib dir for explanation of what's happening here.

	fnums = MRDFITS('fnums.fits')
	TIC
	;bigNumber = 2e7 ; 2e7 produces 1 gb lists
	bigNumber = 2e7*31 ; 2e7 produces 6593 per tile, want at least 2e5 per tile
	u = RANDOMU(seed, bigNumber) ; uniform distribution over 0 to 1
	v = RANDOMU(seed, bigNumber)
	phi = 2.*!dpi*u
	theta = ACOS(2.*v-1.)

	ang2pixClock = TIC('ang2pixClock')
	ANG2PIX_RING, 16, theta, phi, ipring
	TOC, ang2pixClock

	for i=0,N_ELEMENTS(fnums)-1 do begin
		thispix = where(ipring eq fnums[i]) ; indices w/ coords in this tile
		ncoord = n_elements(thispix)
		coordind = lindgen(ncoord)

		glon = phi[thispix[coordind]]*180./!dpi
		glat = (theta[thispix[coordind]]-!dpi/2.)*180./!dpi
		; Transform from galactic healpix to ecliptic observations
		euler, glon, glat, elon, elat, select=6
		euler, glon, glat, ra, dec, select=2

		coordNum = TRANSPOSE([[ipring[thispix]], [glon], [glat], [elon], [elat], [ra], [dec]])

		;writeClock = TIC('writeClock-'+STRTRIM(STRING(fnums[i])))
		;WRITE_CSV, 'coordHPnum'+STRTRIM(STRING(fnums[i]))+'.csv', dat
		;TOC, writeClock

		; writing to CSV takes 20* as long...
		savClock = TIC('savClock-'+STRTRIM(STRING(fnums[i]),2))
		SAVE, coordNum, FILENAME='coordHPnum'+STRTRIM(STRING(fnums[i]),2)+'.sav'
		TOC, savClock

	endfor
	
	TOC
END

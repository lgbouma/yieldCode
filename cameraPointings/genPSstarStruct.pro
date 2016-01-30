;PURPOSE: given input .dat file of tileNumber - psNumber - eclLong - eclLat,
;generate and save a star structure with those attributes ('psStruct.sav'), 
;to be used in starSurvey.pro

PRO genPSstarStruct
	FMT = 'I,I,F,F'
	READCOL, 'psCoords.dat', F=FMT, tileNum, psNum, eLong, eLat

	nStars = N_ELEMENTS(elat)
	cat = REPLICATE({pointingstruct}, nStars)
	
	cat.coord.elon = eLong
	cat.coord.elat = eLat
	cat.coord.healpix_n = tileNum
	cat.tileNum = tileNum
	cat.psNum = psNum

	SAVE, cat, FILENAME='psStruct.sav'
END

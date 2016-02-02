;PURPOSE: given input .dat file of tileNumber - psNumber - eclLong - eclLat,
;for example, "tileCoordsMed.dat"
;generate and save a star structure with those attributes,
;to be used in starSurvey.pro (which generates the # pointings each coordt gets)

PRO genStarStruct, fName
	FMT = 'I,I,F,F'
	;READCOL, 'psCoords.dat', F=FMT, tileNum, psNum, eLong, eLat
	;READCOL, 'tileCoords.dat', F=FMT, tileNum, psNum, eLong, eLat
	READCOL, fName, F=FMT, tileNum, psNum, eLong, eLat

	nStars = N_ELEMENTS(elat)
	cat = REPLICATE({pointingstruct}, nStars)
	
	cat.coord.elon = eLong
	cat.coord.elat = eLat
	cat.coord.healpix_n = tileNum
	cat.tileNum = tileNum
	cat.psNum = psNum

	;SAVE, cat, FILENAME='psStruct.sav' ; 16/01/30: original use was to calc how many pointings PSs got
	;SAVE, cat, FILENAME='tileStruct.sav' ; 16/01/31: how many pointings do nominal tile coords get?
	;SAVE, cat, FILENAME='tileStructMed.sav' ; 16/02/01: tileStruct.sav coords wonky b/c of MEAN vs MEDIAN
	STOP
	inName = STRSPLIT(fname, '.', /EXTRACT)
	print, inName
	outName = inName[0] + 'Struct.sav'	
	print, outName
	SAVE, cat, FILENAME=outName ; 16/02/02: breaking out ext missions
END

pro prepForRead
	dat = read_csv_col('asteroseisStars.csv')
	
	nSG = n_elements(dat.glon)
	PRINT, 'Total ', nSG, ' possible subgiant stars', STRING(10B)
	; relevant IDL routine: ang_2_pix (..)

	phi = dat.glon*!dpi/180. ; convert to rad
	theta = dat.glat*!dpi/180.
	ang2pix_ring, 16, theta, phi, ipring

	; that's actually _it_. this doesn't merit a subroutine. it's going into hacky tile_wrap.

end

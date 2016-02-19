PRO getCoords
	bigNumber = 2e7 ; 2e7 produces 1 gb lists
	u = RANDOMU(seed, bigNumber) ; uniform distribution over 0 to 1
	v = RANDOMU(seed, bigNumber)
	phi = 2.*!dpi*u*(1/24.)
	theta = ACOS(2.*v-1.)  ; uniform over whole sphere. might need to segment...

	STOP
	ANG2PIX_RING, 16, theta, phi, ipring

	;dat = [[theta], [phi], [ipring]]
	dat = TRANSPOSE([[theta], [phi], [ipring]])

	WRITE_CSV, 'coordWithHPnum.csv', dat

END

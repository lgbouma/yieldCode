function ps_sel, tmag, teff, mass, rad, ph_p, $ 
  minrad=minrad, per=per, rn_pix=rn_pix, geom_area=geom_area, npnt=npnt, ps_only=ps_only
;+
;PURPOSE: Select stars that are postage stamps in this healpix tile. 
;		(down to resolution of healpix tiles, at the time of writing). 
;INPUTS:
;	tmag: tess magnitude of targets object (an array of starStructs restored from trilegal files).
;	teff: effective temperature of targets object
;	mass: mass of targets object
;	rad: radius of targets object
;	ph_p: photon fluxes for T=10 vs Teff. Nominally, takes in "ph_fits"
;		which is a 9*4 array of floats, numbers between -50 and 70.
;	minrad: radius needed for number of stars meeting detection criterion
;		to be about 200k. In S+2015, it's 2.25R_e. 
;	per: nominal 20-day period assumed (used to calc transit depth)
;	rn_pix: 15. # of pixels in fiducial aperture (equivalent to saying nominal
;		star brightness is 9.5 I_c mags).
;	geom_area: 
;	npnt: array of (float) numbers 1., 2., ... 12, of length 2908.
;		I.e., this is supposed to be "this tile gets this many pointings."
;RETURNS:
;	sel: stars selected as postage stamps (an array of their indices).
;COMMENTS:
;	Description of method Section6.7 S+2015. As of S+15 process miscounted
;	PS observing time (gave all at least one pointing) (see 'npts.fits')
;	Bug fix: 2016-01-31, Luke Bouma. Implemented in tile_wrapper as "fix" of
;	npnt_fits by skipping entire tile if tile has no postage stamps (& is ps_only).
;   	In this case, ps_sel.pro should never be called.

	nStars = N_ELEMENTS(tmag)		; number of stars (total) for this tile
	if (KEYWORD_SET(minrad)) then minrad=minrad else minrad=2.27 ; diff from Sullivan+ 2015 value.
	if (KEYWORD_SET(per)) then per=per else per=20.0 ; fiducial planet minimum rad
	assert, (npnt ge 0), 'at least need npnt to be well-defined!'
	if ps_only eq 1 then begin
		assert, npnt ne 0, 'NumPointings in ps_sel is 0. This tile should have been skipped.'
	endif

	sz_ph_p = size(ph_p)
	nfilt = sz_ph_p[1]
	ph_filt = dblarr(nfilt, nstars)
	ph_star = dblarr(nstars)
	if (KEYWORD_SET(geom_area)) then geom_area=geom_area else geom_area = 69. ;(74.6) cm^2
	if (KEYWORD_SET(rn_pix)) then rn_pix=rn_pix else rn_pix=10. ; pixels in read noise calcn
	recipteff = 4000./teff
	; compute observed photon number flux thru tess bandpass by sum over PRF(wavelength) [sec 6.1 S+15]
	for j=0, nfilt-1 do begin 
		ph_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff + $
					   ph_p[j,2]*recipteff^2. + ph_p[j,3]*recipteff^3.
	endfor
	ph_filt[where(ph_filt lt 0.0)] = 0.0
	ph_star = 1.0*10.^(-0.4*(tmag-10.))*total(ph_filt, 1) ; observed photon flux for each star
	AU_IN_RSUN = 215.093990942
	REARTH_IN_RSUN = 0.0091705248
	a = mass^(1./3.) * (float(per)/365.25)^(2./3.) ; in AU
	dur = rad * float(per) / (!DPI*a*AU_IN_RSUN) ; duration of nominal transit (days)
	nTransitPerPntg = 1.3 ; avg accounts for on/off time, checked in 'numerics/nTransitPerPntg.py'
	exptime = nTransitPerPntg*float(npnt)*dur*24.*3600 ; exposure time of transit (sec) 
	dep = (REARTH_IN_RSUN * float(minrad) / rad)^2.0
	sig = dep/7.3
	rn = rn_pix*sqrt(4.0*exptime/2.0) ; read noise. Read/stack every 2sec, *4 for CCDs(?)
	minphot = 1.5*(1.+sqrt(1.+4.*sig^2.*rn^2.))/(2.*sig^2.)
	; note minphot/(exptime*geom_area) is the minimum flux required for signal detection
	; also pretty clear here that shouldn't have exptime = 0, but idl doesn't raise div0 errors...
	sel = where(ph_star gt (minphot/(exptime*geom_area)) and teff gt 1500 and teff lt 15000)

	print, 'Selecting ', n_elements(sel),' postage stamps out of ', nstars, ' stars.'
	STOP
	return, sel

END

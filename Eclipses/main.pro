PRO main
  fnums = mrdfits('fnums.fits')	; File containing healpix numbers (skips galactic plane healpix tiles)
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                0, $ ; HEBs
                0  ] ; BTPs
  tile_wrapper, '/home/luke/Dropbox/tessSimulation/trilegal/', fnums, $
	'tileSkipNewNom.fits', $; output file name. Nb: fits files add extensions (they dont overwrite)
  	n_trial=10, $	; number of trials (10 good for reasonable statistics)
	eclass=eclass, $; from above
	ps_only=1, $	; 1=only run postage stamps, 0=run ffis as well
	detmag=0, $	; If you want the sim to return a magnitude-limited catalog, set this to the limit
			;    0=run the TESS model for detection
        pla_err=0, $    ; 0=run with nominal occurrence rates, +1 for upper bounds, -1 for lower
        ;prf_file='psfs/dfrac_t75_f3p31_3.fits' ; ideal PRF file
	prf_file='psfs/dfrac_asbuilt_75c_0f.fits', $ ; built PRF from Deb Woods
	prototypeMode=0	; 1=run fast prototypes, 0=run full simulations.
END

PRO main
  fnums = mrdfits('fnums.fits')	; File containing healpix numbers (skips galactic plane healpix tiles)
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                0, $ ; HEBs
                0  ] ; BTPs
  tile_wrapper, '../../outStarLib/', fnums, $
	'tempOut.fits', $; output file name (w/ fits ext). Nb: fits files dont overwrite
  n_trial=10, $	; number of trials (10 good for reasonable statistics)
	eclass=eclass, $ ; from above
	ps_only=1, $	; 1=only run postage stamps, 0=run ffis as well
	detmag=0, $	; If you want the sim to return a magnitude-limited catalog, set this to the limit
              ; 0=run the TESS model for detection. (LB: I think is V-mag? or TESS mag..)
	pla_err=0, $    ; 0=run with nominal occurrence rates, +1 for upper bounds, -1 for lower
	;prf_file='psfs/dfrac_t75_f3p31_3.fits' ; ideal PRF file
	prf_file='psfs/dfrac_asbuilt_75c_0f.fits', $ ; built PRF from Deb Woods
	prototypeMode=0, $	; 0=full simulations. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)
	fCamCoordPri='../cameraPointings/ns_nominal_camCoord.dat', $ ;where cams pnt for primary mission.
	fCamCoordExt='../cameraPointings/ns_ep_camCoord.dat', $ ;cam coords for ext mission. '' if no ext.
  psPriFile='../../preProcessing/sTIC-selection/primary/sTIC-eq1allCuts-primary.sav', $ ; primary PSs
  psExtFile='../../preProcessing/sTIC-selection/ext-poles/sTIC-eq1allCuts-ext-poles.sav', $ ;ext PSs
            $ ; If primary only, leave as ''
	fTilesCounts='../cameraPointings/ns_nominal_tilesCounts.sav', $ ;pointingStruct w tileNum, coord, nPntg
        ; Clearly, fTilesCounts needs to _match_ fCamCoord (it's generated separately).
	burtCatalog=0 ; are you making catalogs to send to Jenn Burt so she can plan APF RV followup?
END

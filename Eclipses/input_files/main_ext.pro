PRO main_ext
  cd, '../'
  n_trial = fix(getenv('n_trial'))
  prototype_mode = fix(getenv('prototype_mode')) ; fix to integer type
  ext_mission_name = getenv('ext_mission_name') 
  date_sub = getenv('date_sub')
  ; e.g., '160526_shemi_nhemi_nhemi_10t.fits'
  out_fname = date_sub + '_' + ext_mission_name + '_t' + strtrim(string(n_trial),2) + '.fits'
  fcam_coord_ext = getenv('fcam_coord_ext')
  fcam_coord_pri = getenv('fcam_coord_pri')
  ps_ext_file = getenv('ps_ext_file')
  ps_pri_file = getenv('ps_pri_file')
  ffi_ext_file = getenv('ffi_ext_file')
  ffi_pri_file = getenv('ffi_pri_file')

  fnums = mrdfits('fnums.fits')	; File containing healpix numbers (skips galactic plane healpix tiles)
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                0, $ ; HEBs
                0  ] ; BTPs
  tile_wrapper, '../../outStarLib/', fnums, $
  out_fname, $; output file name (w/ fits ext). Nb: fits files dont overwrite
  n_trial=n_trial, $	; number of trials (10 good for reasonable statistics)
  eclass=eclass, $ ; from above
  ps_only=0, $	; 1=only run postage stamps, 0=run ffis as well
  detmag=0, $	; If you want the sim to return a magnitude-limited catalog, set this to the limit
              ; 0=run the TESS model for detection. (LB: I think is V-mag? or TESS mag..)
  pla_err=0, $    ; 0=run with nominal occurrence rates, +1 for upper bounds, -1 for lower
  ;prf_file='psfs/dfrac_t75_f3p31_3.fits' ; ideal PRF file
  prf_file='psfs/dfrac_asbuilt_75c_0f.fits', $ ; built PRF from Deb Woods
  prototypeMode=prototype_mode, $	; 0=full sim. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)
  fCamCoordPri=fcam_coord_pri, $ ;where cams pnt for primary mission.
  fCamCoordExt=fcam_coord_ext, $ ;cam coords for ext mission. '' if no ext.
  psPriFile=ps_pri_file, $ ; primary PSs
  ffiPriFile=ffi_pri_file, $ ; primary FFIs
  psExtFile=ps_ext_file, $ ; If primary only, leave as ''
  ffiExtFile=ffi_ext_file, $ ;
  burtCatalog=0, $ ; are you making catalogs to send to Jenn Burt so she can plan APF RV followup?
  batchJob=1 ; 1 to save output to ../output_files, if running from input_files. Else 0.
END

# Simulating Tess' Exoplanet detectIoNs: STEIN
(everyone gets an anocronym, right?)

This code simulates the population of planets that we expect to detect with TESS.
The most relevant references to understand what it's doing are:
  1. Bouma et al 2016 (hereafter B+16)
  2. Sullivan et al 2015 (hereafter S+15)
  3. The actual comments L.B. wrote up as headers in the code.
  4. `usersmanual.pdf` in the `Doc` folder (dated pre-L.B.'s work).

## Installation overview
The entire code comes in three parts:
`preProcessing`, `yieldCode`, and `postProcessing`.
"Running the yield code" can be done quite quickly once all this is set up. 

For anything to work, you need to obtain `preProcessing` as well as two other separate libraries 
(`outStarLib` and `idlutils`) from Luke Bouma (lgbouma@mit.edu). 

## Installation instructions
0. Be on a mac or linux system. Have a functioning version of IDL (I'm running
    on IDL v8.5).
1. Make a new folder somewhere.
2. Download `preProcessing`, `outStarLib`, and `idlutils` from the dropbox
  link Luke sent you. Put them in your new folder.
  `idlutils` is a library of useful scripts, slightly customized from what Goddard
  puts out.
  `outStarLib` is 60Gb of the stellar data for the bright, medium, and dim
  catalogs S+15 Sec3.1 talks about.
  `preProcessing` is scripts that take us from `outStarLib` to which stars 
  we care about. There is a repository with the version history for all the
  preProcessing scripts as well, but it doesn't have the data you care about, 
  which are the _output_ of those scripts:
  selected stars for which we simulate transits & their observations at 2
  & 30 minutes in `yieldCode`.
  These are saved in `preProcessing/sTIC-selection` and have names like `shemi_nhemi_*_4M*/*.sav`
  where the different words label different spacecraft pointings.
3. Inside your new folder, type:

  ```
  git clone https://github.com/lgbouma/yieldCode
  git clone https://github.com/lgbouma/postProcessing
  ```

  So your directory structure is now:

  ```
  ├── idlutils
  ├── outStarLib
  ├── postProcessing
  ├── preProcessing
  ├── yieldCode
  ```

  Then, do:

  ```
  cd yieldCode
  git checkout ps_ffi_joint
  cd ../postProcessing
  git checkout ps_ffi_joint
  ```

  The `ps_ffi_joint` branch can do both PSs and FFIs in a timely manner (which is
  probably what you want. If you only want PSs, you can stay on master branch, but
  it's a bit dated [see versioning logs]).

4. You now have all the files. Let's make them work.
  To do so, add the `idlutils` directory to your IDL path.
  For instance, use [these instructions](http://slugidl.pbworks.com/w/page/28913708/Adding%20Programs%20to%20Your%20IDL%20Path).
  As a check, you want to be able to run idl from the command line, and then type:

  ```
  IDL> assert, 1
  ```

  and get

  ```
  % Compiled module: ASSERT.
  ```

  returned, instead of an error about not knowing what `assert` means.
5. You should now be able to run the yieldCode!
  Do it by calling `Eclipses/main.pro`, which is a wrapper to `tile_wrapper.pro`, which
  is a loop over all the healpix tiles for all the stars in the sky for which you'll be simulating
  transits and observations.
  (If you want to speed things up, `split_for` would parallelize this, but I haven't implemented it).
  After editing `main.pro` to your desired input parameters (see below), do something like:
  ```
  cat main.pro > my_log_name.out && idl -e main &>> my_log_name.out & 
  ```
  which copies your input file to a log file then executes the main `yieldCode` loop
  in a background process while appending the text output to the same log file.

  Your desired input parameters (i.e. the `main.pro` file) for a first check might be:
  ```
  PRO main
  fnums = mrdfits('fnums.fits') ; File containing healpix numbers (skips galactic plane healpix tiles)
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                0, $ ; HEBs
                0  ] ; BTPs
  tile_wrapper, '../../outStarLib/', fnums, $
  'temp_test.fits', $; output file name (w/ fits ext). Nb: fits files dont overwrite
  n_trial=1, $  ; number of trials (20 good for reasonable statistics)
  eclass=eclass, $ ; from above
  ps_only=0, $  ; 1=only run postage stamps (depricated?), 0=run ffis as well
  detmag=0, $ ; If you want the sim to return a magnitude-limited catalog, set this to the limit
              ; 0=run the TESS model for detection. (LB: I think is V-mag? or TESS mag..)
  pla_err=0, $    ; 0=run with nominal occurrence rates, +1 for upper bounds, -1 for lower
  ;prf_file='psfs/dfrac_t75_f3p31_3.fits' ; ideal PRF file
  prf_file='psfs/dfrac_asbuilt_75c_0f.fits', $ ; built PRF from Deb Woods
  prototypeMode=2, $  ; 0=full simulations. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)
  fCamCoordPri='../cameraPointings/shemi_nhemi_orbits.dat', $ ;where cams pnt for primary mission.
  fCamCoordExt='../cameraPointings/nhemi_orbits.dat', $ ;cam coords for ext mission. '' if no ext.
  psPriFile='../../preProcessing/sTIC-selection/shemi_nhemi_nhemi_4Morb/shemi_nhemi_200k.sav', $ ; primary PSs
  ffiPriFile='../../preProcessing/sTIC-selection/shemi_nhemi_nhemi_4Morb/shemi_nhemi_38M.sav', $ ; primary FFIs (pseudo-FFIs; drawn from 3.8million)
  psExtFile='../../preProcessing/sTIC-selection/shemi_nhemi_nhemi_4Morb/shemi_nhemi_nhemi_200k.sav', $
            $ ; If primary only, leave as ''
  ffiExtFile='../../preProcessing/sTIC-selection/shemi_nhemi_nhemi_4Morb/shemi_nhemi_nhemi_38M.sav', $
            $ ; If primary only, leave as ''
  burtCatalog=0, $ ; are you making catalogs to send to Jenn Burt so she can plan APF RV followup?
  batchJob=0 ; 1 to save output to ../output_files, if running from input_files. Else 0.
END
  ```
  which should run a baby-version of the simulation in less than a minute.
  If this works, congratulations! You've simulated looking at 0.3% of the sky with TESS, for
  a `shemi_nhemi_nhemi` mission (2 year primary, one year staying in the northern hemisphere).
  You should have two output files, `temp_test-pri.fits` and `temp_list-ext.fits`.
  These are the output files from the primary and extended missions with all the information
  about the simulated planets you detected (note that the way this
  works, if you just want primary missions without any kind of extended mission junk,
  you only need to pay attention to `temp_test-pri.fits`).

  If you want to look at the entire sky, change `prototypeMode`.
  If you want something statistical, up `n_trials`, say to around 20.
6. You're now ready to see what you detected. Head over the `postProcessing` to see
  that README file.

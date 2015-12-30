PRO requirements 

  psfsize = 1.0
  psfstr = '_1p'+strtrim(round(10.*(psfsize-1.0)),2)
  ;psfstr = '_dither' + psfstr
  filen = '../Sims/ss_r200.sav' 
  frac_file = '../ExpTimeCalc/frac24'+psfstr+'.fits' 
  fov = 24.0
  seg = 13
  geomarea = 74.6
  readnoise=10.0
  thresh = 7.0
  tranmin = 2.0
  n_trial = 3
  REARTH_IN_RSUN = 0.0091705248

  
  fname = 'req_r200_'+strtrim(seg,2)+'x'+strtrim(round(fov),2)+'.fits'
  period_boundary = [1.0, 2.0, 3.42, 5.85, 10.0, 17.1, 29.2, 50.0, $
			85.5, 146.2]
  mag_bins = findgen(46)/3. + 5.
  logsig_bins = findgen(31)/10. - 5.
  sig_bins  = 10.^logsig_bins
  n_mag = n_elements(mag_bins)-1
  n_sig = n_elements(sig_bins)-1

  se_magsig = dblarr(n_mag, n_sig, n_trial)
  et_magsig = dblarr(n_mag, n_sig, n_trial)
  mn_magsig = dblarr(n_mag, n_sig, n_trial)
  hz_magsig = dblarr(n_mag, n_sig, n_trial)


  period_boundary = [1.0, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0]; 418.0]
  p8 = 3
  radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0] ; Fressin gives 22 as upper limit but we lower it here to 15
  planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']

  rate_fressin = dblarr(10,5) ; period bin, radius bin
  rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
  rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
  rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
  rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
  rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
  rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
  rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
  rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
  rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
  rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
  
  rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
  rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
  rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
  rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
  rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
  rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
  rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
  rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
  rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
  rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
 ;rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
 rate_fressin = rate_fressin/100.

  restore, filen
  n_stars = n_elements(star)
  print, 'Number of stars: ', n_stars
  randomp, mn_rad, -2,  n_stars, range_x=[2.0, 4.0], seed=seed
  randomp, se_rad, -2, n_stars, range_x=[1.25, 2.0], seed=seed
  randomp, et_rad, 0,  n_stars, range_x=[0.8, 1.25], seed=seed

  hz = n_elements(star[0].planet.p)-1
  n_per  = 10
  n_rad  = 5
  nsbefore=0.0
  bwtnfrac=0.0
  print, n_per, n_mag
  
  for ii=0, n_trial-1 do begin
    sp_name = '../Sims/sp_'+strtrim(ii,2)+'.sav'
    spo_name = '../Sims/spo'+'.sav'
    rad_planets, struct=star, infile=filen ;nfil=filen, outfile=sp_name
    req_observe, struct=star, filen='rad', geomarea=geomarea, fov=fov, $ ;infil=sp_name,outfile=spo_name
    	readnoise=readnoise, thresh=thresh, tranmin=tranmin, frac_file=frac_file
    for jj=0, n_per-1 do begin
      for mm=0, n_mag-1 do begin
        ;print, ii, jj, mm
        thismag = where((star.mag.i ge mag_bins[mm]) and $
			(star.mag.i lt mag_bins[mm+1]) and $
			(star.planet[jj].ntra_obs ge 2))
        thismag_hz = where((star.mag.i ge mag_bins[mm]) and $
			(star.mag.i lt mag_bins[mm+1]) and $
			(star.planet[hz].ntra_obs ge 2) and $
			(star.planet[hz].p ge period_boundary[jj]) and $
			(star.planet[hz].p lt period_boundary[jj+1]))
        if (thismag[0] ne -1) then begin
          print, 'imag: ', mag_bins[mm], ' period: ', jj, ' stars: ', n_elements(thismag)
          exptime = double(star[thismag].planet[jj].ntra_obs) * $
       			 star[thismag].planet[jj].dur * 24.0
    	
	  ; Earths      
  	  dep = (REARTH_IN_RSUN * et_rad[thismag] / star[thismag].r)^2.
          sig = dep/(7.*(1. + star[thismag].dil))
          prec = sig * sqrt(exptime)
          for nn=0, n_sig-1 do begin
	    thissig = where((prec ge sig_bins[nn]) and $
			  (prec lt sig_bins[nn+1]))
            if (thissig[0] ne -1) then begin
	      et_magsig[mm,nn,ii] = et_magsig[mm,nn,ii]+rate_fressin[jj,0]*n_elements(thissig)
	    endif
          endfor

	  ; Super-Earths
  	  dep = (REARTH_IN_RSUN * se_rad[thismag] / star[thismag].r)^2.
          sig = dep/(7.*(1. + star[thismag].dil))
          prec = sig * sqrt(exptime)
          for nn=0, n_sig-1 do begin
	    thissig = where((prec ge sig_bins[nn]) and $
			  (prec lt sig_bins[nn+1]))
            if (thissig[0] ne -1) then begin
	      se_magsig[mm,nn,ii] = se_magsig[mm,nn,ii]+rate_fressin[jj,1]*n_elements(thissig)
	    endif
          endfor
 
	  ; Mini-Neptunes
  	  dep = (REARTH_IN_RSUN * mn_rad[thismag] / star[thismag].r)^2.
          sig = dep/(7.*(1. + star[thismag].dil))
          prec = sig * sqrt(exptime)
          for nn=0, n_sig-1 do begin
	    thissig = where((prec ge sig_bins[nn]) and $
			  (prec lt sig_bins[nn+1]))
            if (thissig[0] ne -1) then begin
	      mn_magsig[mm,nn,ii] = mn_magsig[mm,nn,ii]+rate_fressin[jj,2]*n_elements(thissig)
	    endif
          endfor
        endif

	if (thismag_hz[0] ne -1) then begin
	
	  ; HZ Super-Earths
  	  dep = (REARTH_IN_RSUN * se_rad[thismag_hz] / star[thismag_hz].r)^2.
          sig = dep/(7.*(1. + star[thismag_hz].dil))
          prec = sig * sqrt(exptime)
          for nn=0, n_sig-1 do begin
	    thissig = where((prec ge sig_bins[nn]) and $
			  (prec lt sig_bins[nn+1]))
            if (thissig[0] ne -1) then begin
	      hz_magsig[mm,nn,ii] = hz_magsig[mm,nn,ii]+rate_fressin[jj,1]*n_elements(thissig)
	    endif
          endfor
	  ; HZ Earths
  	  dep = (REARTH_IN_RSUN * et_rad[thismag_hz] / star[thismag_hz].r)^2.
          sig = dep/(7.*(1. + star[thismag_hz].dil))
          prec = sig * sqrt(exptime)
          for nn=0, n_sig-1 do begin
	    thissig = where((prec ge sig_bins[nn]) and $
			  (prec lt sig_bins[nn+1]))
            if (thissig[0] ne -1) then begin
	      hz_magsig[mm,nn,ii] = hz_magsig[mm,nn,ii]+rate_fressin[jj,0]*n_elements(thissig)
	    endif
          endfor
        endif
      endfor
    endfor
  endfor
  mwrfits,mean(mn_magsig, dimension=3),'mn_'+fname
  mwrfits,mean(se_magsig, dimension=3),'se_'+fname
  mwrfits,mean(et_magsig, dimension=3),'et_'+fname
  mwrfits,mean(hz_magsig, dimension=3),'hz_'+fname
end

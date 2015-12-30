PRO report, nstar3d=nstar3d, nplanet3d=nplanet3d

  psfsize = 1.0
  ;psfstr = '_4700k'
  psfstr = '_1p'+strtrim(round(10.*(psfsize-1.0)),2)
  ;psfstr = '_asbuilt' + psfstr
  filen = '../Sims/ss_det_13x24.sav' 
  frac_file = '../ExpTimeCalc/frac24'+psfstr+'.fits' 
  fov = 24.0
  seg = 13
  geomarea = 67.5 ;74.6
  readnoise=10.0
  thresh = 7.0
  tranmin = 2.0
  n_trial = 10
  
  magbins = dindgen(22)/2.0+6.75 
  fname = '../Report/rep_0_'+strtrim(seg,2)+'x'+strtrim(round(fov),2)+'_'+strtrim(round(geomarea),2)+psfstr
  period_boundary = [1.0, 2.0, 3.42, 5.85, 10.0, 17.1, 29.2, 50.0, $
			85.5, 146.2]
  elatbins = dindgen(10)*10.
;  set_plot,'ps'
;  device,filen='sim_'+filen+'.eps',/encapsulated,xsize=6.5,ysize=9.5,/inches,$
;           /color,bits_per_pixel=8,/isolatin,/helvetica
;  !p.font=0
;  !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.5
  
;  plot_survey, filen=survey_name

  ; Sky map illustrating stars and scan strategy.  
  ; Number of searchable stars for planets of different sizes.  
  ; Expected number of planets.  
  ; App mag of 1st, 10th, 100th brightest star for each type of planet.  
  ; HZ planets and brightnesses, star types.

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
     ;rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
     rate_fressin = rate_fressin/100.

  restore, filen
  nstars = n_elements(star)
  print, 'Number of stars: ', nstars
 
 ; ysize=9.5, !p.charsize=1.5 for full-page graphics 
  set_plot,'ps'
     device,filen=fname+'_1.eps',/encapsulated,xsize=7.0,ysize=3.2,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.charsize=0.8

  !P.Multi=[0,2,1] 
  colors=['black', 'navyblue', 'royalblue', 'aquamarine', $
          'green', 'lawngreen', 'yellow', 'yellow', 'orange', $
         'orange', 'orangered', 'orangered','red','red', 'red', 'red']
  nplot=50000
  plotme = floor(double(nstars)*randomu(seed,nplot))
  color_index = star[plotme].npointings
  ;q = where(color_index gt seg)
  ;if (q[0] ne -1) then color_index[q] = 13
  
  ;aitoff, star[plotme].coord.elon, star[plotme].coord.elat, x, y
  x = cos(star[plotme].coord.elon*!DPI/180.)*(90.-abs(star[plotme].coord.elat))
  y = sin(star[plotme].coord.elon*!DPI/180.)*(90.-abs(star[plotme].coord.elat))

  ;LoadCT, 27 
  cgplot, x, y, psym=16, symsize=0.1, /nodata, $
        xtit='Scan Polar Projection' 
        ;xtit='ecliptic longitude', $
        ;ytit='ecliptic latitude'

  for i=0,seg+1 do begin
    q=where(color_index eq i)
    ;c = floor(float(i)*255./float(seg))
    ;print, n_elements(q), c
    if q[0] ne -1 then begin
      ;LoadCT, 27 
      cgoplot, x[q], y[q], psym=16, color=colors[i], symsize=0.1
    endif
  endfor
 
  barx = indgen(16)
  bary = dblarr(16)
  ;colors = replicate('gray', 16)
  for ii=0,seg+1 do begin
    bary[ii] = n_elements(where(star.npointings eq ii))
  endfor
  gd = where(bary gt 1)
  names = string(barx, format='(I3)')

  cgbarplot, alog10(bary[gd]), barcoords=barx[gd], barnames=names[gd], colors=colors[gd], $
	xtitle='Number of Pointings', ytitle='log(Number of Stars)'
 ;, ytickformat='(E4.1)'
  device,/close
  nocoverage = 100.*double(bary[0])/double(nstars)
  fullcoverage = 100.*double(n_elements(where(star.npointings ge seg)))/double(nstars)
  
  hz = n_elements(star[0].planet.p)-1
  n_per  = 10
  n_rad  = 5
  nsbefore=0.0
  bwtnfrac=0.0
  n_mag = n_elements(magbins)-1
  n_elat = n_elements(elatbins)-1
;  nstar3b = dblarr(n_trial, n_per+1, n_rad)
  nstar3d = dblarr(n_trial, n_per, n_rad)
  nstar2d = dblarr(n_per, n_rad)
  nstar2derr = dblarr(n_per, n_rad)
  nplanet3d = dblarr(n_trial, n_per, n_rad)
  nplanet2d = dblarr(n_per, n_rad)
  nplanet2derr = dblarr(n_per, n_rad)
  nplanethz=dblarr(n_trial)
  npexthz=dblarr(n_trial)
  npjwstexthz=dblarr(n_trial)
  npjwsthz=dblarr(n_trial)
  mags2 = dblarr(n_trial, n_mag)
  mags4 = dblarr(n_trial, n_mag)
  teff2hz = dblarr(n_trial)
  teff2hzs = dblarr(n_trial) 
  teff2all = dblarr(n_trial) 
  teff2alls = dblarr(n_trial) 
  elat2 = dblarr(n_trial, n_elat)
  satmag = dblarr(n_trial)

  s_small_hz = dblarr(n_trial)
  s_jwst_hz =  dblarr(n_trial)
  s_jwst2_hz = dblarr(n_trial)
  
  for ii=0, n_trial-1 do begin
    sp_name = '../Sims/sp_'+strtrim(ii,2)+'.sav'
    spo_name = '../Sims/spo'+'.sav'
    rad_planets, struct=star, infile=filen ;nfil=filen, outfile=sp_name
    rad_observe, struct=star, filen='rad', geomarea=geomarea, fov=fov, $ ;infil=sp_name,outfile=spo_name
    	readnoise=readnoise, thresh=thresh, tranmin=tranmin, frac_file=frac_file
    ;restore, spo_name

    ; Things to report:
    ; Saturation of host stars
    satmag[ii] = max(star[where(star.sat)].mag.i)
    
    ; Teff of HZ hosts
    hz2 = where((star.planet[hz].r le 2) and (star.planet[hz].det))
    teff2hz[ii] = mean(star[hz2].teff)
    teff2hzs[ii] = stddev(star[hz2].teff)
    ; Teff of 8-day hosts
    all2 = where((star.planet[p8].r le 2) and (star.planet[p8].det gt 0))
    teff2all[ii] = mean(star[all2].teff)
    teff2alls[ii] = stddev(star[all2].teff)
 
    ; # stars with small HZ planets
    s_small_hz[ii] = total((star.planet[hz].r le 2.0) and $
			   star.planet[hz].det)
  
    ; # stars with HZ planets in continuous viewing zone
    s_jwst_hz[ii]  = total((abs(star.coord.elat) gt 85.0) and $
			      (star.planet[hz].r le 2.0) and $
			      (star.planet[hz].det))
    s_jwst2_hz[ii] = total((abs(star.coord.elat) gt 75.0) and $
			      (star.planet[hz].r le 2.0) and $
			      (star.planet[hz].det))
   
    ; Magnitudes of small HZ hosts 
    for mm=0,n_mag-1 do begin
	nmag = where((star.planet[hz].r lt 4) and (star.planet[hz].det) and $
		     (star.mag.i lt magbins[mm+1]) and (star.mag.i ge magbins[mm]))
        if (nmag[0] ne -1) then mags4[ii,mm] = n_elements(nmag)
	nmag = where((star.planet[hz].r lt 2) and (star.planet[hz].det) and $
		     (star.mag.i lt magbins[mm+1]) and (star.mag.i ge magbins[mm]))
        if (nmag[0] ne -1) then mags2[ii,mm] = n_elements(nmag)
    endfor    

    ;Ecliptic latitude of small HZ hosts
    for nn=0,n_elat-1 do begin
	nelat = where((star.planet[hz].r le 2.0) and (star.planet[hz].det) and $
		     (abs(star.coord.elat) lt elatbins[nn+1]) and $
		     (abs(star.coord.elat) ge elatbins[nn]))
        if (nelat[0] ne -1) then elat2[ii,nn] = n_elements(nelat)
    endfor

    print, 'Trial: ',ii,' Small-planet hosts: ',s_small_hz[ii], $
           ' Small-planet hosts above 75 deg Lat: ', s_jwst2_hz[ii], $
           ' Total hosts: ',total(star.planet[0].det)
      
  for jj=0, n_per-1 do begin
    dlogp = alog(period_boundary[jj+1])-alog(period_boundary[jj])   
    for kk=0, n_rad-1 do begin
      dinvr = 1./radius_boundary[kk]-1./radius_boundary[kk+1]   
	  ; Small HZ planets anywhere     
	  nhzbefore =     ((star.planet[hz].r le radius_boundary[kk]) and $
	  		   (star.planet[hz].r le 2.0) and $
			   (star.planet[hz].p ge period_boundary[jj]) and $
			   (star.planet[hz].p lt period_boundary[jj+1]) and $
			   (star.planet[hz].det))
	  nhzbtwn =  where((star.planet[hz].r le radius_boundary[kk+1]) and $
			   (star.planet[hz].r gt radius_boundary[kk]) and $
			   (star.planet[hz].r le 2.0) and $
			   (star.planet[hz].p ge period_boundary[jj]) and $
			   (star.planet[hz].p lt period_boundary[jj+1]) and $
			   (star.planet[hz].det))
          nplanethz[ii] = nplanethz[ii] + rate_fressin[jj,kk] * total(nhzbefore)
          if (nhzbtwn[0] ne -1) then begin
	    nplanethz[ii] = nplanethz[ii] + rate_fressin[jj,kk]* $ 
		 total(1./star[nhzbtwn].planet[hz].r - 1./radius_boundary[kk+1])/dinvr
          endif
   	  ; Small HZ planets above 75 deg
	  nhzbefore =      ((star.planet[hz].r le radius_boundary[kk]) and $
	  		   (star.planet[hz].r le 2.0) and $
			   (star.planet[hz].p ge period_boundary[jj]) and $
			   (star.planet[hz].p lt period_boundary[jj+1]) and $
			   (abs(star.coord.elat) gt 75) and $
			   (star.planet[hz].det))
	  nhzbtwn =  where((star.planet[hz].r le radius_boundary[kk+1]) and $
			   (star.planet[hz].r gt radius_boundary[kk]) and $
			   (star.planet[hz].r le 2.0) and $
			   (star.planet[hz].p ge period_boundary[jj]) and $
			   (star.planet[hz].p lt period_boundary[jj+1]) and $
			   (abs(star.coord.elat) gt 75) and $
			   (star.planet[hz].det))
          npjwsthz[ii] = npjwsthz[ii] + rate_fressin[jj,kk] * total(nhzbefore)
          if (nhzbtwn[0] ne -1) then begin
	    npjwsthz[ii] = npjwsthz[ii] + rate_fressin[jj,kk]* $ 
		 total(1./star[nhzbtwn].planet[hz].r - 1./radius_boundary[kk+1])/dinvr
          endif
	  ; Small extended HZ planets anywhere     
	  allhz =    where((star.planet[jj].r le radius_boundary[kk+1]) and $
	  		   (star.planet[jj].r le 2.0) and $
			   (star.p_hz_out gt period_boundary[jj]) and $
			   (star.p_hz_in  lt period_boundary[jj+1]) and $
			   (star.planet[jj].det))
          if (allhz[0] ne -1) then begin
	    p1 = star[allhz].p_hz_in
            p2 = star[allhz].p_hz_out
            r1 = star[allhz].planet[jj].r
	    r2 = radius_boundary[kk+1]
            p1rep = where(p1 lt period_boundary[jj])
            if (p1rep[0] ne -1) then p1[p1rep] = period_boundary[jj]
            p2rep = where(p2 gt period_boundary[jj+1])
            if (p2rep[0] ne -1) then p2[p2rep] = period_boundary[jj+1]
            r1rep = where(r1 lt radius_boundary[kk])
            if (r1rep[0] ne -1) then r1[r1rep] = radius_boundary[kk]
            npexthz[ii] = npexthz[ii] + rate_fressin[jj,kk] * $
	      total((alog(p2)-alog(p1))*(1./r1-1./r2))/(dlogp*dinvr)
          endif
	  ; Small extended HZ planets anywhere     
	  allhz =    where((star.planet[jj].r le radius_boundary[kk+1]) and $
	  		   (star.planet[jj].r le 2.0) and $
			   (star.p_hz_out gt period_boundary[jj]) and $
			   (star.p_hz_in  lt period_boundary[jj+1]) and $
			   (abs(star.coord.elat) gt 75) and $
			   (star.planet[jj].det))
          if (allhz[0] ne -1) then begin
	    p1 = star[allhz].p_hz_in
            p2 = star[allhz].p_hz_out
            r1 = star[allhz].planet[jj].r
	    r2 = radius_boundary[kk+1]
            p1rep = where(p1 lt period_boundary[jj])
            if (p1rep[0] ne -1) then p1[p1rep] = period_boundary[jj]
            p2rep = where(p2 gt period_boundary[jj+1])
            if (p2rep[0] ne -1) then p2[p2rep] = period_boundary[jj+1]
            r1rep = where(r1 lt radius_boundary[kk])
            if (r1rep[0] ne -1) then r1[r1rep] = radius_boundary[kk]
            npjwstexthz[ii] = npjwstexthz[ii] + rate_fressin[jj,kk] * $
	      total((alog(p2)-alog(p1))*(1./r1-1./r2))/(dlogp*dinvr)
          endif
        ;ns1 = total((star.planet[jj].r lt radius_boundary[kk]) and $
	;	    (star.planet[jj].det))
        sbefore = where((star.planet[jj].r lt radius_boundary[kk]) and $
                 	 (star.planet[jj].det))
        sbtwn =  where((star.planet[jj].r ge radius_boundary[kk]) and $
		 	(star.planet[jj].r lt radius_boundary[kk+1]) and $
                 	(star.planet[jj].det))
        if (sbefore[0] ne -1) then $
	  nsbefore=double(n_elements(sbefore)) $
	  else nsbefore=0.0
        if (sbtwn[0] ne -1) then $
          btwnfrac = (1./star[sbtwn].planet[jj].r - 1./radius_boundary[kk+1])/dinvr $
          else btwnfrac = 0.0
        nstar3d[ii,jj,kk] = nsbefore + total(btwnfrac)
        nplanet3d[ii,jj,kk] = rate_fressin[jj,kk]*nstar3d[ii,jj,kk]
          
      endfor
      
      ; Average over period bins
     ; if (jj lt 1) then begin
     ;   nstar3d[ii,jj,*] = nstar3b[ii,jj,*]
     ; endif else begin
     ; if (jj gt 0) then begin
     ;   nstar3d[ii,jj-1,*] = (nstar3b[ii,jj-1,*]+nstar3b[ii,jj,*])/2.0
     ;   nplanet3d[ii,jj-1,*] = rate_fressin[jj-1,*]*nstar3d[ii,jj-1,*]
     ; endif
    endfor
  endfor
  nstar2d = round(mean(round(nstar3d), dimension=1))
  nstar2derr = ceil(stddev(round(nstar3d), dimension=1))
  nplanet2d = round(mean((nplanet3d), dimension=1))
  per_tot = total(nplanet3d,2)
  earths = per_tot[*,0]
  speths = per_tot[*,1]
  smneps = per_tot[*,2]
    
  nplanet2derr = ceil(stddev((nplanet3d), dimension=1))
  s_hz_mean = mean(s_small_hz)
  s_hz_err = stddev(s_small_hz)
  s_jwst_mean = mean(s_jwst2_hz)
  s_jwst_err = stddev(s_jwst2_hz)
  p_exthz_mean = mean(npexthz)
  p_exthz_err = stddev(npexthz)
  p_jwst_exthz_mean = mean(npjwstexthz)
  p_jwst_exthz_err = stddev(npjwstexthz)
  p_hz_mean = mean(nplanethz)
  p_hz_err = stddev(nplanethz)
  p_jwst_mean = mean(npjwsthz)
  p_jwst_err = stddev(npjwsthz)
  earths_mean = mean(earths,dimension=1)
  earths_err = stddev(earths,dimension=1)
  speths_mean = mean(speths,dimension=1)
  speths_err = stddev(speths,dimension=1)
  smneps_mean = mean(smneps,dimension=1)
  smneps_err = stddev(smneps,dimension=1)
  ;print, 'Stars:'
  ;print, nstar2d
  ;print, 'Planets:'
  ;print, nplanet2d
  ;print, 'Planet Errors:'
  ;print, nplanet2derr
  print, 'Earths: ', earths_mean, ' +/- ', earths_err
  print, 'Super-Earths: ', speths_mean, ' +/- ', speths_err
  print, 'Sub-Neptunes: ', smneps_mean, ' +/- ', smneps_err
  openw, lun, fname + '.tex', /get_lun
  fmt = '(A30, A5, F5.1, A5)'
  printf, lun, f=fmt, 'Field of view [deg]: ', ' & ', fov, ' \\'
  printf, lun, f=fmt, 'Number of Segments: ', ' & ', seg, ' \\'
  printf, lun, f=fmt, 'Geometric Area [cm$^2$]: ', ' & ', geomarea, ' \\'
  printf, lun, f=fmt, 'PSF Scaling from Baseline: ', ' & ', psfsize, ' \\'
  printf, lun, f=fmt, 'Read Noise [e-]: ', ' & ', readnoise, ' \\'
  printf, lun, f=fmt, 'Minimum SNR: ', ' & ', thresh, ' \\'
  printf, lun, f=fmt, 'Saturation Magnitude: ', ' & ', median(satmag), ' \\'
  fmt = '(A30, A5, I8, A5)'
  printf, lun, f=fmt, 'Minimum number of transits: ', ' & ', tranmin, ' \\'
  printf, lun, f=fmt, 'Number of stars: ', ' & ', nstars, ' \\'
  printf, lun, f=fmt, 'Number of trials: ', ' & ', n_trial, ' \\'
  fmt = '(A30, A5, A30, A5)'
  printf, lun, f=fmt, 'Stars Never Observed', ' & ', 'Stars Always Observed', ' \\'
  fmt = '(F5.1, A5, F5.1, A5)'
  printf, lun, f=fmt, nocoverage, '$\%$ & ', fullcoverage, '$\%$ \\'
  fmt = '(A75, I5, A5, I5, A5)'
  ;printf,lun,f=fmt,'dwell_time = ', DWELL_TIME
  ;printf,lun,f=fmt,'geom_area = ', geom_area
  ;printf,lun,f=fmt,'pix_scale = ', pix_scale
  ;printf,lun,f=fmt,'e_pix_ro = ', e_pix_ro
  ;printf,lun,f=fmt,'sys_limit = ', sys_limit
  ;printf,lun,f=fmt,'snr_min = ', SNR_MIN
  ;printf,lun,f=fmt,'ntra_obs_min = ', ntra_obs_min
  ;fmt = '(A25, A40)'
  ;printf,lun,f=fmt,'survey file = ', survey_file

  printf, lun, f=fmt, '$\leq 2R_{\oplus}$ HZ hosts: & ', s_hz_mean, '$\pm$', s_hz_err, '\\'
  printf, lun, f=fmt, '$\leq 2R_{\oplus}$ HZ hosts above 75$^{\circ}$ Lat.:  & ', $
	s_jwst_mean, '$\pm$', s_jwst_err, '\\'
  fmt = '(A75, F8.1, A10, F8.2, A10)'
  printf, lun, f=fmt, '$\leq 2R_{\oplus}$  HZ planets: & ', p_hz_mean, '$\pm$', p_hz_err, '\\'
  printf, lun, f=fmt, '$\leq 2R_{\oplus}$  HZ planets above 75$^{\circ}$ Lat.: & ', p_jwst_mean, '$\pm$', p_jwst_err, '\\'
  printf, lun, f=fmt, 'Mean T$_{\rm eff}$ of all $\leq 2R_{\oplus}$ hosts: & ', $
	round(mean(teff2all)), '$\pm$', round(mean(teff2alls)), '\\'
  printf, lun, f=fmt, 'Mean T$_{\rm eff}$ of $\leq 2R_{\oplus}$ HZ hosts: & ', $
	round(mean(teff2hz)), '$\pm$', round(mean(teff2hzs)), '\\'
  ;printf, 'Teff of <4 Re hosts: ', mean(teff4), '$\pm$', stddev(teff4)
  fmt = '(A20, I5, A5, I5, A5, I6, A5, I6, A5, I7, A5, I7, A5)'
  printf, lun, f=fmt, 'Stars & ', $
	total(nstar2d[*,0]), '$\pm$', ceil(sqrt(total(nstar2derr[*,0]^2.0))), ' & ', $ 
	total(nstar2d[*,1]), '$\pm$', ceil(sqrt(total(nstar2derr[*,1]^2.0))), ' & ', $ 
	total(nstar2d[*,2]), '$\pm$', ceil(sqrt(total(nstar2derr[*,2]^2.0))), ' \\ '
  fmt = '(A20, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5, F8.1, A5, F8.2, A5)'
  printf, lun, f=fmt, 'Planets & ', $
	earths_mean,  '$\pm$', earths_err, ' & ', $ 
	speths_mean, '$\pm$', speths_err, ' & ', $ 
	smneps_mean, '$\pm$', smneps_err, ' & ', $
	p_exthz_mean, '$\pm$', p_exthz_err, ' & ', $
	p_jwst_exthz_mean, '$\pm$', p_jwst_exthz_err, ' & ', $
	p_hz_mean, '$\pm$', p_hz_err, ' & ', $
        p_jwst_mean, '$\pm$', p_jwst_err, ' \\'
  close,lun
  free_lun, lun  

  ;if (keyword_set(nplanet2d)) then nplanet2d=nplanet2d
  ;if (keyword_set(nstar2d)) then nstar2d=nstar2d


  set_plot,'ps'
     device,filen=fname+'_2.eps',/encapsulated,xsize=7.0,ysize=3.2,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.charsize=0.8

  !P.Multi=[0,2,1] 
  ; More plotting
  p = [0.2,0.2,0.9,0.8]
  cgLoadCT, 0, /reverse
  axfmt = {Xticks:9,XTickname:string(round(period_boundary), format='(I3)'), $
           Yticks:5, YTickname:string(radius_boundary, format='(F5.1)')}
  ns2d=alog10(double(round(nstar2d)))*((255.-90.)/5.)+90.
  lo = where((ns2d lt 90) or ~finite(ns2d))
  if (lo[0] ne -1) then ns2d[lo] = 0
  hi = where(ns2d gt 255)
  if (hi[0] ne -1) then ns2d[hi] = 255
  cgImage, ns2d, position=p, xrange=[0,9],yrange=[0,5], $
  ;      minvalue=90, $
  ;      minvalue=0.0, maxvalue=5,$
  ;      bottom=0.0, top=5.0, $ 
	axkeywords=axfmt, /axes, $
	xtitle='Period [Days]', $
        ytitle='Radius [R!Dearth!N]', title='Stars'
  ;cgColorBar,  position=[p[0], p[3]+0.05, p[2], p[3]+0.07], $
  cgColorBar,  position=[p[0], p[3]+0.1, p[2], p[3]+0.15], $
	range=[1,100000], xlog=1, bottom=90
   
 ; nstarr = nstar2d[4,*]
  ;nstarp = mean(nstar3b[*,*,3],dimension=1)

  ;cgPlot, radius_boundary, nstarr, /xlog, /ylog, $
  ;	yrange=[10.0,ceil(max(nstarr))], $
  ;   	xrange=[1.0,10.0], $
  ;	ytitle='Observable Stars', $
  ;	xtitle='R!Dmin!N [R!Dearth!N]'  
  ;cgPlot, period_boundary, nstarp, /xlog, /ylog, $  
  ;	yrange=[10.0,ceil(max(nstarp))], $
  ;     	xrange=[1.0,150.0], $
  ;	ytitle='Observable Stars', $
  ;	xtitle='Period [Days]'  

  p = [0.2,0.2,0.9,0.8]
  cgLoadCT, 0, /reverse 
  np2d=alog10(double(nplanet2d))*((255.-90.)/3.)+90.
  lo = where((np2d lt 90) or ~finite(ns2d))
  if (lo[0] ne -1) then np2d[lo] = 0
  hi = where(np2d gt 255)
  if (hi[0] ne -1) then np2d[hi] = 255
  cgImage, np2d, position=p, xrange=[0,9],yrange=[0,5],$
       ; minvalue=90, $
       ; minvalue=0.0, maxvalue=3,$ 
       ; bottom=0.0, top=3.0, $
        axkeywords=axfmt, /axes, $
	xtitle='Period [Days]', $
	ytitle='Radius [R!Dearth!N]', title='Planets'
  cgColorBar, position=[p[0], p[3]+0.1, p[2], p[3]+0.15], $
	range=[1,1000], xlog=1, bottom=90

  device,/close
  set_plot,'ps'
     device,filen=fname+'_3.eps',/encapsulated,xsize=7.0,ysize=3.2,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.charsize=0.8

  !P.Multi=[0,2,1] 
  elatx = round((elatbins[0:n_elat-1]+(elatbins[1:n_elat]))/2.0)
  names = string(elatx, format='(I3)')
  ;nullnames = indgen(floor(n_mag/2.0))*2+1
  ;names[nullnames] = ' '
 
  elaty = mean(elat2, dimension=1)
  ;gd = where(magy gt 0)
  colors = replicate('black', n_elements(elatx))
  cgbarplot,elaty, barcoords=elatx, barnames=names, colors=colors, $
	xtitle='Ecliptic Latitude [deg]', ytitle='Observable R!Dp!N<2R!Dearth!N HZ Hosts'
 
  gd = where(bary gt 1)
  names = string(barx, format='(I3)')
  colors = replicate('black', n_elements(barx))
  cgbarplot, alog10(bary[gd]), barcoords=barx[gd], barnames=names[gd], colors=colors[gd], $
	xtitle='Number of Pointings', ytitle='log(Number of Stars)'
  
;  magx = (magbins[0:n_mag-1]+(magbins[1:n_mag]))/2.0	
;  names = string(magx, format='(I3)')
;  nullnames = indgen(floor(n_mag/2.0))*2+1
;  names[nullnames] = ' '
  
;  magy = mean(mags2, dimension=1)
;  colors = replicate('black', n_elements(magx))
;  cgbarplot,magy, barcoords=magx, barnames=names, colors=colors, $
;	xtitle='I Magnitude', ytitle='Observable R!Dp!N<2R!Dearth!N HZ Hosts'

  device,/close
  
  save, filen=spo_name, star

  set_plot,'x'
  !p.font=-1
  !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
end

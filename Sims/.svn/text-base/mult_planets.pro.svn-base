pro mult_planets, sstruct=sstruct, pstruct=pstruct, infile=infile, outfile=outfile, verbose=verbose, $
	csr=csr, fressin=fressin, petigura=petigura, dressing=dressing, req=req, nodep=nodep, plotname=plotname

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 0.6395 ; in m/s

  if (keyword_set(sstruct)) then star = sstruct else restore, infile

  nstars = n_elements(star)

  if(keyword_set(dressing)) then begin
    hotstars = where(star.teff ge 4000.)
    nhotstars = total(star.teff ge 4000.)
    coolstars = where(star.teff lt 4000.)
    ncoolstars = total(star.teff lt 4000.)
  endif else begin
    hotstars = lindgen(nstars)
    nhotstars = nstars
    coolstars = -1
    ncoolstars = 0
  endelse

  template_planet = {$
                    n: 0, $   ; planets per star (set to one for now)
                    r: 0.0, $ ; can eventually replace by dblarr(5), $
 		    m: 0.0, $ ; mass
                    p: 0.0, $ ; Period (days)
                    a: 0.0, $ ; Semimajor axis (AU)
                    s: 0.0, $ ; incident flux in units of Sun-->Earth flux
                    b: 0.0, $ ; Impact parameter (0-1)
                    k: 0.0, $ ; rv amplitude
                    tra: 0, $ ; Transit? (boolean)
                    dep: 0.0, $ ; Transit depth (0-1)
                    dur: 0.0, $ ; Transit duration (days)
                    ntra_obs: 0, $ ; Number of transits observed
                    det: 0, $  ; Detected?
                    snr: 0.0, $ ; SNR in phase-folded lightcurve
                    snrtran: 0.0, $ ; SNR per transit
                    durpar: 0.0, $    ; duration of in+engress (days)
                    snrgress: 0.0, $   ; SNR of ingress/egress
 		    hostid: 0L $
                    }



  if (keyword_set(verbose)) then verbose=1 else verbose=0     

  if (keyword_set(csr)) then begin

; Previous CSR assumptions
     nplanets=nstars
     planet = replicate(template_planet, nstars)
     randomp, period_ran, -1.0, nstars, range_x = [2., 500.], seed=seed
     planet.p = period_ran
     randomp, radius_ran, -2.0, nstars, range_x = [0.5, 15.0], seed=seed
     planet.r = radius_ran
     planet.hostid = lindgen(nstars)

  endif else if (keyword_set(fressin)) then begin
     period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
     radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0]
     planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']
     
     rate_fressin = dblarr(11,5) ; period bin, radius bin
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
     rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
     rate_fressin = rate_fressin/100.
 
     rate_dressing = rate_fressin ;just for now
     
     nplanets = total(round(rate_fressin * nhotstars)) + $
		total(round(rate_dressing * ncoolstars))
     ; Pre-allocate for speed
     planet = replicate(template_planet, total(nplanets))
     idx0 = 0L
     for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           binplanets = round(rate_fressin[i,j] * nhotstars)
           if (binplanets gt 0) then begin
	      tmp_planet = replicate(template_planet, binplanets)
	      tmp_planet.hostid = hotstars[floor(double(nhotstars)*randomu(seed, binplanets))]
              if(j eq 0) then radpow = 0.0 else radpow = -1.7
              randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
              randomp, radii, radpow, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              tmp_planet.r = radii
              tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
              ;print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
              planet[idx] = tmp_planet
              delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
     if (keyword_set(dressing)) then begin
       for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           binplanets = round(rate_dressing[i,j] * ncoolstars)
           if (binplanets gt 0) then begin
	      tmp_planet = replicate(template_planet, binplanets)
	      tmp_planet.hostid = coolstars[floor(double(ncoolstars)*randomu(seed, binplanets))]
              if(j eq 0) then radpow = 0.0 else radpow = -1.7
              randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
              randomp, radii, radpow, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              tmp_planet.r = radii
              tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
              ;print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
              planet[idx] = tmp_planet
              delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
     endif
  endif else if (keyword_set(petigura)) then begin
    radius_boundary = [1.0, 2.0, 4.0, 8.0, 16.0]
    period_boundary = [0.8, 2.0, 3.4, 6.25, 12.5, 25.0, 50.0, 100.0, 200.0, 400.0]
    rate_petigura = dblarr(9,4) ; period bin, radius bin
    rate_petigura[0,*] = [0.17, 0.035, 0.4, 0.015]
    rate_petigura[1,*] = [0.74, 0.18, 0.006, 0.067]
    rate_petigura[2,*] = [1.49, 0.73, 0.11, 0.17]
    rate_petigura[3,*] = [4.9, 3.5, 0.3, 0.2]
    rate_petigura[4,*] = [6.6, 6.1, 0.8, 0.2]
    rate_petigura[5,*] = [7.7, 7.0, 0.4, 0.6]
    rate_petigura[6,*] = [5.8, 7.5, 1.3, 0.6]
    rate_petigura[7,*] = [3.2, 6.2, 2.0, 1.1]
    rate_petigura[8,*] = [0.0, 5.0, 1.6, 1.3]
    rate_petigura = rate_petigura/100.0
    nplanets = total(round(rate_petigura * nstars))
    ; Pre-allocate for speed
    planet = replicate(template_planet, total(nplanets))
    idx0 = 0L
    for i=8,0,-1 do begin
        for j=3,0,-1 do begin
           binplanets = round(rate_petigura[i,j] * nstars)
           if (binplanets gt 0) then begin
	      tmp_planet = replicate(template_planet, binplanets)
	      tmp_planet.hostid = floor(double(nstars)*randomu(seed, binplanets))
  	      if (keyword_set(dressing)) then begin
		periods = fltarr(binplanets) + 10^((alog10(period_boundary[i])+alog10(period_boundary[i+1]))/2.0)
		radii   = fltarr(binplanets) + 10^((alog10(radius_boundary[j])+alog10(radius_boundary[j+1]))/2.0)
              endif else begin
                randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
                randomp, radii, -2.0, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              endelse 
              tmp_planet.r = radii
              tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
              print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
              planet[idx] = tmp_planet
              delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
  endif else if(keyword_set(req)) then begin 
    planet = replicate(template_planet, nstars)
    planet.hostid = lindgen(nstars)
    planet.r = 2.0
    planet.p = float(req)
  endif

; Star properties to be randomized
  star.dx = floor(10.*randomu(seed, nstars))
  star.dy = floor(10.*randomu(seed, nstars))
; Random orbital orientation
  if (keyword_set(req)) then begin
    star.cosi = 0.0 ;star[pla].r/(AU_IN_RSUN * star[pla].planet.a) * (-1.0 + 2.0*randomu(seed, nplanets));
  endif else begin
    star.cosi = -1.0 + 2.0*randomu(seed, nstars)
  endelse

; Work out orbital distance and impact parameter
  allid = planet.hostid
  planet.a = (star[allid].m)^(1./3.) * (planet.p/365.25)^(2./3.); in AU
  planet.s = (star[allid].r)^2.0 * (star[allid].teff/5777.0)^4. / (planet.a)^2. ; indicent flux wrt sun-earth value

  planet.b = (planet.a*AU_IN_RSUN / star[allid].r) * star[allid].cosi; assumes circular orbit
 
  planet.m = 10.55*4.*!dpi/3.*(planet.r/3.9)^3. * $ ; assume MgSiO_3 composition (Seager 2007)
     (1.0 + (1.0 - 3./5.*0.541)*(2./3.*!dpi*(planet.r/3.9)^2.)^0.541)

; RV amplitude
  planet.k = RV_AMP*planet.p^(-1./3.) * planet.m * $ 
	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 

	;2.0*!dpi*sqrt(1.0-star[pla].cosi^2.)*star[pla].planet.a * AU_DAY_IN_CM_S * planet_mass /  $
	;	(star[pla].planet.p * star[pla].m * MSUN_IN_MEARTH)

; Work out transit properties
  tra = where(abs(planet.b) lt 1.0)
  traid = planet[tra].hostid
  planet[tra].tra = 1
  planet[tra].dep = (REARTH_IN_RSUN * planet[tra].r / star[traid].r )^2.0
  planet[tra].dur = star[traid].r * planet[tra].p * sqrt(1.-(planet[tra].b)^2.) / (!PI*planet[tra].a*AU_IN_RSUN)
  planet[tra].durpar = planet[tra].r * $
	REARTH_IN_RSUN * planet[tra].p / $
        sqrt(1.-(planet[tra].b)^2.) / $
        (!PI*planet[tra].a*AU_IN_RSUN)
  
  print, 'Created ', n_elements(tra), ' transiting planets out of ', nplanets, ' total.'

  if (keyword_set(outfile)) then save,filen=outfile, planet

  pstruct=planet

  if (verbose) then begin 
; Report breakdown of different types of planets

  print, ' ' & print, 'Planet breakdown:'

  readcol, 'planet_definitions.txt', /silent, comment='#', f='A,D,D,D,D', $
           planet_type, p_min, p_max, r_min, r_max
           
  fmt = '(A20, A9, I3, A3, I3, A9, F5.2, A3, F5.2, A4, I9, F7.2, A1)'
  for i=0,n_elements(planet_type)-1 do begin
     q = where(planet.p ge p_min[i] and planet.p lt p_max[i] and $
               planet.r gt r_min[i] and planet.r lt r_max[i])
     printf,-1, f=fmt, planet_type[i], ' (Period ', p_min[i], ' - ', p_max[i], ', Radius ', r_min[i], ' - ', r_max[i], $
            ') = ', n_elements(q), 100.*float(n_elements(q))/float(n_elements(planet)), '%'
  endfor
  
; Report breakdown of different types of TRANSITING planets

  print, ' ' & print, 'Transiting Planet breakdown:'

  readcol, 'planet_definitions.txt', /silent, comment='#', f='A,D,D,D,D', $
           planet_type, p_min, p_max, r_mn, r_max
           
  fmt = '(A20, A9, I3, A3, I3, A9, F5.2, A3, F5.2, A4, I9, F7.2, A1)'
  for i=0,n_elements(planet_type)-1 do begin
     q = where(planet.p ge p_min[i] and planet.p lt p_max[i] and $
               planet.r gt r_min[i] and planet.r lt r_max[i] and planet.tra eq 1)
     printf,-1, f=fmt, planet_type[i], ' (Period ', p_min[i], ' - ', p_max[i], ', Radius ', r_min[i], ' - ', r_max[i], $
            ') = ', n_elements(q), 100.*float(n_elements(q))/float(n_elements(planet)), '%'
  endfor

  !p.multi=[0,2,2]
  plotsym,0,/fill
  !p.charsize=2

  plotind = floor(double(nplanets)*randomu(seed, round(double(nplanets)/100.)))
  print, 'Plotting ', n_elements(plotind), ' points'
  plot, planet[plotind].p, planet[plotind].r, psym=3, /ynozero, $
        title='All planets', $
        /xlog, xsty=1, xra=[0.5,100], $
      ;  ysty=1, yra=[0.5,4.5], $
        xtit='Period [d]', ytit='Planet radius'

  plot, planet[tra].p, planet[tra].r, psym=3, /ynozero, $
        title='All transiting planets', $
        /xlog, xsty=1, xra=[0.5,100], $
      ;  ysty=1, yra=[0.5,4.5], $
        xtit='Period [d]', ytit='Planet radius'

  plot, planet[tra].p, 24.0*planet[tra].dur, psym=3, $
        title='All Transiting Planets', $
        xtit='Period [days]', /xlog, $
        ytit='Transit Duration [hr]', /ynozero

  plot, planet[tra].p, planet[tra].dep, psym=3, $
        title='All Transiting Planets', $
        xtit='Period [days]', /xlog, $
        /ylog, /ynozero, $
         ytit='Transit Depth'

  endif


  if (keyword_set(plotname)) then begin
    ; Plot the period and radius distributions 
    set_plot, 'ps'
    device,filen=plotname,/encapsulated,xsize=7.0,ysize=7.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
    !p.font=0
    !p.charsize=0.8
    !p.multi=[0,2,2]
    plotsym,0,/fill

  plotind = floor(double(nplanets)*randomu(seed, round(double(nplanets)/10.)))
  print, 'Plotting ', n_elements(plotind), ' points'
  plot, planet[plotind].p, planet[plotind].r, psym=3, /ynozero, $
        title='10% of all planets', $
        /xlog, /ylog, xsty=1, xra=[0.5,500], $
        ysty=1, yra=[0.5,30], $
        xtit='Period [d]', ytit='Planet radius'

  cgHistoPlot, alog10(planet.r), nbins=80, xtitle='log(Planet Radius)', ytitle='', title='Radius Distribution'

  cgHistoPlot, alog10(planet.p), nbins=80, xtitle='log(Period [days])', ytitle='', title='Period Distribution'

  plot, planet[tra].p, planet[tra].r, psym=3, /ynozero, $
        title='All transiting planets', $
        /xlog, /ylog, xsty=1, xra=[0.5,500], $
        ysty=1, yra=[0.5,30], $
        xtit='Period [d]', ytit='Planet radius'

  device,/close

  set_plot,'x'
  !p.font=-1
  !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1

  endif
end

pro new_planets, struct=struct, infile=infile, outfile=outfile, verbose=verbose, $
	csr=csr, fressin=fressin, petigura=petigura, dressing=dressing, req=req, nodep=nodep

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 0.6395 ; in m/s

  if (keyword_set(infile)) then infile=infile else infile = 'ss.sav'
  restore, infile

  nstars = n_elements(star)

  nhotstars = total(star.teff ge 4000.)
  ncoolstars = total(star.teff lt 4000.)

  if (keyword_set(outfile)) then fname=outfile else fname='sp.sav'
  if (keyword_set(verbose)) then verbose=1 else verbose=0     

 if (keyword_set(csr)) then begin

; Previous CSR assumptions

     star.planet.n = 1
     randomp, period_ran, -1.0, nstars, range_x = [2., 500.], seed=seed
     star.planet.p = period_ran
     randomp, radius_ran, -2.0, nstars, range_x = [0.5, 15.0], seed=seed
     star.planet.r = radius_ran

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
     
     for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           nplanets = round(rate_fressin[i,j] * nstars)
           if (nplanets gt 0) then begin
              indices = round(double(nstars-1)*randomu(seed, nplanets))
  	      if (keyword_set(dressing)) then begin
		periods = fltarr(nplanets) + 10^((alog10(period_boundary[i])+alog10(period_boundary[i+1]))/2.0)
		radii   = fltarr(nplanets) + 10^((alog10(radius_boundary[j])+alog10(radius_boundary[j+1]))/2.0)
              endif else begin
                randomp, periods, -1.0, nplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
                randomp, radii, -2.0, nplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              endelse 
	      star[indices].planet.n = 1
              star[indices].planet.r = radii
              star[indices].planet.p = periods
           endif
        endfor
     endfor
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
      for i=8,0,-1 do begin
        for j=3,0,-1 do begin
          nplanets = round(rate_petigura[i,j] * nstars)
          if (nplanets gt 0) then begin
            indices = round(double(nstars-1)*randomu(seed, nplanets))
            randomp, periods, -1.0, nplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
            randomp, radii, 0.0, nplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
            star[indices].planet.n = 1
            star[indices].planet.r = radii
            star[indices].planet.p = periods
          endif
        endfor
     endfor
  endif else if(keyword_set(req)) then begin   
    star.planet.n = 1
    star.planet.r = 2.0
    star.planet.p = float(req)
  endif

     

  pla = where(star.planet.n gt 0)
  nplanets = n_elements(pla)

  ; Random pixel offset
  star.dx = floor(10.*randomu(seed, nstars))
  star.dy = floor(10.*randomu(seed, nstars))
  ; HZ
  a_hz_in = star.r * (star.teff/5777.0)^2. / sqrt(1.5)
  a_hz_out= star.r * (star.teff/5777.0)^2. / sqrt(0.5)
  star.p_hz_in =  365.25 * sqrt(a_hz_in^3.  * star.m)
  star.p_hz_out=  365.25 * sqrt(a_hz_out^3. * star.m)


; Work out orbital distance and impact parameter

  star[pla].planet.a = (star[pla].m)^(1./3.) * (star[pla].planet.p/365.25)^(2./3.); in AU
  star[pla].planet.s = (star[pla].r)^2.0 * (star[pla].teff/5777.0)^4. / (star[pla].planet.a)^2. ; indicent flux wrt sun-earth value

; Random orbital orientation
  if (keyword_set(req)) then begin
	star[pla].cosi = 0.0 ;star[pla].r/(AU_IN_RSUN * star[pla].planet.a) * (-1.0 + 2.0*randomu(seed, nplanets));
  endif else begin
    star[pla].cosi = -1.0 + 2.0*randomu(seed, nplanets)
  endelse

  star[pla].planet.b = (star[pla].planet.a*AU_IN_RSUN / star[pla].r) * star[pla].cosi; assumes circular orbit
 
  planet_mass = 10.55*4.*!dpi/3.*(star[pla].planet.r/3.9)^3. * $ ; assume MgSiO_3 composition (Seager 2007)
     (1.0 + (1.0 - 3./5.*0.541)*(2./3.*!dpi*(star[pla].planet.r/3.9)^2.)^0.541)

; RV amplitude
  star[pla].planet.k = RV_AMP*star[pla].planet.p^(-1./3.) * planet_mass * $ 
	sqrt(1.0-star[pla].cosi^2.) * (star[pla].m)^(-2./3.) 

	;2.0*!dpi*sqrt(1.0-star[pla].cosi^2.)*star[pla].planet.a * AU_DAY_IN_CM_S * planet_mass /  $
	;	(star[pla].planet.p * star[pla].m * MSUN_IN_MEARTH)

; Work out transit properties

  tra = where(star.planet.n gt 0 and abs(star.planet.b) lt 1.0)
  star[tra].planet.tra = 1
  if (keyword_set(nodep)) then star[tra].planet.dep = 0.0 else $
    star[tra].planet.dep = (REARTH_IN_RSUN * star[tra].planet.r / star[tra].r )^2.0
  star[tra].planet.dur = star[tra].r * star[tra].planet.p * sqrt(1.-(star[tra].planet.b)^2.) / (!PI*star[tra].planet.a*AU_IN_RSUN)

  star[tra].planet.durpar = star[tra].planet.r * $
	REARTH_IN_RSUN * star[tra].planet.p / $
        sqrt(1.-(star[tra].planet.b)^2.) / $
        (!PI*star[tra].planet.a*AU_IN_RSUN)

  if keyword_set(outfile) then save,filen=outfile, star

  if (keyword_set(struct)) then struct=star

  if (verbose) then begin 
; Report breakdown of different types of planets

  print, ' ' & print, 'Planet breakdown:'

  readcol, 'planet_definitions.txt', /silent, comment='#', f='A,D,D,D,D', $
           planet_type, p_min, p_max, r_min, r_max
           
  fmt = '(A20, A9, I3, A3, I3, A9, F5.2, A3, F5.2, A4, I9, F7.2, A1)'
  for i=0,n_elements(planet_type)-1 do begin
     q = where(star.planet.p ge p_min[i] and star.planet.p lt p_max[i] and $
               star.planet.r gt r_min[i] and star.planet.r lt r_max[i])
     printf,-1, f=fmt, planet_type[i], ' (Period ', p_min[i], ' - ', p_max[i], ', Radius ', r_min[i], ' - ', r_max[i], $
            ') = ', n_elements(q), 100.*float(n_elements(q))/float(n_elements(star)), '%'
  endfor
  
; Report breakdown of different types of TRANSITING planets

  print, ' ' & print, 'Transiting Planet breakdown:'

  readcol, 'planet_definitions.txt', /silent, comment='#', f='A,D,D,D,D', $
           planet_type, p_min, p_max, r_mn, r_max
           
  fmt = '(A20, A9, I3, A3, I3, A9, F5.2, A3, F5.2, A4, I9, F7.2, A1)'
  for i=0,n_elements(planet_type)-1 do begin
     q = where(star.planet.p ge p_min[i] and star.planet.p lt p_max[i] and $
               star.planet.r gt r_min[i] and star.planet.r lt r_max[i] and star.planet.tra eq 1)
     printf,-1, f=fmt, planet_type[i], ' (Period ', p_min[i], ' - ', p_max[i], ', Radius ', r_min[i], ' - ', r_max[i], $
            ') = ', n_elements(q), 100.*float(n_elements(q))/float(n_elements(star)), '%'
  endfor

; Plot the period and radius distributions

  !p.multi=[0,2,2]
  plotsym,0,/fill
  !p.charsize=2

  plot, star.planet.p, star.planet.r, psym=3, /ynozero, $
        title='All planets', $
        /xlog, xsty=1, xra=[0.5,100], $
      ;  ysty=1, yra=[0.5,4.5], $
        xtit='Period [d]', ytit='Planet radius'

  plot, star[tra].planet.p, star[tra].planet.r, psym=3, /ynozero, $
        title='All transiting planets', $
        /xlog, xsty=1, xra=[0.5,100], $
      ;  ysty=1, yra=[0.5,4.5], $
        xtit='Period [d]', ytit='Planet radius'

  plot, star[tra].planet.p, 24.0*star[tra].planet.dur, psym=3, $
        title='All Transiting Planets', $
        xtit='Period [days]', /xlog, $
        ytit='Transit Duration [hr]', /ynozero

  plot, star[tra].planet.p, star[tra].planet.dep, psym=3, $
        title='All Transiting Planets', $
        xtit='Period [days]', /xlog, $
        /ylog, /ynozero, $
         ytit='Transit Depth'

  endif

end

pro planets, csr=csr, fressin=fressin, outfile=outfile

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii

  restore, 'spws.sav'

  nstars = n_elements(star)
  if (keyword_set(outfile)) then fname=outfile else fname='sp.sav'
  if (keyword_set(csr)) then begin
     
; Previous CSR assumptions

     star.planet.n = 1
     randomp, period_ran, -1.0, nstars, range_x = [2., 500.], seed=seed
     star.planet.p = period_ran
     randomp, radius_ran, -2.0, nstars, range_x = [0.5, 15.0], seed=seed
     star.planet.r = radius_ran

  endif else if (keyword_set(fressin)) then begin

     period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
     radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 15.0] ; Fressin gives 22 as upper limit but we lower it here to 15
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
     
     for i=0,10 do begin
        for j=0,4 do begin
           nplanets = round(rate_fressin[i,j] * nstars)
           if (nplanets gt 0) then begin
              indices = round(double(nstars-1)*randomu(seed, nplanets))
              randomp, periods, -1.0, nplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
              randomp, radii, -2.0, nplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              star[indices].planet.n = 1
              star[indices].planet.r = radii
              star[indices].planet.p = periods
           endif
        endfor
     endfor
     
  endif else begin

     np = 20
     logp = -0.5 + 2.5*dindgen(np)/double(np-1)
     dlogp = logp[1]-logp[0]
     p = 10^logp

     r1 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
     r2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
     dr = 0.5
     nr = n_elements(r1)

     occ = dblarr(np,nr)
  
     for j=0,nr-1 do begin
        for i=0,np-1 do begin
           logp1 = logp[i]-0.5*dlogp
           logp2 = logp[i]+0.5*dlogp
           p1 = 10^logp1
           p2 = 10^logp2
           occ = int_planet_occurrence(p1,p2,r1[j],r2[j]) * (0.5/dr) * (0.1/dlogp)
           nplanets = round(occ * nstars)
           if (nplanets gt 0) then begin
              indices = round(double(nstars-1)*randomu(seed, nplanets))
              randomp, periods, -1.0, nplanets, range_x = [p1, p2], seed=seed
              randomp, radii, -1.0, nplanets, range_x = [r1[j], r2[j]], seed=seed
              star[indices].planet.n = 1
              star[indices].planet.r = radii
              star[indices].planet.p = periods
           endif
        endfor
     endfor

  endelse

; Random orbital orientation

  pla = where(star.planet.n gt 0)
  nplanets = n_elements(pla)

  ;star[pla].planet.cosi = -1.0 + 2.0*randomu(seed, nplanets)

; Work out orbital distance and impact parameter

  star[pla].planet.a = (star[pla].m)^(1./3.) * (star[pla].planet.p/365.25)^(2./3.); in AU
  star[pla].planet.s = (star[pla].r)^2.0 * (star[pla].teff/5777.0)^4. / (star[pla].planet.a)^2. ; indicent flux wrt sun-earth value
  star[pla].planet.b = (star[pla].planet.a*AU_IN_RSUN / star[pla].r) * star[pla].planet.cosi; assumes circular orbit

; Work out transit properties

  tra = where(star.planet.n gt 0 and abs(star.planet.b) lt 1.0)
  star[tra].planet.tra = 1
  star[tra].planet.dep = (REARTH_IN_RSUN * star[tra].planet.r / star[tra].r )^2.0
  star[tra].planet.dur = star[tra].r * star[tra].planet.p * sqrt(1.-(star[tra].planet.b)^2.) / (!PI*star[tra].planet.a*AU_IN_RSUN)

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
           planet_type, p_min, p_max, r_min, r_max
           
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

  save,filen=fname,star

end

pro koiplanets, outfile=outfile, verbose=verbose

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  
  koidat = mrdfits('koi.fits')
  restore, 's_pws.sav'
  
  n_rad = 5
  n_per = 11
  
  psf_x = mrdfits('../xPSF.fits')
  psf_y = mrdfits('../yPSF.fits')
  psf_f = mrdfits('../OnAxisPSF.fits')



  n_stars = n_elements(star)
  if (keyword_set(outfile)) then fname=outfile else fname='spp_koi.sav'

  if (keyword_set(verbose)) then v=1 else v=0
;  if (keyword_set(csr)) then begin
     
; Previous CSR assumptions

;     star.planet.n = 1
;     randomp, period_ran, -1.0, nstars, range_x = [2., 500.], seed=seed
;     star.planet.p = period_ran
;     randomp, radius_ran, -2.0, nstars, range_x = [0.5, 15.0], seed=seed
;     star.planet.r = radius_ran

;  endif else if (keyword_set(fressin)) then begin

  period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
  radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 15.0] ; Fressin gives 22 as upper limit but we lower it here to 15
  planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']
  period_mean = (period_boundary[0:n_per-1] + period_boundary[1:n_per])/2.0
  radius_mean = (radius_boundary[0:n_rad-1] + radius_boundary[1:n_rad])/2.0
     
;  rate_fressin = dblarr(11,5) ; period bin, radius bin
;  rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
;  rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
;  rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
;  rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
;  rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
;  rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
;  rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
;  rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
;  rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
;  rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
;  rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
;  rate_fressin = rate_fressin/100.
     
  ; Period-dependent quantities
  for ii=0,n_per-1 do begin
    star.planets.p[ii] = period_mean[ii]
    star.planets.a[ii] = (star.m)^(1./3.) * (star.planets.p[ii]/365.25)^(2./3.)	; in AU
    star.planets.s[ii] = (star.r)^2.0 * (star.teff/5777.0)^4. / (star.planets.a[ii])^2. ; indicent flux wrt sun-earth value
    star.planets.b[ii] = (star.planets.a[ii]*AU_IN_RSUN / star.r) * star.cosi; assumes circular orbit  
    star.planets.tra[ii] = (abs(star.planets.b[ii]) lt 1.0)	
    star.planets.dur[ii] = star.r * star.planets.p[ii] * sqrt(1.-(star.planets.b[ii])^2.) / (!DPI*star.planets.a[ii]*AU_IN_RSUN)
  endfor
  if (v) then print, 'Done assigning period-dependent quantities'
  
  ; Radius-dependent quantities
  for jj=0,n_rad-1 do begin  	
    star.planets.r[jj]   = radius_mean[jj]
    star.planets.dep[jj] = (REARTH_IN_RSUN * star.planets.r[jj] / star.r )^2.0
  endfor
  if (v) then print, 'Done assigning radius-dependent quantities'

  save,filen=fname,star

  for ii=0, n_stars-1 do begin
    if(total(star[ii].planets.tra) gt 0) then begin
      calc_opt_npix_offset, star[ii].mag.i, 3600., npix, snr, $
        psf_f = psf_f, psf_x = psf_x, psf_y = psf_y, $
        dxdy=[star[ii].dx,star[ii].dy], $
        ;saturation=SATURATION, $
        elon=star[ii].coord.elon, $
        elat=star[ii].coord.elat, $
        teff=star[ii].teff, fluxfrac=frac
      star[ii].npix = npix
      star[ii].frac = frac
      star[ii].snr  = snr
      if(v) then print, 'Star ', ii, ' NPIX = ', npix
    endif
  endfor

  save,filen=fname,star

end

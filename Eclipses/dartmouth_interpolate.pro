PRO dartmouth_interpolate, starstruct, mass_in, age_in, feh_in, $ ;inputs
    rad=rad, ic=ic, teff=teff, v=v, rc=rc, j=j, h=h, ks=ks ; outputs
 ages = [0.250, 0.299, 0.349, 0.400, 0.450, 0.500, 0.550, 0.599, 0.650, 0.699, $
         0.750, 0.800, 0.849, 0.900, 0.949, 1.000, 1.250, 1.500, 1.750, 2.000, $ 
         2.250, 2.500, 2.750, 3.000, 3.250, 3.500, 3.750, 4.000, 4.250, 4.500, $
	 4.750, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 8.500, 9.000, $
         9.500, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00, 14.50, 15.00]
 nage = n_elements(ages)
 fehs = [-1.0, -0.5, 0.0, 0.2]
 nfe = n_elements(fehs)
 masses = (findgen(223) + 20.)/200.
 nmass = n_elements(masses)
 
 ai = (interpol(findgen(nage), ages, age_in))
 fi = (interpol(findgen(nfe),  fehs, feh_in))
 mi = (interpol(findgen(nmass), masses, mass_in))
 ; error checking
 toohigh = where(ai gt nage-1)
 if (toohigh[0] ne -1) then ai[toohigh] = nage-1

 toohigh = where(fi gt nfe-1)
 if (toohigh[0] ne -1) then fi[toohigh] = nfe-1
 
 toohigh = where(mi gt nmass-1)
 if (toohigh[0] ne -1) then mi[toohigh] = nmass-1
 
 toolow = where(ai lt 0)
 if (toolow[0] ne -1) then ai[toolow] = 0
 
 toolow = where(fi lt 0)
 if (toolow[0] ne -1) then fi[toolow] = 0
 
 toolow = where(mi lt 0)
 if (toolow[0] ne -1) then mi[toolow] = 0
 
; mi = lonarr(n_elements(mass_in))
; for ii=0,n_elements(mass_in)-1 do begin
;   masses = starstruct[*, fi[ii], ai[ii]].m
;   mi[ii] = interpol(findgen(nmass), masses, mass_in[ii])
; end
; massinds = rebin(findgen(nmass), nmass, nfe, nage)
; mi = interpolate(massinds, mass_in, fi, ai)

; toohigh = where(mi gt nmass-1)
; if (toohigh[0] ne -1) then mi[toohigh] = nmass-1
 
; toolow = where(mi lt 0)
; if (toolow[0] ne -1) then mi[toolow] = 0
 
   d_age = starstruct.age 
   d_feh = starstruct.feh
     d_m = starstruct.m 
  d_teff = starstruct.teff 
   d_rad = starstruct.rad 
  d_logL = starstruct.logL
     d_u = starstruct.u 
     d_b = starstruct.b 
     d_v = starstruct.v 
     d_r = starstruct.r 
    d_ic = starstruct.ic
     d_j = starstruct.j 
     d_h = starstruct.h  
    d_ks = starstruct.ks
  
  rad = interpolate(d_rad, mi, fi, ai) 
  teff = interpolate(d_teff, mi, fi, ai) 
  v = interpolate(d_v, mi, fi, ai) 
  rc = interpolate(d_r, mi, fi, ai) 
  ic = interpolate(d_ic, mi, fi, ai) 
  j = interpolate(d_j, mi, fi, ai) 
  h = interpolate(d_h, mi, fi, ai) 
  ks = interpolate(d_ks, mi, fi, ai)


END	


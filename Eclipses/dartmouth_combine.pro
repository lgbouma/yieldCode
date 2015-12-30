PRO dartmouth_combine, fpath, starstruct=starstruct, maxm=maxm, minm=minm
 ages = [0.250, 0.299, 0.349, 0.400, 0.450, 0.500, 0.550, 0.599, 0.650, 0.699, $
         0.750, 0.800, 0.849, 0.900, 0.949, 1.000, 1.250, 1.500, 1.750, 2.000, $ 
         2.250, 2.500, 2.750, 3.000, 3.250, 3.500, 3.750, 4.000, 4.250, 4.500, $
	 4.750, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 8.500, 9.000, $
         9.500, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00, 14.50, 15.00]
 nage = n_elements(ages)
 masses = (findgen(223) + 20.)/200.
 nmass = n_elements(masses)
 feh = [-1.0, -0.5, 0.0, 0.2]
 nfe = n_elements(feh)
; nmass = 265L
 starstruct = replicate({dartstruct}, nfe*nage*nmass)
 starstruct = reform(starstruct, nmass, nfe, nage) 
 maxm = fltarr(nfe,nage)
 minm = fltarr(nfe,nage)
 for fi = 0,nfe-1 do begin
   if (feh[fi] lt 0) then fes='m' else fes='p' 
   for ai = 0,nage-1 do begin     
     fstub =  fpath +'a'+string(ages[ai]*1E3, format='(I05)')+'feh'+fes+string(abs(feh[fi]*10),'(I02)')+'*'
     fname = file_search(fstub)
     if ((n_elements(fname) gt 1) or (fname eq '')) then print, 'wrong file name from '+fstub else begin
       ;  print, 'Found '+fname+' from '+fstub
       ;  #EEP   M/Mo    LogTeff  LogG   LogL/Lo U       B       V       R       I       J       H       Ks      Kp      D51
       readcol, fname, eep, mass, logT, logG, logL, U, B, V, R, Ic, J, H, Ks, Kp, D51
       
 ;      starstruct[*,fi,ai].age = ages[ai]
 ;      starstruct[*,fi,ai].feh = feh[fi]
 ;      starstruct[*,fi,ai].m = mass[0:nmass-1]
 ;      starstruct[*,fi,ai].teff = 10.^logT[0:nmass-1]
 ;      starstruct[*,fi,ai].rad = sqrt(mass[0:nmass-1])/sqrt((10.^logG[0:nmass-1])/27542.3)
 ;      starstruct[*,fi,ai].logL = logL[0:nmass-1]
 ;      starstruct[*,fi,ai].u  = U[0:nmass-1]
 ;      starstruct[*,fi,ai].b  = B[0:nmass-1]
 ;      starstruct[*,fi,ai].v  = V[0:nmass-1]
 ;      starstruct[*,fi,ai].r  = R[0:nmass-1]
 ;      starstruct[*,fi,ai].ic = Ic[0:nmass-1]
 ;      starstruct[*,fi,ai].j  = J[0:nmass-1]
 ;      starstruct[*,fi,ai].h  = H[0:nmass-1]
 ;      starstruct[*,fi,ai].ks = Ks[0:nmass-1]
       rad = sqrt(mass)/sqrt((10.^logG)/27542.3)      
;       print, 'FE=', feh[fi], ' A= ', ages[ai], ' Min R=', min(rad), ' Max R=', max(rad) 
;       if (min(rad) gt 0.35) then stop
       starstruct[*,fi,ai].age = ages[ai]
       starstruct[*,fi,ai].feh = feh[fi]
       starstruct[*,fi,ai].m = masses
       starstruct[*,fi,ai].rad = interpol(rad, mass, masses)
       starstruct[*,fi,ai].rad = starstruct[*,fi,ai].rad < max(rad)
       starstruct[*,fi,ai].rad = starstruct[*,fi,ai].rad > 0.08
       starstruct[*,fi,ai].logL = interpol(logL, mass, masses)
       starstruct[*,fi,ai].logL = starstruct[*,fi,ai].logL < max(logL)
;       starstruct[*,fi,ai].logL = starstruct[*,fi,ai].logL > min(logL)
       starstruct[*,fi,ai].teff = 10.^(interpol(logT, mass, masses))
       starstruct[*,fi,ai].teff = starstruct[*,fi,ai].teff < 10.^(max(logT))
;       starstruct[*,fi,ai].teff = starstruct[*,fi,ai].teff > 10.^(min(logT))
       tc = where(starstruct[*,fi,ai].teff lt 2000.)
       minT = 5777.*10.^(0.25*starstruct[tc,fi,ai].logL - $
			0.5*alog10(starstruct[tc,fi,ai].rad))
       starstruct[tc, fi, ai].teff = minT

       starstruct[*,fi,ai].u  = interpol(U, mass, masses)
       starstruct[*,fi,ai].u = starstruct[*,fi,ai].u < 30. ;max(U)
       starstruct[*,fi,ai].u = starstruct[*,fi,ai].u > min(U)
       starstruct[*,fi,ai].b  = interpol(B, mass, masses)
       starstruct[*,fi,ai].b = starstruct[*,fi,ai].b < 30. ;max(B)
       starstruct[*,fi,ai].b = starstruct[*,fi,ai].b > min(B)
       starstruct[*,fi,ai].v  = interpol(V, mass, masses)
       starstruct[*,fi,ai].v = starstruct[*,fi,ai].v < 30. ;max(V)
       starstruct[*,fi,ai].v = starstruct[*,fi,ai].v > min(V)
       starstruct[*,fi,ai].r  = interpol(R, mass, masses)
       starstruct[*,fi,ai].r = starstruct[*,fi,ai].r < 30. ;max(R)
       starstruct[*,fi,ai].r = starstruct[*,fi,ai].r > min(R)
       starstruct[*,fi,ai].ic = interpol(Ic, mass, masses)
       starstruct[*,fi,ai].ic = starstruct[*,fi,ai].ic < 30. ;max(Ic)
       starstruct[*,fi,ai].ic = starstruct[*,fi,ai].ic > min(Ic)
       starstruct[*,fi,ai].j  = interpol(J, mass, masses)
       starstruct[*,fi,ai].j = starstruct[*,fi,ai].j < 30. ;max(J)
       starstruct[*,fi,ai].j = starstruct[*,fi,ai].j > min(J)
       starstruct[*,fi,ai].h = interpol(H, mass, masses)
       starstruct[*,fi,ai].h = starstruct[*,fi,ai].h < 30. ;max(H)
       starstruct[*,fi,ai].h = starstruct[*,fi,ai].h > min(H)
       starstruct[*,fi,ai].ks = interpol(Ks, mass, masses)
       starstruct[*,fi,ai].ks = starstruct[*,fi,ai].ks < 30. ;max(Ks)
       starstruct[*,fi,ai].ks = starstruct[*,fi,ai].ks > min(Ks)
      
       maxm[fi,ai] = max(mass)
       minm[fi,ai] = min(mass)

     end ; if file found
   end ; age loop
 end ; feh loop
END	


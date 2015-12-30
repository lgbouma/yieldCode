PRO maglim, fstub, fname, dmax=dmax, icmax=icmax, kmax=kmax, homies=homies
  fnames = file_search(fstub)
  dpmult=100
;  restore, 'dartmouth_grid.sav'
;  tt = mrdfits('tic_teff.fits');
;  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
;  for ii=0, numfil-1 do begin
;    fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, dmax=dmax, icmax=icmax, homies=homies, dbl=1
;    numstar[ii] = nstar
;  end
;  print, numfil, ' files contain ', total(numstar), ' stars within ', 10.^(dmax/5.+1.), ' pc.'
;  print, numfil, ' files contain ', total(numstar), ' stars brighter than Ic=', icmax
  allstar = replicate({starstruct}, 1E8)
  idx0 = 0L
  for ii=0, numfil-1 do begin
      restore, fnames[ii]
      if (keyword_set(dmax)) then gd = where(star.mag.dm le dmax)
      if (keyword_set(icmax)) then gd = where(star.mag.ic le icmax)
      if (keyword_set(kmax)) then gd = where(star.mag.k le kmax)
      if (keyword_set(homies)) then gd = where(star.spl)
      if (gd[0] ne -1) then begin
        numstar = n_elements(gd)
        idx = idx0+lindgen(numstar)
        allstar[idx] = star[gd]
        idx0 = idx0+numstar
      end    
      thisfn = repstr(fnames[ii], 'hp', 'bk')
      restore, thisfn
      if (keyword_set(dmax)) then gd = where(star.mag.dm le dmax)
      if (keyword_set(icmax)) then gd = where(star.mag.ic le icmax)
      if (keyword_set(kmax)) then gd = where(star.mag.k le kmax)
      if (keyword_set(homies)) then gd = where(star.spl)
      if (gd[0] ne -1) then begin
        numstar = n_elements(gd)
        for jj=0,dpmult-1 do begin
          idx = idx0+lindgen(numstar)
          allstar[idx] = star[gd]
          idx0 = idx0+numstar
        end
      end   
     print, fnames[ii] 
  end
  allstar = allstar[0:(idx0-1),*]
  save, allstar, filen=fname
END

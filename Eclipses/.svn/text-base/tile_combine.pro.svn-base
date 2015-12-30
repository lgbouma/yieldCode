PRO tile_combine, fstub
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = intarr(numfil)
  for ii=0, numfil-1 do begin
    h = headfits(fnames[ii])
    numstar[ii] = fxpar(h, 'NAXIS1')
  end
  print, numfil, ' files contain ', total(numstar), ' stars.'
  numcol = fxpar(h, 'NAXIS2')
  newtbl = fltarr(total(numstar), numcol)
  idx0 = 0L
  for ii=0, numfil-1 do begin
    dat = mrdfits(fnames[ii])
    idx = idx0+indgen(numstar[ii])
    newtbl[idx,*] = dat
    idx0 = idx0+numstar[ii]
  end
  newname = repstr(fstub, '*', '_all')
  mwrfits, newtbl, newname
END

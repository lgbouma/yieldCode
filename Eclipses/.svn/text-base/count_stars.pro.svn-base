PRO count_stars, fstub
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = intarr(numfil)
  for ii=0, numfil-1 do begin
    h = headfits(fnames[ii])
    numstar[ii] = fxpar(h, 'NAXIS1')
  end
  print, numfil, ' files contain ', total(numstar), ' stars.'
END

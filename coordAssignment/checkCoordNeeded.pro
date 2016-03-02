pro checkCoordNeeded
  ; Find if 2e5 unique coordinates per tile is not enough

  fNums = mrdfits('../Eclipses/fnums.fits')
  nTile = N_ELEMENTS(fNums)
  triPath = '../../trilegal/'
  nElArr = []
  for i=0, nTile-1 do begin
    fname = triPath + 'hp' + string(fNums[i], format='(I04)')+'.sav'
    print, fname
    restore, fname
    nElArr = [nElArr, n_elements(star)]
    ;delvarx, star
    print, max(nElArr)
  endfor 

end

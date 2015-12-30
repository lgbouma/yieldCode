PRO coord_ps, fpre, fnums
  ; Assign physical coordinates and do postage stamp selection
  ; Gather the .sav files
  numfil = n_elements(fnums)
  numstar = intarr(numfil)
  u = randomu(seed, 4D8)
  v = randomu(seed, 4D8)
  phi = 2.*!dpi*u
  theta = acos(2.*v-1.)
  ang2pix_ring, 16, theta, phi, ipring
  ;cgHistoPlot, ipring
  for ii=0, numfil-1 do begin
    fname = fpre+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    nstar = n_elements(star[where(star.sec ne 1)])
    thispix = where(ipring eq fnums[ii])
    ncoord = n_elements(thispix)
    if (nstar gt ncoord) then begin
      print, "Need ", nstar, " coords but got ", ncoord, " on tile ", ii
    endif else begin
      print, "Filling in tile ", ii
      glon = phi[thispix[0:(nstar-1)]]*180./!dpi
      glat = (theta[thispix[0:(nstar-1)]]-!dpi)*180./!dpi
      euler, glon, glat, elon, elat, select=5
      star[where(star.sec ne 1)].coord.elon = elon
      star[where(star.sec ne 1)].coord.elat = elat
      star[where(star.sec ne 1)].coord.elon = glon
      star[where(star.sec ne 1)].coord.elat = glat
      star[where(star.sec ne 1)].coord.healpix_n = fnums[ii]
    end
  end
END

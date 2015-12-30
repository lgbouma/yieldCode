PRO test_planets
  frac = mrdfits('psfs/dfrac24_105_f3p33.fits')
  ph_p = 'ph_T_filt.fits'
  stars = replicate({starstruct}, 2D6)
  stars[0:(1D6-1)].teff = 6000.
  stars[1D6:(2D6-1)].teff = 3000.
  stars[0:(1D6-1)].m = 1.0
  stars[1D6:(2D6-1)].m = 0.1
  stars[0:(1D6-1)].r = 1.0
  stars[1D6:(2D6-1)].r = 0.1
  n = add_planets(stars, eclip, frac, ph_p, dressing=1)
  rp = eclip.r2
  per = eclip.p
  teff = stars[eclip.hostid].teff
  mwrfits, [[rp],[per],[teff]], 'test_planets.fits'
END


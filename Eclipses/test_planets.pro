PRO test_planets
  frac = mrdfits('psfs/dfrac_t75_f3p31_3.fits')
  ph_p = mrdfits('ph_T_filt.fits')
  readcol, 'tband.csv', lam, t
  tband = [[lam],[t]]
  stars = replicate({starstruct}, 2D5)
  stars[0:(1D5-1)].teff = 6000.
  stars[1D5:(2D5-1)].teff = 3000.
  stars[0:(1D5-1)].m = 1.0
  stars[1D5:(2D5-1)].m = 0.1
  stars[0:(1D5-1)].r = 1.0
  stars[1D5:(2D5-1)].r = 0.1
  n = add_planets(stars, eclip, frac, ph_p, tband, dressing=1)
  rp = eclip.r2
  per = eclip.p
  pr = eclip.pr
  mult = eclip.mult
  teff = stars[eclip.hostid].teff
  hid = eclip.hostid
  mwrfits, [[rp],[per],[teff],[mult],[pr],[hid]], 'test_add_planets.fits'
END
sfasdf

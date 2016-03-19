PRO dilute_eclipse_img, eclip, bkgnds, frac, ph_p, dx, dy, dilvec, $
	aspix=aspix, sq_deg=sq_deg, radmax=radmax
  ; 16/03/19 (LB): what's done here is correct and smart. It's a great way to avoid huge catalogs.
  ; Offset 0,0 maps to pixel 7.9,7.9
  ; Offset 9,9 maps to pixel 7,7 (exactl)
  ; So increasing dx, dy -> decreasing pixel value

  ; How many filters in the prf?
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  imgsize = sqrt(sz_frac[4])
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(sq_deg)) then sq_deg=sq_deg else sq_deg=1.0
  if (keyword_set(radmax)) then radmax=radmax else radmax=8
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels
  radas = sqrt(sq_deg/!dpi)*3600.
  radpix = radas/aspix
  print, 'Considering stars out to ', radas, ' arcsec or ', radpix, ' pixels'
  nbk = n_elements(bkgnds)
  
  neclip = n_elements(eclip)
  fov = eclip.coord.fov_ind
  dilvec = dblarr(imgsize*imgsize/4, neclip)

  for ii=0, neclip-1 do begin
    ; Generate some random radii
    randomp, r, 1., nbk, range_x=[0., radpix]
    ; How many of these are primaries and fall within radmax?
    gd = where((r lt radmax) and (bkgnds.sec ne 1)) ; note these bkgnds are deeps on first call
    if ((ii mod 1000.) eq 0) then print, "On eclipse ", ii, " of ", neclip, " with ", n_elements(gd), " stars"
    if (gd[0] ne -1) then begin
     ; r[gd] = 0.0 ; test purposes only
      bkteff = bkgnds[gd].teff
      bkmagt = bkgnds[gd].mag.tsys
      bkrad  = r[gd]
      bktheta = 2.*!dpi*randomu(seed, n_elements(gd))
      bkx = bkrad*cos(bktheta)
      bky = bkrad*sin(bktheta)

      ; Compute the fluxes in each sub-filter for each star
      recipteff = 4000./bkteff
      ph_filt = dblarr(nfilt, n_elements(gd))
      for jj=0, nfilt-1 do begin
        ph_filt[jj,*] = ph_p[jj,0] + ph_p[jj,1]*recipteff + $
           ph_p[jj,2]*recipteff^2. + ph_p[jj,3]*recipteff^3.
      end
      ph_filt[where(ph_filt lt 0.0)] = 0.0
      
      ; Stack up the starz
      dilpix = dblarr(imgsize/2, imgsize/2)
      for kk=0, n_elements(gd)-1 do begin
        dydbl = floor(dy[ii]/10.)
        dyrem = dy[ii]-10*dydbl
        pixdx = floor(bkx[kk] + dx[ii]/10. + 0.1)
        pixdy = floor(bky[kk] +  dyrem/10. + 0.1)
        indx  = round(10*(bkx[kk] - pixdx + dx[ii]/10.))
        indy  = round(10*(bky[kk] - pixdy +  dyrem/10.)) + 10*dydbl
        xsel = indgen(imgsize/2) + imgsize/4 + pixdx
        ysel = indgen(imgsize/2) + imgsize/4 + pixdy 
        thisprf = reform(frac[indx,indy,fov[ii],*,0:(nfilt-1)])
        thisimg = reform(thisprf#ph_filt[*,kk], imgsize, imgsize)
        thisimgx  = thisimg[xsel,*]
        thisimgxy = thisimgx[*,ysel]
        dilpix[*,*] = dilpix[*,*] + $
		10^(-0.4*(bkmagt[kk]-10.))*thisimgxy
      end ; loop over stars
      ;stop
      dilvec[*,ii] = reform(dilpix, imgsize*imgsize/4) 
    end ; if gd
  end ; loop over eclipses
END

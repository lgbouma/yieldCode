PRO starcat, fstub, fname, dmax=dmax, icmax=icmax, kmax=kmax, rmax=rmax, homies=homies, ps=ps
  fnames = file_search(fstub)
;  stop
;  restore, 'dartmouth_grid.sav'
;  tt = mrdfits('tic_teff.fits');
;  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
;  numstar = lonarr(numfil)
;  for ii=0, numfil-1 do begin
;    fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, dmax=dmax, icmax=icmax, homies=homies, dbl=1
;    numstar[ii] = nstar
;  end
;  print, numfil, ' files contain ', total(numstar), ' stars within ', 10.^(dmax/5.+1.), ' pc.'
;  print, numfil, ' files contain ', total(numstar), ' stars brighter than Ic=', icmax
;  nustar = replicate({starstruct}, 1E8)i
  ph_fits = mrdfits('ph_T_filt.fits') ; photon fluxes for T=10 vs Teff
  npnt_fits = mrdfits('npnt.fits')

  idx0 = 0L
  nparam=21
  nustar = fltarr(1e7, nparam)
  for ii=0, numfil-1 do begin
    ;if (numstar[ii] gt 0) then begin
      ;thisfn = repstr(fnames[ii], '.fits', '.sav')
      restore, fnames[ii]
      gd = []
      if (keyword_set(icmax)) then gd = where(star.mag.ic le icmax) $
      else if (keyword_set(kmax)) then gd = where(star.mag.k le kmax) $
      else if (keyword_set(rmax)) then gd = where(star.r le rmax and star.mag.ic le icmax) $
      else if (keyword_set(dmax)) then gd = where(star.mag.dm le dmax) $
      else if (keyword_set(homies)) then gd = where(star.spl) $
      else if (keyword_set(ps)) then begin
        pri = where(star.pri lt 2)
	selpri = ps_sel(star[pri].mag.t, star[pri].teff, star[pri].m, star[pri].r, ph_fits, $
                        rn_pix=15., npnt=npnt_fits[ii])
        gd = pri[selpri]
      endif
      ;else gd = where(star.mag.ic gt 0.0)
      if (gd[0] ne -1) then begin
        m = star[gd].m
        rad = star[gd].r
        teff = star[gd].teff
        v = star[gd].mag.v
        r = star[gd].mag.r
        ic = star[gd].mag.ic
        z = star[gd].mag.z
        j = star[gd].mag.j
        h = star[gd].mag.h
        k = star[gd].mag.k
        dm = star[gd].mag.dm
        av = star[gd].mag.dm
        ps = star[gd].pri + 2*star[gd].sec
        mv = star[gd].mag.mv
        mic = star[gd].mag.mic
        mj = star[gd].mag.mj
        icsys = star[gd].mag.icsys
        jsys = star[gd].mag.jsys
        mvsys = star[gd].mag.mvsys
        micsys = star[gd].mag.micsys
        mjsys = star[gd].mag.mjsys

        poop = [[m],[rad],[teff],[v],[r],[ic],[z],[j],[h],[k],$
	[dm],[av],[ps],[mv],[mic],[mj],[icsys],[jsys],[mvsys],[micsys],[mjsys]]

        idx = idx0+lindgen(n_elements(gd))
        nustar[idx,*] = poop
        idx0 = idx0+n_elements(gd)
        print, fnames[ii], ' on index ', idx0
      end
  end
  nustar = nustar[lindgen(idx0-1),*]
  mwrfits, nustar, fname
END

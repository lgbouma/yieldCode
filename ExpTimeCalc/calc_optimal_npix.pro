pro calc_optimal_npix,ps=ps

  !p.charsize=2

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='optimal_npix.eps',/encapsulated,xsize=6.5,ysize=9,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=2
     syms=0.5

  endif

  !p.multi=[0,1,3]
  
  n = 1000
  r = 5.0*dindgen(n)/double(n) ; pixels
  npix = !PI*r^2.

; following is from Roland, a rough description of the PSF
; encircled energy fraction
;
  fwhm = 1.2
  frac = 1.0 - (fwhm/(2.*!PI)) / ((fwhm/2.)^2 + r^2.)
  plot, r, frac, $
        xtit='aperture radius [pixels]', xra=[0.5,5], xsty=1, $
        ytit='fraction of flux enclosed', yra=[0.5,1.1], ysty=1

; now step through apparent magnitude and optimize aperture size

  n_mags = 100
  imag_min = 6.0
  imag_max = 15.0

  imag = imag_min + (imag_max-imag_min)*dindgen(n_mags)/double(n_mags-1)
  r_opt = 0.*imag

  for i=0,n_mags-1 do begin

     snr = 0.*r

     for j=0,n-1 do begin
        calc_noise, imag[i], 3600., noise, $
                    npix_aper = npix[j], frac_aper = frac[j]
        snr[j] = 1./noise
     endfor
  
     snr_opt = max(snr, j_opt)
     r_opt[i] = r(j_opt)

  endfor

  plot, imag, r_opt, $
        xtit = 'I mag', xra=[6,15], xsty=1, $
        ytit = 'optimal aperture radius'

  npix = !PI*r_opt^2.

  plot, imag, npix, /ylog, $
        xtit = 'I mag', xra=[6,15], xsty=1, $
        ytit = 'optimal number of pixels', yra = [0.5, 50], ysty=1

  fitme = where(npix gt 2.0)
  x = imag[fitme] & y = alog10(npix[fitme])
  result = poly_fit(x,y,1.0,yfit=yfit)
  print, 'Fit to imag: ', result[0], result[1]
  print, 'Fit to imag-10: ', 10.*result[1]+result[0], result[1]
  oplot, x, 10.^yfit, color=fsc_color('Red')

  npix_opt = optimal_npix(imag)
  oplot, imag, npix_opt, color=fsc_color('Cadet Blue')
  
  x0=6.5
  y0=2
  r=0.5
  xyouts, x0,y0,'polynomial fit',color=fsc_color('Red'), charsize=0.75*!P.charsize
  xyouts, x0,y0*r,'model',color=fsc_color('Cadet Blue'), charsize=0.75*!P.charsize

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end


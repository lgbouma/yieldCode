pro plot_noise,ps=ps

  !p.charsize=2

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='noise.eps',/encapsulated,xsize=6.5,ysize=5.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  !p.multi=0
  plotsym,0,/fill

  c = fsc_color(['Red','Orange','Forest Green','Navy','Violet'])

  n=100
  imag = 3.0 + 17.0*dindgen(n)/double(n-1)
  noise = 0.*imag
  noise_star = 0.*imag
  noise_sky = 0.*imag
  noise_ro = 0.*imag
  noise_sys = 0.*imag

  exptime = 3600.0
  sys_limit = 60.0

  calc_noise, imag, exptime, noise, $
                 sys_limit=sys_limit, $
                 noise_star=noise_star, $
                 noise_sky=noise_sky, $
                 noise_ro=noise_ro, $
                 noise_sys=noise_sys

  plot, imag, noise, /nodata, $
        xtit='Apparent I mag', xra=[3,20], xsty=1,$
        ytit='Noise in one-hour exposure', yra=[1E-6,1E-1], /ylog, ysty=1
  
  oplot, imag, noise_star, color=c[0]
  oplot, imag, noise_sky, color=c[1]
  oplot, imag, noise_ro, color=c[2]
  oplot, imag, noise_sys, color=c[3]
  oplot, imag, noise, color=c[4]
  
; now repeat but without optimal apertures

  calc_noise, imag, exptime, noise, $
              sys_limit=sys_limit, $
              npix_aper = 16.0, $
              noise_star=noise_star, $
              noise_sky=noise_sky, $
              noise_ro=noise_ro, $
              noise_sys=noise_sys

  oplot, imag, noise_star, color=c[0], linestyle=1
  oplot, imag, noise_sky, color=c[1], linestyle=1
  oplot, imag, noise_ro, color=c[2], linestyle=1
  oplot, imag, noise_sys, color=c[3], linestyle=1
  oplot, imag, noise, color=c[4], linestyle=1

  x0 =4
  y0 =5D-2
  r = 0.5
  xyouts, x0, y0*[1,r,r^2,r^3,r^4], ['Star','Sky','Readout','Sys','Total'], color=c

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

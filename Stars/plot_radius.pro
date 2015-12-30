pro plot_radius,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='radius.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica,decomposed=1
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  restore,filen='../dartmouth_models.sav'

  !p.multi=0
  
  q = where(age ge 1.0 and age le 8.0)

  plot, rs[q], jcvmag[q], /nodata, $
        xtit='stellar radius [R!Dsun!N]', xra=[0.1,1.0], xsty=1, $
        ytit='abs V or J mag', yra=[11,3], ysty=1
  oplot, rs[q], jcvmag[q], psym=3, color=fsc_color('Cadet Blue')
  oplot, rs[q], jmag[q], psym=3, color=fsc_color('Medium Gray')
  
  q = where(met ge -0.2 and met le 0.2 and age ge 1.0 and age le 8.0 and $
           rs ge 0.1 and rs le 1.0)
  oplot, rs[q], jmag[q], psym=3, color=fsc_color('Dark Gray')

  x = rs[q] & y = jmag[q]
  poly_jmag_given_radius = poly_fit(x, y, 5, yfit=yfit)
  x_calc = 0.1 + 0.9*dindgen(100)/99.
  y_calc = poly(x_calc, poly_jmag_given_radius)
  oplot, x_calc, y_calc, color=fsc_color('Red')

  x0 = 0.15
  y0 = 4.0
  dy = 0.6
  xyouts, x0, y0, alignment=0, 'Dartmouth V band', color=fsc_color('Cadet Blue')
  xyouts, x0, y0+dy, alignment=0, 'Dartmouth J band', color=fsc_color('Medium Gray')
  xyouts, x0, y0+1.7*dy, alignment=0, $
          '   (-0.2 < [Fe/H] < 0.2)', color=fsc_color('Dark Gray'),$
          charsize=0.75*!p.charsize
  xyouts, x0, y0+3.*dy, alignment=0, 'polynomial fit', color=fsc_color('Red')
  
  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

  save, filen='poly_jr.sav', poly_jmag_given_radius

end

pro plot_mass,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='mass.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica,decomposed=1
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  restore,filen='../dartmouth_models.sav'

  !p.multi=0
  
  q = where(age ge 1.0 and age le 8.0)

  plot, mass[q], jcvmag[q], /nodata, $
        xtit='stellar mass [M!Dsun!N]', xra=[0.1,1.0], xsty=1, $
        ytit='abs V or J mag', yra=[11,3], ysty=1
  oplot, mass[q], jcvmag[q], psym=3, color=fsc_color('Cadet Blue')
  oplot, mass[q], jmag[q], psym=3, color=fsc_color('Medium Gray')
  
  q = where(met ge -0.2 and met le 0.2 and age ge 1.0 and age le 8.0)
  oplot, mass[q], jmag[q], psym=3, color=fsc_color('Dark Gray')

  mj_delfosse = 6.0 + 4.0*dindgen(100)/99.
  m_delfosse = delfosse_mass_given_mj(mj_delfosse)
  oplot, m_delfosse, mj_delfosse, color=fsc_color('Red')

  x = mass[q] & y = jmag[q]
  poly_jmag_given_mass = poly_fit(x,y,5.0)
  x_calc = 0.1 + 0.9*dindgen(100)/99.
  y_calc = poly(x_calc,poly_jmag_given_mass)
;  oplot, x_calc, y_calc, color=fsc_color('Red')

  x0 = 0.15
  y0 = 4.0
  dy = 0.5
  xyouts, x0, y0, alignment=0,       'Dartmouth V band', color=fsc_color('Cadet Blue')
  xyouts, x0, y0+dy, alignment=0,    'Dartmouth J band', color=fsc_color('Medium Gray')
  xyouts, x0, y0+1.7*dy, alignment=0, '   (-0.2 < [Fe/H] < 0.2)', color=fsc_color('Dark Gray'),$
          charsize=0.75*!p.charsize
  xyouts, x0, y0+3.*dy,  alignment=0,'Delfosse et al. (2000)', color=fsc_color('Red')
;  xyouts, x0, y0+3.*dy, alignment=1, 'polynomial fit', color=fsc_color('Red')
  
  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

pro plot_dartmouth_teff,ps=ps

  !p.charsize=2

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='dartmouth_teff.eps',/encapsulated,xsize=6.5,ysize=5.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  common DATA, m, r

  restore,filen='../dartmouth_models.sav'

  !p.multi=0
  plotsym,0,/fill
  
  q = where(age ge 1.0 and age le 8.0)

  plot, rs[q], teff[q], /nodata, $
        xtit='radius [M!Dsun!N]', xsty=1, xra=[0.05, 1.1], $
        ytit='T!Deff!N [K]', ysty=1, yra=[3000,6200]

  oplot, rs[q], teff[q], psym=3, color=fsc_color('Gray')

  q = where(rs ge 0. and rs le 1.2 and $
            met ge -0.2 and met le 0.2 and age ge 1.0 and age le 8.0)
  r=rs[q] & t = teff[q]

  oplot, r, t, psym=3; , color=fsc_color('Dark Gray')

  poly_teff_given_radius = poly_fit(r, t, 9, yfit=t_calc)
  q = sort(r)
  oplot, r[q], t_calc[q], color=fsc_color('Red')
  
  x0 = 0.1
  y0 = 6000
  dy = 120
  xyouts, x0, y0, 'Dartmouth models', color=fsc_color('Gray')
  xyouts, x0, y0-dy, ' (age 1-8 Gyr, -0.2 < [Fe/H]< 0.2)', charsize=0.75*!p.charsize
  xyouts, x0, y0-2.2*dy, 'polynomial fit', color=fsc_color('Red')

  save, filen='poly_teff_given_radius.sav', poly_teff_given_radius

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

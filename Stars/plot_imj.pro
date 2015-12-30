pro plot_imj,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='imj.eps',/encapsulated,xsize=6.5,ysize=6.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  readcol, 'imj.txt', f='D,D', comment='#', /silent, $
           j, imj

  q = where(j le 11.5)
  j=j[q] & imj=imj[q]
  i = j+imj

  restore, filen='../dartmouth_models.sav'
  q = where(met ge -0.2 and met le 0.2 and age le 8.0)
  imj_d = jcimag[q]-jmag[q]
  j_d = jmag[q]
  i_d = jcimag[q]
  
  !p.multi=[0,1,2]

  plot, imj, j, psym=8, yra=[12,2], ysty=1, xtit='I-J', ytit='abs J mag'
  oplot, imj_d, j_d, color=fsc_color('Gray'), psym=3
  oplot, imj, j, psym=8, color=fsc_color('Blue')
  poly_imj_given_j = poly_fit(j, imj, 5, yfit=imj_fit)
;  oplot, imj_fit, j, color=fsc_color('Red')

  j_calc = 13.0*dindgen(100)/99.
  imj_calc = poly(j_calc,poly_imj_given_j)
  oplot, imj_calc, j_calc, color=fsc_color('Red'), linestyle=1
  ok = where(j_calc ge 3.0 and j_calc le 11)
  oplot, imj_calc[ok], j_calc[ok], color=fsc_color('Red')
  
  xyouts, 2.9, 3.0, 'I.N. Reid web site', alignment=1, color=fsc_color('Blue')
  xyouts, 2.9, 4.0, 'Polynomial fit', alignment=1, color=fsc_color('Red')
  xyouts, 2.9, 5.0, 'Dartmouth models', alignment=1, color=fsc_color('Gray')

  print, 'Residuals = ', imj-imj_fit
  print, 'stddev Residual = ', stddev(imj-imj_fit)

  plot, imj, i, psym=8, yra=[14.5,2.5], ysty=1, $
        xtit='I-J', ytit='abs I mag'
  oplot, imj_d, i_d, color=fsc_color('Gray'), psym=3
  oplot, imj, i, psym=8, color=fsc_color('Blue')
  poly_imj_given_i = poly_fit(i, imj, 5.0, yfit=imj_fit)
;  oplot, imj_fit, i, color=fsc_color('Red')

  i_calc = 15.0*dindgen(100)/99.
  imj_calc = poly(i_calc,poly_imj_given_i)
  oplot, imj_calc, i_calc, color=fsc_color('Red'), linestyle=1
  ok = where(i_calc ge 3.5 and i_calc le 13.75)
  oplot, imj_calc[ok], i_calc[ok], color=fsc_color('Red')

  xyouts, 2.9, 4.0, 'I.N. Reid web site', alignment=1, color=fsc_color('Blue')
  xyouts, 2.9, 5.2, 'Polynomial fit', alignment=1, color=fsc_color('Red')
  xyouts, 2.9, 6.4, 'Dartmouth models', alignment=1, color=fsc_color('Gray')

  print, 'Residuals = ', imj-imj_fit
  print, 'stddev Residual = ', stddev(imj-imj_fit)

  save, filen='poly_ij.sav', poly_imj_given_i, poly_imj_given_j

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

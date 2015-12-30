pro plot_vmi,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='vmi.eps',/encapsulated,xsize=6.5,ysize=7.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  readcol, 'rh_vlf.txt', f='D,D,D', comment='#', /silent, $
           v, phi, vmi

  q = where(vmi gt 0 and v le 16.0)
  v=v[q] & vmi=vmi[q]
  i = v - vmi

  restore, filen='../dartmouth_models.sav'
  q = where(met ge -0.2 and met le 0.2 and age le 8)
  vmi_d = jcvmag[q]-jcimag[q]
  v_d = jcvmag[q]
  i_d = jcimag[q]
  
  !p.multi=[0,1,2]

  plot, vmi, v, psym=8, yra=[16,3], ysty=1, xtit='V-I', ytit='abs V mag'
  oplot, vmi_d, v_d, color=fsc_color('Gray'), psym=3
  oplot, vmi, v, psym=8, color=fsc_color('Blue')

  poly_vmi_given_v = poly_fit(v, vmi, 6.0, yfit=vmi_fit)

  y = 3.0 + 13.0*dindgen(100)/99.
  x = poly(y,poly_vmi_given_v)
  oplot, x, y, color=fsc_color('Red')

  xyouts, 3.8, 4.5, 'I.N. Reid web site', alignment=1, color=fsc_color('Blue')
  xyouts, 3.8, 5.5, 'Polynomial fit', alignment=1, color=fsc_color('Red')
  xyouts, 3.8, 6.5, 'Dartmouth models', alignment=1, color=fsc_color('Gray')

  print, 'Residuals = ', vmi-vmi_fit
  print, 'stddev Residual = ', stddev(vmi-vmi_fit)

;;;

  plot, vmi, i, psym=8, yra=[13,2], ysty=1, xtit='V-I', ytit='abs I mag'
  oplot, vmi_d, i_d, color=fsc_color('Gray'), psym=3
  oplot, vmi, i, psym=8, color=fsc_color('Blue')

  poly_vmi_given_i = poly_fit(i, vmi, 6.0, yfit=vmi_fit)

  y = 13.0*dindgen(100)/99.
  x = poly(y,poly_vmi_given_i)
  oplot, x, y, color=fsc_color('Red'), linestyle=1
  q = where(y ge 3 and y le 12.)
  oplot, x[q], y[q], color=fsc_color('Red')

  xyouts, 3.8, 3.0, 'I.N. Reid web site', alignment=1, color=fsc_color('Blue')
  xyouts, 3.8, 4.0, 'Polynomial fit', alignment=1, color=fsc_color('Red')
  xyouts, 3.8, 5.0, 'Dartmouth models', alignment=1, color=fsc_color('Gray')

  print, 'Residuals = ', vmi-vmi_fit
  print, 'stddev Residual = ', stddev(vmi-vmi_fit)

  save, filen='poly_vi.sav', poly_vmi_given_v, poly_vmi_given_i

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

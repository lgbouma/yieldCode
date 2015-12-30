pro plot_dartmouth

  restore,filen='../dartmouth_models.sav'

  !p.charsize=2
  !p.multi=[0,1,2]
  plotsym,0,/fill
  
  plot, rs, jcvmag, psym=3, $
        xra=[0,1.2], xsty=1, $
        yra=[15,3], ysty=1
  oplot, rs, jmag, psym=3, color=fsc_color('Red')
  
  q = where(met ge -0.6 and met le 0.2 and rs le 1.0)
  oplot, rs[q], jcvmag[q], psym=3, color=fsc_color('Green')
  oplot, rs[q], jmag[q], psym=3, color=fsc_color('Cyan')

  x = rs[q] & y = jmag[q]
  result = poly_fit(x,y,5.0)
  print, transpose(result)
  x_calc = 0.1 + 0.9*dindgen(100)/99.
  y_calc = poly(x_calc,result)
  oplot, x_calc, y_calc, color=fsc_color('Yellow')

;;;;

  plot, mass, jcvmag, psym=3, $
        xra=[0,1.2], xsty=1, $
        yra=[15,3], ysty=1
  oplot, rs, jmag, psym=3, color=fsc_color('Red')
  
  q = where(met ge -0.2 and met le 0.2 and rs le 1.0)
  oplot, mass[q], jcvmag[q], psym=3, color=fsc_color('Green')
  oplot, mass[q], jmag[q], psym=3, color=fsc_color('Cyan')

  x = mass[q] & y = jmag[q]
  result = poly_fit(x,y,5.0)
  print, transpose(result)
  x_calc = 0.1 + 0.9*dindgen(100)/99.
  y_calc = poly(x_calc,result)
  oplot, x_calc, y_calc, color=fsc_color('Yellow')

  mj_delfosse = 6.0 + 4.0*dindgen(100)/99.
  m_delfosse = delfosse_mass_given_mj(mj_delfosse)
  oplot, m_delfosse, mj_delfosse, color=fsc_color('Orange')

end

pro plot_dartmouth_mr,ps=ps

  !p.charsize=2

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='dartmouth_mr.eps',/encapsulated,xsize=6.5,ysize=5.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  common DATA, m, r

  restore,filen='../dartmouth_models.sav'

  !p.multi=0
  plotsym,0,/fill
  
  plot, mass, rs, /nodata, $
        xtit='mass [M!Dsun!N]', xsty=1, /xlog, xra=[0.1, 1.2], $
        ytit='radius [R!Dsun!N]', ysty=1, /ylog, yra=[0.1, 1.2]

  oplot, mass, rs, psym=3, color=fsc_color('Gray')

  q = where(mass ge 0.05 and mass le 1.1 and $
            rs ge 0.05 and rs le 1.1 and $
            met ge -0.2 and met le 0.2 and age ge 1.0 and age le 8.0)
  m = mass[q] & r=rs[q]

  oplot, m, r, psym=3; , color=fsc_color('Dark Gray')

;  p_guess = [0.37, 0.35, 0.81, 0.81, 0.97, 1.27]
;  p_scale = [0.1,   0.1,   0.1,  0.1, 0.1, 0.1]
;  p = amoeba(1E-8, function_name='chisq', p0=p_guess, scale=p_scale)
;  print, p, chisq(p)
  
  m_calc = 0.01 + 1.2*dindgen(100)/99.
  r_calc = radius_given_mass(m_calc)
;  oplot, m_calc, r_calc, color=fsc_color('Red')

  r_calc = 0.01 + 1.2*dindgen(100)/99.
  m_calc = mass_given_radius(r_calc)
  oplot, m_calc, r_calc, color=fsc_color('Red')
  
  x0 = 0.12
  y0 = 0.9
  r = 0.9
  xyouts, x0, y0, 'Dartmouth models', color=fsc_color('Gray')
  xyouts, x0, y0*r, '   (age 1-8 Gyr, -0.2 < [Fe/H]< 0.2)', charsize=0.75*!p.charsize
  xyouts, x0, y0*r^3, 'broken power law', color=fsc_color('Red')


  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

function chisq, p

  common DATA, m, r

  r_calc = radius_given_mass(m, p=p)
  ans = total( (r-r_calc)^2.0 )
  
  return, ans
  
end


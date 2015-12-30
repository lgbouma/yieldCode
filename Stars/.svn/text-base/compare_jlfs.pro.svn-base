pro compare_jlfs,ps=ps

  !p.multi=0
  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='compare_jlfs.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  common DATA, x, y

  readcol, 'jlfs.txt', f='D,D,D,D,D,D,D,D,D', comment='#', /silent, $
           mj, rh, rca07, cruz8_sys, cruz8_tot, cruz20_tot, covey, boch_sys, boch_tot

  vol = 4./3.*!PI*(10.0^3.)

  rh *= vol*2.0               ; the 2.0 is because the tabulations are in terms of 0.5-mag intervals
  rca07 *= vol*2.0
  cruz8_sys *= vol*2.0
  cruz8_tot *= vol*2.0
  cruz20_tot *= vol*2.0
  covey *= vol*2.0
  boch_sys *= vol*2.0
  boch_tot *= vol*2.0

  c = fsc_color(['Forest Green','Cadet Blue','Navy','Blue Violet','Red','Dark Gray'])

  plot, mj, rh, psym=10, /nodata, $
        xtit='abs J mag', xra=[-1,15], xsty=1, $
        ytit='(stars mag!E-1!N) within 10 pc'

  oplot, mj, cruz8_tot, psym=10, color=c[0]
  oplot, mj, cruz8_sys, psym=10, color=c[0], linestyle=1
  oplot, mj, cruz20_tot, psym=10, color=c[1]
  oplot, mj, covey, psym=10, color=c[2]
  oplot, mj, boch_tot, psym=10, color=c[3]
  oplot, mj, boch_sys, psym=10, color=c[3], linestyle=1
  oplot, mj, rh, psym=10, thick=8, color=c[5]
  
  x0 = -0.25
  y0 = 90 & dy = 7
  xyouts, x0, y0, 'Reid & Hawley (2005)', color=c[5]
  xyouts, x0, y0-dy, 'Cruz et al. (2007) 8 pc', color=c[0]
  xyouts, x0, y0-2.*dy, 'Cruz et al. (2007) 20 pc', color=c[1]
  xyouts, x0, y0-3.*dy, 'Covey et al. (2008)', color=c[2]
  xyouts, x0, y0-4.*dy, 'Bochanski et al. (2010)', color=c[3]
  xyouts, x0, y0-5.*dy, 'Analytic model', color=c[4]
  
  q = where(mj ge 3.0 and mj le 11.5)

  x = mj[q]
  y = rh[q]

;  p_guess = [84, 7.25, 3.2, 1.8, 1.5, 0.6]
;  p_scale = [5,  0.25, 0.5, 0.4, 0.3, 0.1]

  p_guess = [80, 5.0, 0.6, 3.0, 0.4]
  p_scale = [5,  0.5, 0.2, 0.5, 0.4]

  p = amoeba(1E-9, function_name='chisq', p0=p_guess, scale=p_scale)

  print, p
  
  x_calc = 15.0*dindgen(1000)/999.
  y_calc = jlf_model(x_calc, p=p)
  oplot, x_calc, y_calc, color=c[4], linestyle=2
  q = where(x_calc ge 3.0 and x_calc le 11.5)
  oplot, x_calc[q], y_calc[q], color=c[4]

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

function chisq, p

  common DATA, x, y

  jlf_calc = jlf_model(x, p=p)
  ans = total( (jlf_calc-y)^2. )
  
  return, ans
  
end

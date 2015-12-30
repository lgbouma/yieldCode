pro compare_zheng,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='compare_zheng.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica,decomposed=1
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif


  readcol, 'zheng_ilf.txt', f='D,D,D', comment='#', /silent, $
           mi, phi, dphi

  phi = phi/1d3
  dphi = dphi/1d3
  vol = 4./3.*!PI*(10.0^3.)
  n = phi*vol
  dn = dphi*vol

  !p.multi=0
  plotsym,0,/fill
  ploterror, mi, n, dn, /nohat, psym=8, $
             xtit='abs I mag', xra=[3,15], xsty=1, $
             ytit='(stars mag!E-1!N) within 10 pc', yra=[0,110], ysty=1

  restore, filen='poly_ij.sav'

  mi_mod = 4.0 + 10.0*dindgen(30)/29.
  n_mod = 0.0*mi_mod

  for i=0,n_elements(mi_mod)-1 do begin
     
     mi1 = mi_mod[i]-0.5
     mi2 = mi_mod[i]+0.5
     imj1 = poly(mi1,poly_imj_given_i)
     imj2 = poly(mi2,poly_imj_given_i)
     mj1 = mi1 - imj1
     mj2 = mi2 - imj2

     n_mod[i] = qromb('jlf_model', mj1, mj2, k=6)
;     print, mi1, mi2, mj1, mj2, n_mod[i]

  endfor

  oplot, mi_mod, n_mod, color=fsc_color('Red')

  xyouts, 3.5, 100, alignment=0, 'Zheng et al. (2004)'
  xyouts, 3.5, 92, alignment=0, 'Model', color=fsc_color('Red')

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

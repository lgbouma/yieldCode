pro compare_vband,ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='compare_vband.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica,decomposed=1
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  readcol, 'rh_vlf.txt', f='D,D,D', comment='#', mv_rh, phi_rh, junk
  readcol, 'recons_vlf.txt', f='D,D', comment='#', mv_recons, n_recons
  readcol, 'zheng_vlf.txt', f='D,D,D', mv_z, phi_z, dphi_z, /silent

  vol = 4./3.*!PI*(10.0^3.)
  n_rh = phi_rh*vol
  n_z = phi_z*1d-3 * vol

  !p.multi=0

  plot, mv_rh, n_rh, psym=10, $
        xtit='abs V mag', xra=[0,21], xsty=1, $
        ytit='(stars mag!E-1!N) within 10 pc', yra=[0,70], ysty=1
  oplot, mv_recons, n_recons, psym=10, color=fsc_color('Cadet Blue')
  oplot, mv_z, n_z, psym=10, color=fsc_color('Orange')

  n_mod = 0.0*mv_rh
  restore,filen='poly_vi.sav'
  restore,filen='poly_ij.sav'

  mv_mod = 4. + 12.0*dindgen(30)/29.0
  n_mod = 0.*mv_mod

  for i=0,n_elements(mv_mod)-1 do begin

     mv1 = mv_mod[i]-0.5
     mv2 = mv_mod[i]+0.5
     vmi1 = poly(mv1,poly_vmi_given_v)
     vmi2 = poly(mv2,poly_vmi_given_v)
     mi1 = mv1-vmi1
     mi2 = mv2-vmi2
     imj1 = poly(mi1,poly_imj_given_i)
     imj2 = poly(mi2,poly_imj_given_i)
     mj1 = mi1-imj1
     mj2 = mi2-imj2

     n_mod[i] = qromb('jlf_model', mj1, mj2, k=6)
;     print, mv1, mv2, mi1, mi2, mj1, mj2, n_mod[i]

  endfor
  oplot, mv_mod, n_mod, color=fsc_color('Red')

  x0 = 1
  y0 = 65
  dy = 4

  xyouts, x0, y0, alignment=0, 'Reid & Hawley (2005)'
  xyouts, x0, y0-dy, alignment=0, 'RECONS 2009.0', color=fsc_color('Cadet Blue')
  xyouts, x0, y0-2*dy, alignment=0, 'Zheng et al. (2001)', color=fsc_color('Orange')
  xyouts, x0, y0-3*dy, alignment=0, 'Model', color=fsc_color('Red')

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

pro compare_vlfs,ps=ps

  !p.multi=0
  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='compare_vlfs.eps',/encapsulated,xsize=6.5,ysize=5.0,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.25
     syms=0.5

  endif

  readcol, 'recons_vlf.txt', f='D,D', mv_recons, n_recons, /silent
  readcol, 'rgh_vlf.txt', f='D,D,D', mv_rgh, phi_1, phi_tot, /silent
  readcol, 'rh_vlf.txt', f='D,D,D', mv_rh, phi_rh, vmi_rh, /silent
  readcol, 'zheng_vlf.txt', f='D,D,D', mv_z, phi_z, dphi_z, /silent

  vol = 4./3.*!PI*(10.0^3.)
  n_rgh = phi_tot/1D4 * vol
  n_rh = phi_rh * vol
  n_z = phi_z*1d-3 * vol
 
  plot, mv_rgh, n_rgh, psym=10, $
        xtit='abs V mag', xra=[-1,23], xsty=1, $
        ytit='(stars mag!E-1!N) within 10 pc', yra=[0,70], ysty=1

  oplot, mv_recons, n_recons, psym=10, color=fsc_color('Cadet Blue')  
  oplot, mv_rh, n_rh, psym=10, color=fsc_color('Forest Green')  
  oplot, mv_z, n_z, psym=10, color=fsc_color('Orange')  

  x0 = 0
  y0 = 65
  dy = 5
  xyouts, x0, y0, 'Reid, Gizis, & Hawley (2002)'
  xyouts, x0, y0-dy, 'RECONS LF 2009.0', color=fsc_color('Cadet Blue')
  xyouts, x0, y0-2*dy, 'Reid & Hawley (2005)', color=fsc_color('Forest Green')
  xyouts, x0, y0-3*dy, 'Zheng et al. (2001)', color=fsc_color('Orange')

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

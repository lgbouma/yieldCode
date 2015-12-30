pro count_stars, ps=ps

  restore,'star_properties.sav'

  !p.charsize=3
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='count_stars.eps',/encapsulated,xsize=6.5,ysize=4,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1
     syms=0.5

  endif

  vol = (4./3.)*!PI*(10.0^3.)
  phi = dndr*(r[2]-r[1])/0.1/vol; stars per cubic parsec per radius bin
  h = 300.0

  x0 = 0.9
  y0 = 80
  q = 0.4

  ilim = [10,11,12,13,14]
  c = fsc_color(['Red','Orange','Green Yellow','Green','Blue'])
  
  for i=0,n_elements(ilim)-1 do begin

     dmax = 10.0*10.0^(0.2*(ilim[i]-imag))
     nstars = number_of_stars_expz(dmax, phi, 300.0) 

     if (i eq 0) then begin
        plot, r, nstars, /nodata, $
              xtit='stellar radius [R!Dsun!N]', xra=[0.1,1], xsty=1, $
              ytit='no. stars with I < I!Dlim!N, per radius bin', yra=[1,3e5], ysty=1, /ylog 
     endif

     oplot, r, nstars, psym=10, color=c[i]
     print, ilim[i], total(nstars)
     
     label = 'I <' + strcompress(ilim[i]) + ':  N!Dtot!N = ' + strcompress(round(total(nstars)))
     xyouts, alignment=1, x0, y0*q^i, label, color=c[i]

  endfor

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end


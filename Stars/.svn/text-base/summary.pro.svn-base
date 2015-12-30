pro summary, ps=ps

  !p.charsize=2
  plotsym,0,/fill
  syms=1

  if (keyword_set(ps)) then begin
  mwrfits, newtbl, newname


     set_plot,'ps'
     device,filen='summary.eps',/encapsulated,xsize=6.5,ysize=6.5,/inches,$
           /color,bits_per_pixel=8,/isolatin,/helvetica,decomposed=1
     !p.font=0
     !p.thick=6 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=1.
     syms=0.5

  endif

  rmin = 0.1
  rmax = 1.0
  nr = 100
  r = rmin + (rmax-rmin)*dindgen(nr)/double(nr-1)
  dr = r[1]-r[0]

; load up and plot the various point relationships

  restore, filen='poly_ij.sav'
  restore, filen='poly_vi.sav'
  restore, filen='poly_jr.sav'
  restore, filen='poly_teff_given_radius.sav'

  m = mass_given_radius(r)
  teff = poly(r,poly_teff_given_radius)
  jmag = poly(r,poly_jmag_given_radius)
  imj = poly(jmag,poly_imj_given_j)
  imag = jmag+imj
  vmi = poly(imag,poly_vmi_given_i)
  vmag = imag+vmi
  
  rra = [0.08,1.0]
  
  !p.multi=[0,2,2]

  plot, r, m, $
        xtit='radius [R!DSun!N]', xra=rra, xsty=1, $
        ytit='mass [M!DSun!N]', yra=[0.05,1.1], ysty=1

  plot, r, teff, $
        xtit='radius [R!DSun!N]', xra=rra, xsty=1, $
        ytit='effective temperature [K]', yra=[3000,6000], ysty=1

  plot, r, jmag, /nodata, $
        xtit='radius [R!DSun!N]', xra=rra, xsty=1, $
        ytit='abs mag', yra=[17,3], ysty=1
  oplot, r, jmag, color=fsc_color('Red')
  oplot, r, imag, color=fsc_color('Orange')
  oplot, r, vmag, color=fsc_color('Blue')

  x0 = 0.15
  y0 = 5
  dy = 1
  xyouts, x0, y0, 'J band', color=fsc_color('Red')
  xyouts, x0, y0+dy, 'I band', color=fsc_color('Orange')
  xyouts, x0, y0+2*dy, 'V band', color=fsc_color('Blue')

  dndr = 0.*r

  for i=0,nr-1 do begin
     
     r1 = r[i]-0.5*dr
     r2 = r[i]+0.5*dr
     jmag1 = poly(r2,poly_jmag_given_radius)
     jmag2 = poly(r1,poly_jmag_given_radius)

     dndr[i] = qromb('jlf_model', jmag1, jmag2, k=6) * (0.1/dr)
;     print, r1, r2, j1, j2, dndr[i]

  endfor

  plot, r, dndr, $
        xtit='radius [R!Dsun!N]', xra=rra, xsty=1, $
        ytit='dn/dR, stars (0.1 R!Dsun!N)!E-1!N within 10 pc'; , yra=[0,90], ysty=1

  save, filen='star_properties_pws.sav', r, m, teff, jmag, imag, vmag, dndr

  if (keyword_set(ps)) then begin
     
     device,/close
     
     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1
     
  endif

end

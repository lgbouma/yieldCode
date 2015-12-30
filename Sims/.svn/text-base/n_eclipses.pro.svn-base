function n_eclipses, period, duration
; if the transits have period P, and I observe for a continuous
; interval D, how many transits do I observe?

  n = 0.0*period
  e = period*randomu(seed,n_elements(period))

  observed = where(e lt duration, complement=missed)

  if (missed[0] ne -1) then n[missed] = 0.0

  if (observed[0] ne -1) then begin
     n[observed] = 1.0 + floor( (duration[observed] - e[observed])/period[observed] )
  endif

  return, n

end

pro test_n_eclipses

  m = 100
  period = 1.0 + fltarr(m)
  duration = 0.0 + 10.0*findgen(m)/float(m)

  n = n_eclipses(period, duration)

  !p.charsize=2
  plotsym,0,/fill
  plot, duration, n, psym=3, yra=[-1,11], xra=[-1,11], xsty=1, ysty=1

end

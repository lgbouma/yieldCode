function n_eclip, $
	period, obstime, npoint, offday, $ 
	periblank=periblank, apoblank=apoblank
; For a transit period with npoint observations of duration obstime,
; how many transits do I observe?
; Also, allow for an amount of blanking time at apogee (earth/moon) and 
;  perigee (downlink).
; Must also input the phase offset of the transit (in days)

  duration = npoint*obstime
 
  if (keyword_set(apoblank)) then apoblank=apoblank else apoblank=0.0
  if (keyword_set(periblank)) then periblank=periblank else periblank=0.0

  n = 0.0*period

  observed = where(offday lt duration, complement=missed)

  if (missed[0] ne -1) then n[missed] = 0.0

  if (observed[0] ne -1) then begin
    n[observed] = 1.0 + floor( (duration[observed] - offday[observed])/period[observed] )
    n0 = n
    for m=0,max(n)-1 do begin
      trantime = (offday+m*period) mod obstime
      ; Does the mth transit fall into the blanking zone?
      apoblanked =  where((trantime gt (obstime-periblank-apoblank)/2.0) and $
                          (trantime lt (obstime-periblank+apoblank)/2.0) and $
                          (m lt n0))
      periblanked = where((trantime gt (obstime-periblank)) and $
                          (trantime lt obstime) and $
                          (m lt n0))
      if(apoblanked[0] ne -1)  then n[apoblanked] = n[apoblanked]-1
      if(periblanked[0] ne -1) then n[periblanked] = n[periblanked]-1
    endfor
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

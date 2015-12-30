function planet_occurrence, logp, r
; gives you dN/(dr*dlogp) as a function of logp and r
; p is in days, r in earth radii
; output is number per factor-of-ten in p and unit interval in r

  p = [ 0.305163   ,   1.02878  ,   0.374655 ,  -0.0189619  ,  0.0598877]

  n_calc = 0.*r

  for i=0,n_elements(r)-1 do begin
     c = p[0] + p[3]*r[i]^2
     logp0 = p[1] + p[4]*r[i]^2
     dlogp = p[2]
     if (logp le logp0) then n_calc[i] = c*exp( -0.5*( (logp-logp0)/dlogp)^2. ) else n_calc[i]=c
  endfor

  return, n_calc

end

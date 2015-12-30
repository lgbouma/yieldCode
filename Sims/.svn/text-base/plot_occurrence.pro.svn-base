pro plot_occurrence

  np = 100

  logp = -0.5 + 2.5*dindgen(np)/double(np-1)
  dlogp = logp[1]-logp[0]
  p = 10^logp

  r1 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
  r2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
  dr = 0.5
  nr = n_elements(r1)

  !p.multi=0

  plot, p, 0.*p, /nodata, $
        xra=[0.5,100], xsty=1, /xlog, $
        yra=[0.0001,0.015], ysty=1;, /ylog

  occ = dblarr(np,nr)
  
  for j=0,nr-1 do begin
     for i=0,np-1 do begin
        logp1 = logp[i]-0.5*dlogp
        logp2 = logp[i]+0.5*dlogp
        p1 = 10^logp1
        p2 = 10^logp2
        occ[i,j] = int_planet_occurrence(p1,p2,r1[j],r2[j]) * (0.5/dr) * (0.1/dlogp)
     endfor
     oplot, p, occ[*,j]
  endfor
  
end

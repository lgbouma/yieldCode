function radius_given_mass, m, p=p

  if (keyword_set(p)) then begin
     p=p
  endif else begin
     p = [ 0.370370  ,   0.356669 ,    0.813026   ,  0.813049  ,   0.969737  ,    1.26946]
  endelse

  m1 = p[0]
  r1 = p[1]
  
  m2 = p[2]

  a1 = p[3]
  a2 = p[4]
  a3 = p[5]

  r = r1*(m/m1)^a1

  q = where(m gt m1 and m le m2)
  if (q[0] ne -1) then r[q] = r1*(m[q]/m1)^a2

  q = where(m gt m2)
  if (q[0] ne -1) then begin
     r2 = r1*(m2/m1)^a2
     r[q] = r2*(m[q]/m2)^a3
  endif

  return, r

end

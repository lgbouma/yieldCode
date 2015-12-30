pro assign_xyz, n, dmax, h, x, y, z, d
; assign xyz coordinates satisfying d<dmax

  x = dmax * 2.0*(randomu(seed,n)-0.5)
  y = dmax * 2.0*(randomu(seed,n)-0.5)
  z = h * randomu(seed,n,gamma=1)
  d = sqrt(x^2. + y^2. + z^2.)
  q = where(d gt dmax)

  n_iter = 0
  while (q[0] ne -1) do begin
     n_iter = n_iter + 1
     m = n_elements(q)
     x[q] = dmax * 2.0*(randomu(seed,m)-0.5)
     y[q] = dmax * 2.0*(randomu(seed,m)-0.5)
     z[q] = h * randomu(seed,m,gamma=1)
     d[q] = sqrt(x[q]^2. + y[q]^2. + z[q]^2.)
     q = where(d gt dmax)
  endwhile

  q = randomu(seed,n)
  lo = where(q le 0.5, complement=hi)
  z[lo] = -z[lo]

end

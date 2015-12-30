;gridmaker.pro (makes a 2D array of x values and y values, given size vector).

; e.g.
; gridmaker,numx,numy,dx,dy,ix,iy

pro gridmaker,numx,numy,dx,dy,ix,iy,type=type

n_x = (numx-1) / dx + 1
n_y = (numy-1) / dy + 1

ix=intarr(n_x,n_y)
iy=intarr(n_x,n_y)
xline=indgen(n_x)*dx
yline=indgen(n_y)*dy
wx=replicate(1,n_x)
wy=replicate(1,n_y)
ix(*,*)=xline#wy
iy(*,*)=wx#yline

if n_elements(type) eq 0 then return

if type eq 'float' then begin
  ix = float(ix)
  iy = float(iy)
endif

if type eq 'double' then begin
  ix = double(ix)
  iy = double(iy)
endif

return
end

PRO calc_cr_noise, sout, counts=counts
  imglen = 8
  imgnum = 256L
  nblock = 15
  blocksz = 2L^(indgen(nblock))
  counts = dblarr(imgnum*imgnum, imglen*imglen)
  stds = dblarr(nblock, imglen*imglen)
  frac = mrdfits('bigfrac24_105_f3p33.fits')
  cr = mrdfits('cosmics_5p0persecond_120seconds_0000_c1.fits')
  mask2d = intarr(16,16)
  mid = indgen(imglen)+imglen/2
  mask1 = mask2d
  mask2 = mask2d
  mask1[*,mid] = 1
  mask2[mid,*] = 1
  mask2d = mask1*mask2
  mask1d = reform(mask2d, 16*16)

  thisfrac = frac[0,0,1,*,3]
  midfrac = frac[where(mask1d)]
  sind = reverse(sort(midfrac))
  ;sind = sind[0:(imglen*imglen)-1]
  for ii=0,imgnum-1 do begin
    for jj=0, imgnum-1 do begin
      xind = indgen(imglen) + imglen*ii
      yind = indgen(imglen) + imglen*jj
      ;print, ii, jj, median(xind), median(yind)
      crx = cr[*,xind]
      crxy = crx[yind,*]
      counts[imgnum*ii+jj, *] = total(crxy[sind], /cumulative)
    end
  end
  ; cut into blocks and find std and mean vs. co-adding time
  for kk=0,nblock-1 do begin
    c3d = reform(counts, blocksz[kk], imgnum*imgnum/blocksz[kk], 64)
    m2d = mean(c3d, dimension=1)
    s2d = stddev(m2d, dimension=1)
    stds[kk,*] = s2d
  end
 ; s = stddev(counts, dimension=1)
   c2d = reform(counts, 512, 65536/512, 64)
   sout = stddev(c2d, dimension=2)
END

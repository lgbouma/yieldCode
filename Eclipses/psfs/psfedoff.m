function frac1d = psfedoff(dat, dx, dy)
tot = trapz(trapz(dat));
%dat = load(fname);
npix=12;
%pixsc=60; 
pixsc=101;
%16 for 60, 
%24 for 250, 
%124 for 225, 
%276 for 187, 
%184 for 168
% 324 for 140
% 74 for 190
% 212 for 225
frac2d = zeros(npix, npix);
%r2d = zeros(npix);
%121 for Deb files, 122 for Kris files
 for ii=1:npix
     for jj=1:npix
         thisx = (((ii-1)*pixsc):(ii*pixsc))+1+round(dx*pixsc/10);
         thisy = (((jj-1)*pixsc):(jj*pixsc))+1+round(dy*pixsc/10);
         %[min(thisx) max(thisx) min(thisy) max(thisy)]
         frac2d(ii,jj) = trapz(trapz(dat(thisx,thisy)));
         %r2d(ii,jj) = sqrt((mean(thisx)-xcen)^2+(mean(thisy)-ycen)^2);
     end
 end
 imagesc(frac2d);
 drawnow;
 
 frac1d = reshape(frac2d,1,npix*npix)/tot;
 %r1d = reshape(r2d,1,npix*npix);
 
 %frac = cumsum(sort(frac1d, 'descend'))/tot;
 %frac = frac1d/tot;
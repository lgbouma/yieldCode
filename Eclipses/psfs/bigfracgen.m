dat = load('TESS_wSi_Sept23model_ASBUILT_75C_foc0.01.mat');
fout = 'dfrac_asbuilt_75c+1f.fits'; % file name contains temperature and back focus
dither = true; % file name is dfrac if true, frac if not

angs = [0 5 12.5 15.9];
pixs = angs/24;
numfilt = 9;
nrot = 2;
nang = 4;
pixsc = 101;

onedat = dat.PRF_lambda(1,1).PSFimage;
imsize = size(onedat);
imlen = imsize(1);
[xx, yy] = meshgrid(1:imlen, 1:imlen);


% Time step is 0.2 s, so need 600 points for 2-minute stack
njit = 600;
% Errors are in pixels, so mult by pixsc for sub-pixel sampling
if dither
    err = load('atterr.dat');
    adx = err(5000:(5000+njit-1), 2)*pixsc;
    ady = err(5000:(5000+njit-1), 3)*pixsc;
    adx = adx-median(adx);
    ady = ady-median(ady);
end

npix=144;
%masterfrac = zeros(10,10,4,npix);
bigfrac = zeros(20,10,nang,npix,numfilt);
%bigrad = zeros(20,10,nang,npix,numfilt);
for rr=1:nrot
  for mm=1:nang
    for nn=1:numfilt    
      frac = zeros(10,10,npix);
      thisdat = dat.PRF_lambda(mm,nn).PSFimage;
      shiftdat = zeros(size(thisdat));
      if dither
        for jj=1:njit
         if (~mod(jj,100)) display(num2str(jj)); end
         shiftdat = shiftdat + ...
             interp2(xx, yy, thisdat, xx+adx(jj), yy+ady(jj), 'linear', 0);       
        end
      else
        shiftdat = thisdat;
      end
      if (rr==2 || mm==4)
          shiftdat = imrotate(shiftdat, 45, 'bilinear', 'crop');
      end
      %[indx indy] = centroid(thisdat);
      for ii=0:9
        for jj=0:9
            frac = psfedoff(shiftdat, ii, jj);
            %bigrad(ii*rr,jj,mm,:,nn) = r1d;
            bigfrac(((ii+1)+(rr-1)*10),jj+1,mm,:,nn) = frac;    
        end
      end
    end
  end
end

fitswrite(bigfrac, fout)

%{
break;

plotx = 10;
ploty = 10;
figure;
plot(squeeze(frac_onaxis(plotx,ploty,:)), 'k');
hold on
plot(squeeze(frac_halfedge(plotx,ploty,:)), 'b');
plot(squeeze(frac_edge(plotx,ploty,:)), 'g');
plot(squeeze(frac_corner(plotx,ploty,:)), 'r');

figure;
semilogx(squeeze(mean(mean(frac_onaxis))), 'k');
hold on
semilogx(squeeze(mean(mean(frac_halfedge))), 'b');
semilogx(squeeze(mean(mean(frac_edge))), 'g');
semilogx(squeeze(mean(mean(frac_corner))), 'r');
xlabel('Number of Pixels in Aperture');
ylabel('Fraction of Enclosed Flux');
legend('On Axis', 'Half Edge', 'Edge', 'Corner');
%}

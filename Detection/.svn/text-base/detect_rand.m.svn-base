% 1-minute cadence to data, so all times in minutes
% Campaign duration:
D = round(27.34*24*60);
% Minimum Period:
Pmin = round(4*60);
% Maximum Period
Pmax = D/2;
% Period grid
%dp = (1:10000)/2571.6;
%dp = ones(10000,1);
%P = round(cumsum(dp)+Pmin);
P=Pmin:Pmax;
NP = length(P);
% Transit width for circular orbit:
W = (1.3*60)*(P/(24*60)).^(1/3);
%Wmin = floor(78*(4.6/24)^(1/3));
%Wmax = ceil(78*(27.34)^(1/3));
%W = Wmin:Wmax;
%NW = length(W);
% Number of transits per trial period:
NT = floor(D./P);
% Number of widths
NW = floor(P./W);
% Re-do width
W = floor(P./NW);
% Period grid
%P = round(NT.*(W/78).^3)./NT;
%Number of stars
Nstar = 5e3;
% Mean of transit trials:
TSM = zeros(1,sum(NW));
tsmind = cumsum(NW);
% Transit sub-sampling
NSUB = 2;
% Save threshold
NSIG = 4;
% Number of averaged points
L = W.*NT;
% STD in phase-folded light curve
sig_sum = 1./sqrt(L);
% Events out
hisigs = [];

for jj=1:Nstar
  % Normally-distributed timeseries:
  TS = randn(D,1);
  for kk=1:NSUB
    for ii=1:NP   
      % circ-shift for sub-sampling
      TS = circshift(TS, round((kk-1)*W(ii)/NSUB));
      % Truncate the lightcurve to a multiple of the period
      TST = TS(1:(W(ii)*NT(ii)*NW(ii)));
      % Phase-fold and Period-fold the light curve
      TSP = reshape(TST,W(ii),NT(ii),NW(ii));
      % Take mean at each possible transit
      thistsm = squeeze(mean(squeeze(mean(TSP))))/sig_sum(ii);
      if (ii==1)
        TSM(1:tsmind(ii)) = thistsm;
      else
        TSM(tsmind(ii-1)+1:tsmind(ii)) = thistsm;
      end
    end
    events = (TSM>NSIG);
    nevents = sum(events);
    display([num2str(nevents) ' events on trial ' num2str(jj)]);
    if (nevents>0)
      hisigs = [hisigs TSM(events)];
    end 
  end
end

save 'hisigs.mat' hisigs
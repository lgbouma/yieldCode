% 2-minute cadence to data, so all times in minutes
% Campaign duration:
etime=2;
D = round(355*24*60/etime); % days * hours * exposures/hr

% Transit width for circular orbit:
%W = (1.3*60)*(P/(24*60)).^(1/3);
Wmin = floor(78*(4.6/24)^(1/3)/etime);
Wmax = ceil(78*(355/2)^(1/3)/etime);
W = Wmin:Wmax;

% Number of transit widths:
NW = floor(D./W);

% Period in exposures
P = 24*60*(etime*W/78).^3/etime;

% Period in width range
Pwmin = floor(P./W);
Pwmax = ceil(P./W);

% Number of guaranteed transits
NTmin = ceil(NW./Pwmin);
NTmax = ceil(NW./Pwmax);

%Number of stars
Nstar = 6e5;
% Mean of transit trials:
TSM = zeros(1,sum(NW));
tsmind = cumsum(NW);
% Transit sub-sampling
NSUB = 2;
% Save threshold
NSIG = 4;
% Number of averaged points
%L = W.*NT;
% STD in phase-folded light curve
%sig_sum = 1./sqrt(L);
% Events out
sigs = 4:0.1:8;
num = 1:25;

for p=1:25
  hisighist = zeros(length(sigs),1);
  for jj=1:Nstar
    % Normally-distributed timeseries:
    TS = randn(D,1);
    for kk=1:NSUB
      for ii=1:length(W)  
        % circ-shift for sub-sampling
        if kk>1
          os = round((kk-1)*W(ii)/NSUB);
          % Truncate the lightcurve to a multiple of the widths
          TST = [TS((os+1):(W(ii)*NW(ii))); TS(1:os)];
        else
          TST = TS(1:(W(ii)*NW(ii)));
        end
        % Width-fold the lightcurve
        % Take mean over each possible transit
        TSW = mean(reshape(TST,NW(ii),W(ii)),2);
        cts = ones(size(TSW));
        
        break;
      
        % Zero-pad
        TSWmin = [TSW; zeros(NTmin(ii)*Pwmin(ii)-length(TSW),1)];
        ctsmin = [cts; zeros(NTmin(ii)*Pwmin(ii)-length(TSW),1)];
        TSWmax = [TSW; zeros(NTmax(ii)*Pwmax(ii)-length(TSW),1)];
        ctsmax = [cts; zeros(NTmax(ii)*Pwmax(ii)-length(TSW),1)];
        csmin = sum(reshape(ctsmin, Pwmin(ii), NTmin(ii)),2);
        csmax = sum(reshape(ctsmax, Pwmax(ii), NTmax(ii)),2);
      
        % Modified mean
        TSPmin = sum(reshape(TSWmin, Pwmin(ii), NTmin(ii)),2)./csmin;          
        TSPmax = sum(reshape(TSWmax, Pwmax(ii), NTmax(ii)),2)./csmax;
      
        % Take the ntransit>1 cases 
        TSMmin = TSPmin(csmin>1).*sqrt(W(ii)*csmin(csmin>1));
        TSMmax = TSPmax(csmax>1).*sqrt(W(ii)*csmax(csmax>1));
        hisighist = hisighist + histc(TSMmin, sigs);
        hisighist = hisighist + histc(TSMmax, sigs);
%       minevent = (TSMmin>NSIG);
%       maxevent = (TSMmax>NSIG);
%       nminevents = sum(minevent);
%       nmaxevents = sum(maxevent);
%       %display([num2str(nevents) ' events on trial ' num2str(jj)]);
%       if (nminevents>0)
%         hisigs = [hisigs TSMmin(minevent)'];
%         %display([num2str(nminevents) ' events on trial ' num2str(jj)]);
%       end 
%       if (nmaxevents>0)
%         hisigs = [hisigs TSMmax(maxevent)'];
%       end  
      end
    end
    if ~mod(jj,100) 
        display(['Trial ' num2str(jj)]);
    end
  end
  fname = ['sighist' num2str(num(p)) '.mat'];
  save(fname, 'hisighist');
end
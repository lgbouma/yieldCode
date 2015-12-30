parfor p=1:2
% 2-minute cadence to data, so all times in minutes
% Campaign duration:
etime=2; % in minutes

% Timeseries length
D = round(178*2*24*60/etime); % days * hours * exposures/hr
C = ones(D,1);

% Transit width for circular orbit:
%W = (1.3*60)*(P/(24*60)).^(1/3);
%Wmax1 = floor(78*(4.6/24)^(1/3)/etime);
%Wmax2 = ceil(78*(355/2)^(1/3)/etime);
%Wmax = Wmax1:Wmax2;
PmaxDays = 178; %355ish/2
Wmax = round(13*60*(PmaxDays/365)^(1/3)*0.5^(-1/3)/etime); %1228 mins/etime

% Max. Period in exposures
Pmax = PmaxDays*24*60/etime;

% Min. width corresponding to b~1
Wmin = 2;
Wstep = 2;
W = Wmin:Wstep:Wmax;

% Min. Period in exposures (Roche limit)
Proach = 1.7*60/etime; 
Pmin   = 365*24*60/etime*(W/(13*60/etime)).^3*0.5;
Pmin(Pmin < Proach) = Proach;

% Number of widths
NW = floor(D./W);

% Periods measured in widths
PminW = ceil(Pmin./W);
PmaxW = floor(Pmax./W);

% Number of transits
NTmin = floor(NW./PminW);
NTmax = floor(NW./PmaxW);

% Transit sub-sampling
NSUB = 2;
% Save threshold
NSIG = 4;

%Number of stars
%Nstar = 6e5;
Nstar = 1e5;

% Mean of transit trials:
TSM = zeros(1,sum(NW));
tsmind = cumsum(NW);

% Number of averaged points
%L = W.*NT;
% STD in phase-folded light curve
%sig_sum = 1./sqrt(L);
% Events out
sigs = 4:0.1:8;
num = 1:25;

  hisighist = zeros(length(sigs),1);
  losighist = zeros(length(sigs),1);
  for jj=1:Nstar
    % Normally-distributed timeseries:
    TS = randn(D,1);    
    for ii=1:length(W)
      for kk=1:2        
        % Truncate the lightcurve to a multiple of the widths
        if kk>1
          % But first, circ-shift for sub-sampling
          os = W(ii)/Wstep;          
          TST = [TS((os+1):(W(ii)*NW(ii))); TS(1:os)];
        else
          TST = TS(1:(W(ii)*NW(ii)));
        end
        CT = C(1:(W(ii)*NW(ii)));
        % Width-fold the lightcurve,
        % Take mean over each possible transit
        TSW = mean(reshape(TST,NW(ii),W(ii)),2);
        CW = sum(reshape(CT, NW(ii),W(ii)),2);
      
        PW = PminW(ii):PmaxW(ii);
        NT = floor(NW(ii)./PW);
        for pp=1:length(PW)
            
          % Truncate again
          TSWT = TSW(1:(PW(pp)*NT(pp)));  
          CWT = CW(1:(PW(pp)*NT(pp)));
          
          % Re-shape
          TSP = mean(reshape(TSWT, PW(pp), NT(pp)), 2);
          CWP = sum(reshape(CWT, PW(pp), NT(pp)), 2);
         
          hisighist = hisighist + histc(TSP.*sqrt(CWP), sigs);
          losighist = losighist + histc(abs(TSP.*sqrt(CWP)), sigs);
          
        end
 
      end
    end
    if ~mod(jj,100) 
        display(['Trial ' num2str(jj)]);
    end
  end
  fname = ['bsighist' num2str(num(p)) '.mat'];
  parsave(fname, hisighist, losighist);
end
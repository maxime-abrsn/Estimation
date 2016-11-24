
function [avgBeta, avgCorr, avgErm, mpr, volMkt, adsRep, trDates] = ...
    fEmpData(eDATA, rf, erm, ADSindex, dates, nBins, start, stop, rolling, rollw, rollstep);

% some constants
T = size(eDATA,1);
N = size(eDATA,2);

% regression to get beta with missing values
% get size
if rolling == 0
    aa = [find(diff(round(dates./10000))==1)+1; length(dates)+1]; % get starting position of following year
       % the "find" gives final aa value which is the first date of the last year available in the sample
       % the "length" add the final date value + 1 (ie. aa-1 will give the last available date of each year)
    numberYear = length(aa);
end

tmp = size(eDATA,2);
if rolling == 0
    limit = numberYear;
elseif rolling == 1
    limit = floor((T-rollw)/rollstep)+1;
else
    disp('Please specify the type of window length');
    stop
end
  
beta = zeros(limit,tmp);
betaSort = zeros(limit,tmp);

correlations = zeros(tmp,tmp,limit);

voltot = nan(limit,tmp);  % total volatility
volsys = nan(limit,tmp);
volidio = nan(limit,tmp);

avgErm = nan(limit,1);  % average expected market return    
adsRep = nan(limit,1);  % to store average ads for the window
mpr = nan(limit,1);     % to store market price of risk (market portfolio)
volMkt = nan(limit,1);  % total volatility of market portfolio

counter = 0;
for i=1:limit     
    if rolling == 0
        trDates(i,1) = dates(aa(i)-1); % last date of that window
        if i==1
            er = eDATA(1:aa(i)-1,:);
            x  = erm(1:aa(i)-1,:); 
            rfi = rf(1:aa(i)-1,:); 
            adsRep(i) = nanmean(ADSindex(1:aa(i)-1,2));
        else
            er = eDATA(aa(i-1):aa(i)-1,:);
            x = erm(aa(i-1):aa(i)-1,:);
            rfi = rf(aa(i-1):aa(i)-1,:);
            adsRep(i) = nanmean(ADSindex(aa(i-1):aa(i)-1,2));
        end
    elseif rolling == 1
        trDates(i,1) = dates(rollstep*(i-1)+rollw);   % last date of that window
        er = eDATA(rollstep*(i-1)+1:rollstep*(i-1)+rollw,:);
        x = erm(rollstep*(i-1)+1:rollstep*(i-1)+rollw,:);
        rfi = rf(rollstep*(i-1)+1:rollstep*(i-1)+rollw,:);
        adsRep(i) = nanmean(ADSindex(rollstep*(i-1)+1:rollstep*(i-1)+rollw,2));
    end 
    
    % OLS regression
    
    [nobs, nvar] = size(x);   
    x = [ones(nobs,1) x];    % include constant in regression
    nvar = nvar + 1;
    
    if nobs < 10000
        [q, r] = qr(x,0);
        xpxi = (r'*r)\eye(nvar);
    else 
        xpxi = (x'*x)\eye(nvar);
    end;
    mbeta = xpxi*(x'*er);    
    beta(i,:) = mbeta(2,:);   % constant first, then E(Rm - Rf)
    
    % sort according to beta
    tmp = sortrows([beta(i,:)' er'], 1);
    erSort = tmp(:,2:end)';
    betaSort(i,:) = tmp(:,1)';        
      
    % pairwise (quarterly, yearly) correlations using daily estimates
    correlations(:,:,i) = corr(erSort);
    clear tmp
    
    % volatility using daily estimate
    voltot(i,:) = nanstd(erSort);
    volsys(i,:) = betaSort(i,:) .* nanstd(x(:,2));
    volidio(i,:) = sqrt(voltot(i,:).^2 - volsys(i,:).^2);
    
    % average market return and market price of risk
    avgErm(i) = nanmean(x(:,2));  % first column is just constant
    mpr(i) = (avgErm(i) - mean(rfi))/var(x(:,2));
    
    % volatility of market portfolio
    volMkt(i) = nanstd(x(:,2));
end


% calculate average pairwise correlation within each bin 
for i = 1:limit
    % remove auto-correlations and stack triangular values
    tm = triu(correlations(:,:,i),1);
    tm = tm.';
    value = tm(find(tm));
    avgCorr(i,1) = nanmean(value); 

    for j=1:nBins
        tmp = triu(correlations(start(j):stop(j),start(j):stop(j),i),1);
        tmp = tmp.';
        values = tmp(find(tmp));
        avgCorr(i,j+1) = nanmean(values);     
    end
end

clear tm tmp tmp2 tmp3 value values values2 values3

% average beta for each bin

avgBetaAll = nanmean(betaSort,2);
tmp = sqrt(252);
avgVoltotAll = tmp* nanmean(voltot,2);
avgVolsysAll = tmp *nanmean(volsys,2);
avgVolidioAll = tmp* nanmean(volidio,2);
for j=1:nBins
    avgBeta(:,j) = nanmean(betaSort(:,start(j):stop(j)),2);
    avgVoltot(:,j) = tmp * nanmean(voltot(:,start(j):stop(j)),2);  % annualized
    avgVolsys(:,j) = tmp * nanmean(volsys(:,start(j):stop(j)),2);
    avgVolidio(:,j) = tmp * nanmean(volidio(:,start(j):stop(j)),2);
end
volMkt = tmp * volMkt;  



end
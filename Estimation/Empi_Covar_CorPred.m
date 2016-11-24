clear all
close all
clc

%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose sample (raw DAILY returns (DATA, TxN), dates(Tx1), NOT excess return)
% should start before 19630102
% should end on 20141231 (but it is working well is end date is earlier)

%%% OPTIONS %%%%
portfolio = 4; %(0) for 48FF, (1) for 100BtM (2) for beta decile, (3) for 10BtM, (4) for top 500 cap 20beta

startDate =  19710104; %19710104 (middle: 19931207) %19630102 (middle: 19900122); % cut the sample to interesting times
endEstDate = 20141231;  % end the estimation window

rolling = 0; % rolling window or not, 1 or 0
rollw = 216; % window size, 1 year
rollstep = 1;  % step size for rolling, half year

saveGraph = 0;
saveRes =0 ;
saveStatics = 0 ;

origin_Pads=0 ; % 0 for linear mapping ADS -> P, 1 for non-linear mapping (see line 299 -> .. )


checkreg = 1 ;  % to perform regression
constant = 1; % regression with contant (1) or not (0)
averageBeta = 1; % averageBeta to compute th. corr (1), mapping for every stock (0), for 2 nearest stock to average beta (2)




%==========================================================================================
% Choose data

if portfolio == 0
    load('48FF_daily.mat')
elseif portfolio ==1
    load('100BTM_Size_daily.mat')
elseif portfolio == 2
    load('beta_decile_daily.mat') % average beta below 1 at some point, possible?
elseif portfolio == 3
    load('10BTM_Size_daily_corrected.mat') % average beta below 1 at some point, possible?
    %load('10BTM_Size_daily.mat') % average beta below 1 at some point, possible?
elseif portfolio ==4
    load('top500_betaYr_20BetaBin_daily.mat')
    DATA = DATA.*100; % to make it as the others 4 datasets
end

 
%============================================================================
% Divide data into different bins


% Choose number of beta bins to divide data into
%--------------------------------------------------------
% French data are in percent, change the value if data are already divided
cent = 100;
if (portfolio == 2 | portfolio == 3) 
    % position of interest to generate corr dynamic
    pos1 = 1;
    pos2 = 2;
    % number of bins
    nBins = 2;
    gap = 0;% gap between bins
    
else 
    % position of interest to generate corr dynamic
    pos1 = 2;
    pos2 = 3;
    % number of bins
    nBins = 4;
    gap = 0;%10;% gap between bins
end

if gap*(nBins-1) >= size(DATA,2)
    disp('The Gap is too big')
    stop
end

% cut the sample in n bins of equivalent size
%-------------------------------------------------
n = floor((size(DATA,2)-gap*(nBins-1))/nBins); % multiple of 3 are better to avoid discarding data

for i=1:nBins
    if i == 1
        start(i) = 1;
        stop(i) = n;
    elseif i == nBins
        start(i) = n*(i-1) + 1 + (i-1) * gap;
        stop(i) = size(DATA,2);
    else
        start(i) = n*(i-1) + 1 + (i-1)* gap;
        stop(i) = n*i + (i-1)* gap;
    end
end

%========================================================================
% load factors (riskfree rate, mkt - rf, ADS index)

% range 19260701 to 20141231
load('FFerm_rf.mat')

% ADS
% range 19600301 to 20141227
load ADSindex.mat

%=======================================================================
% find dataset 

% The end data date = (btw DATA, factor and ADS) that has lowest date

if ADSindex(end,1) <= dates(end) && ADSindex(end,1) <= FFdates(end)
    finalDate = ADSindex(end,1);
elseif dates(end) <= FFdates(end);
    finalDate = dates(end);
else
    finalDate = FFdates(end);
end
if(finalDate < endEstDate)
    disp('Estimation window is too long, above data-end date');
    
end
    

% cut below startDate and above finalDate
tmp = find(dates<startDate | dates > finalDate);
DATA(tmp,:) = [];
dates(tmp,:) = [];
clear tmp
tmp = find(FFdates < startDate | FFdates > finalDate);
FF3facDaily(tmp,:) = [];
FFdates(tmp,:) = [];
clear tmp
tmp = find(ADSindex(:,1) < startDate | ADSindex(:,1) > finalDate);
ADSindex(tmp,:) = [];
clear tmp


% ADS has more entries due to weekends
tmp = find(ismember(ADSindex(:,1),dates));
ADSindex = ADSindex(tmp,:);


% replace missing values
tmp = find(DATA == -99.99 | DATA ==-999);
DATA(tmp) = NaN;
clear tmp 

% divide by 100 to get decimal returns
eDATA = DATA/cent;


% select some data of interest
rf = FF3facDaily(:,5)/cent;
erm = FF3facDaily(:,2)/cent;   % E(Rm - Rf)

% some constants
T = size(DATA,1);
N = size(DATA,2);


%============================================================================
% Divide data into beta bins, and calculate avr pairwise correl within bins


% regression to get beta with missing values
% get size
if rolling == 0
    aa = [find(diff(round(dates./10000))==1)+1; length(dates)+1]; % get starting position of following year
       % the "find" gives final aa value which is the first date of the last year available in the sample
       % the "length" add the final date value + 1 (ie. aa-1 will give the last available date of each year)
    numberYear = length(aa);
end

tmp = size(DATA,2);
if rolling == 0
    beta = zeros(numberYear,tmp);
    betaSort = zeros(numberYear,tmp);
    quartExRet = zeros(numberYear,tmp);
    correlations = zeros(tmp,tmp,numberYear);
    limit = numberYear;
    voltot = nan(limit,tmp);  % total volatility
    volsys = nan(limit,tmp);
    volidio = nan(limit,tmp);
elseif rolling == 1
    limit = floor((T-rollw)/rollstep)+1;
    beta = zeros(limit,tmp);
    betaSort = zeros(limit,tmp);
    quartExRet = zeros(limit,tmp);
    correlations = zeros(tmp,tmp,limit);
    voltot = nan(limit,tmp);
    volsys = nan(limit,tmp);
    volidio = nan(limit,tmp);
else
    disp('Please specify the type of window length');
    stop
end
adsRep = nan(limit,1);  % to store average ads for the window

counter = 0;
for i=1:limit     
    if rolling == 0
        trDates(i,1) = dates(aa(i)-1); % last date of that window
        if i==1
            er = eDATA(1:aa(i)-1,:);
            x  = erm(1:aa(i)-1,:); 
            adsRep(i) = nanmean(ADSindex(1:aa(i)-1,2));
        else
            er = eDATA(aa(i-1):aa(i)-1,:);
            x = erm(aa(i-1):aa(i)-1,:);
            adsRep(i) = nanmean(ADSindex(aa(i-1):aa(i)-1,2));
        end
    elseif rolling == 1
        trDates(i,1) = dates(rollstep*(i-1)+rollw);   % last date of that window
        er = eDATA(rollstep*(i-1)+1:rollstep*(i-1)+rollw,:);
        x = erm(rollstep*(i-1)+1:rollstep*(i-1)+rollw,:);
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

% average beta for each bin (all in first columns, then each bin)
avgBeta = nanmean(betaSort,2);
tmp = sqrt(252);
avgVoltot = tmp* nanmean(voltot,2);
avgVolsys = tmp *nanmean(volsys,2);
avgVolidio = tmp* nanmean(volidio,2);
for j=1:nBins
    avgBeta(:,j+1) = nanmean(betaSort(:,start(j):stop(j)),2);
    avgVoltot(:,j+1) = tmp * nanmean(voltot(:,start(j):stop(j)),2);  % annualized
    avgVolsys(:,j+1) = tmp * nanmean(volsys(:,start(j):stop(j)),2);
    avgVolidio(:,j+1) = tmp * nanmean(volidio(:,start(j):stop(j)),2);
end




%=================================================================
% extreme market conditions

ads1 = prctile(adsRep,10);
ads2 = prctile(adsRep,90);
idxnorm = (adsRep>ads1 & adsRep<ads2); % normal time
idxex = (adsRep<ads1 | adsRep>ads2);  % extreme value
idxexB = (adsRep<ads1);  % bad time
idxexG = (adsRep>ads2);  % good tim


statEmpi.beta = nanmean(avgBeta(:,1)); % avg beta
statEmpi.betaLow = nanmean(avgBeta(:,2)); % avg low beta
statEmpi.betaHigh = nanmean(avgBeta(:,end)); % avg high beta

statEmpi.corr = nanmean(avgCorr(:,1)); %avg corr sample
statEmpi.corrLow = nanmean(avgCorr(:,2)); %avg corr low beta sample
statEmpi.corrHigh = nanmean(avgCorr(:,end)); %avg corr high beta sample
statEmpi.corrstd = nanstd(avgCorr(:,1)); 
statEmpi.corrstdL = nanstd(avgCorr(:,2)); 
statEmpi.corrstdH = nanstd(avgCorr(:,end));

statEmpi.correx = nanmean(avgCorr(idxex,1)); 
statEmpi.correxL = nanmean(avgCorr(idxex,2)); 
statEmpi.correxH = nanmean(avgCorr(idxex,end)); 
statEmpi.correxB = nanmean(avgCorr(idxexB,1)); 
statEmpi.correxBL = nanmean(avgCorr(idxexB,2)); 
statEmpi.correxBH = nanmean(avgCorr(idxexB,end)); 
statEmpi.correxG = nanmean(avgCorr(idxexG,1)); 
statEmpi.correxGL = nanmean(avgCorr(idxexG,2)); 
statEmpi.correxGH = nanmean(avgCorr(idxexG,end)); 
statEmpi.corrnorm = nanmean(avgCorr(idxnorm,1)); 
statEmpi.corrnormL = nanmean(avgCorr(idxnorm,2)); 
statEmpi.corrnormH = nanmean(avgCorr(idxnorm,end)); 

statEmpi.corrSex = nanstd(avgCorr(idxex,1)); 
statEmpi.corrSexL = nanstd(avgCorr(idxex,2)); 
statEmpi.corrSexH = nanstd(avgCorr(idxex,end)); 
statEmpi.corrSexB = nanstd(avgCorr(idxexB,1)); 
statEmpi.corrSexBL = nanstd(avgCorr(idxexB,2)); 
statEmpi.corrSexBH = nanstd(avgCorr(idxexB,end)); 
statEmpi.corrSexG = nanstd(avgCorr(idxexG,1)); 
statEmpi.corrSexGL = nanstd(avgCorr(idxexG,2)); 
statEmpi.corrSexGH = nanstd(avgCorr(idxexG,end)); 
statEmpi.corrSnorm = nanstd(avgCorr(idxnorm,1)); 
statEmpi.corrSnormL = nanstd(avgCorr(idxnorm,2)); 
statEmpi.corrSnormH = nanstd(avgCorr(idxnorm,end)); 

statEmpi.voltot = nanmean(avgVoltot(:,1));
statEmpi.voltotLow = nanmean(avgVoltot(:,2));
statEmpi.voltotHigh = nanmean(avgVoltot(:,end));
statEmpi.voltotex = nanmean(avgVoltot(idxex,1));
statEmpi.voltotexL = nanmean(avgVoltot(idxex,2));
statEmpi.voltotexH = nanmean(avgVoltot(idxex,end));
statEmpi.voltotexB = nanmean(avgVoltot(idxexB,1));
statEmpi.voltotexBL = nanmean(avgVoltot(idxexB,2));
statEmpi.voltotexBH = nanmean(avgVoltot(idxexB,end));
statEmpi.voltotexG = nanmean(avgVoltot(idxexG,1));
statEmpi.voltotexGL = nanmean(avgVoltot(idxexG,2));
statEmpi.voltotexGH = nanmean(avgVoltot(idxexG,end));
statEmpi.voltotnorm = nanmean(avgVoltot(idxnorm,1));
statEmpi.voltotnormL = nanmean(avgVoltot(idxnorm,2));
statEmpi.voltotnormH = nanmean(avgVoltot(idxnorm,end));

statEmpi.volsys = nanmean(avgVolsys(:,1));
statEmpi.volsysLow = nanmean(avgVolsys(:,2));
statEmpi.volsysHigh = nanmean(avgVolsys(:,end));
statEmpi.volsysex = nanmean(avgVolsys(idxex,1));
statEmpi.volsysexL = nanmean(avgVolsys(idxex,2));
statEmpi.volsysexH = nanmean(avgVolsys(idxex,end));
statEmpi.volsysexB = nanmean(avgVolsys(idxexB,1));
statEmpi.volsysexBL = nanmean(avgVolsys(idxexB,2));
statEmpi.volsysexBH = nanmean(avgVolsys(idxexB,end));
statEmpi.volsysexG = nanmean(avgVolsys(idxexG,1));
statEmpi.volsysexGL = nanmean(avgVolsys(idxexG,2));
statEmpi.volsysexGH = nanmean(avgVolsys(idxexG,end));
statEmpi.volsysnorm = nanmean(avgVolsys(idxnorm,1));
statEmpi.volsysnormL = nanmean(avgVolsys(idxnorm,2));
statEmpi.volsysnormH = nanmean(avgVolsys(idxnorm,end));

statEmpi.volidio = nanmean(avgVolidio(:,1));
statEmpi.volidioLow = nanmean(avgVolidio(:,2));
statEmpi.volidioHigh = nanmean(avgVolidio(:,end));
statEmpi.volidioex = nanmean(avgVolidio(idxex,1));
statEmpi.volidioexL = nanmean(avgVolidio(idxex,2));
statEmpi.volidioexH = nanmean(avgVolidio(idxex,end));
statEmpi.volidioexB = nanmean(avgVolidio(idxexB,1));
statEmpi.volidioexBL = nanmean(avgVolidio(idxexB,2));
statEmpi.volidioexBH = nanmean(avgVolidio(idxexB,end));
statEmpi.volidioexG = nanmean(avgVolidio(idxexG,1));
statEmpi.volidioexGL = nanmean(avgVolidio(idxexG,2));
statEmpi.volidioexGH = nanmean(avgVolidio(idxexG,end));
statEmpi.volidionorm = nanmean(avgVolidio(idxnorm,1));
statEmpi.volidionormL = nanmean(avgVolidio(idxnorm,2));
statEmpi.volidionormH = nanmean(avgVolidio(idxnorm,end));

tmp = [statEmpi.correxB statEmpi.corrnorm statEmpi.correxG  ...
        statEmpi.correxBL statEmpi.corrnormL statEmpi.correxGL  ...
        statEmpi.correxBH statEmpi.corrnormH statEmpi.correxGH;  ...
       statEmpi.corrSexB statEmpi.corrSnorm statEmpi.corrSexG  ...
        statEmpi.corrSexBL statEmpi.corrSnormL statEmpi.corrSexGL  ...
        statEmpi.corrSexBH statEmpi.corrSnormH statEmpi.corrSexGH];
statEmpi.corrmat = tmp;

tmp = [statEmpi.voltotexB statEmpi.voltotnorm statEmpi.voltotexG  ...
        statEmpi.voltotexBL statEmpi.voltotnormL statEmpi.voltotexGL  ...
        statEmpi.voltotexBH statEmpi.voltotnormH statEmpi.voltotexGH;  ...
       statEmpi.volsysexB statEmpi.volsysnorm statEmpi.volsysexG  ...
        statEmpi.volsysexBL statEmpi.volsysnormL statEmpi.volsysexGL  ...
        statEmpi.volsysexBH statEmpi.volsysnormH statEmpi.volsysexGH;  ...
       statEmpi.volidioexB statEmpi.volidionorm statEmpi.volidioexG  ...
        statEmpi.volidioexBL statEmpi.volidionormL statEmpi.volidioexGL  ...
        statEmpi.volidioexBH statEmpi.volidionormH statEmpi.volidioexGH];
statEmpi.volmat = tmp;

tmp = [statEmpi.corr statEmpi.corrnorm statEmpi.correx  ...
        statEmpi.corrLow statEmpi.corrnormL statEmpi.correxL  ...
        statEmpi.corrHigh statEmpi.corrnormH statEmpi.correxH;  ...
       statEmpi.corrstd statEmpi.corrSnorm statEmpi.corrSex  ...
        statEmpi.corrstdL statEmpi.corrSnormL statEmpi.corrSexL  ...
        statEmpi.corrstdH statEmpi.corrSnormH statEmpi.corrSexH; ...
       statEmpi.voltot statEmpi.voltotnorm statEmpi.voltotex  ...
        statEmpi.voltotLow statEmpi.voltotnormL statEmpi.voltotexL  ...
        statEmpi.voltotHigh statEmpi.voltotnormH statEmpi.voltotexH;  ...
       statEmpi.volsys statEmpi.volsysnorm statEmpi.volsysex  ...
        statEmpi.volsysLow statEmpi.volsysnormL statEmpi.volsysexL  ...
        statEmpi.volsysHigh statEmpi.volsysnormH statEmpi.volsysexH;  ...
       statEmpi.volidio statEmpi.volidionorm statEmpi.volidioex  ...
        statEmpi.volidioLow statEmpi.volidionormL statEmpi.volidioexL  ...
        statEmpi.volidioHigh statEmpi.volidionormH statEmpi.volidioexH];
statEmpi.covarmat = tmp;

disp('============================================================')
disp(strcat('portfolio ',{' '},num2str(portfolio)));
disp('return, vol (direct cal), correl, correl dispersion')
[252*nanmean(nanmean(eDATA))  sqrt(252)*nanmean(nanstd(eDATA))  statEmpi.corr   statEmpi.corrstd]
%
disp('----------------------------------------')
disp('Big column: All, low beta, high beta')
disp('smaller column: bad time, normal time, good time');
disp('row: correl, correl dispersion');
statEmpi.corrmat
disp('row: voltotal, volsys, volidio');
statEmpi.volmat
%
disp('-----------------------------------------');
disp('return, voltot, correl, correl dispersion')
[252*nanmean(nanmean(eDATA))  statEmpi.voltot  statEmpi.corr   statEmpi.corrstd]
disp('Big column: All, low beta, high beta')
disp('smaller column:all time, normal time, extreme time');
disp('row: correl, correl dispersion, voltotal, volsys, volidio');
statEmpi.covarmat
% 
disp('----------------------------------------');
disp('proportion of systematic volatility to total volatility (ie. v_s^2/v_t^2)');
(statEmpi.covarmat(4,:)).^2 ./ (statEmpi.covarmat(3,:).^2)





%{
disp('row: correl, cor dispersion, cor extreme, cor normal, voltot, vol ex, vol norm, volsys, volsys ex, volsys norm, volidio, volidio ex, voldio norm');
[statEmpi.corr  statEmpi.corrLow statEmpi.corrHigh; ...
 statEmpi.corrstd  statEmpi.corrstdL statEmpi.corrstdH; ...
 statEmpi.correx  statEmpi.correxL     statEmpi.correxH; ...
 statEmpi.corrnorm  statEmpi.corrnormL   statEmpi.corrnormH; ...
 statEmpi.voltot  statEmpi.voltotLow statEmpi.voltotHigh; ...
 statEmpi.voltotex  statEmpi.voltotexL     statEmpi.voltotexH; ...
 statEmpi.voltotnorm  statEmpi.voltotnormL   statEmpi.voltotnormH; ...
 statEmpi.volsys  statEmpi.volsysLow statEmpi.volsysHigh; ...
 statEmpi.volsysex  statEmpi.volsysexL     statEmpi.volsysexH; ...
 statEmpi.volsysnorm  statEmpi.volsysnormL   statEmpi.volsysnormH; ...
 statEmpi.volidio  statEmpi.volidioLow statEmpi.volidioHigh; ...
 statEmpi.volidioex  statEmpi.volidioexL     statEmpi.volidioexH; ...
 statEmpi.volidionorm  statEmpi.volidionormL   statEmpi.volidionormH]
 %}

if saveStatics == 1
    path3 = ['./Results/statEmpi',num2str(portfolio),'.mat'];
    save(path3, 'statEmpi')
end 




%========================================================================
% infer p_t from ADS index

path = ['./results/gmmout',num2str(portfolio),'Estim.mat'];
load(path);

lambda_12 = gout.b(9);
lambda_21 = gout.b(10);
        
pHs = lambda_21 / (lambda_12 + lambda_21);
pLs = 1- pHs;

tmp2 = adsRep;

if origin_Pads==0

%disp('verify pads to be inside 0-1')
for i=1:size(tmp2,1)
    if tmp2(i,1) < 0
        pads(i,1) = pHs + tmp2(i)/((abs(min(tmp2)))/pHs);
    elseif tmp2(i,1) > 0
        pads(i,1) = pHs + tmp2(i)/((abs(max(tmp2)))/pLs);
    else
        pads(i,1) = pHs;
    end
end
clear tmp2

else
  
disp('verify pads to be inside 0-1')
for i=1:size(tmp2,1)
    if tmp2(i,1) < 0
        pads(i,1) = pHs + (1+(tmp2(i)-min(tmp2)).^1.2/3)*tmp2(i)/((abs(min(tmp2)))/pHs);
    elseif tmp2(i,1) > 0
        pads(i,1) = pHs + (1+(max(tmp2)-tmp2(i)).^1.2/3)*tmp2(i)/((abs(max(tmp2)))/pLs);
    else
        pads(i,1) = pHs;
    end
end
clear tmp2
end



%=======================================================================
% Choose estimation sample and out of sample data


% Choose estimation window and out of sample window
tmp2 = find(trDates<=endEstDate,1,'last');  % the last element that satisfies the condition
avgCorrIn = avgCorr(1:tmp2,:);
padsIn = pads(1:tmp2,:);
avgCorrOut = avgCorr(tmp2+1:end,:);
padsOut = pads(tmp2+1:end,:);
clear tmp2


%========================================================================
% Regression within sample for dependence between correlation over time
   % beta coefficient = b1 + b2*p + b3*p^2 (p = probability good state)

betav = nan(4,nBins+1);   % 4 coeffs: constant, b1, b2, b3
sebetav = nan(4,nBins+1);
r2adjv = nan(1,nBins+1);

Tin = size(avgCorrIn,1)-1;
rhv1_1 = repmat(padsIn(2:end),1,nBins+1) .*avgCorrIn(1:end-1,:);
rhv2_1 = repmat(padsIn(2:end).^2,1,nBins+1) .*avgCorrIn(1:end-1,:);

% rhv1 = repmat(padsIn(1:end-1),1,nBins+1) .*avgCorrIn(1:end-1,:);
% rhv2 = repmat(padsIn(1:end-1).^2,1,nBins+1) .*avgCorrIn(1:end-1,:);

padsIn2=padsIn-mean(padsIn);

rhv1_2= repmat(padsIn2(2:end),1,nBins+1) .*avgCorrIn(1:end-1,:);
rhv2_2 = repmat(padsIn2(2:end).^2,1,nBins+1) .*avgCorrIn(1:end-1,:);

for j = 1:2;  % a lazy way of using old codes, of course not elegant!
    if j==1;
        rhv1 = rhv1_1;
        rhv2 = rhv2_1;
        tmp = 'pads';
    else
        rhv1 = rhv1_2;
        rhv2 = rhv2_2;
        tmp = 'demeaned pads';
    end;


% regression with gmm corrected standard error
lags = 18; % lag for stderr calculation 
weight = 1; % newey-west weighting
for i = 1:(nBins+1)
    lhv = avgCorrIn(2:end,i);
    rhv = [ones(Tin,1) avgCorrIn(1:end-1,i) rhv1(:,i) rhv2(:,i)];
    [bv,sebv,R2v,R2vadj,v,F] = olsgmm(lhv,rhv,lags,weight);
    betav(:,i) = bv;
    sebetav(:,i) = sebv;
    r2adjv(i) = R2vadj;
end;

if j==1;
    disp('------------------------------------------------------')
    %portfolio
    %endEstDate
    rolling
    if(rolling==1)
        rollw
        rollstep
    end
    disp('regression of p_t on p_{t-1} using pads')
    disp('Rho_t = c + b*Rho_{t-1} + e, where b = b0 + b1*p + b2*p^2');
    disp('columns = all, very low beta, low beta, medium beta, high beta; row = c, b0, b1, b2');
else
    disp('-----------------------------------------------------')
    disp('regression of p_t on p_{t-1} using demeaned pads');
end
betav
tbetav_gmm = betav ./ sebetav
r2adjv

%{
% regression with White standard error
for i = 1:(nBins+1)
    lhv = avgCorrIn(2:end,i);
    rhv = [ones(Tin,1) avgCorrIn(1:end-1,i) rhv1(:,i) rhv2(:,i)];
    [bhat,vcov]=olswhite(lhv,rhv);
    betav(:,i) = bhat;
    sebetav(:,i) = sqrt(diag(vcov));
end;

tbetav_white = betav ./ sebetav
%}


% graph link between p and coefficient 
bp = repmat(betav(2,:),Tin,1) + repmat(betav(3,:),Tin,1) .* repmat(padsIn(2:end,:),1,nBins+1) ...
    + repmat(betav(3,:),Tin,1) .* repmat((padsIn(2:end,:).^2),1,nBins+1); 


figure('Position', [100, 100, 1049, 895]);
subplot(2,1,1)
plot(padsIn(2:end,:),bp(:,1),'.')
strtitle = ['coeffcient b as a function of ',tmp,' in regression rho_t = c + b*rho_{t-1} + e'];
title(strtitle)
legend('all stocks')
subplot(2,1,2)
plot(padsIn(2:end,:),bp(:,2:end),'.')
title('coeffcient b as a function of pabs for each beta bin')
if nBins == 2
    legend('Low beta','High beta')
else
    legend('Very low beta', 'Low beta', 'Medium beta', 'High beta')
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%







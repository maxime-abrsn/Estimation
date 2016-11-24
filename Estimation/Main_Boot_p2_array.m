
%{
%====================================================================
Main_ModelEstim.m

Author: Thuy Duong To
version: Apr 2016

Purpose: to estimate the model of 
        "Pairwise correlation dynamics and incomplete information"

NOTE: DESIGNED FOR CLUSTER ARRAY JOB (BOOTSTRAP)
    - run bootstrap multiple time, each with different seeds
    - read in the number of run and record different output file

Program:
- use 5 different datasets
- use GMM estimation:
 + 10 parameters: 
  discount rate: rho
  risk aversion: R (gamma in paper)
  growth rate high state and low state: theta H/L
  intensity, high to low and low to high: lambda 12/ 21
  aggregate cashflow volatility: sigmaA
  spread of mean reversion: alpha
  long term mean: beta
  indiosyncratic volatility: delta_i
    (note that aggregate cashflow X is normalized to 100?)
 + Constraints: only use correlation
  cross correlation, at different economic condition (recesstion/expansion),
     and at different beta value (low/high)   -> (4 constraints)
  + instrumental variables: ones, lagged mkt volatility, lagged mkt sharpe ratio
- Setting
  + rolling window option
  + use the average beta to copute the theoretical correlation / covar
- Bootstrap: Use Hall book for correcting the stderr of estimate
  + overlapping scheme
  + independent processes  



%==========================================================================
%}

clear all
close all
clc
tic;  % start measuring program running time

addpath(genpath(pwd));

%dbstop if error;     % only for debug


%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose sample (raw DAILY returns (DATA, TxN), dates(Tx1), NOT excess return)
% should start before 19630102
% should end on 20141231 (but it is working well is end date is earlier)

%%% OPTIONS %%%%
portfolio = 2; %(0) for 48FF, (1) for 100BtM (2) for beta decile, (3) for 10BtM, (4) for top 500 cap 20beta

startDate =  19710104; %19710104 (middle: 19931207) %19630102 (middle: 19900122); % cut the sample to interesting times
endEstDate = 20001231;  % end the estimation window

flag_rolling = 0;
rollw = 252;  % window size, 252 days = 1 year
rollstep = 21;  % rolling step, 21 days = 1 month

flag_origin_Pads=0;  % 0 for linear mapping ads to P, 1 for nonlinear mapping
flag_averageBeta = 1; % 1 = use average beta to compute theoretical correlation, 
                        % 0 for mapping of every stock, 2 for nearest stock to average beta
                      % NOTE: current version only use 1

flag_gmm = 0;   % run gmm esitmation (1), or not (0)
flag_graph = 1;  
flag_reg = 0;   % regression empirical correl on estimated correl
flag_save = 0;  % save gmm and bootstrap results


flag_bootstrap = 0;
nBoot = 100;



%==========================================================================================
% output file
diary off;
%prt1 = ['foutput',num2str(portfolio),'.txt'];
%diary(prt1);

%=========================================================================================
% only for Lunix cluster array job
run_i = str2num(getenv('PBS_ARRAYID'));  
    % very strange results if read in as string and keep as string
    % after convert, behave well.

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

if portfolio == 2 
    % position of interest to generate corr dynamic
    pos1 = 1;
    pos2 = 2;
    % number of bins
    nBins = 2;
    gap = 0;% gap between bins
    % range for Y
    rangeY = [-1 11];%[-2 8];

    % Preferences
    params0.R = 3; % Risk aversion, 4
    params0.rho=0.05; %0.05

    % Aggregate shock
    params0.thetaH=0.08; %0.08
    params0.thetaL=-0.02; %-0.02
    params0.deltaA=0.10; %0.1

    % Id shock
    params0.alpha=0.5; % 0.5
    params0.beta=5; % 1
    params0.deltaI=0.8; % 0.8
    
    % Markov chain
    params0.lambda_12 = 0.23;
    params0.lambda_21 = 1.15;
    
elseif portfolio == 3;
    % position of interest to generate corr dynamic
    pos1 = 1;
    pos2 = 2;
    % number of bins
    nBins = 2;
    gap = 0;% gap between bins
    % range for Y
    rangeY = [-1 11];%[-2 8];

    % Preferences
    params0.R = 3; % Risk aversion, 4
    params0.rho=0.05; %0.05

    % Aggregate shock
    params0.thetaH=0.08; %0.08
    params0.thetaL=-0.02; %-0.02
    params0.deltaA=0.10; %0.1

    % Id shock
    params0.alpha=0.7; % 0.5
    params0.beta=5; % 1
    params0.deltaI=0.8; % 0.8
        
    % Markov chain
    params0.lambda_12 = 0.23;
    params0.lambda_21 = 1.15;
    
elseif portfolio ==4;

    % position of interest to generate corr dynamic
    pos1 = 2;%2;
    pos2 = 3;%3;
    % number of bins
    nBins = 4;%4;
    gap = 0;%10;% gap between bins
    % range for Y
    rangeY = [-1 11];%[-2 8];

    % Preferences
    params0.R = 3; % Risk aversion, 4
    params0.rho=0.05; %0.05

    % Aggregate shock
    params0.thetaH=0.08; %0.08, INcrease difference to get more volatily results
    params0.thetaL=-0.02; %-0.02
    params0.deltaA=0.10; %0.1 % decrease to get more amplitude
    
    % Id shock
    params0.alpha=0.5; % 0.5 
    params0.beta=5; % 1
    params0.deltaI=0.8; % 0.8 
    
    % Markov chain
    params0.lambda_12 = 0.23;
    params0.lambda_21 = 1.15;
    
else 
    % position of interest to generate corr dynamic
    pos1 = 2;%2;
    pos2 = 3;%3;
    % number of bins
    nBins = 4;%4;
    gap = 0;%10;% gap between bins
    % range for Y
    rangeY = [-1 11];%[-2 8];

    % Preferences
    params0.R = 3; % Risk aversion, 4
    params0.rho=0.05; %0.05

    % Aggregate shock
    params0.thetaH=0.08; %0.08
    params0.thetaL=-0.02; %-0.02
    params0.deltaA=0.10; %0.1

    % Id shock
    params0.alpha=0.25; % 0.5
    params0.beta=5; % 1
    params0.deltaI=0.8; % 0.8
    
    % Markov chain
    params0.lambda_12 = 0.23;
    params0.lambda_21 = 1.15;
    
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
TT = size(DATA,1);
N = size(DATA,2);

% final data:
% eDATA, rf, erm, ADSindex, dates



%============================================================================
% Divide data into beta bins, and calculate avr pairwise correl within bins

[avgBeta, avgCorr, avgErm, mpr, volMkt, adsRep, trDates] = ... 
    fEmpData(eDATA, rf, erm, ADSindex, dates, nBins, start, stop, flag_rolling, rollw, rollstep);


%=================================================================================================
% GMM estimation
%=================================================================================================


gmmopt.W0 = 'Z';			% Initial weighting matrix, identity
%gmmopt.W='fweight';				% Subsequent wtg matrix optimal, inverse Spectral Density from gmmS
gmmopt.S= 'H';				% type of spectral density matrix: W,NW, G, H, AM, P
gmmopt.Strim = 2;              % 2 = demean Z'e
gmmopt.Slast = 2;           % 2 = recalculate S at final para and update W
gmmopt.prt = 1;				% Control printing, 0 = none, [1] = screen, else file
                                % there is an erro in 0 printing option.
gmmopt.plot = 0;            % [1], 0

gmmopt.infoz.hess = 'dfp';    % available method: dfp, bfgs, gn, marq, sd

prt1 = ['portfolio ',num2str(portfolio),'; rolling = ',num2str(flag_rolling),' (',num2str(rollw),'/',num2str(rollstep),')'];
prt2 = ['optim setup: ',num2str(gmmopt.W0),', ',num2str(gmmopt.S),', ',num2str(gmmopt.infoz.hess)];
disp(prt1);
disp(prt2);

gmmopt.vname = strvcat('R','rho','thetaH','thetaL','deltaA','alpha','beta','deltaI','lambda12','lambda21');
                % variable names
para = [params0.R; params0.rho; params0.thetaH; params0.thetaL; params0.deltaA; ...
        params0.alpha; params0.beta; params0.deltaI; params0.lambda_12; params0.lambda_21];
                % starting values
                
% observed, independent and instrument variables
X = avgErm(2:end);  % average expected market return
%Z = [ones(limit-1,1)   252*10*avgErm(1:end-1)    mpr(1:end-1)];  
    % for the market return: first annualize, then rescale so that Z is close to 1
    % 2 instruments: mkt return, Sharpe ratio
limit = size(avgErm,1);   
Z = [ones(limit-1,1)   volMkt(1:end-1)    mpr(1:end-1)];  
y = ones(limit-1,1);  % dummy variables

% other parameter needed for calculating moments
gmmopt.infoz.avgCorr = avgCorr(2:end,:);   % use the first lag for instrument
%gmmopt.infoz.avgVolsys = avgVolsys(2:end,:);
%gmmopt.infoz.avgVolidio = avgVolidio(2:end,:);
gmmopt.infoz.mpr = mpr(2:end,:);
gmmopt.infoz.volMkt = volMkt(2:end,:);
gmmopt.infoz.avgBeta = avgBeta(2:end,:);

gmmopt.infoz.rangeY = rangeY;
gmmopt.infoz.ADSindex = ADSindex;
gmmopt.infoz.trDates = trDates(2:end,:);
gmmopt.infoz.adsRep = adsRep(2:end,:);
gmmopt.infoz.pos = [pos1 pos2];
gmmopt.infoz.origin_Pads = flag_origin_Pads;
gmmopt.infoz.portfolio = portfolio;

gmmopt.infoz.flag_graph = 0;
gmmopt.infoz.flag_reg = 0;

% moment conditions
gmmopt.infoz.momt='fmoments_Corr';		% Moment Conditions

% --- Estimate the model with gmm() -------------------------------------
if flag_gmm == 1;
    gout=gmm_4Debug(para,gmmopt,y,X,Z);
else % if model already estimated (previous run), read output in
    path = ['./results/gmmout',num2str(portfolio),'Estim.mat'];
    load(path);
    %load('./results/gmmout2Estim.mat');
end;

% Checking the matching between empirical and  model values: 
    % graph and % regression
%gmmopt.infoz.flag_graph = flag_graph;
%gmmopt.infoz.flag_reg = flag_reg;


[mhat,ehat] = fmoments_Corr(gout.b,gmmopt.infoz,[],y,X,Z);
  % size mhat: 3 instruments x 4 basis moments
  % size ehat: limit (new T) x 4



%================================================================
% Moving block bootstrap

if flag_bootstrap ==1;
    
    disp('bootstrap gmm');
    
    rng('shuffle');  % seeds based on current time 
                    % used to run parallel bootstrap (ensure the bootstrap samples are different)
        
    % optimal bootsrap length
    l = opt_block_length_REV_dec07(erm);
    l = floor(mean(l));        
    
    % mT (centre moment for overlapping scheme)
    
    %k = 0;   % no dependence structure
    T = limit - 1; 
    w = nan(T,1);
    w(1:l-1) = (1:1:l-1)  / l;
    w(l:T-l+1) = 1;
    w(T-l+2 : T) = (T - (T-l+2:1:T) +1 ) / l;
    
    wf = repmat(w(1:T),1,size(ehat,2)) .* ehat(1:T,:);   
    mT = sum(wf) / (T-l+1);  % size = 1 x size(e,2)
    
    % gmm input for all bootstrap samples  
    
    gmmopt.prt = 1;  % no printing option (but still need to set infoz.prt to 1)
        % (to total no printing to work: need to change setting in all files. I have only changed
        %   setting in gmm.m file [written in gmm_4Debug.m])
    gmmopt.infoz.prt = 1;   % printing for optimization routine
        % for no printing, need to change step2 called by minz
    
    gmmopt.infoz.flag_graph = 0;
    gmmopt.infoz.flag_reg = 0;
    
    gmmopt.infoz.mT = mT;  % centre moment calculated at theta hat
    gmmopt.infoz.momthat = 'fmoments_Corr';
    gmmopt.infoz.thetahat = gout.b;
    gmmopt.infoz.Shat = gout.S;
    gmmopt.infoz.l = l;  % block length
    
    gmmopt.infoz.momt='fmoments_Corr_boot';		% Moment Conditions (recentered)
    
    % Moving block bootstrap (no dependence)
    
    boot_parai = nan(nBoot,size(para,1));
    boot_Ji = nan(nBoot,1);
      
    
    for bi = 1:nBoot;
        
        prt1 = ['bootstrap number ',num2str(bi)];
        disp(prt1);

        % ceil(T/m) Uniform random numbers over 1...T-l+1
        u = ceil((TT-l+1)*rand(ceil(TT/l),1));
        u = bsxfun(@plus,u,0:l-1)';
        % Transform to col vector, and remove excess
        u = u(:);
        u = u(1:TT);
        
        % Data sample simulation
        eDATA_bi = eDATA(u,:);
        rf_bi = rf(u);
        erm_bi = erm(u);
        ADSindex_bi = ADSindex(u,:);
        
        % Moment using bootstrap data
                
        [avgBeta_bi, avgCorr_bi, avgErm_bi, mpr_bi, volMkt_bi, adsRep_bi, trDates_bi] = ...
            fEmpData(eDATA_bi, rf_bi, erm_bi, ADSindex_bi, dates, nBins, start, stop, flag_rolling, rollw, rollstep);
        
        
        % observed, independent and instrument variables
        Xi = avgErm_bi(2:end);  % average expected market return
        %Z = [ones(limit-1,1)   252*10*avgErm(1:end-1)    mpr(1:end-1)];
        % for the market return: first annualize, then rescale so that Z is close to 1
        % 2 instruments: mkt return, Sharpe ratio
        Zi = [ones(limit-1,1)   volMkt_bi(1:end-1)    mpr_bi(1:end-1)];
        yi = ones(limit-1,1);  % dummy variables
        
        % other parameter needed for calculating moments
        gmmopt.infoz.avgCorr = avgCorr_bi(2:end,:);   % use the first lag for instrument
        %gmmopt.infoz.avgVolsys = avgVolsys(2:end,:);
        %gmmopt.infoz.avgVolidio = avgVolidio(2:end,:);
        gmmopt.infoz.mpr = mpr_bi(2:end,:);
        gmmopt.infoz.volMkt = volMkt_bi(2:end,:);
        gmmopt.infoz.avgBeta = avgBeta_bi(2:end,:);
                
        gmmopt.infoz.ADSindex = ADSindex_bi;
        gmmopt.infoz.trDates_bi = trDates(2:end,:);  % need to keep original one if sample chosen according to year
        gmmopt.infoz.adsRep = adsRep_bi(2:end,:);
        
        % all of the following not changed with bootstrap
        %gmmopt.infoz.rangeY = rangeY;  
        %gmmopt.infoz.pos = [pos1 pos2];
        %gmmopt.infoz.origin_Pads = flag_origin_Pads;
        %gmmopt.infoz.portfolio = portfolio;
        %gmmopt.infoz.flag_graph = 0;
        %gmmopt.infoz.prt = 0;   % printing for optimization routine        
        
        % special input for bootstrap
    
        %{
                
        % First step GMM estimator 
        %-----------------------------------------         
        
        gmmopt.infoz.momt='fmoments_Corr_boot';		% Moment Conditions (recentered)
        
        gout_bi=gmm_4Debug(para,gmmopt,yi,Xi,Zi);     % estimation
        
        
        % Second step GMM estimator
        %------------------------------------------
        
        % recalcualte centre moment at first step theta estimation
        [mhat,ehat] = fmoments_Corr_boot(gout_bi.b,gmmopt.infoz,[],y,X,Z);
        wf_bi = repmat(w(1:T),1,size(ehat,2)) .* ehat(1:T,:);   
        mT_bi = sum(wf_bi) / (T-l+1);
        %}        
  
        gout_bi=gmm_bootstrap(para,gmmopt,yi,Xi,Zi);     % estimation, correct J and tstat for bootstrap procedure
                    
        boot_parai(bi,:) = gout_bi.b';
        boot_Ji(bi) = gout_bi.Jbooti;
        boot_atstat(bi,:) = gout_bi.atstat';  % recentre the distribution
    end
    
       
    
    % the t-test
    tmp = (repmat(abs(gout.t'),nBoot,1) < boot_atstat);
    boot_p = sum(tmp) / nBoot; % proportion of distn on the right tail
    boot_p = boot_p';
    tmp = (repmat(gout.J,nBoot,1) < boot_Ji);
    boot_Jp = sum(tmp) / nBoot;
    
    % reprint the output of original GMM (because no priting option doesn't work for now)
    prt1 = ['portfolio ',num2str(portfolio),'; rolling = ',num2str(flag_rolling),' (',num2str(rollw),'/',num2str(rollstep),')'];
    prt2 = ['optim setup: ',num2str(gmmopt.W0),', ',num2str(gmmopt.S),', ',num2str(gmmopt.infoz.hess)];
    disp(prt1);
    disp(prt2);
    if isfield(gmmopt.infoz,'ftol')			% This checks if
        gout.ftol = gmmopt.infoz.ftol;		% moments = 0 for
    else						% Just-identified model
        gout.ftol = 1e-7;
    end
    gout.hprt = 0; gout.eprt = 0;
    if gmmopt.prt ~= 0
        gout.eprt = 1;
        prt_gmm(gout,gmmopt.vname,gmmopt.prt)
    end
    gout = rmfield(gout,'ftol');
    
    % print the output of bootstrap
    prt1 = ['Bootstrap p value (proportion of test statistics >= critical value), nBoot = ',num2str(nBoot)];
    disp(prt1);
    [gout.b boot_p]
    disp('p value for J-stat');
    boot_Jp

end

%===============================================================================
% final analysis
%===============================================================================

if flag_save == 1

% save estimation results
    
    gout.gmmoptW0 = gmmopt.W0;
    gout.gmmoptS = gmmopt.S;
    gout.gmmoptInfozHess = gmmopt.infoz.hess;  
        
    if flag_bootstrap == 1
        gout.nBoot = nBoot;
        gout.bootp = boot_p;
        gout.bootJp = boot_Jp;
        gout.boot_parai = boot_parai;
        gout.boot_Ji = boot_Ji;
        gout.boot_atstat = boot_atstat;
        gout.boot_optlen = l;
    end
        
    %path = ['./results/gmmout',num2str(portfolio),'Estim.mat'];  
    %path = ['./results/gmmout',num2str(portfolio),'_part1.mat'];  % temporary, if we run many instances of bootstrap
     
    % for Linux array job submit
    path = ['./results/gmmout',num2str(portfolio),'_part',num2str(run_i),'.mat'];
     
     save(path,'gout');
end

%=========================================================================
% Graph for empirical correlation 

if flag_graph == 1;
    
vym = round(trDates./100);  % take the year and month
vym = round(trDates./10000);  % take the year (choose either option)
Tg = length(vym);
ng = 5;   % need to choose an appropriate step length for the graphs (if rolling = 0 then ng=1 means each year)
            % X tick jumps by ng values 

figure('Position', [100, 100, 1049, 895]);
subplot(3,1,1)
plot(avgCorr(:,1))
title('average correlation')
legend('all stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
subplot(3,1,2)
plot([avgCorr(:,1+pos1) avgCorr(:,1+pos2)])
title('average correlation per bin')
    legend('Low beta','High beta')
%legend('low','mid','high')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
subplot(3,1,3)
plot([avgBeta(:,pos1) avgBeta(:,pos2)])
title('average beta per bin')
    legend('Low beta','High beta')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
      path1 = ['./Graphs/EmpCorrBeta_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/EmpCorrBeta_',num2str(portfolio),'.fig'];
      savefig(path3);
      saveas(gcf,path1);
    
%{
figure('Position', [100, 100, 1049, 395]);
subplot(1,2,1)
    plot(adsRep,avgCorr(:,1),'o')
    legend('All');
    title('Average correlation');
subplot(1,2,2)
    plot(adsRep,avgCorr(:,1+pos1),'*')
    hold on
    plot(adsRep,avgCorr(:,1+pos2),'o')
    legend('Low beta', 'High beta')

%if saveGraph == 1
      path1 = ['./Graphs/EmpCorrAds_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/EmpCorrAds_',num2str(portfolio),'.fig'];
      savefig(path3);
      saveas(gcf,path1);
      %saveTightFigure(gcf,path1);
%end
%}
      
      
%=======================================================================
% Graph of correlation against ads

%{
figure('Position', [100, 100, 1049, 895]);
subplot(3,2,1)
    plot(adsRep,avgVoltotAll,'o')
    legend('All');
    title('Total volatility')
subplot(3,2,2)
    plot(adsRep,avgVoltot(:,pos1),'*')
    hold on
    plot(adsRep,avgVoltot(:,pos2),'o')
    legend('Low beta', 'High beta')
subplot(3,2,3)
    plot(adsRep,avgVolsysAll,'o')
    legend('All');
    title('Systematic volatility')
subplot(3,2,4)
    plot(adsRep,avgVolsys(:,pos1),'*')
    hold on
    plot(adsRep,avgVolsys(:,pos2),'o')
    legend('Low beta', 'High beta')
subplot(3,2,5)
    plot(adsRep,avgVolidioAll,'o')
    legend('All');
    title('Idiosyncratic volatility')
subplot(3,2,6)
    plot(adsRep,avgVolidio(:,pos1),'*')
    hold on
    plot(adsRep,avgVolidio(:,pos2),'o')
    legend('Low beta', 'High beta')

%if saveGraph == 1
      path1 = ['./Graphs/EmpVolAds_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/EmpVolAds_',num2str(portfolio),'.fig'];
      savefig(path3);
      saveas(gcf,path1);
      %saveTightFigure(gcf,path1);
%end
%}



%========================================================================
% Check the matching between empirical and  model values
% graph and regression

gmmopt.infoz.flag_graph = flag_graph;
gmmopt.infoz.flag_reg = flag_reg;

% other parameter needed for calculating moments
gmmopt.infoz.avgCorr = avgCorr(2:end,:);   % use the first lag for instrument
gmmopt.infoz.mpr = mpr(2:end,:);
gmmopt.infoz.volMkt = volMkt(2:end,:);
gmmopt.infoz.avgBeta = avgBeta(2:end,:);
gmmopt.infoz.ADSindex = ADSindex;
gmmopt.infoz.trDates = trDates(2:end,:);
gmmopt.infoz.adsRep = adsRep(2:end,:);

[m,e] = fmoments_Corr(gout.b,gmmopt.infoz,[],y,X,Z);

end

%===========================================================
% program running time

disp('running time');
code_runtime=toc;
time_to_run={'hours','minutes','seconds';...
floor(code_runtime/3600),floor(mod(code_runtime,3600)/60),mod(code_runtime,60)}

diary off;

%==============================================================

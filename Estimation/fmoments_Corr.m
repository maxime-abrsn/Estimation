
%=======================================================================
%  moment functions for our estimation
%=======================================================================

function [m,e] = fmoments_Corr(para,infoz,stat,y,X,Z);

params.R = para(1); 
params.rho=para(2);  
params.thetaH=para(3);  
params.thetaL=para(4);  
params.deltaA=para(5); 
params.alpha=para(6);  
params.beta=para(7);  
params.deltaI=para(8);  
params.lambda_12=para(9);  
params.lambda_21=para(10); 
    
avgCorr = infoz.avgCorr;   % use the first lag for instrument
%avgVolsys = infoz.avgVolsys;
%avgVolidio = infoz.avgVolidio;
mpr = infoz.mpr;
volMkt = infoz.volMkt;
avgBeta = infoz.avgBeta;

rangeY = infoz.rangeY;
ADSindex = infoz.ADSindex;
trDates = infoz.trDates;
origin_Pads = infoz.origin_Pads;
tmp = infoz.pos;
pos1 = tmp(1);
pos2 = tmp(2);
clear tmp;

flag_graph = infoz.flag_graph;
flag_reg = infoz.flag_reg;
portfolio = infoz.portfolio;

% Calculating theoretical covar matrix

[B, Mbeta, P ,Y, crossCorrTh, volTotalTh, volsysTh, volidioTh, volMktTh] = gmm_theoric_covar(rangeY,params,flag_graph,portfolio);

%{
% Method 1: Use last ADS in the rolling window

% keep only ADS dates of interest
tmp = find(ismember(ADSindex(:,1),trDates));
if length(tmp) == length(trDates)
    %disp('last date of ADS matches last date of sample')   % this is the good case
else
    tmp = [tmp; length(ADSindex)];
    disp('ADS index length is shorter than sample')
end
    
tmp2 = ADSindex(tmp,2);
datesADS = ADSindex(tmp,1);
%}

%  Method 2: use the representative ADS of the rolling window
tmp2 = infoz.adsRep;

% infer p_t from ADS index

% Mapping may not be appropriate, not necessarily symmetric
% change that....

pHs = params.lambda_21 / (params.lambda_12 + params.lambda_21);
pLs = 1- pHs;

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
clear tmp tmp2

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
clear tmp tmp2
end

% use the Average Beta to calculate theoretical correlation / covar

    lowB = avgBeta(:,pos1);
    highB = avgBeta(:,pos2);
    tmpl = (lowB-min(lowB))/(max(highB)-min(lowB));  % create a relativity index from min(lowB) to max(highB)
    tmph = (highB-min(lowB))/(max(highB)-min(lowB));

    clear j

    % find nearest position of ADS in P
    C = abs(bsxfun(@minus,P',pads'));
    [~,idx] = min(C(:,1:size(C,2)));
    
    
    
    
    % find nearest position of Mbeta in P&Y
    for i=1:length(idx)
        clear C
        tmp = Mbeta(idx(i),:);   %  note Mbeta: row i = Pi, column j = Yj
        newBL = tmpl(i)*(max(tmp)-min(tmp))+min(tmp);  % transform to relativity index
        newBH = tmph(i)*(max(tmp)-min(tmp))+min(tmp);
                
        % transform real beta into theoretical
        if newBL <= min(tmp) && idx(i) < length(P)   
            idx2(i) = 1;
        elseif newBL >= max(tmp)
            idx2(i) = length(Y);
        else
            C = abs(bsxfun(@minus,tmp,newBL));
            [~,idx2(i)] = min(C(:,1:size(C,2)));    %  low beta position
        end
        if newBH <= min(tmp) && idx(i) < length(P)
            idx3(i) = 1;
        elseif newBH >= max(tmp)
            idx3(i) = length(Y);
        else
            C = abs(bsxfun(@minus,tmp,newBH));
            [~,idx3(i)] = min(C(:,1:size(C,2)));    % high beta position
        end
    end
    clear tmpl tmph tmp

    % covariance dynamic
    for i=1:length(idx)     % we already take out the first lag for instrumental variable when passing the variable in
        simCorr(i,1) = B(idx2(i),idx(i));   % note that for B: row i = Yi, column j = Pj
        simCorr(i,2) = B(idx3(i),idx(i));
        simVoltotal(i,1) = volTotalTh(idx(i),idx2(i));  % note that for vol matrix: row i = Pi, column j = Yj
        simVoltotal(i,2) = volTotalTh(idx(i),idx3(i));
        simVolsys(i,1) = volsysTh(idx(i),idx2(i));
        simVolsys(i,2) = volsysTh(idx(i),idx3(i));
        simVolidio(i,1) = volidioTh(idx(i),idx2(i));
        simVolidio(i,2) = volidioTh(idx(i),idx3(i));
        simVolMkt(i) = volMktTh(idx(i));  % not dependent on Y
    end
    simVolMkt = simVolMkt';   % transform to column vector
    
    % high and low states
    hidx = (pads>0.5);
    lidx = (pads<= 0.5);
    T = length(hidx);
    e1 = zeros(T,1);
    e2 = zeros(T,1);
    e3 = zeros(T,1);
    e4 = zeros(T,1);
    
    

    % moments matching
    
    e1(hidx) = simCorr(hidx,1) - avgCorr(hidx,1+pos1);  % correlation, low beta, high state
    e2(hidx) = simCorr(hidx,2) - avgCorr(hidx,1+pos2);  % correl, high beta, high state
    e3(lidx) = simCorr(lidx,1) - avgCorr(lidx,1+pos1);  % correlation, low beta, low state
    e4(lidx) = simCorr(lidx,2) - avgCorr(lidx,1+pos2);  % correl, high beta, low state
        
    e = [e1 e2 e3 e4];

    m = reshape(Z'*e/T,[],1);
    
      
    %    [1 2] - [1 2 3];
     
     
%=====================================================================================
% Final analysis (graphs)
%=====================================================================================

if(flag_graph == 1);
    
vym = round(trDates./100);  % take the year and month
vym = round(trDates./10000); % take the year only (choose either, depending on what we want)
Tg = length(vym);
ng = 5;   % need to choose an appropriate step length for the graphs (if rolling = 0 then ng=1 means each year)
            % X tick jumps by ng values


figure('Position', [100, 100, 1049, 895]);
subplot(3,1,1)
plot(pads)
%legend('ADS probability','Location','NorthWest')
title('ADS probability')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
subplot(3,1,2)
plot(simCorr(:,1))
hold on
plot(avgCorr(:,1+pos1),'r--')
hold off
title('Correlation for Low beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
%legend('simulated','empirical','Location','NorthWest')
subplot(3,1,3)
plot(simCorr(:,2))
hold on
plot(avgCorr(:,1+pos2),'r--')
hold off
title('Correlation for High beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
%legend('simulated','empirical','Location','NorthWest')
%if saveGraph == 1
      path1 = ['./Graphs/Corr_SimVsEmpi_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/Corr_SimVsEmpi_',num2str(portfolio),'.fig'];

      savefig(path3);  
      saveas(gcf,path1);
      %saveTightFigure(gcf,path3);  % doesn't seem to work with figures containing multiple subplot
%end

  % save result
     mCor.trDates = trDates;
     mCor.pads = pads;
     mCor.simCor = [simCorr(:,1) simCorr(:,2)];
     mCor.avgCor = [avgCorr(:,1+pos1) avgCorr(:,1+pos2)];
     path2 = ['./results/mCor',num2str(portfolio),'.mat'];
     save(path2,'mCor')

%{
figure('Position', [100, 100, 1049, 895]);
subplot(3,1,1)
plot(pads)
%legend('ADS probability','Location','NorthWest')
title('ADS probability')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
subplot(3,1,2)
plot(simVolsys(:,1))
hold on
plot(avgVolsys(:,pos1),'r--')
hold off
title('Systematic volatility, low beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
subplot(3,1,3)
plot(simVolsys(:,2))
hold on
plot(avgVolsys(:,pos2),'r--')
hold off
title('Systematic volatility, high beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
%if saveGraph == 1
      path1 = ['./Graphs/VolSys_SimVsEmpi_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/VolSys_SimVsEmpi_',num2str(portfolio),'.fig'];
      savefig(path3);
      saveas(gcf,path1);
      %saveTightFigure(gcf,path1);
%end


figure('Position', [100, 100, 1049, 895]);
subplot(3,1,1)
plot(pads)
%legend('ADS probability','Location','NorthWest')
title('ADS probability')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
subplot(3,1,2)
plot(simVolidio(:,1)); %,'m:')
hold on
plot(avgVolidio(:,pos1), 'k--');  %,'k-.')
hold off
title('Idiosyncratic volatility, Low beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
%legend('simulated','empirical','Location','NorthWest')
subplot(3,1,3)
plot(simVolidio(:,2)); %,'m:')
hold on
plot(avgVolidio(:,pos2), 'k--');  %,'k-.')
hold off
title('Idiosyncratic volatility, High beta stocks')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
%legend('simulated','empirical','Location','NorthWest')
%if saveGraph == 1
      path1 = ['./Graphs/VolIdio_SimVsEmpi_',num2str(portfolio),'.png'];
      path3 = ['./Graphs/VolIdio_SimVsEmpi_',num2str(portfolio),'.fig'];
      savefig(path3);
      saveas(gcf,path1);
      %saveTightFigure(gcf,path1);
%end

%}
end


%=====================================================================================
% Final analysis (regression)
%=====================================================================================

if(flag_reg == 1);
    
    % OLS regression to get R2  (ols(y,x)), low and high beta)
    reg.Low = ols(avgCorr(:,1+pos1),[ones(size(simCorr(:,1))) simCorr(:,1)]);
    reg.High = ols(avgCorr(:,1+pos2),[ones(size(simCorr(:,2))) simCorr(:,2)]);
    
    % only subset, low and high state (measured by P)
    
    reg.LowBLowP = ols(avgCorr(lidx,1+pos1), [ones(size(simCorr(lidx,1))) simCorr(lidx,1)]);
    reg.HighBLowP = ols(avgCorr(lidx,1+pos2), [ones(size(simCorr(lidx,1))) simCorr(lidx,2)]);
    reg.LowBHighP = ols(avgCorr(hidx,1+pos1), [ones(size(simCorr(hidx,1))) simCorr(hidx,1)]);
    reg.HighBHighP = ols(avgCorr(hidx,1+pos2), [ones(size(simCorr(hidx,1))) simCorr(hidx,2)]);
    
     % test H0: g1 != 1
     reg.Low.tdiff1 = (reg.Low.beta(end)-1)/reg.Low.bstd(end);
     reg.LowBLowP.tdiff1 = (reg.LowBLowP.beta(end)-1)/reg.LowBLowP.bstd(end);
     reg.LowBHighP.tdiff1 = (reg.LowBHighP.beta(end)-1)/reg.LowBHighP.bstd(end);
     reg.High.tdiff1 = (reg.High.beta(end)-1)/reg.High.bstd(end);
     reg.HighBLowP.tdiff1 = (reg.HighBLowP.beta(end)-1)/reg.HighBLowP.bstd(end);
     reg.HighBHighP.tdiff1 = (reg.HighBHighP.beta(end)-1)/reg.HighBHighP.bstd(end);
          
     % print result
     disp('Regression of model Corr on empirical corr');
     disp('R-square / coeff / tstat');
     tmp = [reg.Low.rsqr  reg.Low.beta(end)  reg.Low.tdiff1;
         reg.LowBLowP.rsqr  reg.LowBLowP.beta(end)  reg.LowBLowP.tdiff1;
         reg.LowBHighP.rsqr  reg.LowBHighP.beta(end)  reg.LowBHighP.tdiff1;
         reg.High.rsqr    reg.High.beta(end)   reg.High.tdiff1;
         reg.HighBLowP.rsqr    reg.HighBLowP.beta(end)   reg.HighBLowP.tdiff1;
         reg.HighBHighP.rsqr    reg.HighBHighP.beta(end)   reg.HighBHighP.tdiff1];
     disp([['Low B All P  '; 'Low B Low P  '; 'Low B HighP  '; 'HighB All P  '; 'HighB Low P  '; 'HighB HighP  '] num2str(tmp)]);        
     
     % save result
     path2 = ['./results/reg',num2str(portfolio),'.mat'];
     save(path2,'reg')
     
end
 
end

%=======================================================================
%  moment functions for bootstrap estimation
% recentered the momment
%=======================================================================

function [m,e] = fmoments_Corr_boot(para,infoz,stat,y,X,Z);

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

mT = infoz.mT;    % centered moment
   

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
    
    % recentre the  moment
    e = e - repmat(mT,T,1);    
    
    m = reshape(Z'*e/T,[],1);
    
      
     
end
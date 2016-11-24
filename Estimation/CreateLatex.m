

% project: pairwise correlation
% create table for latex and calculate few interesting values from estimation results

clear all;
close all;

portfolio = [4; 0; 1; 3];
nBoot = [10400; 11200; 11000; 11000];

%gout.vname = strvcat('R','rho','thetaH','thetaL','deltaA','alpha','beta','deltaI','lambda12','lambda21');
vname = {'$\gamma$'; '$\rho$'; '$\theta_H$'; '$\theta_L$'; ...
            '$\delta^A$'; '$\alpha$'; '$\eta$'; '$\delta^i$'; ...
            '$\lambda_{12}$'; '$\lambda_{21}$'};
    

nport = size(portfolio,1);
npar = size(vname,1);

b = nan(npar,nport);
pb =  nan(npar,nport);
J = nan(nport,1);
pJ = nan(nport,1);

for i = 1:nport;
    fname = ['./results/gmmout',num2str(portfolio(i)),'_',num2str(nBoot(i)),'bootstrap.mat'];
    %fname = ['gmmout',num2str(portfolio(i)),'_',num2str(nBoot(i)),'bootstrap.mat'];
    load(fname);
    b(:,i) = gout.b;
    pb(:,i) = gout.bootp_10k;
    J(i) = gout.J;
    pJ(i) = gout.bootJp_10k;
end

bs = cell(size(b)); % string for significant level
 for i=1:npar;
     for j=1:nport;
         if pb(i,j) <= 0.01;
             bs{i,j} = '***';
         elseif pb(i,j) <= 0.05;
             bs{i,j} = '**';
         elseif pb(i,j) <=0.1;
             bs{i,j} = '*';
         end
     end
 end
 Js = cell(size(J)); % string for significant level
 for i = 1:nport;
     if pJ(i) <= 0.01
         Js{i} = '***';
     elseif pJ(i) <= 0.05;
         Js{i} = '**';
     elseif pJ(i) <=0.1;
         Js{i} = '*';
     end
 end

disp('half life of idio CF');
log(2)./b(6,:)    % ln(2) / alpha

disp('model implied equity risk premium - mkt sharpe ratio');
b(1,:) .* b(5,:)   % risk aversion * systematic CF vol

disp('holding time in high state');
1./b(9,:)    % 1 / lambda12
12./b(9,:)   %  measured in months
mean(12./b(9,:))

disp('holding time in low state');
1./b(10,:)    % 1 / lambda21
12./b(10,:)
mean(12./b(10,:))



disp('======================================================');
disp(' parameter estimates ');
disp('======================================================');

fprintf('\n \n');
disp('---------------------------------------------');
disp('\begin{table}');
disp('\centering');
disp('\begin{tabular}{lcccc}');
disp('\hline');
disp('Parameter/ & Dataset 1 & Dataset 2 & Dataset 3 & Dataset 4 \\');
disp(' [p-value] & S\&P500 & 48industry & 100BTM & 10beta-BTM \\');
disp('\hline');
for i = 1: npar ;
    fprintf('%s',vname{i});
    for j = 1:nport;
        fprintf(' & %4.4f', b(i,j));
        fprintf('%s',bs{i,j});  % fprintf doesn't work with cell so have to do loop one by one
    end
    fprintf(' \\\\ \n');  % to place \\ at the end of the line, then start a new line
    fprintf('%s', '\vspace{0.25cm}');
    for j = 1:nport;
        fprintf(' & [%4.4f]',pb(i,j));
    end
    fprintf(' \\\\ \n');
end
disp('\hline');
fprintf('%s','J-stat');
for i = 1:nport
    fprintf(' & %4.2f', J(i));
end
fprintf(' \\\\ \n');
fprintf('%s', '\vspace{0.25cm}');
fprintf('%s','J-stat p-value');
for i = 1:nport
    fprintf(' & [%4.4f]', pJ(i));
end
fprintf(' \\\\ \n');
disp('\hline');
disp('\end{tabular}');
disp('\end{table}');



%=============================================================
% Regression result

b = nan(nport,4);
t = nan(nport,4);
df = nan(1,4);
rsqr = nan(nport,4);


for i = 1:nport;
    fname = ['./results/reg',num2str(portfolio(i)),'.mat'];
    load(fname);
    b(i,1) = reg.LowBLowP.beta(2);
    b(i,2) = reg.LowBHighP.beta(2);
    b(i,3) = reg.HighBLowP.beta(2);
    b(i,4) = reg.HighBHighP.beta(2);
    t(i,1) = reg.LowBLowP.tdiff1;
    t(i,2) = reg.LowBHighP.tdiff1;
    t(i,3) = reg.HighBLowP.tdiff1;
    t(i,4) = reg.HighBHighP.tdiff1;
    rsqr(i,1) = reg.LowBLowP.rsqr;
    rsqr(i,2) = reg.LowBHighP.rsqr;
    rsqr(i,3) = reg.HighBLowP.rsqr;
    rsqr(i,4) = reg.HighBHighP.rsqr;
    df(i,1) = reg.LowBLowP.nobs - 2;
    df(i,2) = reg.LowBHighP.nobs - 2;
    df(i,3) = reg.HighBLowP.nobs - 2;
    df(i,4) = reg.HighBHighP.nobs - 2;
end
pt = tcdf(-abs(t),df)*2;

ts = cell(size(t)); % string for significant level
 for i=1:nport;
     for j=1:4; 
         if pt(i,j) <= 0.01;
             ts{i,j} = '***';
         elseif pt(i,j) <= 0.05;
             ts{i,j} = '**';
         elseif pt(i,j) <=0.1;
             ts{i,j} = '*';
         end
     end
 end

 pname = {'Dataset 1 (S\&P500)'
        'Dataset 2 (48industry)'
        'Dataset 3 (100BTM)'
        'Dataset 4 (100beta-BTM)'};
    
    
disp('======================================================');
disp(' regression result ');
disp('======================================================');

fprintf('\n \n');
disp('\begin{tabular}{lcccc}');
disp('\toprule');
disp(' & Low beta, & Low beta,  & High beta, & High beta, \\');
disp(' & Low P & High P & Low P & High P \\');
disp('\hline');

for i = 1: nport ;
    fprintf('%s',pname{i});
    for j = 1:4;
        fprintf(' & %4.4f', b(i,j));
        fprintf('%s',ts{i,j});  % fprintf doesn't work with cell so have to do loop one by one
    end
    fprintf(' \\\\ \n');  % to place \\ at the end of the line, then start a new line
    fprintf('%s', '\vspace{0.25cm}');
    for j = 1:4;
        fprintf(' & (%4.4f)',t(i,j));
    end
    fprintf(' \\\\ \n');
end
disp('\bottomrule');
disp('\end{tabular}');
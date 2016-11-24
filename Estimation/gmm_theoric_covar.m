function [B, Mbeta, P ,Y, crossCorr1, Sigmatot, Sigmaag, Sigmaid, SigmaMkt] = gmm_theoric_covar(rangeY,params,finalFlag,portfolio)



x = 200; % cut sample due to arbitrary starting values


% power on deltaI*Y^
puissance = 0;

% Markov chain
lambda_12 = params.lambda_12;
lambda_21 = params.lambda_21;
staprobaG=lambda_21/(lambda_12+lambda_21);

% Preferences
R = params.R; % Risk aversion, 4
gamma=1-R; % 
rho=params.rho; %0.05

% Aggregate shock
thetaH=params.thetaH; %0.08
thetaL=params.thetaL; %-0.02
deltaA=params.deltaA; %0.1
X(1)=100;

% Id shock
alpha=params.alpha; % 0.5
beta=params.beta; % 1
deltaI=params.deltaI; % 0.8

% density of pt
% simulation parameters
N = 1000; % number of trajectories   
T = 30; % upper limit
h = 1e-01; % step size
M = T/h; % number of steps
    

%risk adjusted process
thH=gamma*thetaH+0.5*gamma*(gamma-1)*deltaA^2;
thL=gamma*thetaL+0.5*gamma*(gamma-1)*deltaA^2;

%constants
a = (rho+alpha+lambda_12+lambda_21-thH ) / ( (rho+alpha)*(rho+alpha+lambda_12+lambda_21) - thH*(rho+alpha+lambda_21-thL)-(rho+alpha+lambda_12)*thL );
b = (thH-thL ) / ( (rho+alpha)*(rho+alpha+lambda_12+lambda_21) - thH*(rho+alpha+lambda_21-thL)-(rho+alpha+lambda_12)*thL );
abar = (rho+lambda_12+lambda_21-thH ) / ( (rho)*(rho+lambda_12+lambda_21) - thH*(rho+lambda_21-thL)-(rho+lambda_12)*thL );
bbar = (thH-thL ) / ( (rho)*(rho+lambda_12+lambda_21) - thH*(rho+lambda_21-thL)-(rho+lambda_12)*thL );

%loop on Y_i and p_t
ymin=rangeY(1);
ymax=rangeY(2);


Y=linspace(ymin,ymax,100);
P=linspace(0,1,100);
% 
%  clear Y
%  Y = beta*ones(1,100);

for i=1:length(P)
    for j=1:length(Y)
        Sigmaag20(i,j) = deltaA + (((Y(j)-beta)*b+beta*bbar)*P(i)*(1-P(i))*(thetaH-thetaL)/deltaA)/((Y(j)-beta)*(a+b*P(i))+beta*(abar+bbar*P(i)));
    end
end
for i=1:length(P)
    Sigmaid(i,:)    = ( Y + beta .* ( (abar+bbar*P(i))- (a+b*P(i)) )./(a+b*P(i)) ).^(-1)*deltaI.*(Y.^puissance);
    Sigmaag(i,:)    = deltaA + ( (Y-beta)*b + beta*bbar )./( (Y-beta).*(a+b*P(i)) + beta.*(abar+bbar*P(i)) ) .* P(i)*(1-P(i)) .* (thetaH-thetaL)/deltaA;
    
    if i == 25
        j = 25;
        Sigmaag2    = deltaA + ( (Y-beta)*b + beta*bbar )./( (Y-beta).*(a+b*P(i)) + beta.*(abar+bbar*P(i)) ) .* P(i)*(1-P(i)) .* (thetaH-thetaL)/deltaA;
        sigIa       = deltaA + ( (Y-beta)*b + beta*bbar )./( (Y-beta).*(a+b*P(j)) + beta.*(abar+bbar*P(j)) ) .* P(j)*(1-P(j)) .* (thetaH-thetaL)/deltaA;
    end
    Sigmatot(i,:)   = sqrt(Sigmaid(i,:).^2 + Sigmaag(i,:).^2);
    
    SigmaMkt(i) = deltaA + ( beta*bbar )./( beta.*(abar+bbar*P(i)) ) .* P(i)*(1-P(i)) .* (thetaH-thetaL)/deltaA ;
    Mbeta(i,:)      = Sigmaag(i,:) ./ SigmaMkt(i);
end

% NOTE: Struture for sigma and beta 
% row i = probability i
% column j = Y_j


%{
figure('Position', [1000, 1000, 800, 600]);
mesh(Y,P,Mbeta)
title('Beta')
xlabel('Id. Shock')
ylabel('Proba.')
path1 = ['./Graphs/Theoretical_betas.png'];
path3 = ['./Graphs/Theoretical_betas.fig'];
saveTightFigure(gcf,path3);
saveas(gcf,path1);
%}


%pairwise correlations

%vary Y(s) hold P in [0.05, 0.5, 0.95]
%low Y = low Beta, high PER

for i=1:length(P) % all pt
    for j=1:length(Y) % All Yi, Yj
        crossCorr1(:,j,i) = (Sigmaag(i,j).*Sigmaag(i,:))./(Sigmatot(i,j).*Sigmatot(i,:));
    end

end
    
% NOTE: struture for crossCorr1:
% each page (last index) is for each probability, ie, (:,:,i) is for p_i
% in each page, crossCorr between different Y (rows Y_i, column Y_j)
        

% Effect of probability on cross-correlation for a given idiosyncratic shock
% low corresponds to beta (reverse for PER)


[Ma,~,Na]=size(crossCorr1);%# A is your matrix
indx=cumsum([1:(Ma+1):Ma^2; Ma^2.*ones(Na-1,Ma)]); %#diagonal indices
B=crossCorr1(indx');

% NOTE: structure for B
% each column is for one probability, ie. Column i for p_i
% each row is for Y, ie. correl between Y_i and Y_i


%==========================================================================
% GRAPH

if(finalFlag ==1);
       
    
    %=================================================  
    % The effect of Y and P on various measures
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,Mbeta)
    zlabel('Beta')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_betas_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_betas_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,Sigmatot)
    zlabel('Total volatility')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_SigmaTot_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_SigmaTot_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,Sigmaid)
    zlabel('Idiosyncratic volatility')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_SigmaId_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_SigmaId_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,Sigmaag)
    zlabel('Systematic volatility')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_SigmaAg_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_SigmaAg_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);

    
    % effect of Y and P on PER and Stock price
    
    for i=1:length(P)
        PER(i,:)        = ( (Y-beta).*(a+b*P(i)) + beta.*(abar+bbar*P(i)) ) ./ Y;
        StockPrice(i,:) = X.*( (Y-beta).*(a+b*P(i)) + beta.*(abar+bbar*P(i)) );
    end
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,PER)
    zlabel('Price-earning ratio')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_PER_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_PER_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure('Position', [1000, 1000, 800, 600]);
    mesh(Y,P,StockPrice)
    zlabel('Price')
    xlabel('Id. Shock')
    ylabel('Probability of high state')
    path1 = ['./Graphs/Theory_YP_Price_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_YP_Price_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    
    %=================================================================
    % Correlation - dependency on Y, P, beta
    
    figure('Position', [1000, 0, 800, 600]);
    mesh(Y,P,B')
    %title('Relation between Beta and Correlation for all states of the economy')
    set(gca,'XTick',[Y(1),Y(50),Y(100)])
    set(gca,'XTickLabel',{'Low/low beta','Medium/medium beta','High/high beta'})
    xlabel('Id. shock')
    ylabel('Probability of high state')
    zlabel('Cross-correlation')
    path1 = ['./Graphs/Theory_Corr_BetaP_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_Corr_BetaP_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure('Position', [1000, 0, 800, 600]);
    j = [1;50;100];
    for i = 1:3
        plot(P,B(j(i),:));
        hold on
    end
    hold off
    legend('Low beta/ Low beta','Medium beta/ Medium beta', 'High beta/ High beta');
    xlabel('Probability');
    ylabel('Cross-correlation');
    path1 = ['./Graphs/Theory_Corr_P',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_Corr_P',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);    
    
    % Effect of probability on cross-correlation for a given idiosyncratic shock
    % low corresponds to beta (reverse for PER)
    
    lowlowcorr = reshape(crossCorr1(1,1,:),100,1);
    midmidcorr = reshape(crossCorr1(50,50,:),100,1);
    highhighcorr = reshape(crossCorr1(100,100,:),100,1);
    
    figure
    plot (P,lowlowcorr,P,midmidcorr,P,highhighcorr) 
    legend('Low beta/ Low beta','Medium beta/ Medium beta', 'High beta/ High beta');
    xlabel('Probability');
    ylabel('Cross-correlation');
    path1 = ['./Graphs/Theory_Corr_P_v2_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_Corr_P_v2_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);    
    
    % link between Beta and P/E for different values of P when Y moves 

    figure
    plot(Mbeta(5,1:end-10),PER(5,1:end-10))
    title('P=0.05')
    xlabel('Beta')
    ylabel('PER')
    path1 = ['./Graphs/Theory_BetaP_lowP_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_BetaP_lowP_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure
    plot(Mbeta(50,1:end-10),PER(50,1:end-10))
    title('P=0.50')
    xlabel('Beta')
    ylabel('PER')
    path1 = ['./Graphs/Theory_BetaP_medP_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_BetaP_medP_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    figure
    plot(Mbeta(95,1:end-10),PER(95,1:end-10))
    title('P=0.95')
    xlabel('Beta')
    ylabel('PER')
    path1 = ['./Graphs/Theory_BetaP_highP_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_BetaP_highP_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);

    % Correlation dynamic
    %===============================
    
    % volatility from pt
    fp = zeros(100,100,100);
    fpdf = zeros(100,100,100);
    fy = zeros(100,100,100);
    fydy = zeros(100,100,100);
    
    fpdpVerif3 = (diff(highhighcorr,1)./0.0101).*P(2:end)'.*(1-P(2:end)').*(thetaH-thetaL)./deltaA;
    fpdpVerif1 = (diff(lowlowcorr,1)./0.0101).*P(2:end)'.*(1-P(2:end)').*(thetaH-thetaL)./deltaA;
    fpdpVerif2 = (diff(midmidcorr,1)./0.0101).*P(2:end)'.*(1-P(2:end)').*(thetaH-thetaL)./deltaA;
    
    figure
    title('fp dynamic')
    plot(P(2:end)',fpdpVerif1,P(2:end)',fpdpVerif2,P(2:end)',fpdpVerif3)
    legend('Low Beta / Low Beta','Mid Beta / Mid Beta','High Beta / High Beta')
    path1 = ['./Graphs/Theory_fp_',num2str(portfolio),'.png'];
    path3 = ['./Graphs/Theory_fp_',num2str(portfolio),'.fig'];
    saveTightFigure(gcf,path3);
    saveas(gcf,path1);
    
    % 2 methods should give the same results, but not the case. dY? drift?
    
    j = 50;
    
    test2 = diff(Sigmaag20(j,:)./(1/length(Y)));
    figure
    title('fsigA dynamic')
    plot(Y(2:end)',test2)
    
    test3 = 1.*((abar*b-a*bbar)*beta.*P(j).*(1-P(j)).*(thetaH-thetaL)./deltaA)./((Y.*(a+b*P(j)+beta*(abar-a+P(j)*(bbar-b)))).^2);
    figure
    title('fsigA dynamic derive')
    plot(Y',test3)

    
end

%[1 2] - [1 2 3]

end


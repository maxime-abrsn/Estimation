

% create graph

clear all;
close all;

portfolio = [4; 0; 1; 3];
nport = size(portfolio,1);

pname = {'Dataset 1 (S&P500): low beta'
        'Dataset 1 (S&P500): high beta'
        'Dataset 2 (48industry): low beta'
        'Dataset 2 (48industry): high beta'
        'Dataset 3 (100BTM): low beta'
        'Dataset 3 (100BTM): high beta'
        'Dataset 4 (100beta-BTM): low beta'
        'Dataset 4 (100beta-BTM): high beta'};
    
load ./results/mCor0.mat;
Tg = size(mCor.simCor,1);
simCorL = nan(Tg,nport);
simCorH  = nan(Tg,nport);
avgCorL = nan(Tg,nport);
avgCorH = nan(Tg,nport);

for i = 1:nport;
    clear mCor;
    fname = ['./results/mCor',num2str(portfolio(i)),'.mat'];
    load(fname);
    Tg = size(mCor.simCor,1);  % has to do this because portfolio 4 has 1 less obs than other portfolios
    simCorL(1:Tg,i) = mCor.simCor(:,1);
    simCorH(1:Tg,i) = mCor.simCor(:,2);
    avgCorL(1:Tg,i) = mCor.avgCor(:,1);
    avgCorH(1:Tg,i) = mCor.avgCor(:,2);
end


vym = round(mCor.trDates./100);  % take the year and month
vym = round(mCor.trDates./10000); % take the year only (choose either, depending on what we want)
Tg = length(vym);
ng = 5;   % need to choose an appropriate step length for the graphs (if rolling = 0 then ng=1 means each year)
            % X tick jumps by ng values


figure('Position', [100, 100, 1049, 895]);
subplot(5,2,1)
plot(mCor.pads)
%legend('ADS probability','Location','NorthWest')
title('ADS probability')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
    
subplot(5,2,2)
plot(mCor.pads)
%legend('ADS probability','Location','NorthWest')
title('ADS probability')
    set(gca,'xlim',[0 Tg]);
    set(gca,'XTick',[1:ng:Tg]);
    set(gca,'XTickLabel',[vym(1:ng:end)]);
    set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
    set(gca,'XTick',[-2:ng:Tg+2]);  
    set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);

for i = 1:nport;
    subplot(5,2,2*i+1)
        plot(simCorL(:,i))
        hold on
        plot(avgCorL(:,i),'r--')
        hold off
        title(pname(2*i-1))
        set(gca,'xlim',[0 Tg]);
        set(gca,'XTick',[1:ng:Tg]);
        set(gca,'XTickLabel',[vym(1:ng:end)]);
        set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
        set(gca,'XTick',[-2:ng:Tg+2]);
        set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
        %legend('simulated','empirical','Location','NorthWest')
        
    subplot(5,2,2*i+2)
        plot(simCorH(:,i))
        hold on
        plot(avgCorH(:,i),'r--')
        hold off
        title(pname(2*i))
        set(gca,'xlim',[0 Tg]);
        set(gca,'XTick',[1:ng:Tg]);
        set(gca,'XTickLabel',[vym(1:ng:end)]);
        set(gca,'xlim',[-2 Tg+2]); % temporary set up, for no rolling case
        set(gca,'XTick',[-2:ng:Tg+2]);
        set(gca,'XTickLabel',[vym(1:ng:end)-2 ; 2015]);
        
end

%legend('simulated','empirical','Location','NorthWest')
%if saveGraph == 1
      path1 = ['./Graphs/Corr_SimVsEmpi_Paper.png'];
      path3 = ['./Graphs/Corr_SimVsEmpi_Paper.fig'];
      saveTightFigure(gcf,path3);
      saveas(gcf,path1);
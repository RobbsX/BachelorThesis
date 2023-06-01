% Make histograms of layer 1 and 2 using fitrnet

clearvars;
clc;
set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultLineLineWidth', 2);

%% Choose to use Model or not
useModel = 1;


%% Files & User
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Bachelorarbeit"));

%% import summarycell
load('MACDsummarycell.mat') % MACDsummarycell
load('D7summarycell.mat') % D7summarycell
load('BB1summarycell.mat') % BB1summarycell
load('BB2summarycell.mat') % BB2summarycell
load('FitrnetStdsummarycell.mat')
load('FitrnetTSLTPsummarycell.mat')
load('FitnetStdsummarycell.mat')
load('FitnetTSLTPsummarycell.mat')


% Name=1, PnL=2, avgReturnPA=4, PnL BnH=5, avgRetPA BnH=7,
% MaxDD=8,TradesNum=9, Exposure=11, RoMaD=12, SR=13,slopeLR=14

%% Delete not used shares, PnL = 0 
MACDsummarycell = MACDsummarycell(cellfun(@(a)~isempty(a)&&a~=0,MACDsummarycell(:,2)),:);
D7summarycell = D7summarycell(cellfun(@(a)~isempty(a)&&a~=0,D7summarycell(:,2)),:);
BB1summarycell = BB1summarycell(cellfun(@(a)~isempty(a)&&a~=0,BB1summarycell(:,2)),:);
BB2summarycell = BB2summarycell(cellfun(@(a)~isempty(a)&&a~=0,BB2summarycell(:,2)),:);
FitrnetStdsummarycell = FitrnetStdsummarycell(cellfun(@(a)~isempty(a)&&a~=0,FitrnetStdsummarycell(:,2)),:);
FitrnetTSLTPsummarycell = FitrnetTSLTPsummarycell(cellfun(@(a)~isempty(a)&&a~=0,FitrnetTSLTPsummarycell(:,2)),:);
FitnetStdsummarycell = FitnetStdsummarycell(cellfun(@(a)~isempty(a)&&a~=0,FitnetStdsummarycell(:,2)),:);
FitnetTSLTPsummarycell = FitnetTSLTPsummarycell(cellfun(@(a)~isempty(a)&&a~=0,FitnetTSLTPsummarycell(:,2)),:);

%% convert to double for boxplot
MACDsummarydouble = cell2mat(MACDsummarycell(:,2:14));
D7summarydouble = cell2mat(D7summarycell(:,2:14));
BB1summarydouble = cell2mat(BB1summarycell(:,2:14));
BB2summarydouble = cell2mat(BB2summarycell(:,2:14));
FitrnetStdsummarydouble = cell2mat(FitrnetStdsummarycell(:,2:14));
FitrnetTSLTPsummarydouble = cell2mat(FitrnetTSLTPsummarycell(:,2:14));
FitnetStdsummarydouble = cell2mat(FitnetStdsummarycell(:,2:14));
FitnetTSLTPsummarydouble = cell2mat(FitnetTSLTPsummarycell(:,2:14));

%% Strategy beats benchmark %
BMACD = length(MACDsummarydouble(MACDsummarydouble(:,3)>MACDsummarydouble(:,6))==1)/size(MACDsummarydouble,1)*100; %#ok<*SZARLOG>
BD7 = length(D7summarydouble(D7summarydouble(:,3)>D7summarydouble(:,6))==1)/size(D7summarydouble,1)*100; 
BBB1 = length(BB1summarydouble(BB1summarydouble(:,3)>BB1summarydouble(:,6))==1)/size(BB1summarydouble,1)*100; 
BBB2 = length(BB2summarydouble(BB2summarydouble(:,3)>BB2summarydouble(:,6))==1)/size(BB2summarydouble,1)*100; 
BFitrStd = length(FitrnetStdsummarydouble(FitrnetStdsummarydouble(:,3)>FitrnetStdsummarydouble(:,6))==1)/size(FitrnetStdsummarydouble,1)*100;
BFitrTSLTP = length(FitrnetTSLTPsummarydouble(FitrnetTSLTPsummarydouble(:,3)>FitrnetTSLTPsummarydouble(:,6))==1)/size(FitrnetTSLTPsummarydouble,1)*100;
BFitStd = length(FitnetStdsummarydouble(FitnetStdsummarydouble(:,3)>FitnetStdsummarydouble(:,6))==1)/size(FitnetStdsummarydouble,1)*100; 
BFitTSLTP = length(FitnetTSLTPsummarydouble(FitnetTSLTPsummarydouble(:,3)>FitnetTSLTPsummarydouble(:,6))==1)/size(FitnetTSLTPsummarydouble,1)*100; 


%% Create table of summarycell
MACDsummarytable = cell2table(MACDsummarycell);
MACDsummarytable.Properties.VariableNames = {'Share','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 
D7summarytable = cell2table(D7summarycell);
D7summarytable.Properties.VariableNames = {'Share','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 
BB1summarytable = cell2table(BB1summarycell);
BB1summarytable.Properties.VariableNames = {'Share','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 
BB2summarytable = cell2table(BB2summarycell);
BB2summarytable.Properties.VariableNames = {'Share','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 

%% Make summarylines
% MACD
MACDsummaryline = cell(1, 14); MACDsummaryline{1,1} = 'MACD';
for ss = 2:14
    MACDsummaryline{1,ss} = median(cell2mat(MACDsummarycell(:,ss)));
end
% Double7
D7summaryline = cell(1, 14); D7summaryline{1,1} = 'Double7';
for ss = 2:14
    D7summaryline{1,ss} = median(cell2mat(D7summarycell(:,ss)));
end
% BB1
BB1summaryline = cell(1, 14); BB1summaryline{1,1} = 'BB1';
for ss = 2:14
    BB1summaryline{1,ss} = median(cell2mat(BB1summarycell(:,ss)));
end
% BB2
BB2summaryline = cell(1, 14); BB2summaryline{1,1} = 'BB2';
for ss = 2:14
    BB2summaryline{1,ss} = median(cell2mat(BB2summarycell(:,ss)));
end
% FitrnetStd
FitrnetStdsummaryline = cell(1, 14); FitrnetStdsummaryline{1,1} = 'FitrnetStd';
for ss = 2:14
    FitrnetStdsummaryline{1,ss} = median(cell2mat(FitrnetStdsummarycell(:,ss)));
end
% FitrnetTSLTP
FitrnetTSLTPsummaryline = cell(1, 14); FitrnetTSLTPsummaryline{1,1} = 'FitrnetTSLTP';
for ss = 2:14
    FitrnetTSLTPsummaryline{1,ss} = median(cell2mat(FitrnetTSLTPsummarycell(:,ss)));
end
% FitnetStd
FitnetStdsummaryline = cell(1, 14); FitnetStdsummaryline{1,1} = 'FitnetStd';
for ss = 2:14
    FitnetStdsummaryline{1,ss} = median(cell2mat(FitnetStdsummarycell(:,ss)));
end
% FitnetTSLTP
FitnetTSLTPsummaryline = cell(1, 14); FitnetTSLTPsummaryline{1,1} = 'FitnetTSLTP';
for ss = 2:14
    FitnetTSLTPsummaryline{1,ss} = median(cell2mat(FitnetTSLTPsummarycell(:,ss)));
end


%% Create table of strategies summarylines
strategycomparisonStrategies = [MACDsummaryline;D7summaryline;BB1summaryline;BB2summaryline];
strategycomparisonStrategies = cell2table(strategycomparisonStrategies);
strategycomparisonStrategies.Properties.VariableNames = {'Strategy','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 

%% Create table of NN summarylines
%strategycomparisonNN = [FitrnetStdsummaryline;FitrnetTSLTPsummaryline;FitnetStdsummaryline;FitnetTSLTPsummaryline];
strategycomparisonNN = [FitrnetStdsummaryline;FitrnetTSLTPsummaryline];
strategycomparisonNN = cell2table(strategycomparisonNN);
strategycomparisonNN.Properties.VariableNames = {'Strategy','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 

%% Create table of all summarylines
strategycomparison = [MACDsummaryline;D7summaryline;BB1summaryline;BB2summaryline;...
    FitrnetStdsummaryline;FitrnetTSLTPsummaryline;FitnetStdsummaryline;FitnetTSLTPsummaryline];
strategycomparison = cell2table(strategycomparison);
strategycomparison.Properties.VariableNames = {'Strategy','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 


% writetable(strategycomparisonStrategies, 'comparisonStrategies.xlsx')
% writetable(strategycomparisonNN, 'comparisonNN.xlsx')
% writetable(strategycomparison, 'comparison.xlsx')


%% Cut to smallest summarydouble
cutStrategies = min([size(MACDsummarydouble(:,1),1) size(D7summarydouble(:,1),1)...
    size(BB1summarydouble(:,1),1) size(BB2summarydouble(:,1),1)]);
% cutNN = min([size(FitrnetStdsummarydouble(:,1),1) size(FitrnetTSLTPsummarydouble(:,1),1)...
%     size(FitnetStdsummarydouble(:,1),1) size(FitnetTSLTPsummarydouble(:,1),1)]);
cutNN = min([size(FitrnetStdsummarydouble(:,1),1) size(FitrnetTSLTPsummarydouble(:,1),1)]);
cutAll = min([size(MACDsummarydouble(:,1),1) size(D7summarydouble(:,1),1)...
    size(BB1summarydouble(:,1),1) size(BB2summarydouble(:,1),1)...
    size(FitrnetStdsummarydouble(:,1),1) size(FitnetTSLTPsummarydouble(:,1),1) ]);

%% BoxPlot avg Return p.a. for Strategies
avgRetPAStrategies = [MACDsummarydouble(1:cutStrategies,3),D7summarydouble(1:cutStrategies,3),...
    BB1summarydouble(1:cutStrategies,3),BB2summarydouble(1:cutStrategies,3) BB2summarydouble(1:cutStrategies,6)];
colorsStrategies = {'#2D42C2','#405A00','#932055','#33486B','#74356E'};

% MACD=#2D42C2 D7=#405A00 BB1=#932055 BB2=#33486B benchmark=#74356E
figure;
ax1 = axes();
hold(ax1);
for i=1:size(avgRetPAStrategies,2)
    boxchart(i*ones(size(avgRetPAStrategies(:,i))),avgRetPAStrategies(:,i),'BoxFaceColor',colorsStrategies{i})
end
yline(0,'-k')
xlabel('Technical Strategy')
ylabel('Average Return per anuum')
title('Average Return p.a. of Technical Strategies')
legend('MACD','Double7','Bollinger Band 1','Bollinger Band 2','Buy and Hold')
hold off;


%% Plot avg Return p.a. for NN
% avgRetPANN = [FitrnetStdsummarydouble(1:cutNN,3),FitrnetTSLTPsummarydouble(1:cutNN,3),...
%     FitnetStdsummarydouble(1:cutNN,3),FitnetTSLTPsummarydouble(1:cutNN,3) FitrnetTSLTPsummarydouble(1:cutNN,6)];
% colorsNN = {'#405A00','#10930E','#2D42C2','#0068FF','#74356E'};
avgRetPANN = [FitrnetStdsummarydouble(1:cutNN,3),FitrnetTSLTPsummarydouble(1:cutNN,3),...
     FitrnetTSLTPsummarydouble(1:cutNN,6)];
colorsNN = {'#405A00','#10930E','#74356E'};

% MACD=#2D42C2 D7=#405A00 BB1=#932055 BB2=#33486B benchmark=#74356E
figure;
ax2 = axes();
hold(ax2);
for i=1:size(avgRetPANN,2)
    boxchart(i*ones(size(avgRetPANN(:,i))),avgRetPANN(:,i),'BoxFaceColor',colorsNN{i})
end
yline(0,'-k')
xlabel('Neural Network Strategy')
ylabel('Average Return per anuum')
title('Average Return p.a. of Neural Network Strategies')
%legend('Fitrnet Std','Fitrnet TSLTP','Fitnet Std','Fitnet TSLTP','Buy and Hold')
legend('Fitrnet Std','Fitrnet TSLTP','Buy and Hold')
hold off


%% Plot avg Return p.a. for Strategies and NN
avgRetPA = [MACDsummarydouble(1:cutAll,3),D7summarydouble(1:cutAll,3),...
    BB1summarydouble(1:cutAll,3),BB2summarydouble(1:cutAll,3) ...
    FitrnetStdsummarydouble(1:cutAll,3),FitrnetTSLTPsummarydouble(1:cutAll,3),...
    FitrnetTSLTPsummarydouble(1:cutAll,6)];
colors = {'#2D42C2','#405A00','#932055','#33486B',...
    '#405A00','#10930E','#74356E'};


xLabels = {'MACD','Double7','Bollinger Band 1','Bollinger Band 2',...
    'Fitrnet Std','Fitrnet TSLTP','Buy and Hold'}; 
figure;
b1 = boxchart( avgRetPA(:,[1,2,3,4,5,6,7]) ); % 
yline(0,'-k')
xlabel('Strategy')
ylabel('Average Return per anuum')
title('Average Return p.a. of Technical Strategies and Neural Network Strategies')
set(gca, 'XTickLabels', xLabels) 
hold off



%% Box Plot other verisons (not used)
%xLabelsStrategies = {'MACD','Double7','Bollinger Band 1','Bollinger Band 2','Buy and Hold'}; 
% BPavgRetplot = figure;
% b1 = boxplot(avgRetPA(:,[1,2,3,4,5]),'ColorGroup',colors,'symbol',''); % 
% yline(0,'-k')
% xlabel('Strategy')
% ylabel('Average Return per anuum')
% title('Average Return p.a. Boxplot')
% set(gca, 'XTickLabels', xLabels) 
% 
% 
% MACD=#2D42C2 D7=#405A00 BB1=#932055 BB2=#33486B benchmark=#74356E
% figure;
% ax2 = axes();
% hold(ax2);
% for i=1:size(avgRetPA,2)
%     boxchart(i*ones(size(avgRetPA(:,i))),avgRetPA(:,i),'BoxFaceColor',colors{i})
% end
% yline(0,'-k')
% xlabel('Strategy')
% ylabel('Average Return per anuum')
% title('Average Return p.a. of Technical Strategies and Neural Network Strategies')
% legend('MACD','Double7','Bollinger Band 1','Bollinger Band 2',...
%     'Fitrnet Std','Fitrnet TSLTP','Fitnet Std','Fitnet TSLTP','Buy and Hold')
% hold off



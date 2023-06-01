%% Bollinger Band Strategy
clear
clc
set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultLineLineWidth', 2);

%%%% Data Import %%%%%%%%%%%%%%%
timeRange = 'Daily';
%timeRange = 'Weekly';
%timeRange = 'Monthly';

addpath(genpath("/Users/robinkulha/Documents/MATLAB/Bachelorarbeit"));
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Mat Data"));

files = dir(fullfile('/Users/robinkulha/Documents/MATLAB/','Mat Data','*',timeRange,'*.mat*'));
N = length(files) ;   % total number of files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assign cell for summarytable
BB2summarycell = cell(length(N),14);
%load('BB2summarytable.mat')

tic

% loop through each file 
for i = 1:160 % ACHTUNG: 1:N 
    thisfile = files(i).name;

    % Using Mat files 
    Data = load(fullfile(files(i).folder, thisfile));
    if strcmp(timeRange, 'Daily')
        Data = Data.data;
    else
        Data = Data.(['data', timeRange]);
    end
    %Data(2) = Open, Data(3) = High, Data(4) = Low, Data(5) = Close
    
    % Make Data coherent, so that there are no gaps
    Data = cutData(Data, timeRange);
    
    %% Prepare WFA
    
    nPeriod = 100;
    nISTPeriod = 5;
    
    if size(Data(:,1),1) < 2000
        disp('Too small Data size. This set will be skipped');
        filename = strrep(thisfile(1:end-5-length(timeRange)),"_"," ");
        BB2summarycell{i, 1} = filename;
        BB2summarycell{i, 2} = 0;
        BB2summarycell{i, 3} = 0;
        BB2summarycell{i, 4} = 0;
        BB2summarycell{i, 5} = 0;
        BB2summarycell{i, 6} = 0;
        BB2summarycell{i, 7} = 0;
        BB2summarycell{i, 8} = 0;
        BB2summarycell{i, 9} = 0;
        BB2summarycell{i, 10} = 0;
        BB2summarycell{i, 11} = 0;
        BB2summarycell{i, 12} = 0;
        BB2summarycell{i, 13} = 0;
        BB2summarycell{i, 14} = 0;
        %summarycell{i, 15} = slopeMC;
        %summarycell{i, 16} = probPosRes;
        continue
    end

%     % Adapt nPeriod considering Data length
%     ratioPeriod = 15;
%     if size(Data(:,1),1)/nPeriod < ratioPeriod
%         ratioISTPeriod = nISTPeriod / nPeriod;
%         nPeriod = round(size(Data(:,1),1)/ratioPeriod, 0);
%         nISTPeriod = round(nPeriod * ratioISTPeriod,0);
%         % Omit too small data and set results to NaN
%         if nISTPeriod == 0 && nPeriod > 1
%             nISTPeriod = 1;
%         elseif nISTPeriod == 0 || nPeriod < 2
%             disp('Too small Data size. This set will be skipped');
%             continue
%         end
%     end
    
    % Get a table of FLA dates - IS and OS - from start to end
    wfaPeriod = myFLA(Data, timeRange, nPeriod, nISTPeriod);
    
    %% Assigning zeros
    Data_OS_all = [];
    tradesAll_OS_all = [];
    position_all = [];
    tradesAllClean_OS_all = [];
    %paramzeros = zeros(height(wfaPeriod),4);
    
        %% Forward Looking Analysis
    parfor w = 1:height(wfaPeriod)
        %% Get data fitting to the current walk
        [Data_IS, Data_OS] = wfaData(Data, wfaPeriod, w);
        
        %% IS Training
        
        % prepare optimisation
        RiskManagement = 'conservative';
        plotparam = 0;
        BB = 2;
        fitData = {RiskManagement, plotparam, BB};
        
        % Use datasz for checking if chromosome is within the constriants
        datasz = length(Data_IS);
        
        % Optimise parameters
        [~,periodBB,nStd,SL,TP] = optGenBB(Data_IS,fitData,datasz);
        
        % Get IS Values
        [EC_IS,lastEC_IS,~,tradesAllClean_IS,~] = sBB(Data_IS,periodBB,nStd,SL,TP,RiskManagement,plotparam,BB);
        
        % save used parameters
        paramarray = [periodBB,nStd,SL,TP];
        %paramzeros(w,1:end) = paramarray;
        
        
        %% OS Testing 
        % Add IS Data to evaluate EMAs before original Data_OS needs ddata
        cut = periodBB;
        Data_OS_work = [Data_IS(end-cut+2:end,:); Data_OS];
        % end-cut+2, because if eg. cut = 302, the 302 period is calculated
        
        [~,endvalueEC2_OS,tradesAll_OS,tradesAllClean_OS,position_OS] ...
            = sBB(Data_OS_work,periodBB,nStd,SL,TP,RiskManagement,plotparam,BB);
        
        
        % Filling arrays
        Data_OS_all = [Data_OS_all; Data_OS]; %#ok<*AGROW>
        tradesAll_OS_all = [tradesAll_OS_all; tradesAll_OS];
        position_all = [position_all; position_OS];
        tradesAllClean_OS_all = [tradesAllClean_OS_all; tradesAllClean_OS];
        
        disp(w);
    end
    
    
    %% Prepare data
    % Cumulate sum for EC_OS_all
    EC_OS_all = cumsum(tradesAll_OS_all(:,2));
    
    
    %% Variables for Summarytable %% 
    
    % Profit and Loss
    PnL = EC_OS_all(end);
    % Total Return [%]
    totRet = (PnL / Data_OS_all(1,5))*100;
    % Annual Return [%]
    years = length(unique(str2num(datestr(tradesAll_OS_all(:,1),'yyyy')))); %#ok<ST2NM>
    avgPercPA = totRet/years;
    
    % Profit of Buy and Hold strategy
    benchmark = Data_OS_all(end,5) - Data_OS_all(1,5);
    % Total Return benchmark [%]
    totRetBench = (benchmark / Data_OS_all(1,5))*100; 
    % Annual Return benchmark [%]
    avgPercPABench = totRetBench/years;
    
    % Max Drawdown
    [MaxDD, MaxDDIndex] = MAXDRAWDOWN(EC_OS_all);
    % Number of Trades
    TradesNum = length(tradesAllClean_OS_all);
    % PL per Trade
    PLpTrade = PnL/TradesNum;
    % Ratio of being invested the Market (=Exposure)
    exposure = sum(sum(position_all~=0))/length(position_all);
    % Return over Max Drawdown (RoMaD) avgPercPA/MaxDD
    RoMaD = avgPercPA/MaxDD;
    % Sharp Ratio (Return - Risk Free Return)/Std
    SR = (PnL - 0)/std(EC_OS_all);
    
    % Linear Regression
    linLen = (1:length(EC_OS_all));
    [p,S] = polyfit(linLen,EC_OS_all, 1); slopeLR = p(1);
    [f, delta] = polyval(p,(1:length(EC_OS_all)),S);
    % Monte Carlo Simulation
    if tradesAllClean_OS_all(:,2) ~= 0
        [slopeMC, probPosRes] = MonteCarloSim(tradesAllClean_OS_all(:,2), 0);
    else
        slopeMC = NaN;
        probPosRes = NaN;
    end
    
    
   %% Save results of share i to summarycell
   filename = strrep(thisfile(1:end-5-length(timeRange)),"_"," ");
   BB2summarycell{i, 1} = filename;
   BB2summarycell{i, 2} = PnL;
   BB2summarycell{i, 3} = totRet;
   BB2summarycell{i, 4} = avgPercPA;
   BB2summarycell{i, 5} = benchmark;
   BB2summarycell{i, 6} = totRetBench;
   BB2summarycell{i, 7} = avgPercPABench;
   BB2summarycell{i, 8} = TradesNum;
   BB2summarycell{i, 9} = MaxDD;
   BB2summarycell{i, 10} = PLpTrade;
   BB2summarycell{i, 11} = exposure;
   BB2summarycell{i, 12} = SR;
   BB2summarycell{i, 13} = RoMaD;
   BB2summarycell{i, 14} = slopeLR;
   %summarycell{i, 15} = slopeMC;
   %summarycell{i, 16} = probPosRes;
   
   plots = 1;
   if plots == 1 && i == 1
       %=========== Plot OS and benchmark =====================
       % Prepare Close data
       Close = Data_OS_all(1:end,5) - Data_OS_all(1,5); % 0 because EC_OS_all also has zero
       
       figEC_OS = figure;
       plot(EC_OS_all,'Color','#932055'); % BB1 EC:#932055 BB2 EC:#33486B
       hold on;
       plot(Close,'Color','#74356E'); % benchmark
       if MaxDD ~= 0
           plot(MaxDDIndex(1):MaxDDIndex(2), ...
               EC_OS_all(MaxDDIndex(1):MaxDDIndex(2)),'Color','#A2142F');
       end
       hold off;
       
       title("Bollinger Band Simple EC and benchmark of " + filename);
       yline(0,'k-');
       xlabel("Time");
       ylabel("Profit (absolute)");
       legend('EC_B_B_1','benchmark','MaxDD')
       %=========================================================
   end
   disp("Share Number " + i);
end

timepassed = toc;

%% Delete not used shares, either PnL = 0 or slopeMC = NaN
BB2summarycell = BB2summarycell(cellfun(@(a)~isempty(a)&&a>0,BB2summarycell(:,2)),:);
BB2summarycell = BB2summarycell(cellfun(@(a)~isempty(a)&&a>0,BB2summarycell(:,14)),:);

%% Create Summarytable of all shares
summarytable = cell2table(BB2summarycell);
summarytable.Properties.VariableNames = {'Name','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 

%% Create one line of the strategy result

summaryline = cell(1, 14);
summaryline{1,1} = 'BB Simple';
for ss = 2:14
    summaryline{1,ss} = mean(cell2mat(BB2summarycell(:,ss)));
end

summarytable = table({filename},PnL,totRet,avgPercPA,benchmark,totRetBench,avgPercPABench,MaxDD,TradesNum,PLpTrade,...
    exposure,RoMaD,SR,slopeLR); 
summarytable.Properties.VariableNames = {'Name','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 


%writetable(summarytable, 'MACDISOS.xlsx')

% to concatenate more summarylines 
% greatsummary = [summarylineMACD;summarylineXXX]

%% compare BB Simple mit hist (EMA = 20, nStd = 2)

EChistBB1 = EChistBB1(length(EChistBB1)-length(Close)+1:length(EChistBB1));
figcompare = figure;
plot(Close,'Color','#74356E'); % benchmark
hold on;
plot(EC_OS_all,'Color','#33486B'); % BB1 EC:#932055 BB2 EC:#33486B
plot(EChistBB1,'Color','#2D42C2');
hold off;

title("Bollinger Band Complex of " + filename);
yline(0,'k-');
xlabel("Time");
ylabel("Profit (absolute)");
legend('benchmark','EC_B_B_2','EC_h_i_s_t')




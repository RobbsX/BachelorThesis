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


% loop through each file 
for i = 1 % ACHTUNG: 1:N 
    thisfile = files(i).name ;

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
    
    TP = 5; SL = 2;
    periodHL=7;
    periodEMA=200;
    [EC,endvalueEC2,tradesAll,tradesAllClean,position] = sDouble7(Data,periodHL,periodEMA,SL,TP,'conservative',0);
    
    % Cut Data just like in sMACD
    cut = periodEMA;
    Data = Data(cut:end,1:5);
    
    
    % Profit and Loss 
    PnL = EC(end);
    % Total Return [%]
    totRet = sum(tradesAllClean(:,3));
    % Annual Return [%]
    %[yearPerc,avgPercPA] = annualReturn(tradesAll);
    
    % Profit of Buy and Hold strategy
    benchmark = Data(end,5) - Data(1,5);
    % Total Return benchmark [%]
    totRetBench = Data(end,5)/Data(1,5)-1;
    % Annual Return benchmark [%]
    %avgPercPABench = totRetBench/size(yearPerc(:,1),1);
    
    % Max Drawdown
    [MaxDD, MaxDDIndex] = MAXDRAWDOWN(EC);
    % Number of Trades
    TradesNum = sum(tradesAllClean(:,2)~=0);
    % PL per Trade
    PLpTrade = PnL/TradesNum;
    % Ratio of being invested the Market (=Exposure)
    exposure = sum(position~=0)/length(position); % Updated 31.8.
    % Return over Max Drawdown (RoMaD) avgPercPA/MaxDD
    %RoMaD = avgPercPA/MaxDD;
    % Sharp Ratio (Return - Risk Free Return)/Std
    SR = (PnL - 0)/std(EC);
    
    % Linear Regression
    linLen = (1:length(EC));
    [p,S] = polyfit(linLen,EC, 1); slopeLR = p(1);
    [f, delta] = polyval(p,(1:length(EC)),S);
    % Monte Carlo Simulation
    if and(tradesAllClean(:,2) ~= 0, size(tradesAllClean(:,2),1)>1)
        [slopeMC, probPosRes] = MonteCarloSim(tradesAllClean(:,2), 0);
    else
        slopeMC = NaN;
        probPosRes = NaN;
    end

end



% Plot EC
filename = split(thisfile, "_"); filename = filename{1};
figure;
plot([0;EC],'Color','#2D42C2');
hold on;
plot(Data(1:end,5) - Data(1,5),'Color','#74356E');
if MaxDD ~= 0
    plot(MaxDDIndex(1):MaxDDIndex(2), EC(MaxDDIndex(1):MaxDDIndex(2)),'Color','#A2142F');
end

title("Equity Curve of " + filename);
yline(0,'k-');
xlabel("Time");
ylabel("Profit (absolute)");
legend('EC','Benchmark')
hold off;


% For more Shares --> lists for summary table

summarytable = table({filename},PnL,totRet,avgPercPA,benchmark,totRetBench,avgPercPABench,MaxDD,TradesNum,PLpTrade,...
    exposure,RoMaD,SR,slopeLR); 
summarytable.Properties.VariableNames = {'Name','PnL','totRet','avgReturnPA','PnL BnH','totRet BnH','avgRetPA BnH','MaxDD','TradesNum','PnLpTrade',...
    'Exposure','RoMaD','Sharp Ratio','slopeLR'}; 




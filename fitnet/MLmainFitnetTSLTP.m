clearvars;
clc;
warning ('off','all'); 
set(0,'DefaultFigureWindowStyle','docked')

%% Choose to use Model or not
useModel = 1;


%% Files & User
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Bachelorarbeit"));
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Mat Data"));

timeRange = 'Daily';
%timeRange = 'Weekly';
%timeRange = 'Monthly';
files = dir(fullfile('/Users/robinkulha/Documents/MATLAB/','Mat Data','*',timeRange,'*.mat*'));
N = length(files) ;   % total number of files 

%% Assign cell for summarytable
summarycell = cell(length(N),14);

%% Prepare struct
% MScollectionFitnet = struct; % just used for the first time
load('MScollectionFitnet.mat') % Real name? 
if useModel == 0
    istep = length(fieldnames(MScollectionFitnet))+1:N;
else
    istep = 1:length(fieldnames(MScollectionFitnet));
end
MS = fieldnames(MScollectionFitnet);

for iallocate = istep
    
i = str2double(MS{iallocate}(end));

thisfile = files(i).name ;
filename = strrep(thisfile(1:end-5-length(timeRange)),"_"," ");
% Using Mat files
Data = load(fullfile(files(i).folder, thisfile)); 
if strcmp(timeRange, 'Daily')
    Data = Data.data;
else
    Data = Data.(['data', timeRange]);
end
%Data(2) = Open, Data(3) = High, Data(4) = Low, Data(5) = Close

%% Prepare Data
Data = cutData(Data, timeRange);
Data = Data(:, 1:5);

%% Define parameters for training
xnodes = 20; % input nodes, CE best result at 6
ynodes = 1; % output nodes, must be 1
predictOffset = 5; % x days from the day the prediction used the last data point
% Could be dependent on xnodes (eg. 5*xnodes) and optimized.
% The period can't be very long because the market changes and the
% prediction would be very bad. (18.5.21)


nPeriod = 100; % Used: 100 & 5. Test: 5 & 2.
nISTPeriod = 5;

%% Prepare and/or load the ModelSet
if useModel == 0
    ModelSet = cell(nPeriod-nISTPeriod ,8);
else
    currentMS = MS{iallocate};
    ModelSet = MScollectionFitnet.(currentMS);  
end


%% Get a table of FLA dates - IS and OS - from start to end
wfaPeriod = myFLA(Data, timeRange, nPeriod, nISTPeriod);


%% Prepare lists
% Data lists
Ylist = [];
Xlist = [];
Dlist = [];
Data_OS_list = [];
% Result lists
positionlist = [];
profitOSlist = [];
tradesAllCleanlist = [];
tradesAlllist = [];
lenOptSLTP = 10; SLTPcollection = cell(lenOptSLTP, 4);



%% Forward Looking Analysis
tic
for w = 1:height(wfaPeriod)
    
    % Get data fitting to the current walk
    [Data_IS, Data_OS] = wfaData(Data, wfaPeriod, w, predictOffset+xnodes);
    
    %% In Sample %%
    
    predstepIS = 1:predictOffset:size(Data_IS(:,1),1) -predictOffset-xnodes+1;
    
    % Set zeros
    X_IS = zeros(length(predstepIS), xnodes);
    
    % Split Data
    for tt = 1:length(predstepIS) % 1:size(Data_IS(:,1),1) -predictOffset-xnodes
        X_IS(tt,:) = Data_IS((1:xnodes) +predstepIS(tt)-1, 5)'; % train data (last day of training is 'today' or 'current day')
    end
    D_IS = Data_IS(xnodes+predictOffset :predictOffset: size(Data_IS(:,1),1), 5); % Control
    
    % Normalisation x_n = (x-min)/(max-min)
    minXIS = min(X_IS(:,1)); %minX = min(X_IS(X_IS(:,1)~=0));
    maxXIS = max(X_IS(:,1));
    XnormIS = (X_IS-minXIS)/(maxXIS-minXIS) *(0.9999-0.0001)+0.0001;
    DnormIS = (D_IS-min(D_IS(D_IS(:,1)~=0)))/(max(D_IS(D_IS(:,1)~=0))-min(D_IS(D_IS(:,1)~=0))) *(0.9999-0.0001)+0.0001;
    
    check = 0; c = 0; maxitaration = 2;
    if useModel == 0 % use to test whole code quickly or to save and continue process
        while check == 0
            % Train model
            [net, ~] = MLoptFitnet(X_IS, D_IS);
            Y_IS = net(X_IS');
            %Y_IS = YnormIS * (maxXIS-minXIS) + minXIS;
            stdMdl = std(Y_IS)/mean(X_IS(:,end)); % check real values' std
            c = c + 1; disp("c = " + c);
            if stdMdl > 0.03 % above 1 probably is a decent std
                check = 1;
            elseif c >= maxitaration
                disp("For " + w + " no trades.");
                break % no trades will be made
            end
        end
    else % Use ModelSet
        net = ModelSet{w,1};
    end
    
    
    
    %% Out of Sample %%
    
    predstepOS = 1:predictOffset:size(Data_OS(:,1),1)-predictOffset-xnodes+1;
    % Set zeros
    X_OS = zeros(length(predstepOS), xnodes);
    % X_OS = zeros(size(Data_OS(:,1),1)-xnodes-predictOffset, xnodes);
    
    
    % Setting testData for evaluation
    for tt = 1:length(predstepOS)
        X_OS(tt,:) = Data_OS((1:xnodes) +predstepOS(tt)-1, 5)'; % train data (last day of training is 'today' or 'current day')
    end
    D_OS = Data_OS(xnodes+predictOffset :predictOffset: size(Data_OS(:,1),1), 5);
    
    % Normalisation x_n = (x-min)/(max-min)
    minXOS = min(X_OS(X_OS(:,1)~=0));
    maxXOS = max(X_OS(X_OS(:,1)~=0));
    XnormOS = (X_OS-minXOS)/(maxXOS-minXOS) *(0.9999-0.0001)+0.0001;
    DnormOS = (D_OS-min(D_OS(D_OS(:,1)~=0)))/(max(D_OS(D_OS(:,1)~=0))-min(D_OS(D_OS(:,1)~=0))) *(0.9999-0.0001)+0.0001;
    
    
    
    %%% Make predictions of IS and OS
    Y_IS = net(X_IS');  Y_OS = net(X_OS');
    Y_IS = Y_IS';       Y_OS = Y_OS';
    
    
    
    % RSME
    rsme_IS = sqrt(mean((Y_IS-D_IS).^2));
    rmse_OS = sqrt(mean((Y_OS-D_OS).^2));
    
    % Correlation of IS
    correlation_IS = corrcoef(D_IS,Y_IS); correlation_IS = correlation_IS(2);
    
    % Correlation of OS
    correlation_OS = corrcoef(D_OS,Y_OS); correlation_OS = correlation_OS(2);
    
    % RSQ OS
%     sse = sum((D_OS - Y_OS).^2);
%     sst = sum((D_OS - mean(D_OS)).^2);
%     RSQOS = 1-sse/sst;
    
    %% Strategy
    % To set strategy
    Data_OSnew = Data_OS(predictOffset+xnodes :predictOffset: end,:);
    vola = std(X_IS(:,1:end))/mean(X_IS(:,1:end)); % No OS data, because std() is only one const. value, instead of an array! --> IS data!
    
    %%% Get Trades to save in ModelSet
    position = zeros(size(Y_OS,1),1);
    RiskManagement = 'conservative';
    if c < maxitaration
        if strcmp(RiskManagement, 'conservative')
            position( and(Y_OS>X_OS(:,end), Y_OS>X_OS(:,end)*(1+vola*(1/3))) ) =  1;
            position( and(Y_OS<X_OS(:,end), Y_OS<X_OS(:,end)*(1-vola*(1/3))) ) = -1;
        else
            position( Y_OS>X_OS(:,end) ) =  1;
            position( Y_OS<X_OS(:,end) ) = -1;
        end
    end
    
    if w <= lenOptSLTP %w <= 0 
        profitOS = 0; % because not used
        %% Optimise SL and TP in x periodes of OS data %%
        if w == 1
            penalty = 0.02; % inizialise penalty
        end
        % SLTPparams = [xSL,xTP,penalty];
        fitData = {position,predictOffset,xnodes,RiskManagement,0,penalty};
        SLTPparams = optGenTSLTP(Data_OS,fitData);
        
        % store results
        for xx = 1:3
            SLTPcollection{w,xx} = SLTPparams{xx};
        end

        
        % use mean of optSLTP collection
        if w == lenOptSLTP
            avgCollection = 1;
            if avgCollection == 1 % use mean of the  results
                xSL = mean(cell2mat(SLTPcollection(:,1)));
                xTP = mean(cell2mat(SLTPcollection(:,2)));
            else % use only the best result
                choice = SLTPcollection(cell2mat(SLTPcollection(:,end))==max(cell2mat(SLTPcollection(:,end))),:);
                xSL = choice{1};
                xTP = choice{2};
            end
            penalty = 0; % will not be used after W == lenOptSLTP
            SLTPparams = {xSL,xTP,penalty};
        end
    elseif w > lenOptSLTP %else%
        %% Use optimal SL und TP for next OS data 
        Data_OSoffset = Data_OS(predictOffset+xnodes :predictOffset: end,:);
        [EC,profitOS,tradesAll,tradesAllClean] = ...
            evalTradesOffestTSLTP(Data_OS,position,SLTPparams,predictOffset,xnodes,RiskManagement,0);
        
        % Fill lists
        % EClist = [EClist; EC];
        positionlist = [positionlist; position];
        profitOSlist = [profitOSlist; profitOS];
        tradesAllCleanlist = [tradesAllCleanlist; tradesAllClean];
        tradesAlllist = [tradesAlllist; tradesAll];
        

        % Save the arrays to make sure everything is alright
        Ylist = [Ylist; Y_OS];  %#ok<*AGROW>
        Xlist = [Xlist; X_OS];
        Dlist = [Dlist; D_OS];
        Data_OS_list = [Data_OS_list; Data_OSoffset];
        
        %% Correlation Plot
        corrplot = 1;
        if corrplot == 1 && i==1 && w == 1
            figure;
            subplot(2,1,1)
            plot(D_IS,D_IS,'-r','LineWidth',2)
            hold on;
            plot(D_IS,Y_IS,'bo','MarkerSize',8,'MarkerFaceColor','m')
            title("Correlation of IS: " + correlation_IS);
            
            subplot(2,1,2)
            plot(D_OS,D_OS,'-r','LineWidth',2)
            hold on;
            plot(D_OS,Y_OS,'bo','MarkerSize',8,'MarkerFaceColor','m')
            title("Correlation of OS: " + correlation_OS);
        end
        
        
    end
    
    % Save trainedModels to compare the RMSE
    ModelSet{w,1} = net;
    ModelSet{w,2} = [rsme_IS, rmse_OS];
    ModelSet{w,3} = [correlation_IS, correlation_OS];
    ModelSet{w,4} = [net.layers{1}.size, net.layers{2}.size];
    ModelSet{w,5} = net.layers{1}.transferFcn;
    ModelSet{w,6} = net.trainFcn;
    ModelSet{w,7} = std(Y_OS);
    ModelSet{w,8} = profitOS;
    
    % Estimate time remaining
    %     timePassed = toc;
    disp( w + "/" + height(wfaPeriod) );
    %     timeRemaining = (height(wfaPeriod)-w)*timePassed;
    %     disp( "Estimated time remaining: " + timeRemaining + "s");
end


%% Prepare data and EC
% Delete first values (predictOffset+xnodes) of Data_OS_list
% because there can't be a prediction for them
Data_OS_list = Data_OS_list(predictOffset+xnodes:end,:);
ECplotChoice = 'daily'; % or 'tradesOnly'
if strcmp(ECplotChoice, 'tradesOnly')
    EClist = cumsum(tradesAllCleanlist(:,2));
elseif strcmp(ECplotChoice, 'daily')
    EClist = cumsum(tradesAlllist(:,2));
end
% EClist = cumsum(EClist); % Wrong!


%% Plot Price and Prediction with RSME

% absolut difference
diff = Ylist-Dlist; 
% RSME
rmse = sqrt(mean((diff).^2));

plotting = 1;
if plotting == 1  && i == 1
    figure;
    subplot(2,1,1)
    plot(Dlist,'Color','#2D42C2','LineWidth',1);
    hold on;
    plot(Ylist,'Color','#932055','LineWidth',1);
    hold off;
    xlabel(timeRange);
    ylabel("USD absolut");
    title("Forecast");
    legend(["Price" "Forecast"]);
    
    % Plot error
    subplot(2,1,2)
    stem(diff)
    xlabel(timeRange)
    ylabel("Error")
    title("RMSE = " + rmse)
end 





%% Variables for Summarytable %%

% Profit and Loss
PnL = EClist(end);
% Total Return [%]
totRet = (PnL / Data_OS_list(1,5))*100;
% Annual Return [%]
years = length(unique(str2num(datestr(tradesAlllist(:,1),'yyyy')))); %#ok<ST2NM>
avgPercPA = totRet/years;

% Profit of Buy and Hold strategy
benchmark = Data_OS_list(end,5) - Data_OS_list(1,5);
% Total Return benchmark [%]
totRetBench = (benchmark/Data_OS_list(1,5))*100;
% Annual Return benchmark [%]
avgPercPABench = totRetBench/years;

% Max Drawdown
[MaxDD, MaxDDIndex] = MAXDRAWDOWN(EClist);
% Number of Trades
TradesNum = length(tradesAllCleanlist);
% PL per Trade
PLpTrade = PnL/TradesNum;
% Ratio of being invested the Market (=Exposure)
exposure = sum(positionlist~=0)/length(positionlist); 
% Return over Max Drawdown (RoMaD) avgPercPA/MaxDD
RoMaD = avgPercPA/MaxDD;
% Sharp Ratio (Return - Risk Free Return)/Std
SR = (PnL - 0)/std(EClist);
% Benchmark
benchmark = Data_OS_list(end,2) - Data_OS_list(1,2);

% Linear Regression
linLen = (1:length(EClist));
[p,S] = polyfit(linLen,EClist, 1); slopeLR = p(1);
[f, delta] = polyval(p,(1:length(EClist)),S);
% Monte Carlo Simulation
if and(tradesAllCleanlist(:,2) ~= 0, size(tradesAllCleanlist(:,2),1)>1)
    [slopeMC, probPosRes] = MonteCarloSim(tradesAllCleanlist(:,2), 0);
else
    slopeMC = NaN;
    probPosRes = NaN;
end

%% Create Table of Models for better reading
ModelSetTable = cell2table(ModelSet);
ModelSetTable.Properties.VariableNames = {'Model','RSME','RSQ','Layer Sizes',...
    'Activation','Training Algo','OS Std','OS Profit'};


%% Save results of share i to summarycell
summarycell{i, 1} = filename;
summarycell{i, 2} = PnL;
summarycell{i, 3} = totRet;
summarycell{i, 4} = avgPercPA;
summarycell{i, 5} = benchmark;
summarycell{i, 6} = totRetBench;
summarycell{i, 7} = avgPercPABench;
summarycell{i, 8} = TradesNum;
summarycell{i, 9} = MaxDD;
summarycell{i, 10} = PLpTrade;
summarycell{i, 11} = exposure;
summarycell{i, 12} = SR;
summarycell{i, 13} = RoMaD;
summarycell{i, 14} = slopeLR;
%summarycell{i, 15} = slopeMC;
%summarycell{i, 16} = probPosRes;    

%% Plot all OS restuls after each other
if plotting == 1 && i == 1
    %=========== Plot OS and benchmark =====================
       % Prepare Close data. NOTE: The Data has been 'compressed' using
       % periodOffset
       Close = Data_OS_list(1:end,5) - Data_OS_list(1,5);
       
       figEC_OS = figure;
       plot(EClist,'Color','#0068FF'); % fitnet TSLTP 
       hold on;
       plot(Close,'Color','#74356E'); % benchmark 
       if MaxDD ~= 0
           plot(MaxDDIndex(1):MaxDDIndex(2), ...
               EClist(MaxDDIndex(1):MaxDDIndex(2)),'r-');
       end
       hold off;
       
       title("OS EC and benchmark of " + filename);
       yline(0,'k-');
       xlabel("Time");
       ylabel("Profit (absolut)");
       legend('EC_T_S_L_T_P','benchmark','MaxDD')
       %=========================================================
end

if useModel == 0
    %% Save Model in Struct
    saveName = "ModelSet" + i;
    MScollectionFitnet.(saveName) = ModelSet;
end


disp("Share Number " + i);
end

timepassed = toc;


%% Create Summarytable of all shares
summarytable = cell2table(summarycell);
summarytable.Properties.VariableNames = {'Name','PnL', '%total', '%pa','Benchmark', '%totalB', '%paB','TradesNum',...
    'MaxDD','PLperTrade','Exposure','Sharp Ratio','RoMaD','Slope LR'}; %,'Slope MC','probPosRes'






%% Not used

% % Statistically Significant RSME 
% 
% RSMEs = ModelSet(:,2);
% RSMEs = cell2mat(RSMEs);
% % Linear Regression for RSME
% lmRSME = fitlm(RSMEs(:,1) , RSMEs(:,2));
% 
% 
% % Statistically Significant Correlation
% 
% Corr = ModelSet(:,3);
% Corr = cell2mat(Corr);
% % Linear Regression for RSME
% lmCorr = fitlm(Corr(:,1) , Corr(:,2));
% 
% 
% % Correlation of OS Corr and Std
% if plotting == 2
%     stdd = ModelSet(:,7);
%     stdd = cell2mat(stdd);
%     figure;
%     plot(Corr(:,2),Corr(:,2),'-r','LineWidth',2)
%     hold on;
%     plot(Corr(:,2),stdd,'bo','MarkerSize',8,'MarkerFaceColor','m')
%     xlabel("OS Correlation");
%     ylabel("Standard Deviation");
%     title("Correlation of OS Correlation and Std");
% end



%%% Plots save 17.8.21
% 
% Y_IS_MS250 = predict(ModelSet{w,1}.Regression, XnormIS);
% Y_OS_MS250 = predict(ModelSet{w,1}.Regression, XnormIS);
% YnormIS_MS250 = predict(ModelSet{w,1}.Regression, XnormIS);
% YnormOS_MS250 = predict(ModelSet{w,1}.Regression, XnormIS);
% Y_IS_MS250 = YnormIS_MS250 * (maxXIS-minXIS) + minXIS;
% Y_OS_MS250 = YnormOS_MS250 * (maxXOS-minXOS) + minXOS;
% rsme_IS_MS250 = sqrt(mean((Y_IS-D_IS).^2));
% rmse_OS_MS250 = sqrt(mean((Y_OS-D_OS).^2));
% 
% [Y_IS_NNFlm,Xf,Af] = myNeuralNetworkFunction(X_IS);
% [Y_OS_NNFlm,Xf,Af] = myNeuralNetworkFunction(X_OS);
% rsme_IS_NNF = sqrt(mean((Y_IS_NNFlm-D_IS).^2));
% rmse_OS_NNF = sqrt(mean((Y_OS_NNFlm-D_OS).^2));
% 
% figure;plot(D_IS, 'b-');
% hold on;
% plot(Y_IS, 'g-o');
% plot(Y_IS_NNFlm, 'c-o');
% plot(Y_IS_MS250, 'r-o');
% legend('Close','2Layer','Net Fit','3Layer')
% 
% [Y_IS_Bays,Xf,Af] = myNeuralNetworkFunctionBays(X_IS)
% rsme_IS_Bays = sqrt(mean((Y_IS_Bays-D_IS).^2));
% rsme_OS_Bays = sqrt(mean((Y_OS_Bays-D_IS).^2));
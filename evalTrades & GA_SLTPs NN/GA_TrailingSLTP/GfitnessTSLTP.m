function [P, P_dec] = GfitnessTSLTP(data, P, bit_size, rng, fitData)
% data...........Data: 1=Open, 2=High, 3=Low, 4=Close. Will be used to
% backtest
% P..............Input Population Matrix
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
%
% P..............Output Population Matrix

position = fitData{1};
predictOffset = fitData{2};
xnodes = fitData{3};
RiskManagement = fitData{4};
plotparam = fitData{5};

P_bi = cell(length(P),length(bit_size)+1);
P_dec = cell(length(P),length(bit_size)+1);
fit = cell(length(P),1);

% Delete old fitness values (needed to concatinate RSME using parfor)
P = P(:,1);
P_dec = P_dec(:,1:length(bit_size));
P_bi = P_bi(:,1:length(bit_size));

for i = 1:length(P) 
    
    ende = 0;
    c = 0;
    % split the single chromosome read it
    for j = bit_size
        start = ende + 1;
        ende = ende + j; % fixed. before: bit_size(1,j)
        c = c+1; % count for parameter / bit_size
        P_bi{i,c} = P{i,1}(start:ende);
        P_dec{i,c} = Gdec(P{i,1}(start:ende), rng(1,c), rng(2,c));
    end
    
end
    

for ii = 1:length(P) %parfor
    
    P_oneline = P_dec(ii,:); % due to parfor
    
    % saving the value for calculating the fitness
    xSL = P_oneline{1}; 
    xTP = P_oneline{2}; 
    SLTPparams = {xSL,xTP};
    
    % Calculate the fitness value 
    [~,endvalueEC2] = evalTradesOffestTSLTP(data,position,SLTPparams,predictOffset,xnodes,RiskManagement,plotparam);
    
    % assign fit to the matrix
    fit{ii,1} = endvalueEC2; 
    %P_dec{i,length(bit_size)+1} = validationRMSE;
    %P_bi{i,length(bit_size)+1} = validationRMSE;
end

    % Due to parfor RSMEs addad to P after loop
    P = [P, fit];
    P_dec = [P_dec, fit];
    P_bi = [P_bi, fit];
    
    % sort P using fit
    P = sortrows(P, -2); % - for descending. fit --> profit
    P_dec = sortrows(P_dec, -(length(bit_size)+1));
    P_bi = sortrows(P_bi, -(length(bit_size)+1)); %#ok<NASGU>
    
end



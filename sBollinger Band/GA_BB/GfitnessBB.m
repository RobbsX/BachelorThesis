function [P, P_dec] = GfitnessBB(data, P, bit_size, rng, fitData)
% data...........Data: 1=Open, 2=High, 3=Low, 4=Close. Will be used to
% backtest
% P..............Input Population Matrix
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
%
% P..............Output Population Matrix

RiskManagement = fitData{1};
plotparam = fitData{2};
BB = fitData{3};

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
    
    periodBB = P_oneline{1};
    nStd = P_oneline{2};
    SL = P_oneline{3};
    TP = P_oneline{4};
    
    % Calculate the fitness value 
    [~,endvalueEC2] = sBB(data,periodBB,nStd,SL,TP,RiskManagement,plotparam,BB);

    % assign fit to the matrix
    fit{ii,1} = endvalueEC2; 
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



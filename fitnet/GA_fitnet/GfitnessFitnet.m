function [P, P_dec] = GfitnessFitnet(data, P, bit_size, rng)
% data...........Data: 1=Open, 2=High, 3=Low, 4=Close. Will be used to
% backtest
% P..............Input Population Matrix
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
%
% P..............Output Population Matrix

% Unpack data
trainingData = data{1};
responseData = data{2};

P_bi = cell(length(P),length(bit_size)+1);
P_dec = cell(length(P),length(bit_size)+1);
RSMEs = cell(length(P),1);

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
    
    P_dec{i,1} = round(P_dec{i,1});
    P_dec{i,2} = round(P_dec{i,2});
    
    if P_dec{i,3} == 1
        activation = 'relu';
    elseif P_dec{i,3} == 2
        activation = 'tanh';
    elseif P_dec{i,3} == 3
        activation = 'none';
    elseif P_dec{i,3} == 4
        activation = 'sigmoid';
    end
    P_dec{i,3} = activation;
    
    P_dec{i,4} = round(P_dec{i,4});
    optimizerOpt = ["trainlm", "trainbr", "trainbfg", "trainrp",... 
        "trainscg","traincgb", "traincgf", "traingd"]; %optimizers
    P_dec{i,4} = optimizerOpt(P_dec{i,4});
    
end
    

parfor ii = 1:length(P) %parfor
    
    P_oneline = P_dec(ii,:); % due to parfor
    
    % saving the value for calculating the fitness
    l1 = round(P_oneline{1}); 
    l2 = round(P_oneline{2});
    activation = P_oneline{3};
    optimizer = P_oneline{4};
    params = {l1,l2,activation,optimizer};
    
    % Calculate the fitness value 
    [~, validationRMSE] = MLmodelsFitnet(trainingData, responseData, params);
    
    % assign fit to the matrix
    RSMEs{ii,1} = validationRMSE; % correct order using parfor? probably yes.
    %P_dec{i,length(bit_size)+1} = validationRMSE;
    %P_bi{i,length(bit_size)+1} = validationRMSE;
end

    % Due to parfor RSMEs addad to P after loop
    P = [P, RSMEs];
    P_dec = [P_dec, RSMEs];
    P_bi = [P_bi, RSMEs];
    %P{:,2} = RSMEs;
    %P_dec{:,length(bit_size)+1} = RSMEs;
    %P_bi{:,length(bit_size)+1} = RSMEs; 
    
    % sort P using fit
    P = sortrows(P, 2); % - for descending. RSME: low is better
    P_dec = sortrows(P_dec, length(bit_size)+1);
    P_bi = sortrows(P_bi, length(bit_size)+1); %#ok<NASGU>
    
end



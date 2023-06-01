function SLTPparams = optGenTSLTP(data,fitData)
% data........data: 1=Open, 2=High, 3=Low, 4=Close 
%
% xSL.........
% xTP.........
% ATRstd...... 
% sPeriod.....
%
% SLTPparams = [xSL,xTP,ATRstd,sPeriod,penalty];

%disp("Optimize SLTP");
datasz = size(data(:,1),1);

% define settings for optGen
Popsz = 300; % size of population % use 40
pCross = 0.25; % probability of crossover / preassure of selection. usually around 20-30. high: harder to make changes. low: no/slow converge
pMutate = 0.05;  % probability of mutation per bit
one_or_two_point = 2; % either 1 or 2. stands for one or two point crossover
termination = 1; % if termination 1, the termination condition will be used.

% rng = [0,2; ... % setting the range of each parameter
%        5,3]; % row1: lowest value, row2: highest value
% bit_size = [9,1]; % setting the size of the bits

% for TrailingSLstd
rng = [0.1,0.1; ... % setting the range of each parameter
       7,10]; % row1: lowest value, row2: highest value
bit_size = [7,7]; % setting the size of the bits


P = GgeneratePopTSLTP(Popsz, bit_size, rng, datasz);
[P, P_dec] = GfitnessTSLTP(data, P, bit_size, rng, fitData);

% evolving the population
if sum([P{:,2}]) ~= 0
    for i = 1:20
        P = GrecombineTSLTP(P, pCross, one_or_two_point, bit_size, rng, datasz);
        P = GmutateTSLTP(P, pMutate, bit_size, rng, datasz);
        [P, P_dec] = GfitnessTSLTP(data, P, bit_size, rng, fitData);
        
        disp(i);
        if all([P{1:round(length(P)*0.1),2}] < P{1,2}*1.01) && ...
                all([P{1:round(length(P)*0.1),2}] > P{1,2}*0.99)...
                && termination == 1
            disp("Termination at " + i);
            break % termination condition/Abbruchbedingung
            %makes the code a lot faster but the optimum might not be that good
        end % If the best 1% are equally fit (and not -inf)
    end 
end

% Getting the decimal numbers of the best parameters
xSL = P_dec{1,1};
xTP = P_dec{1,2};


% evaluate panalty = avg profit [%] per day
position = fitData{1};
predictOffset = fitData{2};
xnodes = fitData{3};
RiskManagement = fitData{4};
plotparam = fitData{5};
% penalty = fitData{6}; % only needed in GfitnessSLTP
SLTPparams = {xSL,xTP};

[~,endvalueEC2] = evalTradesOffestTSLTP(data,position,SLTPparams,predictOffset,xnodes,RiskManagement,plotparam);

% Pack results in SLTPparams
SLTPparams = {xSL,xTP,endvalueEC2};


end

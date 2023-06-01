function [endvalueEC2,periodBB,nStd,SL,TP] = optGenBB(data,fitData,datasz)
% data........data: 1=Open, 2=High, 3=Low, 4=Close 
%


% define settings for optGen
Popsz = 500; % size of population % use 40
pCross = 0.25; % probability of crossover / preassure of selection. usually around 20-30. high: harder to make changes. low: no/slow converge
pMutate = 0.05;  % probability of mutation per bit
one_or_two_point = 2; % either 1 or 2. stands for one or two point crossover
termination = 1; % if termination 1, the termination condition will be used.

% periodBB = 1, nStd = 2, SL = 3, TP = 4
rng = [2,0.5,1,1; ... % setting the range of each parameter
       129,3,7,10]; % row1: lowest value, row2: highest value
bit_size = [7,7,9,9]; % setting the size of the bits

P = GgeneratePopBB(Popsz, bit_size, rng, datasz);
[P, P_dec] = GfitnessBB(data, P, bit_size, rng, fitData);

% Insist on getting results
c=0;
while P{1,2} == 0
    P = GgeneratePopBB(Popsz, bit_size, rng, datasz);
    [P, P_dec] = GfitnessBB(data, P, bit_size, rng, fitData);
    c=c+1;
    if c == 5
        endvalueEC2=0;periodBB=50;nStd=1;SL=3;TP=5; % assign any values if there are no results
        return
    end
end

% evolving the population
if sum([P{:,2}]) ~= 0
    for i = 1:100
        P = GrecombineBB(P, pCross, one_or_two_point, bit_size, rng, datasz);
        P = GmutateBB(P, pMutate, bit_size, rng, datasz);
        [P, P_dec] = GfitnessBB(data, P, bit_size, rng, fitData);
        
        %disp(i);
        if all([P{1:round(length(P)*0.1),2}] < P{1,2}*1.01) && ...
                all([P{1:round(length(P)*0.1),2}] > P{1,2}*0.99)...
                && termination == 1
            disp("Termination at " + i);
            break % termination condition/Abbruchbedingung
            %makes the code a lot faster but the optimum might not be that good
        end % If the best 1% are equally fit (and not -inf)
    end 
end

RiskManagement = fitData{1};
plotparam = fitData{2};
BB = fitData{3};

periodBB = P_dec{1,1};
nStd = P_dec{1,2};
SL = P_dec{1,3};
TP = P_dec{1,4};

[~,endvalueEC2] = sBB(data,periodBB,nStd,SL,TP,RiskManagement,plotparam,BB);


end

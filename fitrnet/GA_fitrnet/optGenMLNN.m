function [l1,l2,activation,lambda] = optGenMLNN(Data,datasz)
% Data........Data: 1=Open, 2=High, 3=Low, 4=Close 
%
% l1.........optimal layer 1 size 
% l2.........optimal layer 2 size
% Activation........optimal activation function 
% Lambda........optimal lambda value

%params = {l1, l2, activation, lambda};

% define settings for optGen
Popsz = 40; % size of population % use 40
pCross = 0.25; % probability of crossover / preassure of selection. usually around 20-30. high: harder to make changes. low: no/slow converge
pMutate = 0.05;  % probability of mutation per bit
one_or_two_point = 2; % either 1 or 2. stands for one or two point crossover
termination = 1; % if termination 1, the termination condition will be used.

rng = [1,0,1,0; ... % setting the range of each parameter
    25,25,4,100]; % row1: lowest value, row2: highest value
%emaBitsz = adaptBitsz(max_length); % Adapt Bit size of EMAs
bit_size = [5,5,2,8]; % setting the size of the bits
        

P = GgeneratePopML(Popsz, bit_size, rng, datasz);
[P, P_dec] = GfitnessML(Data, P, bit_size, rng);

% evolving the population
if sum([P{:,2}]) ~= 0
    for i = 1:40
        P = GrecombineML(P, pCross, one_or_two_point, bit_size, rng, datasz);
        P = GmutateML(P, pMutate, bit_size, rng, datasz);
        [P, P_dec] = GfitnessML(Data, P, bit_size, rng);
        
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

% Getting the decimal numbers of the best parameters
l1 = round(P_dec{1,1});
l2 = round(P_dec{1,2});
activation = P_dec{1,3};
lambda = P_dec{1,4}; 

% Create a table to have a better look on the best parameters
%opt = table({'v'; v}, {'period'; period}, {'nstd'; nstd}, {'EMA1'; EMA1},...
%    {'EMA2'; EMA2}, {'EMAsl'; EMAsl}, {'SL'; SL}, {'TP'; TP}); %#ok<NASGU>

end

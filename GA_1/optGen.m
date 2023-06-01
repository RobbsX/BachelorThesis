function [fit,v,period,nstd,EMA1,EMA2,EMAsl,SL,TP] = optGen(Data,max_length,datasz)
% Data........Data: 1=Open, 2=High, 3=Low, 4=Close 
% Popsz.......size of population
% pCross......probability of crossover / preassure of selection. usually around 20-30. high: harder to make changes. low: no/slow converge
% pMutate.....probability of mutation per bit
% one_or_two_point..either 1 or 2. stands for one or two point crossover
% termination.if termination 1, the termination condition will be used.
% bit_size....setting the size of the bits
% rng.........setting the range of each parameter
%
% fit.........optimal value of the profit or loss 
% v...........optimal parameter for the volume of the BB
% period......optimal period used in determining the BB
% EMA1........optimal first EMA used to calculate the MACD
% EMA2........optimal second EMA used to calculate the MACD
% EMAsl.......optimal signal line used for determining the position using MACD
% SL..........optimal Stop Loss
% TP..........optimal Take Profit

% define settings for optGen
Popsz = 300; % size of population
pCross = 0.25; % probability of crossover / preassure of selection. usually around 20-30. high: harder to make changes. low: no/slow converge
pMutate = 0.02;  % probability of mutation per bit
one_or_two_point = 2; % either 1 or 2. stands for one or two point crossover
termination = 1; % if termination 1, the termination condition will be used.

rng = [1,2,1,2,2,2,1,1; ... % setting the range of each parameter
    100,max_length,128,max_length,max_length,max_length,128,128]; % row1: lowest value, row2: highest value
emaBitsz = adaptBitsz(max_length); % Adapt Bit size of EMAs
bit_size = [7,emaBitsz,7,emaBitsz,emaBitsz,emaBitsz,7,7]; % setting the size of the bits
        

P = GgeneratePop(Popsz, bit_size, rng, datasz);
[P, P_dec] = Gfitness(Data, P, bit_size, rng);

% evolving the population
if sum([P{:,2}]) ~= 0
    for i = 1:400
        P = Grecombine(P, pCross, one_or_two_point, bit_size, rng, datasz);
        P = Gmutate(P, pMutate, bit_size, rng, datasz);
        [P, P_dec] = Gfitness(Data, P, bit_size, rng);
        
        if all([P{1:round(length(P)*0.3),2}] < P{1,2}*1.01) && ...
                all([P{1:round(length(P)*0.3),2}] > P{1,2}*0.99)...
                && termination == 1
            break % termination condition/Abbruchbedingung
            %makes the code a lot faster but the optimum might not be that good
        end % If the best 1% are equally fit (and not -inf)
    end 
end

% Getting the decimal numbers of the best parameters
v = P_dec{1,1};
period = round(P_dec{1,2}); % 6 bits
nstd = P_dec{1,3}; % x/40 -> rng from 0 to 3.2
EMA1 = round(P_dec{1,4}); % 6 bits
EMA2 = round(P_dec{1,5}); % 6 bits
EMAsl = round(P_dec{1,6}); % 6 bits
SL = P_dec{1,7}; % 7 bits, x/4 -> 1 to 32
TP = P_dec{1,8}; % 7 bits, x/4 -> 1 to 32
fit = P_dec{1,9};

nstd = nstd/4/10;
SL = SL/4;
TP = TP/4;

% Create a table to have a better look on the best parameters
opt = table({'v'; v}, {'period'; period}, {'nstd'; nstd}, {'EMA1'; EMA1},...
    {'EMA2'; EMA2}, {'EMAsl'; EMAsl}, {'SL'; SL}, {'TP'; TP}); %#ok<NASGU>

end

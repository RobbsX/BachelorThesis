function [P, P_dec] = Gfitness(data, P, bit_size, rng)
% data...........Data: 1=Open, 2=High, 3=Low, 4=Close. Will be used to
% backtest
% P..............Input Population Matrix
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
%
% P..............Output Population Matrix

P_bi = cell(length(P),length(bit_size)+1);
P_dec = cell(length(P),length(bit_size)+1);
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
    
    % saving the value for calculating the fitness
    v = P_dec{i,1};
    period = round(P_dec{i,2}); % 6 bits
    nstd = P_dec{i,3}; % x/40 -> rng from 0 to 3.2
    EMA1 = round(P_dec{i,4}); % 6 bits
    EMA2 = round(P_dec{i,5}); % 6 bits
    EMAsl = round(P_dec{i,6}); % 6 bits
    SL = P_dec{i,7}; % 7 bits, x/4 -> 1 to 32
    TP = P_dec{i,8}; % 7 bits, x/4 -> 1 to 32
    
    nstd = nstd/4/10;
    SL = SL/4;
    TP = TP/4;
    
    % Calculate the fitness value 
    [~, fit] = myECSLTP(data,v,period,nstd,EMA1,EMA2,EMAsl,SL,TP,'conservative',0,1,0);
    
    % assign fit to the matrix
    P{i,2} = fit; 
    P_dec{i,length(bit_size)+1} = fit;
    P_bi{i,length(bit_size)+1} = fit;
end

    % sort P using fit
    P = sortrows(P, -2); % - for descending
    P_dec = sortrows(P_dec, -9);
    P_bi = sortrows(P_bi, -9); %#ok<NASGU>
end



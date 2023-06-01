function [P] = Grecombine(P, pCross, one_or_two_point, bit_size, rng, datasz)
% P...................Population Matrix
% pCross..............probability of crossover / preassure of selection
% one_or_two_point....stands for one or two point crossover
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
% datasz.........length of Data. Only used to check constraints
% 
% P...................Output Population Matrix

% SELECTION 

P_new = cell(length(P),2);

% setting an array of random values
r_seq = rand(length(P),1);

% creating the cumsum of fit values
fit_cumsum = ([P{:,2}] + (abs(P{end,2})))';
fit_cumsum = (fit_cumsum/sum(fit_cumsum(1:end)));
fit_cumsum = cumsum(fit_cumsum);

% random, roulette wheel
for ii = 1:length(P)
    remaining = r_seq(ii,1);
    for jj = 1:length(fit_cumsum)
        remaining = remaining - fit_cumsum(jj);
        if remaining < 0
            P_new{ii,1} = P{jj,1};
            break
        end
    end
end

% CROSSOVER

rCross_seq = rand(length(P),1);
idx = find(rCross_seq<pCross);

% if the length of the random choosen inxdex is odd, randomly...
if mod(length(idx),2) ~= 0
    del_or_add = rand();
    if del_or_add > 0.5
        idx = [idx; (randi([1 500],1,1))]; %#ok<*NASGU> %...add a Chromosome
    else
        idx(randi([1 length(idx)],1,1)) = []; %...delete a Chromosome
    end
end

% choose parents and replace them with their children 
for zz = 1:length(idx)/2
    p = randperm(length(idx), 2);
    p1 = P_new{p(1),1};
    p2 = P_new{p(2),1};
    
    if one_or_two_point == 1
        boolean1 = 0; boolean2 = 0; 
        while boolean1 == 0 || boolean2 == 0
            % One point crossover
            pos = randi([1 length(P_new{1,1})],1,1);
            c1 = [p1(1:pos), p2((pos+1):end)];
            c2 = [p2(1:pos), p1((pos+1):end)];
            boolean1 = GcheckConstr(c1, bit_size, rng, datasz);
            boolean2 = GcheckConstr(c2, bit_size, rng, datasz);
        end
        
        P_new{p(1),1} = c1;
        P_new{p(2),1} = c2;
        
    else % one_or_two_point == 2, or anything else
        boolean1 = 0; boolean2 = 0;
        while boolean1 == 0 || boolean2 == 0
            % Two point crossover
            pos = randi([1 length(P_new{1,1})],1,2);
            pos1 = min(pos); pos2 = max(pos);
            c1 = {[p1(1:pos1), p2((pos1+1):pos2), p1((pos2+1):end)]};
            c2 = {[p2(1:pos1), p1((pos1+1):pos2), p2((pos2+1):end)]};
            boolean1 = GcheckConstr(c1, bit_size, rng, datasz);
            boolean2 = GcheckConstr(c2, bit_size, rng, datasz);
        end
        
        P_new(p(1),1) = c1;
        P_new(p(2),1) = c2;
    end
    
end

% Elitearism
P_new{1,1} = P{1,1};

P = P_new;
end


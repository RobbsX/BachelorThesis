%0 2 4 8 16 32 64 128 256 512
function [P] = GgeneratePopML(Popsz, bit_size, rng, datasz)
% Popsz....Size of the Population - How large will it be?
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
% datasz.........length of Data. Only used to check constraints
% 
% P........Polulation Matrix, Format: Cell, col1: chromosome, col2: 

%bit_size = [5,5,2,8];
bit_sum = sum(bit_size); %34

% build population
P = cell(Popsz,2);
for i = 1:Popsz
    boolean = 0;
    while boolean == 0
        % create a random chromosome
        r = {randi([0 1], 1, bit_sum)};
        boolean = GcheckConstrML(r, bit_size, rng, datasz); %just lambda
    end
    P(i,1) = r;
end

end

% Notes ---------------------
% 500 individuen -> 1 array per individuum
% ein individuum -> 1 array mit total Chromosom length
% Chromosom length =  number of bit of param_1 + num of bit of param_n

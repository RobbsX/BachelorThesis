function [boolean] = GcheckConstrSLTP(chromosome, bit_size, rng, ~)
% chromosome.....chromosome array in a 1x1 cell
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
% datasz.........length of Data. Only used to check constraints
%
% boolean........True=1, Flase=0. True if the chromosome is within the constraints

boolean = 1;
% split the chromosome to read it
param = cell(1, length(bit_size));
ende = 0;
c = 0;
for j = bit_size
    start = ende + 1;
    ende = ende + j; % bit_size(1,j)
    c = c+1; % count for parameter / bit_size
    param{c} = Gdec(chromosome{1}(start:ende), rng(1,c), rng(2,c));
end

% Annahme: xSL < xTP
if param{1} > param{2}
    boolean = 0;
end


end

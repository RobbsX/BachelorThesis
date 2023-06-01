function [boolean] = GcheckConstr(chromosome, bit_size, rng, datasz)
% chromosome.....chromosome array in a 1x1 cell
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
% datasz.........length of Data. Only used to check constraints
%
% boolean........True=1, Flase=0. True if the chromosome is within the constraints

boolean = 0;
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

% check Constraints
% period = 2, EMA1 = 4, EMA2 = 5, EMAsl = 6, SL = 7, TP = 8
if param{6}+1<param{4} && param{4}+1<param{5} && ...
        param{7}+0.25<param{8} && param{2}+1<datasz-param{4} %3+1<19-4.5
    boolean = 1;
end


end
% Constraints:
% read like eg: EMA1 must be smaller than EMA2
% EMA1 < EMA2, EMA1….smaller EMA of MACD, EMA2….bigger EMA of MACD
% EMAsl < EMA1, EMAsl….EMA of the signal line
% SL < TP, TP…Take Profit, SL…Stop Loss
% period < datasz - EMA1
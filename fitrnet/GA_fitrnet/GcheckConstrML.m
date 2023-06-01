function [boolean] = GcheckConstrML(chromosome, bit_size, rng, ~)
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

% check Constraints - boolean = 1 if satisfied
% l1 = 1, l2 = 2, Activation = 3, Lambda = 4
% if param{2}==0
%     if param{3}~=0
%         disp('There was a single layered model!');
%         boolean = 0;
%     end
% end 

% Check Sigmoid --> Lambda == 0
if param{3} == 4
    if param{4}~=0
        boolean = 0;
    end
end

end

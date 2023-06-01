function [num] = Gbit2num(x)
% x......array of 1/0
%
% num....translated array to a integer

bits = length(x);
num = 0;
for i = bits:-1:1
    if x(i) == 1
        num = num + power(2,bits-i);
    end
end

end
function [z] = Gdec(x, rng1, rng2)
% x.....array of bits eg. "1001001"
% rng1..the lower range
% rng2..the upper range
% 
% z.....Output Decimal number which represents the bit array

    z = rng1 + Gbit2num(x)*(rng2-rng1) / (power(2,length(x))-1) ;
    
end
function [P] = GmutateTSLTP(P, pMutate, bit_size, rng, datasz)
% P.........Input Population Matrix
% pMutate...probability of mutation per bit
% bit_size.......length of the bits of each parameter in an array
% rng............2x(length(bitsize)) matrix. upper and lower bound  
% for each parameter
% datasz.........length of Data. Only used to check constraints
% 
% P.........Output Population Matrix

for ii = 1:length(P)
    boolean = 0;
    x = zeros(1,sum(bit_size));
    while boolean == 0
        for jj = 1:length(P{ii,1})
            r = rand();
            if r<pMutate
                x(jj) = ~P{ii,1}(jj);
            else
                x(jj) = P{ii,1}(jj);
            end
        end
        boolean = GcheckConstrTSLTP({x}, bit_size, rng, datasz);
    end
    P{ii,1} = x;
end

end

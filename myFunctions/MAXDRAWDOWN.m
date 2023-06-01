function [MaxDD, MaxDDIndex] = MAXDRAWDOWN(Data)
% Maximum drawdown is defined as the largest drop from a peak to a bottom
% experienced in a certain time period.
% Also works if data is negative. 

[T, N] = size(Data);

MaxDD = zeros(1,N);
MaxDDIndex = ones(2,N);

MaxData = Data(1,:);
MaxIndex = ones(1,N);
for i = 1:T
    MaxData = max(MaxData, Data(i,:));
    q = MaxData == Data(i,:);
    MaxIndex(1,q) = i;
    DD = MaxData - Data(i,:);
    if any(DD > MaxDD)
        p = DD > MaxDD;
        MaxDD(p) = DD(p);
        MaxDDIndex(1,p) = MaxIndex(p);
        MaxDDIndex(2,p) = i;
    end
end

k = MaxDDIndex(1,:) == MaxDDIndex(2,:);
MaxDDIndex(:,k) = NaN;

end

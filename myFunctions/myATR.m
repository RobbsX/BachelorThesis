function [atr] = myATR(Data,period)
% Data............for High, Low & Close
% period..........in days for EMA
%
% ATR.............ATR of the length of Data
% 
% Data(:,1)=Date, Data(:,2)=Open, Data(:,3)=High, Data(:,4)=Low,
% Data(:,5)=Close

observ = length(Data);

% Input data
hi = Data(:,3);
lo = Data(:,4);
cl = Data(:,5);

% True range
h_m_l = hi-lo;                                   % high - low
h_m_c = [0;abs(hi(2:observ)-cl(1:observ-1))];  % abs(high - close)
l_m_c = [0;abs(lo(2:observ)-cl(1:observ-1))];   % abs(low - close)
tr = max([h_m_l,h_m_c,l_m_c],[],2);                 % true range

% Average true range
atr = myEMA(tr,period);

end
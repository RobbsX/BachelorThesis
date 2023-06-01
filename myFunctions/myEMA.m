function EMA = myEMA(data, period, ~)

EMA = movavg(data,'exponential',period);

% old
% if nargin == 2
%     s = 2;
% end
% EMA = zeros(size(data,1),1);
% sf = (s/(1+period));
% 
% for i = period:size(data,1)
%     if i == period
%         EMA(i) = sum(data(1:period))/period;
%     else
%         EMA(i) = data(i)*sf+ EMA(i-1)*(1-sf);
%     end
% end

end
%{
EMA = (value today * (smoothing/(1+days)) 
+ EMA yesterday * (1 - (smoothing/(1+days))
%}
function [nHigh, nLow] = myNdayHighLow(data, period)
% data.....Matrix n*m. data(:,2) = Open, data(:,3) = High, data(:,4) = Low

nHigh = movmax(data(:,2),[period 0]);
nLow = movmax(data(:,3),[period 0]);

end

%% Old
% Define High and Low vectors 
% High = data(:,2);
% Low = data(:,3);
% 
% % Define zeros for the result
% nHigh = zeros(size(data,1),1);
% nLow = zeros(size(data,1),1);
% 
% for ii = period+1:size(data,1)
%     nHigh(ii) = max(High(ii-period:ii-1, 1));
%     nLow(ii) = min(Low(ii-period:ii-1, 1));
% end
% Note: using 'ii-period:ii' it would take today inclusive. eg. today-6days

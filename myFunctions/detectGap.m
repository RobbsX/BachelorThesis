function diff = detectGap(Data, timeRange)
% Data..........Matrix n*m, data(:,1) = dates, data(:,2) = Open, data(:,3) = High, data(:,4) = Low ;
% timeRange.....Data time period in Days, Weeks or Months. Format: String
%
% diff..........Number of days which need to be considered for the WFA

if strcmp(timeRange, 'Weekly')
     t = [5, 9]; 
     normalTime = 7;
elseif strcmp(timeRange, 'Daily')
    t = [1, 5]; 
    normalTime = 1;
elseif strcmp(timeRange, 'Monthly')
    t = [27, 33];
    normalTime = 30;
end

diff = 0;

for ii = 2:size(Data(:,1),1)
    x = Data(ii,1) - Data(ii-1,1);
    if ~(t(1) < x && x < t(2))    
        diff = diff + Data(ii,1)-Data(ii-1,1)-normalTime;
    end 
end
end
function d = cutData(Data, timeRange)
% Data..........Matrix n*m, data(:,1) = dates, data(:,2) = Open, data(:,3) = High, data(:,4) = Low ;
% timeRange.....Data time period in Days, Weeks or Months. Format: String
%
% d.............Data, but 'cutted'. The longest time series which has
% coherent days

if strcmp(timeRange, 'Weekly')
     t = [5, 30]; 
elseif strcmp(timeRange, 'Daily')
    t = [0, 30]; 
elseif strcmp(timeRange, 'Monthly')
    t = [27, 60];
end

c = 1;
range = [];

range(c,1) = 1;
for ii = 2:size(Data(:,1),1)
    x = Data(ii,1) - Data(ii-1,1);
    if ~(t(1) < x && x < t(2))        
        range(c,2) = find(Data(ii,1)==Data(:,1))-1;
        range(c+1,1) = find(Data(ii,1)==Data(:,1));
        c = c + 1;
    end 
end

range(c,2) = find(Data(end,1)==Data(:,1));
range(:,3) = (range(:,2)-range(:,1))/size(Data(:,1),1);
mx = find(max(range(:,3))==range(:,3));
d = Data(range(mx,1):range(mx,2),:);
end


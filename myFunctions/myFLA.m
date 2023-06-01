function wfaPeriod = myFLA(Data, timeRange, nPeriod, nISTPeriod)
% Data.........Closing, Opening, etc. price
% timeRange....Data time period in Days, Weeks or Months. Format: String
% nPeriod......The number of periodes which will be generated
% nISTPeriod...The period which will be the Out of Sample data
% 
% wfaPeriod....A table which contains all dates for each IS and OS time
% frame

% Foward Looking Analysis

% Define Size of Periodes
% if nargin == 1
%     nPeriod = 22;
%     nISTPeriod = 2;
% end

    % Zeros for start- & enddate
    startdate = zeros(nPeriod,1);
    enddate = zeros(nPeriod,1);
    
    startdate1 = Data(1,1);
    enddate1 = Data(size(Data,1));
    
    % Find out small gaps to make the WFA more precise
    diff = detectGap(Data, timeRange);
    
    % Days per period - numdate of Data is used
    dayperperiod = (enddate1 - startdate1 +1 -diff)/nPeriod; 
    
    % Calculate start/end date for each period
    for ii = 1:nPeriod
        % Set start date for each period
        if ii == 1 
            startdate(ii,1) = startdate1;
        else
            startdate(ii,1) = enddate(ii-1,1)+1;
        end

        % Set end date for each period
        if ii == nPeriod 
            enddate(ii,1) = enddate1;
        else
            enddate(ii,1) = (startdate1 - 1) + round(dayperperiod*ii,0);
        end
    end
    
    
    % Summazie into period Table
    period = (1:nPeriod)';
    startdate = datestr(startdate,'yyyy-mm-dd');
    enddate = datestr(enddate,'yyyy-mm-dd');
    periodTable = table(period, startdate,enddate);
    
    % Calculate number of walk
    nWalk = nPeriod - nISTPeriod; 
    
    % Prepare walk-forward period table (wfaPriod)
    walk = (1:nWalk)'; 
    endISTPeriod = walk + nISTPeriod - 1; 
    startIST = datetime(periodTable{walk,2});
    endIST = datetime(periodTable{endISTPeriod,3}); 
    startOST = datetime(periodTable{endISTPeriod + 1,2}); 
    endOST = datetime(periodTable{endISTPeriod + 1,3});
    wfaPeriod = table(walk, startIST, endIST, startOST, endOST); 
    
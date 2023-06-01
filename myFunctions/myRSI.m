function RSI = myRSI(data,R)
% data....close data (similar to MACD)
% R.......The period of R trading days
% 
% RSI.....The RSI value between 0 and 100

RSI = zeros(1,length(data) - R);
goingUp = zeros(1,R);
goingDown = zeros(1,R);



% Use the first R = 14 samples to calculate the first RSI value
% We cannot only use Close(1)
for i = 1:R
    
    % find difference between days
    difference = data(i+1) - data(i);
    
    % if positive change, price goes up, else down
    if difference >= 0
        goingUp(i) = difference;
    else
        goingDown(i) = abs(difference);
    end
end

% take mean of all price increases and decreases
AvgGain = mean(goingUp);
AvgLoss = mean(goingDown);

% if average loss is 0, RSI is 100 else the RSI value is caluclated
if AvgLoss == 0
    RSI(1) = 100;
else
    RS = AvgGain / AvgLoss;
    RSI(1) = 100 - (100/(1+RS));
end
clear goingUp goingDown 


% Calculating the remaining RSI values.
for i = 1+R:length(data)-1
    
    % find difference between days
    difference = data(i+1) - data(i);
    
    % if positive change, price goes up, else down
    if difference >= 0                     
        goingUp = difference;
        goingDown = 0;
    else
        goingDown = abs(difference);
        goingUp = 0;
    end
    
    %%%% Calculate next average gain and loss
    AvgGain = ((AvgGain*(R-1))+goingUp)/R;
    AvgLoss = ((AvgLoss*(R-1))+goingDown)/R;
    
    
    
    % if average loss is 0, RSI is 100 else the RSI value is caluclated
    if AvgLoss == 0
        RSI(i-R) = 100;
    else
        RS = AvgGain / AvgLoss;
        RSI(i+1-R) = 100 - (100/(1+RS)); % calculate ever RSI value for i
    end
end
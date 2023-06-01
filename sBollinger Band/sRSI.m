function [EC,endvalueEC2,tradesAll,tradesAllClean,position] = sRSI(data,v,period,nstd,EMA1,EMA2,EMAsl,SL,TP,RiskManagement,plotparam, BB, RSI)
% Input
% SL ... Stopp integer (%)
% TP ... Take profit integer(%)
% v .... Sensitivity of volatility of the bollinger bands (%) or
% ...... Range of RSI (65:5:85) in case of RSI
% data...Matrix n*m, data(:,1) = dates, data(:,2) = Open, data(:,3) = High, data(:,4) = Low ;
% BB ... Decide weither you want to use the BB Strategy: 0 = No,
% ...... 1 = Yes with simple Strategy, 2 = Yes with complex strategy 
% RSI ...Using the RSI strategy. 0 = No, 1 = Yes

% Output
% EC...............Vector of Equity Curve
% endvalueEC.......Last value of EC
% tradesAll........Full vector of Trades ( length(tradesAll) = length(EC) )
% tradesAllClean...Vector of Trades ( length(tradesAll) = Nr of Trades +1 )
% position.........Vector of positions. Position Long <-> == 1, Position Short <-> == -1
% compTime.........Computing Time ( in sec)

%%% clean data
data(isnan(data(:,5)),:) = [];


data = data(EMA1:end,1:4); % 12 for MACD = value of EMA1
MACD = MACD(EMA1:end,1);
signalline = signalline(EMA1:end,1);
Open = data(:,1);
Close = data(:,4);
v = v/100;

    % RSI Strategy
if RSI == 1
    R = 14;
    RSI = myRSI(data(:,4),R);
    RSI = RSI';
    
    positionRSI = zeros(length(RSI),1);
    positionRSI(RSI>=v) = -1; %zb. 70 / 30
    positionRSI(RSI<(100-v)) = 1;

    % Bollinger Band Strategy
elseif or( BB == 1, BB == 2 )
    [mid,uppr,lowr] = bollinger(Close,period,0,nstd); 
    signboll = zeros(length(data),1);
    
    bolluppr = zeros(length(data),1);
    bolluppr((Open-uppr)>=0)=-1;
    
    bollmid = zeros(length(data),1);
    bollmid((Close-mid)>=0)=1; % 1 = über mid
    bollmid((Close-mid)<0)=-1; % -1 = unter mid
    
    bolllowr = zeros(length(data),1);
    bolllowr((lowr-Open)>=0)=1;
    
    difboll = uppr - lowr;
    mi = min(difboll);
    mx = max(difboll);
    lowvol = mi + ((mx - mi) * v); % v, from 0 to 1

    vol = zeros(length(data),1);
    vol((difboll-lowvol)<0)=1;
    
    
    if BB == 1 % simple BB
        signboll(bolluppr==-1) = -1;
        signboll(bolllowr==1) = 1;
    end 
    
    if BB == 2 % complex BB
        for i = 1:length(bolluppr)
            if vol(i) == 1
                if i <= length(bolluppr) && bolluppr(i) == -1
                    while i <= length(bolluppr) && bollmid(i) == 1 && vol(i) == 1
                        signboll(i) = -1;
                        i = i +1; %#ok<*FXSET>
                    end
                    if i <= length(bolluppr)
                        signboll(i) = -1; % Bc I change position at Close
                    end
                end
                
                if i <= length(bolllowr) && bolllowr(i) == 1
                    while i <= length(bolllowr) && bollmid(i) == -1 && vol(i) == 1
                        signboll(i) = 1;
                        i = i +1;
                    end
                    if i <= length(bolllowr)
                        signboll(i) = 1; % Bc I change position at Close
                    end
                end
            end
        end
    end
    %========== Bollinger Band Test ===========
    %[mid,uppr,lowr] = bollinger(data(:,4));
    %plot(data(:,4));
    %hold on;
    %plot([mid,uppr,lowr]);
    %
    % Uncomment to get a plot of Data + BB
    %==========================================
end


%%% calculate vector position ( 1 or -1 in this strategy )
difference = MACD-signalline;

% Setting seperate positions
position = zeros(length(MACD),1);
position(difference>=0)=1;
position(difference<0)=-1;

if RSI == 1
    position(positionRSI == 1)=1;
    position(positionRSI == -1)=-1;
    
elseif or( BB == 1 , BB == 2) %using simple or complex BB 
    position(signboll == 1)=1;
    position(signboll == -1)=-1;
end




%%% Make trades and get results
[EC,endvalueEC2,tradesAll,tradesAllClean] = evalTrades(data,position,SL,TP,RiskManagement,plotparam);


end



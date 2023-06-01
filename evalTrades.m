function [EC,endvalueEC2,tradesAll,tradesAllClean] = evalTrades(data,position,SL,TP,RiskManagement,plotparam)



%%% Set SL, TP and Open
TP = TP/100;
SL = SL/100;
Open = data(:,2);

%%% Calculate P&L of each Trade
tradesAll = zeros(length(position), 3);
tradesAll(:,1) = data(:,1);
tradeOpen  = Open(1); % opening price for trade calculation
tradeInd  = 1;       %index value for inner for-loop


for ii = 2:length(position)-1
    
    %%%% from long to short
    if or(position(ii) == -1, position(ii) == 0) && ...
            position(ii-1) == 1 && ~and(position(ii) == 0, position(ii-1) == 0)
        dataTrade = data(tradeInd:ii,2:5);
        Trade = 0;
        TradePerc = 0;
        index = 1 ;
        LT = size(dataTrade,1);
        tradeTP = tradeOpen*(1 + TP);
        tradeSL = tradeOpen*(1 - SL);
        while index <= LT
            tmpOpen = dataTrade(index, 1);
            tmpHigh = dataTrade(index, 2);
            tmpLow = dataTrade(index, 3);
            % Open >  StopLoss. StopLoss will be triggered.
            if tmpOpen <= tradeSL
                Trade = tmpOpen-tradeOpen;
                break
                % Open >  TakeProfit. TakeProfit will be triggered.
            elseif  tmpOpen  >=  tradeTP
                Trade = tmpOpen-tradeOpen ;
                break
                % Open ok, but price trigger in the same Day SL and TP.
                % In this case different possibilities :
                %    - random choice between TP and SL ( optionTrades == 0 )
                %    - conservative ( optionTrades == -1 )
                %    - optimistic ( optionTrades == 1)
            elseif and(tmpLow <= tradeSL, tmpHigh >= tradeTP)
                if strcmp(RiskManagement,'neutral')
                    aa = rand;
                    if aa >= 0.5
                        Trade = -SL*tradeOpen;
                    else
                        Trade = TP*tradeOpen;
                    end
                elseif  strcmp(RiskManagement,'conservative')
                    Trade = -SL*tradeOpen;
                elseif strcmp(RiskManagement,'risky')
                    Trade = TP*tradeOpen;
                end
                
                % Open ok, Price trigger the Takeprofit and Low bigger than SL
            elseif and(and(tmpOpen  <  tradeTP,tmpHigh >= tradeTP),tmpLow >=  tradeSL)
                Trade = TP*tradeOpen;
                break;
                
                % Open ok, Price trigger the StopLoss and High lower than TP
            elseif and(and(tmpOpen   > tradeSL ,  tmpLow  <= tradeSL),tmpHigh <= tradeTP)
                Trade =-SL*tradeOpen;
                break
                % Nothing happened, P&L of  actual Trade calculated with next Open
            elseif and(index == LT , Trade == 0)
                Trade =  Open(ii+1)-tradeOpen;
                
            end
            index = index +1;
            
        end
        
        TradePerc = Trade/tradeOpen;
        tradeOpen = Open(ii + 1);
        tradeInd = ii + 1;
        
        %%%% from short to long / Sell short
    elseif or(position(ii) == 1, position(ii) == 0) && ...
            position(ii-1) == -1 && ~and(position(ii) == 0, position(ii-1) == 0)
        dataTrade = data(tradeInd:ii,2:5);
        Trade = 0;
        TradePerc = 0;
        index = 1 ;
        LT = size(dataTrade,1);
        tradeTP = tradeOpen*(1 - TP);
        tradeSL = tradeOpen*(1 + SL);
        while index <= LT
            tmpOpen = dataTrade(index, 1);
            tmpHigh = dataTrade(index, 2);
            tmpLow = dataTrade(index, 3);
            % Open >  StopLoss
            if tmpOpen >= tradeSL
                Trade = tradeOpen - tmpOpen;
                break
                % Open >  TakeProfit
            elseif  tmpOpen  <=  tradeTP
                Trade = tradeOpen - tmpOpen;
                break
                % Open ok, but price trigger in the same Day SL and TP.
                % In this case different possibilities :
                %    - random choice between TP and SL ( optionTrades == 0 )
                %    - conservative ( optionTrades == -1 )
                %    - optimistic ( optionTrades == 1)
            elseif and( tmpLow <= tradeTP,tmpHigh >= tradeSL)
                if strcmp(RiskManagement,'neutral')
                    aa = rand;
                    if aa >= 0.5
                        Trade = -SL*tradeOpen;
                    else
                        Trade = TP*tradeOpen;
                    end
                elseif  strcmp(RiskManagement,'conservative')
                    Trade = -SL*tradeOpen;
                elseif strcmp(RiskManagement,'risky')
                    Trade = TP*tradeOpen;
                end
                
            elseif and(and(tmpOpen  >  tradeTP,tmpHigh <= tradeSL),tmpLow <=  tradeTP)
                Trade = TP*tradeOpen;
                break;
                
            elseif and(and(tmpOpen  <  tradeSL,tmpHigh >= tradeSL),tmpLow >=  tradeTP) %% FEHLER VIDEO ANSCHAUEN
                Trade =-SL*tradeOpen;
                break
                
                % Nothing happened, P&L of  actual Trade calculated with next Open
            elseif and(index == LT , Trade == 0)
                Trade =  tradeOpen - Open(ii+1);
                
            end
            index = index +1;
            
        end
        
        TradePerc = Trade/tradeOpen;
        tradeOpen = Open(ii + 1);
        tradeInd = ii + 1;
        
    else
        Trade = 0;
        TradePerc = 0;
    end
    
    % To open Trades one day after '0' occured in position 
    if position(ii) == 0 && ii ~= length(position)-1
        tradeOpen = Open(ii + 2);
        tradeInd = ii + 2;
    end
    
    tradesAll(ii+1,2) = Trade;
    tradesAll(ii+1,3) = TradePerc;
end

% Delete the last trade if it's on the last day ( length(tradesAll) =
% length(Data) )

tradesAll(end,2)=0; tradesAll(end,3)=0;

aa = find(tradesAll(:,2)~=0); % index where the first trade occurs
if aa ~= 0
    tradesAll(aa(1),2) = 0; % delete first trade
    tradesAll(aa(1),3) = 0; % delete first trade percentage
    tradesAllClean = tradesAll(logical(tradesAll(:,2)),:);
else
    tradesAllClean = [0,0,0];
end

% Calculate EC
EC = cumsum(tradesAll(:,2));
%endvalueEC = EC(end);
EC2 = cumsum(tradesAllClean(:,2));
if isempty(EC2)
    endvalueEC2 = 0;
else
    endvalueEC2 = EC2(end);
end
if plotparam == 2
    benchmark = data(:,5) - data(1,5);
    figure;plot([0;EC],'b-')
    hold on;
    plot(benchmark,'m-');
    title("Equity Curve");
    xlabel("Time");
    ylabel("Profit (absolut)");
    legend('EC','Benchmark')
end



end
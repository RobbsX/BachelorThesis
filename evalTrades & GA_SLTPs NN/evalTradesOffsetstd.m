function [EC,endvalueEC2,tradesAll,tradesAllClean] = evalTradesOffsetstd(data,position,SLTPparams,predictOffset,xnodes,RiskManagement,plotparam)
% This function evaluates the Trades. The SL and TP is calculated using the
% std of the last sPeriod days and a factor. 
% Used by MLmainFitnet.m

if size(data,2) ~= 5
   error("Data must contain Date=1,Open=2,High=3,Low=4,Close=5!");
end

datanew = data(predictOffset+xnodes :predictOffset: end,:); 

% Extract SLTPparams
xSL = SLTPparams{1}/100; % eg. 0,018 
xTP = SLTPparams{2}/100;
sPeriod = SLTPparams{3};

% Calculate the Std using periods
S = movstd(data(:,5),[sPeriod-1 0])./(movavg(data(:,5),'exponential',sPeriod-1));
S = [0;S]; S = S(1:end-1);

%%% Calculate P&L of each Trade
tradesAll = zeros(length(position), 3);
tradesAll(:,1) = datanew(:,1);
% tradeInd  = 1;       %index value for inner for-loop

% To compare ATR with std
% figure; plot(data(:,4), 'b-'); hold on; plot(data(:,4)+atr, 'r-'); plot(data(:,4)-atr, 'r-'); 
    
% predictOffset ~= 1
for ii = 1:length(position)
    %%% For long predict
    if position(ii) == 1
        idx1 = find(data==datanew(ii,1))-predictOffset+1;
        idx2 = find(data==datanew(ii,1));
        dataTrade = data(idx1:idx2, 2:5);
        Strade = S(idx1:idx2, 1);
        tradeOpen  = dataTrade(1,1); % opening price for trade calculation
        Trade = 0;
        index = 1 ;
        LT = size(dataTrade,1);

        tradeSL = dataTrade(:,1).*(1-Strade)*(1-xSL);
        tradeTP = dataTrade(:,1).*(1+Strade)*(1+xTP);
        
        while index <= LT
            tmpOpen = dataTrade(index, 1);
            tmpHigh = dataTrade(index, 2);
            tmpLow = dataTrade(index, 3);
            tmptradeTP = tradeTP(index);
            tmptradeSL = tradeSL(index);
            % Open >  StopLoss. StopLoss will be triggered.
            if tmpOpen <= tmptradeSL
                Trade = tmpOpen-tradeOpen;
                break
                % Open >  TakeProfit. TakeProfit will be triggered.
            elseif  tmpOpen  >=  tmptradeTP
                Trade = tmpOpen-tradeOpen ;
                break
                % Open ok, but price trigger in the same Day SL and TP.
                % In this case different possibilities :
                %    - random choice between TP and SL ( optionTrades == 0 )
                %    - conservative ( optionTrades == -1 )
                %    - optimistic ( optionTrades == 1)
            elseif and(tmpLow <= tmptradeSL, tmpHigh >= tmptradeTP)
                if strcmp(RiskManagement,'neutral')
                    aa = rand;
                    if aa >= 0.5
                        %Trade = -abs(tradeOpen-tmptradeSL);
                        Trade = tmptradeSL-tradeOpen;
                    else
                        %Trade = abs(tradeOpen+tmptradeTP);
                        Trade = tmptradeTP-tradeOpen;
                    end
                elseif  strcmp(RiskManagement,'conservative')
                    Trade = tmptradeSL-tradeOpen;
                elseif strcmp(RiskManagement,'risky')
                    Trade = tmptradeTP-tradeOpen;
                end
                
                % Open ok, Price trigger the Takeprofit and Low bigger than SL
            elseif and(and(tmpOpen  <  tmptradeTP,tmpHigh >= tmptradeTP),tmpLow >= tmptradeSL)
                Trade = tmptradeTP-tradeOpen;
                break;
                
                % Open ok, Price trigger the StopLoss and High lower than TP
            elseif and(and(tmpOpen   > tmptradeSL,tmpLow <= tmptradeSL),tmpHigh <= tmptradeTP)
                Trade = tmptradeSL-tradeOpen;
                break
                % Nothing happened, P&L of  actual Trade calculated with next Open
            elseif and(index == LT , Trade == 0)
                Trade = data(idx2,5)-tradeOpen; % sell at close
                
            end
            index = index +1;
            
        end
        
        TradePerc = Trade/tradeOpen;
%         tradeInd = ii + 1;
        
        %%% For short prediction
    elseif position(ii) == -1
        idx1 = find(data==datanew(ii,1))-predictOffset+1;
        idx2 = find(data==datanew(ii,1));
        dataTrade = data(idx1:idx2, 2:5);
        Strade = S(idx1:idx2, 1);
        tradeOpen  = dataTrade(1,1); % opening price for trade calculation
        Trade = 0;
%         TradePerc = 0;
        index = 1 ;
        LT = size(dataTrade,1);

        tradeSL = dataTrade(:,1).*(1+Strade)*(1+xSL);
        tradeTP = dataTrade(:,1).*(1-Strade)*(1-xTP);
        
        while index <= LT
            tmpOpen = dataTrade(index, 1);
            tmpHigh = dataTrade(index, 2);
            tmpLow = dataTrade(index, 3);
            tmptradeTP = tradeTP(index);
            tmptradeSL = tradeSL(index);
            % Open >  StopLoss
            if tmpOpen >= tmptradeSL %, index ~= 1)
                Trade = tradeOpen - tmpOpen;
                break
                % Open >  TakeProfit
            elseif  tmpOpen  <=  tmptradeTP
                Trade = tradeOpen - tmpOpen;
                break
                % Open ok, but price trigger in the same Day SL and TP.
                % In this case different possibilities :
                %    - random choice between TP and SL ( optionTrades == 0 )
                %    - conservative ( optionTrades == -1 )
                %    - optimistic ( optionTrades == 1)
            elseif and( tmpLow <= tmptradeTP,tmpHigh >= tmptradeSL)
                if strcmp(RiskManagement,'neutral')
                    aa = rand;
                    if aa >= 0.5
                        Trade = tradeOpen-tmptradeSL;
                    else
                        Trade = tradeOpen-tmptradeTP;
                    end
                elseif  strcmp(RiskManagement,'conservative')
                    Trade = tradeOpen-tmptradeSL;
                elseif strcmp(RiskManagement,'risky')
                    Trade = tradeOpen-tmptradeTP;
                end
                
            elseif and(and(tmpOpen  >  tmptradeTP,tmpHigh <= tmptradeSL),tmpLow <=  tmptradeTP)
                Trade = tradeOpen-tmptradeTP;
                break;
                
            elseif and(and(tmpOpen  <  tmptradeSL,tmpHigh >= tmptradeSL),tmpLow >=  tmptradeTP) 
                Trade = tradeOpen-tmptradeSL;
                break
                
                % Nothing happened, P&L of  actual Trade calculated with next Open
            elseif and(index == LT , Trade == 0)
                Trade = tradeOpen - data(idx2,5); % sell at close
            end
            index = index +1;
            
        end
        
        TradePerc = Trade/tradeOpen;
%         tradeInd = ii + 1;
    else
        Trade = 0;
        TradePerc = 0;
    end
    
    tradesAll(ii,2) = Trade;
    tradesAll(ii,3) = TradePerc;
    
end

% Delete the last trade if it's on the last day ( length(tradesAll) =
% length(Data) )

% tradesAll(end,2)=0; % No, bc of offset

aa = find(tradesAll(:,2)~=0); % index where the first trade occurs
if aa ~= 0
    %tradesAll(aa(1),2) = 0; % delete first trade -> not with NN 
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
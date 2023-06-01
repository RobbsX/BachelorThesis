function [EC,endvalueEC2,tradesAll,tradesAllClean,position] = sDouble7(data,periodHL,periodEMA,SL,TP,RiskManagement,plotparam)
% Input
% SL ... Stopp integer (%)
% TP ... Take profit integer(%)
% v .... Sensitivity of volatility of the bollinger bands (%) or
% ...... Range of RSI (65:5:85) in case of RSI
% data...Matrix n*m. data(:,2) = Open, data(:,3) = High, data(:,4) = Low

% Output
% EC...............Vector of Equity Curve
% endvalueEC.......Last value of EC
% tradesAll........Full vector of Trades (length(tradesAll) = length(EC))
% tradesAllClean...matrix of Trades (length(tradesAll) = Nr of Trades +1),
% 1: Date, 2: Trades absolut, 3: Trades percent
% position.........Vector of positions. Position Long <-> == 1, Position Short <-> == -1
% compTime.........Computing Time ( in sec)

%%%%%%%%%%% Define Strategy %%%%%%%%%%%%%

%%% clean data, Define EMA param, calculate EMA and adjust the Data and SL,TP
data(isnan(data(:,5)),:) = [];

%%% Double 7's Strategy Input
EMA = myEMA(data(:,2:5), periodEMA);
[nHigh, nLow] = myNdayHighLow(data(:,2:5), periodHL);

% periodEMA or periodHighLow
cut = max(periodEMA, periodHL);
data = data(cut:end,1:5);

%%% Check if data is empty
if isempty(data)
    EC=0; endvalueEC2=0; tradesAll=0; tradesAllClean=0; position=0;
    return
end

%%% cut EMA, nHigh & nLow as well
EMA = EMA(cut:end);
nHigh = nHigh(cut:end);
nLow = nLow(cut:end);

Open = data(:,2);  
High = data(:,3);
Low = data(:,4);
Close = data(:,5);

position = zeros(size(data,1),1);

for xx = 1:size(data,1)
    %Open Long trade if current Low is lower than n days low & 
    %the market (Close) is above the EMA
    if and( Low(xx)<nLow(xx), Close(xx)>EMA(xx) )
        start = xx; 
        %close the trade if current High is above the n days high
        while and(High(xx) < nHigh(xx), xx<length(High) ) 
            xx = xx+1; %#ok<*FXSET>
        end
        ende = xx-1;
        position(start:ende) = 1;
    
    %Open Short trade if current High is higher than n days high & 
    %the market (Close) is below the EMA
    elseif and( High(xx)>nHigh(xx), Close(xx)<EMA(xx) )
        start = xx;
        %close the trade if current High is above the n days high
        while and( Low(xx)>nLow(xx), xx<length(High) )
            xx = xx+1;
        end
        ende = xx-1;
        position(start:ende) = -1;
    end
end

%%%%%%%%%%%%%%%%% End of Define Strategy %%%%%%%%%%%%%%%%%%%
if plotparam == 2
    %%%% Preparation for Plotting %%%%
    dottsBuy=zeros(length(position),2);
    dottsSell=zeros(length(position),2);
    
    signalBuyLong=zeros(length(position),2);
    signalBuyShort=zeros(length(position),2);
    signalSellLong=zeros(length(position),2);
    signalSellShort=zeros(length(position),2);
    
    for zz = 2:length(position)
        if zz ~= length(position)
            %%% Buy
            if position(zz)==1 && or(position(zz-1)==0,position(zz-1)==-1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyLong(zz,1) = nLow(zz,1); %Y
            end
            if position(zz)==-1 && or(position(zz-1)==0,position(zz-1)==1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyShort(zz,1) = nHigh(zz,1); %Y
            end
            
            %%% Sell
            if or(position(zz)==-1,position(zz)==0) && position(zz-1)==1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellLong(zz,1) = nHigh(zz,1); %Y
            end
            if or(position(zz)==1,position(zz)==0) && position(zz-1)==-1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellShort(zz,1) = nLow(zz,1); %Y
            end
        end
    end
    dottsBuy(dottsBuy==0)=nan;
    dottsSell(dottsSell==0)=nan;
    signalBuyLong(signalBuyLong==0)=nan;
    signalBuyShort(signalBuyShort==0)=nan;
    signalSellLong(signalSellLong==0)=nan;
    signalSellShort(signalSellShort==0)=nan;
    %%%% Preparation End %%%%
    
    %%%%%%%%%% Plot to show & explain the strategy (without SL and TP) %%%%%%%%%
    figure
    part = (1:100); %length(High)
    cndl(High(part,:),Low(part,:),Open(part,:),Close(part,:))
    hold on;
    p1 = plot(nHigh(part,:),'Color','#D539A5','LineWidth',1);
    p2 = plot(nLow(part,:),'Color','#D539A5','LineWidth',1);
    p3 = plot(EMA(part,:),'Color','#0068FF','LineWidth',1);
    
    p4 = plot(dottsBuy(part,1),'go','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p5 = plot(dottsSell(part,1),'ro','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p6 = plot(signalBuyLong(part,1),'g^','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p7 = plot(signalBuyShort(part,1),'gv','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p8 = plot(signalSellLong(part,1),'rv','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p9 = plot(signalSellShort(part,1),'r^','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    
    hold off;
    title("Double7 Strategy");
    xlabel("Time (Daily)");
    ylabel("Close Price");
    legend([p1 p2 p3 p4 p5 p6 p7 p8 p9], "High","Low","EMA","Buy","Sell","Signal Buy Long","Signal Buy Short","Signal Sell Long","Signal Sell Short")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Make trades and get results
[EC,endvalueEC2,tradesAll,tradesAllClean] = evalTrades(data,position,SL,TP,RiskManagement,plotparam);


end
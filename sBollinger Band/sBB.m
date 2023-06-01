function [EC,endvalueEC2,tradesAll,tradesAllClean,position] = sBB(data,periodBB,nStd,SL,TP,RiskManagement,plotparam, BB)
% Input
% SL ... Stopp integer (%)
% TP ... Take profit integer(%)
% v .... Sensitivity of volatility of the bollinger bands (%) 
% data...Matrix n*m, data(:,1) = dates, data(:,2) = Open, data(:,3) = High, data(:,4) = Low ;
% BB ... Decide weither you want to use the BB Strategy: 
% ...... 1 = Yes with simple Strategy, 2 = Yes with complex strategy 

% Output
% EC...............Vector of Equity Curve
% endvalueEC.......Last value of EC
% tradesAll........Full vector of Trades ( length(tradesAll) = length(EC) )
% tradesAllClean...Vector of Trades ( length(tradesAll) = Nr of Trades +1 )
% position.........Vector of positions. Position Long <-> == 1, Position Short <-> == -1
% compTime.........Computing Time ( in sec)

%% clean data
data(isnan(data(:,5)),:) = [];

%% Calculate indicators 
[mid,uppr,lowr] = bollinger(data(:,5),periodBB,0,nStd);

%% Cut data
cut = periodBB;
data = data(cut:end,1:5);
mid = mid(cut:end);
uppr = uppr(cut:end);
lowr = lowr(cut:end);


%% Define Open and Close
Open = data(:,2);  
High = data(:,3);
Low = data(:,4);
Close = data(:,5);

%% Create positions
signboll = zeros(length(data),1);

bolluppr = zeros(length(data),1);
bolluppr((Open-uppr)>=0)=-1;

% make trades until the upper or lower band is reached
bollmid = zeros(length(data),1);
bollmid((Close-mid)>=0)=1; % 1 = über mid
bollmid((Close-mid)<0)=-1; % -1 = unter mid

bolllowr = zeros(length(data),1);
bolllowr((lowr-Open)>=0)=1;


if BB == 11 % simple BB1 - old
    % Setting seperate positions
    position = zeros(length(Close),1);
    position((lowr-Open)>=0)=1;
    position((Open-uppr)>=0)=-1;
end

if BB == 12 % simple BB2 - new
    % Setting seperate positions
    position = zeros(length(Close),1);
    position((lowr-Open)>=0)=-1;
    position((Open-uppr)>=0)=1;
end

if BB == 2 % complex BB
    for i = 1:length(bolluppr)
        if i <= length(bolluppr) && bolluppr(i) == -1
            while i <= length(bolluppr) && bollmid(i) == 1
                signboll(i) = -1;
                i = i +1; %#ok<*FXSET>
            end
            if i <= length(bolluppr)
                signboll(i) = -1; % Bc I change position at Close
            end
        end
        
        if i <= length(bolllowr) && bolllowr(i) == 1
            while i <= length(bolllowr) && bollmid(i) == -1
                signboll(i) = 1;
                i = i +1;
            end
            if i <= length(bolllowr)
                signboll(i) = 1; % Bc I change position at Close
            end
        end
    end
    % Setting seperate positions
    position = zeros(length(Close),1);
%     position(bollmid == 1)=1;
%     position(bollmid == -1)=-1; % new --> to trade until upper/lower bond is reached
    
    position(signboll == 1)=1;
    position(signboll == -1)=-1;
end

%% Prepare for plotting
if plotparam == 2 && BB == 11
    dottsBuy=zeros(length(position),1);
    dottsSell=zeros(length(position),1);
    signalBuyLong=zeros(length(position),1);
    signalBuyShort=zeros(length(position),1);
    signalSellLong=zeros(length(position),1);
    signalSellShort=zeros(length(position),1);
    
    for zz = 2:length(position)
        if zz ~= length(position)
            %%% Buy
            if position(zz)==1 && or(position(zz-1)==0,position(zz-1)==-1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyLong(zz,1) = lowr(zz,1); %Y
            end
            if position(zz)==-1 && or(position(zz-1)==0,position(zz-1)==1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyShort(zz,1) = uppr(zz,1); %Y
            end
            
            %%% Sell
            if or(position(zz)==-1,position(zz)==0) && position(zz-1)==1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellLong(zz,1) = lowr(zz,1); %Y
            end
            if or(position(zz)==1,position(zz)==0) && position(zz-1)==-1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellShort(zz,1) = uppr(zz,1); %Y
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
    
%% Plot BB11
    part = 1:min(length(data(:,1)),100);
    figure;
    cndl(High(part,:),Low(part,:),Open(part,:),Close(part,:))
    hold on;
    p1 = plot([mid,uppr,lowr],'Color','#D539A5','LineWidth',1);
    
    p2 = plot(dottsBuy(part,1),'go','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p3 = plot(dottsSell(part,1),'ro','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p4 = plot(signalBuyLong(part,1),'g^','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p5 = plot(signalBuyShort(part,1),'gv','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p6 = plot(signalSellLong(part,1),'rv','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p7 = plot(signalSellShort(part,1),'r^','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    hold off;
    
    title("BB Simple");
    xlabel("Time (Daily)");
    ylabel("Close Price");
    legend([p1(1) p2 p3 p4 p5 p6 p7], "Bollinger Bands & Mid","Buy","Sell","Signal Buy Long","Signal Buy Short","Signal Sell Long","Signal Sell Short")
end



%% BB 1.2
if plotparam == 2 && BB == 12
    dottsBuy=zeros(length(position),1);
    dottsSell=zeros(length(position),1);
    signalBuyLong=zeros(length(position),1);
    signalBuyShort=zeros(length(position),1);
    signalSellLong=zeros(length(position),1);
    signalSellShort=zeros(length(position),1);
    
    for zz = 2:length(position)
        if zz ~= length(position)
            %%% Buy
            if position(zz)==1 && or(position(zz-1)==0,position(zz-1)==-1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyLong(zz,1) = uppr(zz,1); %Y
            end
            if position(zz)==-1 && or(position(zz-1)==0,position(zz-1)==1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyShort(zz,1) = lowr(zz,1); %Y
            end
            
            %%% Sell
            if or(position(zz)==-1,position(zz)==0) && position(zz-1)==1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellLong(zz,1) = uppr(zz,1); %Y
            end
            if or(position(zz)==1,position(zz)==0) && position(zz-1)==-1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellShort(zz,1) = lowr(zz,1); %Y
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
%% Plot BB12
    part = 1:min(length(data(:,1)),100);
    figure;
    cndl(High(part,:),Low(part,:),Open(part,:),Close(part,:))
    hold on;
    p1 = plot([mid,uppr,lowr],'Color','#D539A5','LineWidth',1);
    
    p2 = plot(dottsBuy(part,1),'go','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p3 = plot(dottsSell(part,1),'ro','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p4 = plot(signalBuyLong(part,1),'g^','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p5 = plot(signalBuyShort(part,1),'gv','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p6 = plot(signalSellLong(part,1),'rv','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p7 = plot(signalSellShort(part,1),'r^','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    hold off;
    
    title("BB Simple 1.2");
    xlabel("Time (Daily)");
    ylabel("Close Price");
    legend([p1(1) p2 p3 p4 p5 p6 p7], "Bollinger Bands & Mid","Buy","Sell","Signal Buy Long","Signal Buy Short","Signal Sell Long","Signal Sell Short")
end


%% BB 2 
if plotparam == 2 && BB == 2
    dottsBuy=zeros(length(position),1);
    dottsSell=zeros(length(position),1);
    signalBuyLong=zeros(length(position),1);
    signalBuyShort=zeros(length(position),1);
    signalSellLong=zeros(length(position),1);
    signalSellShort=zeros(length(position),1);
    
    for zz = 2:length(position)
        if zz ~= length(position)
            %%% Buy
            if position(zz)==1 && or(position(zz-1)==0,position(zz-1)==-1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyLong(zz,1) = lowr(zz,1); %Y
            end
            if position(zz)==-1 && or(position(zz-1)==0,position(zz-1)==1)
                dottsBuy(zz+1,1) = Open(zz+1,1); %Y
                signalBuyShort(zz,1) = uppr(zz,1); %Y
            end
            
            %%% Sell
            if or(position(zz)==-1,position(zz)==0) && position(zz-1)==1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellLong(zz,1) = mid(zz,1); %Y
            end
            if or(position(zz)==1,position(zz)==0) && position(zz-1)==-1
                dottsSell(zz+1,1) = Open(zz+1,1); %Y
                signalSellShort(zz,1) = mid(zz,1); %Y
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
%% Plot BB2
    part = 1:min(length(data(:,1)),100);
    figure;
    cndl(High(part,:),Low(part,:),Open(part,:),Close(part,:))
    hold on;
    p1 = plot([mid,uppr,lowr],'Color','#D539A5','LineWidth',1);
    
    p2 = plot(dottsBuy(part,1),'go','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p3 = plot(dottsSell(part,1),'ro','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p4 = plot(signalBuyLong(part,1),'g^','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p5 = plot(signalBuyShort(part,1),'gv','MarkerSize',10,'MarkerFaceColor',[.3 0.5 0.1]);
    p6 = plot(signalSellLong(part,1),'rv','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    p7 = plot(signalSellShort(part,1),'r^','MarkerSize',10,'MarkerFaceColor',[.7 .1 .2]);
    hold off;
    
    title("BB Complex");
    xlabel("Time (Daily)");
    ylabel("Close Price");
    legend([p1(1) p2 p3 p4 p5 p6 p7], "Bollinger Bands & Mid","Buy","Sell","Signal Buy Long","Signal Buy Short","Signal Sell Long","Signal Sell Short")
end




%%% Make trades and get results
[EC,endvalueEC2,tradesAll,tradesAllClean] = evalTrades(data,position,SL,TP,RiskManagement,plotparam);


end








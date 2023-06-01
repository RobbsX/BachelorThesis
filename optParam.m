function [EC,perioddHL,periodEMA,SL,TP] = optParam(Data) 
% Data......Closing, Opening or ect. price
% BB........BB strategy if 1 or 2
% RSI.......RSI strategy if 1 (if BB and RSI --> RSI will be executed)
%
% ECs.......list of all last values (=P&L) of each equity curve 
% vol.......
% SL_list...list of stop loss used
% TP_list...list of take profit used

% Vectorizing all combinations for better performance?
    vdHL = 3:10; % periodHighLow
    vEMA = 25:25:200; % periodEMA
    vTP = 1:3:16; % TP % 0.25 bis 20
    vSL = 0.25:0.5:4.75; % SL 0.25 bis 20
%     l = length(vdHL)*length(vEMA)*length(vTP)*length(vSL);
%     c = 1;

    % Creating Zeros for performance - parfor possible
    ECs = [];
    perioddHLlist = [];
    periodEMAlist = [];
    TPlist = [];
    SLlist = [];
    for perioddHL = vdHL
        for periodEMA = vEMA
            for TP = vTP
                for SL = vSL
                    [~, lastEC] = sDouble7(Data,perioddHL,periodEMA,SL,TP,'conservative',0);
                    ECs = [ECs; lastEC]; %#ok<*AGROW>
                    perioddHLlist = [perioddHLlist; perioddHL];
                    periodEMAlist = [periodEMAlist; periodEMA];
                    TPlist = [TPlist; TP];
                    SLlist = [SLlist; SL];
                end
            end
        end
    end
    

    merged = horzcat(ECs,perioddHLlist,periodEMAlist,SLlist,TPlist);
    index = find(ECs == max(ECs));
    index = index(1); % falls 2 exakt gleiche Werte in ECs sind
    EC = merged(index,1);
    perioddHL = merged(index,2);
    periodEMA = merged(index, 3);
    SL = merged(index,4);
    TP = merged(index, 5);
    
end

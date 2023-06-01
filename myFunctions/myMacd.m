function [macd, signalline] = myMacd(close, EMA1, EMA2, EMAsl)
% Close........Closing Prive
% EMA1.........Usually 12. Smaller EMA
% EMA2.........Usually 26. Bigger EMA
% EMAsl........Usually 9. EMA for the signal line

% macd.........MACD
% signalline...signalline

% How it should be used as a signal generator:
% If macd is above the signal line, a long signal accurs and vv.

% calculate the MACD
macd = myEMA(close,EMA1,2)-myEMA(close,EMA2,2);
macd(1:EMA2-1) = 0;
% Need to be careful with handling NaN's in the second calculation
idx = isnan(macd);
signalline = [macd(idx); myEMA(macd(~idx),EMAsl,2)];

end
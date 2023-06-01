function [mid,uppr,lowr] = bollinger(cl,varargin)
% middle.....Middle Bollinger Band
% upper......Upper Bollinger Band
% lower......Lower Bollinger Band
%
% cl.........Closing, Opening or etc. price to calculate the BB
% varargin
%   period...the periode of the SMA to calculate nstd
%   weight...weight of the SMA, usually 0.
%   nstd.....standard deviation +/- for the upper and lower band

% Bollinger Bands

% Variable argument input
if isempty(varargin)
    period = 20;
    weight = 0;
    nstd   = 2;
else
    period = varargin{1};
    weight = varargin{2};
    nstd   = varargin{3};
end

% Create output vectors
mid  = nan(size(cl, 1), 1);
uppr = mid;
lowr = mid;

% Create weight vector
wtsvec = ((1:period).^weight) ./ (sum((1:period).^weight));

% Save the original data and remove NaN's from the data to be processed
nnandata = cl(~isnan(cl));

% Calculate middle band moving average using convolution
cmid    = conv(nnandata, wtsvec);
nnanmid = cmid(period:length(nnandata));

% Calculate shift for the upper and lower bands. The shift is a
% moving standard deviation of the data.
mstd = movstd(nnandata,[period-1 0]);
mstd = mstd(period:end);

% Calculate the upper and lower bands
nnanuppr = nnanmid + nstd.*mstd;
nnanlowr = nnanmid - nstd.*mstd;

%disp(['Data: ', num2str(length(cl)), 'period: ', num2str(period)]);
% Return the values
nanVec = nan(period-1,1);
mid(~isnan(cl))  = [nanVec; nnanmid];
uppr(~isnan(cl)) = [nanVec; nnanuppr];
lowr(~isnan(cl)) = [nanVec; nnanlowr];

end
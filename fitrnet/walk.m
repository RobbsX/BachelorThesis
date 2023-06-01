function ModelSet = walk(Data, wfaPeriod, w, predictOffset, xnodes)  % , nPeriod, nISTPeriod

% Get data fitting to the current walk
[Data_IS, ~] = wfaData(Data, wfaPeriod, w, predictOffset+xnodes);

%% In Sample %%

predstepIS = 1:predictOffset:size(Data_IS(:,1),1) -predictOffset-xnodes+1;

% Set zeros
X_IS = zeros(length(predstepIS), xnodes);

% Split Data
for tt = 1:length(predstepIS) % 1:size(Data_IS(:,1),1) -predictOffset-xnodes
    X_IS(tt,:) = Data_IS((1:xnodes) +predstepIS(tt)-1, 5)'; % train data (last day of training is 'today' or 'current day')
end
D_IS = Data_IS(xnodes+predictOffset :predictOffset: size(Data_IS(:,1),1), 5); % Control

% Normalisation x_n = (x-min)/(max-min)
minXIS = min(X_IS(:,1)); %minX = min(X_IS(X_IS(:,1)~=0));
maxXIS = max(X_IS(:,1));
XnormIS = (X_IS-minXIS)/(maxXIS-minXIS) *(0.9999-0.0001)+0.0001;
DnormIS = (D_IS-min(D_IS(D_IS(:,1)~=0)))/(max(D_IS(D_IS(:,1)~=0))-min(D_IS(D_IS(:,1)~=0))) *(0.9999-0.0001)+0.0001;

check = 0; c = 0; maxitaration = 2;
while check == 0
    % Train model
    [trainedModel, ~] = MLopt(XnormIS, DnormIS, 'NN');
    YnormIS = predict(trainedModel.Regression, XnormIS);
    Y_IS = YnormIS * (maxXIS-minXIS) + minXIS;
    stdMdl = std(Y_IS)/mean(X_IS(:,end));
    c = c + 1; disp("c = " + c);
    if stdMdl > 0.03 % above 1 probably is not a decent std!
        check = 1;
    elseif c >= maxitaration
        disp("For " + w + " no better model found.");
        break % no trades will be made
    end
end

%% Save trainedModels to compare the RMSE
ModelSet{1} = trainedModel;


end
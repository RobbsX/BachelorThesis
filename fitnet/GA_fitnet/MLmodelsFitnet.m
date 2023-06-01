function [trainedModel, RSME] = MLmodelsFitnet(trainingData, responseData, params)


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
%inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20'});

%predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20'};
%predictors = inputTable(:, predictorNames);
%response = responseData;

if params{2} == 0
    layers = [params{1}];
else
    layers = [params{1} params{2}];
end

% Neural Net Fit
    % Train a NN with different optimizer
    
    % params{1} = 1:25
    % params{2} = 0:25
    % params{3} = 'relu'  , 'tanh'  , 'none',    'sigmoid'  %activation
    %             'poslin', 'tansig', 'purelin', 'logsig'
    % params{4} = 0:1 % lambda
    % params{5} = 'trainlm', 'trainbr', 'trainbfg', 'trainrp', 'trainscg',...
    % 'traincgb', 'traincgf', 'traincgp', 'trainoss', 'traingdx', ...
    % 'traingdm', 'traingd' %optimizer
    % params = {basisfunction, kernelfunction};
    
    optimizer = params{4};
    net = fitnet(layers, optimizer);
    net = configure(net,trainingData',responseData');
    net.trainParam.showWindow = false;
    
    if strcmp(params{3}, 'relu')
        activationsNNF = 'poslin';
    elseif strcmp(params{3}, 'tanh')
        activationsNNF = 'tansig';
    elseif strcmp(params{3}, 'none')
        activationsNNF = 'purelin';
    elseif strcmp(params{3}, 'sigmoid')
        activationsNNF = 'logsig';
    end
    
    net.layers{1}.transferFcn = activationsNNF;
    if length(layers) == 2
        net.layers{2}.transferFcn = activationsNNF;
    end
    
    net = train(net, trainingData', responseData'); % transposed bc input [I N] and putput [O N]
    Y = net(trainingData');
    RSME = sqrt(perform(net, Y, responseData')); %performance RMSE
    trainedModel = net;


end
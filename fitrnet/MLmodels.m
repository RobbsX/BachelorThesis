function [trainedModel, RSME] = MLmodels(trainingData, responseData, model, params)


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20'};
predictors = inputTable(:, predictorNames);
response = responseData;

if params{2} == 0
    layers = [params{1}];
else
    layers = [params{1} params{2}];
end

if model == "NN"
    % Train a Trilayered Neural Network
    
    % params{1} = 1:25
    % params{2} = 0:25
    % params{3} = 'relu', 'tanh', 'none', 'sigmoid'  %activation
    % params{4} = 0:1 % lambda
    % params = {l1, l2, activation, lambda};
    
    
    mdl = fitrnet(...
        predictors, ...
        response, ...
        'LayerSizes', layers, ...
        'Activations', params{3}, ...
        'Lambda', params{4}, ...
        'IterationLimit', 1000, ...
        'Standardize', true);
    
    % Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
PredictFcn = @(x) predict(mdl, x);
trainedModel.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.Regression = mdl;
%trainedModel.About = 'This struct is a trained model exported from Regression Learner R2021a.';
%trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 20 columns because this model was trained using 20 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Perform cross-validation
partitionedModel = crossval(trainedModel.Regression, 'KFold', 5);

% Compute validation predictions
validationPredictions = kfoldPredict(partitionedModel); %#ok<NASGU>

% Compute validation RMSE
validationRMSE = sqrt(kfoldLoss(partitionedModel, 'LossFun', 'mse'));
RSME = validationRMSE;


elseif model == "NNF" % Neural Net Fit
    % Train a NN with different optimizer
    
    % params{1} = 1:25
    % params{2} = 0:25
    % params{3} = 'relu'  , 'tanh'  , 'none',    'sigmoid'  %activation
    %             'poslin', 'tansig', 'purelin', 'logsig'
    % params{4} = 0:1 % lambda
    % params{5} = 'trainlm', 'trainbr', '' %optimizer
    % params = {basisfunction, kernelfunction};
    
    net = fitnet(layers, 'trainlm');
    net = configure(net,trainingData',response');
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
    net.trainParam.lamda = params{4};
    %optimizer = params{5};
    
    net = train(net, trainingData', response'); % transposed bc input [I N] and putput [O N]
    Y = net(trainingData');
    RSME = sqrt(perform(net, Y, response')); %performance MSE
    trainedModel = net;
end



end
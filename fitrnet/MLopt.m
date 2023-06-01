function [trainedModel, validationRMSE] = MLopt(trainingData, responseData, model)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

Data = {trainingData, responseData};
datasz = size(trainingData(:,1),1);

% GA optimum of GPR und NN
[l1,l2,activation,lambda] = optGenMLNN(Data,datasz);
% [l1,l2,activation,lambda] = optGAVaPS(Data,datasz);
params = {l1, l2, activation, lambda};

% Compare RMSE of test data
[trainedModel, validationRMSE] = MLmodels(trainingData, responseData, model, params);

% save comparision and repeat r=5 times 

end







% Make histograms of layer 1 and 2 using fitrnet

clearvars;
clc;
set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultLineLineWidth', 2);

%% Choose to use Model or not
useModel = 1;


%% Files & User
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Bachelorarbeit"));
% addpath(genpath("/Users/robinkulha/Documents/MATLAB/Mat Data"));

%% Prepare layers
layersAll = [];
funct = 'fitrnet';

if strcmp(funct, 'fitrnet')
    
    %% Prepare / Load struct
    load('MScollectionFitrnet.mat')
    
    %% Loop
    iifiles = fieldnames(MScollectionFitrnet);
    for ii = 1:numel(iifiles)
        ModelSet = MScollectionFitrnet.(iifiles{ii});
        
        layers = zeros(size(ModelSet(:,1),1),2);
        for l = 1:size(cell2mat(ModelSet(:,1)),1)
            if length(ModelSet{l}.Regression.ModelParameters.LayerSizes) == 1
                layers(l,1:2) = [ModelSet{l,1}.Regression.ModelParameters.LayerSizes 0];
            elseif length(ModelSet{l}.Regression.ModelParameters.LayerSizes) == 2
                layers(l,1:2) = [ModelSet{l,1}.Regression.ModelParameters.LayerSizes];
            end
        end
        
        
        
        
        layersAll = [layersAll; layers]; %#ok<AGROW>
        
    end
    
    
    %% Plot
    nbins = 25;
    figure;
    subplot(1,2,1);
    histogram(layersAll(:,1),nbins,'FaceColor','#405A00','Normalization','probability');
    xlabel("Layer 1 Size");
    ylabel("Number of Layers in Size [%]");
    title("Size of 1st Hidden Layer of Fitrnet");
    
    hold on;
    subplot(1,2,2);
    histogram(layersAll(:,2),nbins,'FaceColor','#405A00','Normalization','probability');
    hold off;
    xlabel("Layer 2 Size");
    ylabel("Number of Layers in Size [%]");
    title("Size of 2nd Hidden Layer of Fitrnet");
    
    
else
    %% Prepare / Load struct
    load('MScollectionFitnet.mat')
    
    %% Loop
    iifiles = fieldnames(MScollectionFitnet);
    for ii = 1:numel(iifiles)
        %ii = iifile{1}(end);
        %currentMS = "ModelSetFitr" + i;
        ModelSet = MScollectionFitnet.(iifiles{ii}); 
        
        layers = zeros(size(ModelSet(:,1),1),2);
        for l = 1:size(ModelSet(:,1),1)
            if length(ModelSet{1}.layers)-1 == 1
                layers(l,1:2) = [ModelSet{l}.layers{1}.size 0];
            else
                layers(l,1:2) = [ModelSet{l}.layers{1}.size ModelSet{l}.layers{2}.size];
            end
        end
        
        
        
        
        layersAll = [layersAll; layers]; %#ok<AGROW>
        
    end
    
    
    %% Plot
    nbins = 25;
    figure;
    subplot(1,2,1);
    histogram(layersAll(:,1),nbins,'FaceColor','#2D42C2','Normalization','probability');
    xlabel("Layer 1 Size");
    ylabel("Number of Layers in Size [%]");
    title("Size of 1st Hidden Layer of Fitnet");
    
    hold on;
    subplot(1,2,2);
    histogram(layersAll(:,2),nbins,'FaceColor','#2D42C2','Normalization','probability');
    hold off;
    xlabel("Layer 2 Size");
    ylabel("Number of Layers in Size [%]");
    title("Size of 2nd Hidden Layer of Fitnet");
    
    
    
    
end



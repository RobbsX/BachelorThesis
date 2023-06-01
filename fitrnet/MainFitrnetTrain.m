clearvars;
clc;
warning ('off','all'); 
set(0,'DefaultFigureWindowStyle','docked')

% For optimising / training the model, parfor is used. Because of that,
% some changes were made.

%% Files & User
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Bachelorarbeit"));
addpath(genpath("/Users/robinkulha/Documents/MATLAB/Mat Data"));

timeRange = 'Daily';
%timeRange = 'Weekly';
%timeRange = 'Monthly';
files = dir(fullfile('/Users/robinkulha/Documents/MATLAB/','Mat Data','*',timeRange,'*.mat*'));
files_name = vertcat({files.name})';
files_folder = vertcat({files.folder})';
N = length(files) ;   % total number of files 

%% Prepare / Load struct
% MScollectionFitrnet = struct; % just used for the first time
load('MScollectionFitrnet.mat')
MS = fieldnames(MScollectionFitrnet);

MScollectionFitrnetCell = struct2cell(MScollectionFitrnet);
saveNameCell = MS;

%% For RG from 44 to 110 - not 125
% istep = 44:100;
% istep = 532; % test with Trane Techn. as it is the larges file
istep = 100:106; % 4 shares should take 9h on my macbook

tTotal = tic;
for iallocate = istep
%% Re-assign Struct

%% Get data
thisfile = files_name{iallocate};
filename = strrep(thisfile(1:end-5-length(timeRange)),"_"," ");
% Using Mat files
Data = load(fullfile(files_folder{iallocate}, thisfile)); 
if strcmp(timeRange, 'Daily')
    Data = Data.data;
else
    Data = Data.(['data', timeRange]);
end
%Data(2) = Open, Data(3) = High, Data(4) = Low, Data(5) = Close

%% Prepare Data
Data = cutData(Data, timeRange); % make ddata coherent, without gaps
Data = Data(:, 1:5);


%% Define parameters for training
xnodes = 20; % input nodes, CE best result at 6
ynodes = 1; % output nodes, must be 1
predictOffset = 5; % x days from the day the prediction used the last data point 


nPeriod = 100; % Used: 100 & 5. Test: 5 & 2.
nISTPeriod = 5;

if size(Data(:,1),1) < 3000
    disp('Too small Data size. This set will be skipped');
    %istep = [istep, istep(end)+1];
    continue
end

%% Prepare and/or load the ModelSet
ModelSet = cell(nPeriod-nISTPeriod,1);

%% Get a table of FLA dates - IS and OS - from start to end
wfaPeriod = myFLA(Data, timeRange, nPeriod, nISTPeriod);

tShare = tic;
time_walk = zeros(height(wfaPeriod),1); % save times of walks
%% Forward Looking Analysis
parfor w = 1:height(wfaPeriod)
    tWalk = tic;
    ModelSetWalk = walk(Data, wfaPeriod, w, predictOffset, xnodes); %, nPeriod, nISTPeriod
    % ModelSetWalk = {rand(12000,1000)}; % test run with rand cell
    ModelSet(w) = ModelSetWalk;
    
    time_walk(w) = toc(tWalk);
    disp("Walk "+w+", Time: "+time_walk(w));
end

time_share = toc(tShare);
disp("Share Nr. " + iallocate + ", Time so far: " + time_share);

%% Save results of share i to summarycell - not because of parfor

saveName = char("ModelSetFitr" + iallocate);
saveNameCell{iallocate,1} = saveName;
MScollectionFitrnetCell{iallocate,1} = ModelSet;

disp("Share Number " + iallocate);
end

time_total = toc(tTotal);

%% Delete empty cells
MScollectionFitrnetCell = MScollectionFitrnetCell(cellfun(@isempty,MScollectionFitrnetCell) == 0);
saveNameCell = saveNameCell(cellfun(@isempty,saveNameCell) == 0);

%% Convert to struct
MScollectionFitrnet = cell2struct(MScollectionFitrnetCell(:,1), saveNameCell, 1);




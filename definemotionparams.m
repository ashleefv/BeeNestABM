clear all
sensitivity = 0; 
% if 0, just run the stddev calculations 
% if 1, run the sensitivity and the static stddev calculations for constant
% coefficients
% if 2, run through a range of coefficient values rather than true
% sensitivity

%% Parameter estimation function for bee pesticide project
% Uses the "simulationOuput" function for f and manipulates 
% coefficients x of the model until the error between the model Y predictions
% f(x) for activity, distanceToNestMates, and distanceToBrood (technically all
% the nest structures) and the Ydata is minimized. The initial guess x0 is 
% called coefficientsguess, and the source data are Ydata. 
%
% Input: read through .mat files:
%   colony data set through allDataCol*.mat
%
%   estimatedData: bump transition probabilities calculated directly for
%       various data sets
%
% Output: 
%   coefficients = [lambdaBrood, lambdaFullFood,
%       lambdaEmptyFood,environmentalStimuliWeight];
%% Take data from Colony 1 all bees pre-exposure tracked for 3600 seconds at 2 frames/second
% Initial state of the agents
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
brood = colony.broodPre;
nestData = colony.preNest; %Currently using data from untreated condition
nestData = calculateActivityMatrix(nestData);

tags = colony.orTagTreat; %aka treatmentList
taglist = colony.tags; %Read out colony's list of marked bees
% number of time points
totalTimePoints = 600; % 5min of simulation time % should take about 45 min if totalTimePoints = 7200
vis = 0;
numBees = size(colony.preNest,2);

nestSimulationData = zeros(totalTimePoints,numBees,5); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);

exposure_state = 'pre'; 


LB = [0,0,0,0]; % lower bound for coefficients
UB = [1,1,1,1]; % upper bound for coefficients

%Calculate means for empirical data
[empiricalMeans ~] = calculateSummaryStatistics(nestData,brood);

% Repack the empirical values into a vector for optimization
Ydata = [empiricalMeans.activity; empiricalMeans.distanceToNestmates; ...
    empiricalMeans.distanceToBrood; empiricalMeans.distanceToCenter; ...
    empiricalMeans.porTimeOnNest];

%% define motion parameters
% used from inSilicoExp committed on 10-10-17
% values from optimization == 1 case
%empirically derived parameters combined for all groups pre-exposure
params = readtable('parameterEstimates.csv');

%For setting motion parameters, are these being drawn from the empirical values (specific to
    %individual and colony), or by reshuffling the distribution of values by cohort? (non-empirical)?
    empirical = 1;
    pre = 1;
    motionParamsCurrent = generateMotionParamsObject(params, colonyNumber, taglist, pre, empirical);
    save('motionsParamspre.mat','motionParamsCurrent')
    
    pre = 0;
    motionParamsCurrent = generateMotionParamsObject(params, colonyNumber, taglist, pre, empirical);
    save('motionsParamspost.mat','motionParamsCurrent')

%% Script to run the model for in silico simulation cases
clear all
close all
format long e

%% Inputs
% Nest chamber size
nestMaxX = 0.2; % m
nestMaxY = 0.2; % m

% Number of bees
numBees = 111;

% Number of repeated simulations with the same conditions
numTrials = 2;  

% Fraction of bees treated
fractTreated = 0.5;

% Nest structures
% fraction of domain covered by nest structures, put 
% uniformly distributed spacingNestStructures apart within that central zone
fractNestStructures = 0.5;
% spacing of nest structure objects
spacingNestStructures = 0.01; % m

% Attraction weights
environmentalStimuliWeightUntreated = 0.06675; % untreated, fitted value = 0.06675
environmentalStimuliWeightTreated = 0.0519; % treated, fitted value = 0.0519

% number of time points
totalTimePoints = 600; % 5 min of simulation time at 2 Hz

% view plots: 1 for turning plotting on, 0 for plotting off
vis = 0;

%% Calculate quantities based on the inputs

%% Make list of tags 
% based on numBees and fractTreated where the number of treated bees is 
% rounded down to the nearest integer
orTagTreat = zeros(numBees,1); %aka treatmentList
numBeesTreated = floor(numBees*fractTreated);
orTagTreat(1:numBeesTreated) = 3;
% recalculated fraction treated based on the integer numbers of bees treated and untreated
actualFractTreated = numBeesTreated/numBees; 

%% Initializae simulation output storage matrix
nestSimulationData = zeros(totalTimePoints,numBees,7); 

%% Multiple colony averaged velocity distribution
% Overwrites speed distribution with the estimate obtained by combining
% data across all colonies (instead of a single colony)
% averages specified [ ] colonies pre-exposure velocity distributions
estimatedData.PDF_Dist_Speed_Log = compute_speed_dist_allColonies([1 2 3 4],0); 

%% Multiple colony averaged activity probability distribution
% as in simulationOutputSpatial, only the pre-exposure values are used for
% the initial activities
estimatedData.Activity_Prob_Dist_Pre = compute_Activity_Prob_Dist_allColonies([1 2 3 4],0);

%% Initial positions
% random x & y distribution throughout the nest
estimatedData.X_pos_initial = nestMaxX*rand(numBees,1);
nestSimulationData(1,:,1) = estimatedData.X_pos_initial; % x position
estimatedData.Y_pos_initial =  nestMaxY*rand(numBees,1);
nestSimulationData(1,:,2) = estimatedData.Y_pos_initial; % y position

%% Using values set above, this section properly formats the coefficients. 
% DO NOT CHANGE THE NEXT TWO LINES
paramcase = 4; % simultaneously change the attraction weights for the untreated and treated cohorts
coefficients = [environmentalStimuliWeightUntreated environmentalStimuliWeightTreated];

%% Take nest dimensions into nestSize vector to be passed into simulationOutputSpatial
nestSize = [nestMaxX nestMaxY];

%% Create new brood structure
centerNest = [nestMaxX/2 nestMaxY/2];
nestStructuresCoords = [centerNest.*(1 - fractNestStructures);... % min coords
                        centerNest.*(1 + fractNestStructures)]; % max coords
numX = (nestStructuresCoords(2,1)- nestStructuresCoords(1,1))./spacingNestStructures+1;
numY = (nestStructuresCoords(2,2)- nestStructuresCoords(1,2))./spacingNestStructures+1;
totalnum = round(numX*numY);
brood = zeros(totalnum,3);
brood(:,3) = '1'; % label is brood
x = nestStructuresCoords(1,1); % xmin
y = nestStructuresCoords(1,2); % ymin
for j = 1:numY
    for i = 1:numX
        index = round(i+(j-1)*numX);
        brood(index,1) = x;
        brood(index,2) = y;
        x = x+spacingNestStructures;
        if x > nestStructuresCoords(2,1) % xmax
            x = nestStructuresCoords(1,1); % xmin
        end
    end
    y = y+spacingNestStructures;   
end
%plot(brood(:,1),brood(:,2),'o') % to view the new brood object

%% define motion parameters
% motionParamsCurrent=;

% Option to load empirical estimates
%load(['motionsParams', 'pre', '.mat'],'motionParamsCurrent');

%Generate 
params = readtable('parameterEstimates.csv');
motionParams = generateMotionParamsObject_general(params, orTagTreat);

%% Load colony structure for passing colony specific info to simulationOutputSpatial
colony.brood = brood;
colony.numBees = numBees;
colony.orTagTreat = orTagTreat;

%% Exposure state 
% There is no exposure state needed for further calculations. The
% parameters in this file set up everything needed without using any of the
% values hard-coded in rules.m
exposure_state = '';
%% Run model
for i = 1:numTrials
    nestSimulationData = simulationOutputSpatial(colony, estimatedData, exposure_state, totalTimePoints,motionParams, vis,coefficients,paramcase,nestSize);

    %Post-process a single trial
    [means distributions]=calculateSummaryStatistics(nestSimulationData,brood);
    activity = means.activity;
    distanceToNestmates = means.distanceToNestmates;
    distanceToBrood = means.distanceToBrood;
    distanceToCenter = means.distanceToCenter;
    porTimeOnNest = means.porTimeOnNest;
    interactionRate = means.interactionRate;
    Ycalc = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest; interactionRate];
    YcalcStochastic(:,i) = Ycalc;
end

%% Post-process all the trials
YcalcStochasticMean = mean(YcalcStochastic,2)
clf
figure(1)
hold on
if numTrials>1
    YcalcStochasticStddev = std( YcalcStochastic')'
    errorbar(YcalcStochasticMean,YcalcStochasticStddev,'bo')
else
    plot(YcalcStochasticMean,'o')
end
axis([0 7 0 1.5])
legend(['simulation mean & std dev results for ' num2str(numTrials) ' runs'])
hold off
set(gca,'XTickLabel',{' ','activity','distanceToNestmates','distanceToBrood','distanceToCenter','porTimeOnNest','interactionRate',' '})

function coefficients = optimizerShellSpatial
clear all
close all
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
% number of time points
totalTimePoints = 600; % 5min of simulation time % should take about 45 min if totalTimePoints = 7200
vis = 0;
numBees = size(colony.preNest,2);

nestSimulationData = zeros(totalTimePoints,numBees,6); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);

exposure_state = 'pre'; 

startcoefficientsguess(1) = 0.1; %COHORTenvironmentalStimuliWeight
startcoefficientsguess(2) = 0.1; %COHORTlambdaBrood 
startcoefficientsguess(3) = 0.1; %COHORTlambdaFullFood 
startcoefficientsguess(4) = 0.1; %COHORTlambdaEmptyFood

LB = [0,0,0,0]; % lower bound for coefficients
UB = [1,1,1,1]; % upper bound for coefficients

%Calculate means for empirical data
[empiricalMeans ~] = calculateSummaryStatistics(nestData,brood);

% Repack the empirical values into a vector for optimization
Ydata = [empiricalMeans.activity; empiricalMeans.distanceToNestmates; ...
    empiricalMeans.distanceToBrood; empiricalMeans.distanceToCenter; ...
    empiricalMeans.porTimeOnNest];
fval = 1;
fvalthreshold = 0.012; %estimated by inspection of typical min f(x) with starting guess of coef = 0.1

%% define motion parameters
% used from inSilicoExp_Working committed on 9-22-17
% values from optimization == 1 case
%empirically derived parameters combined for all groups pre-exposure

%On nest parameters
    AIU1 = .01272; %Active to inactive unbumped, group 1
    IAU1 = .0325; %Inactive to active, unbumped, group 1...
    AIB1 = .01003; %Active to inactive, bumped, group 1...
    IAB1 = .3728;
    
    motionParamsOnNest = [AIU1 AIU1 AIU1 AIU1; IAU1 IAU1 IAU1 IAU1; AIB1 AIB1 AIB1 AIB1; IAB1 IAB1 IAB1 IAB1];
    
    %Off nest parameters
    AIU1 = .02959; %Active to inactive unbumped
    IAU1 = .2022; %Inactive to active, unbumped, group 1...
    AIB1 = .01749; %Active to inactive, bumped, group 1...
    IAB1 = .35560;
    
    %For "full" model
    motionParamsOffNest = [AIU1 AIU1 AIU1 AIU1; IAU1 IAU1 IAU1 IAU1; AIB1 AIB1 AIB1 AIB1; IAB1 IAB1 IAB1 IAB1];
    
    %Combine on and off-nest parameters sets
    motionParamsFull = cat(3,motionParamsOnNest, motionParamsOffNest);
    motionParamsCurrent = motionParamsFull;
%% Minimize
options = optimset('Display','iter','PlotFcns',@optimplotfval,'tolfun',2e-3,'tolx',1e-2); %both tolfun & tolx must be satisfied
% tolfun estimated from the natural stochastic standard deviation of the
% output from 100 repeated trials with the same coefficients:
% outputStochasticStddev from localSensitivity
bestcoef = startcoefficientsguess;
bestfval = 1;
numTrials = 100;
store_coefficients = zeros(numTrials,4);
for i = 1:numTrials
    i
    coefficientsguess = startcoefficientsguess;
    fval = 1;
    % tolx set because to have 2 sig figs of accuracy with the coefficient values
%     while fval > fvalthreshold 
        % coefficients = fmincon(@(coefficients) objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, vis),coefficientsguess,'','','','',LB,UB,'',options) % default function tolerance is 1e-6
        [coefficients,fval,exitflag,output] = fminsearch(@(coefficients) objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, motionParamsCurrent,vis),coefficientsguess,options)
        coefficientsguess = coefficients;
%     end
    store_coefficients(i,:) = coefficients;
    store_fval(i) = fval;
    if fval < bestfval
        bestcoef = coefficients;
        bestfval = fval;
    end
end
store_fval
bestfval
fvalstddev = std(store_fval)
bestcoef
[y,j] = min(store_fval)
coefficients(j,:)
meancoef = mean(coefficients,1)
stddevcoef = std(coefficients,0,1)

    function output = objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, motionParamsCurrent,vis)
    	nestSimulationData = simulationOutputSpatial(colony, estimatedData, exposure_state, totalTimePoints,motionParamsCurrent, vis,coefficients);
        [means distributions]=calculateSummaryStatistics(nestSimulationData,brood);
        activity = means.activity;
        distanceToNestmates = means.distanceToNestmates;
        distanceToBrood = means.distanceToBrood;
        distanceToCenter = means.distanceToCenter;
        porTimeOnNest = means.porTimeOnNest;
        Ycalc = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest];
        output = sum(((Ycalc-Ydata)).^2);%normalize to get all on the same order
    end
save('spatialOptimizer.mat')
end


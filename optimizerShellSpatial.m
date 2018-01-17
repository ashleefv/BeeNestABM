function coefficients = optimizerShellSpatial
clear all
close all
format long e
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
exposure_state = 'post'; 
if strcmp(exposure_state,'pre')
    brood = colony.broodPre;
    nestData = colony.preNest; %using data from untreated condition
    numBees = size(colony.preNest,2);

elseif strcmp(exposure_state,'post')
    brood = colony.broodPost;
    nestData = colony.postNest; %using data from treated condition
    numBees = size(colony.postNest,2);
end
nestData = calculateActivityMatrix(nestData);

tags = colony.orTagTreat; %aka treatmentList
taglist = colony.tags; %Read out colony's list of marked bees
% number of time points
totalTimePoints = 600; % 5min of simulation time % should take about 45 min if totalTimePoints = 7200
vis = 0;

nestSimulationData = zeros(totalTimePoints,numBees,5); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);


% set coefficients
%  coefficients = [0.1 0.1 0.1 0.1];
 paramcase = 3; 
% case 0 is where the COHORTenvironmentalStimuliWeight is the only
% parameter to estimate for the pre exposure case parameter estimation;
% includes cases A, B, C, & D
% case 1 is four coefficients for the pre exposure case parameter
% estimation where brood, full and empty food each get a weighting along
% with nestmates
% case 2 is two coefficients for the post exposure case parameter estimation
% case 3 is one coefficient for the post exposure case parameter estimation
%OPTION USED IN RESULTS IN EMPIRICAL PAPER: PARAMCASE 3 CASE B
    % cases A & B: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics, pre-exposure
    % cases C & D: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics EXCEPT interactionRate, pre-exposure 
    % cases A & C: pre-exposure estimate coef 4 and 5 distinctly with
    % paramcase 2
    % cases B & D: post-exposure estimate coef 5 distinctly and set 
    % coef 4 = to control/untreated with paramcase 3
if paramcase == 0
    % case 0; pre-exposure Cases A, B, C, & D 
    startcoefficientsguess(1) = 0.06675; %COHORTenvironmentalStimuliWeight
    LB = [0]; % lower bound for coefficients
    UB = [1]; % upper bound for coefficients
elseif paramcase == 1 
    % case 1 is four coefficients for the pre exposure case
    startcoefficientsguess(1) = 1; %COHORTenvironmentalStimuliWeight
    startcoefficientsguess(2) = 1; %COHORTlambdaBrood 
    startcoefficientsguess(3) = 0; %COHORTlambdaFullFood 
    startcoefficientsguess(4) = 0; %COHORTlambdaEmptyFood
    LB = [0,0,0,0]; % lower bound for coefficients
    UB = [1,1,1,1]; % upper bound for coefficients
elseif paramcase == 2
    % case A & C: pre-exposure estimate
    % coef 4 and 5 distinctly
    startcoefficientsguess(1) =  0.06675*.75; %COHORTenvironmentalStimuliWeight(4)
    startcoefficientsguess(2) =  0.06675*.5; %COHORTenvironmentalStimuliWeight(5)
    LB = [0,0]; % lower bound for coefficients
    UB = [1,1]; % upper bound for coefficients
elseif paramcase == 3
    % case B & D & E: post-exposure estimate
    % coef 5 distinctly and set coef 4 = to control/untreated
    startcoefficientsguess(1) =  0.055; %COHORTenvironmentalStimuliWeight(5) %  
    LB = [0]; % lower bound for coefficients
    UB = [1]; % upper bound for coefficients
end

%%% put into objective function
    % case A & B: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics, pre-exposure
    % case C & D: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics EXCEPT interactionRate, pre-exposure 
    
%Calculate means for empirical data
if paramcase == 3
    % Case E
    [empiricalMeans empiricalDistributions] = calculateSummaryStatistics(nestData(:,tags==3,:),brood);
    % Cases B & D
%     [empiricalMeans empiricalDistributions] = calculateSummaryStatistics(nestData,brood);
else
    [empiricalMeans empiricalDistributions] = calculateSummaryStatistics(nestData,brood);
end
% Repack the empirical values into a vector for optimization
% cases A & B
Ydata = [empiricalMeans.activity; empiricalMeans.distanceToNestmates; ...
    empiricalMeans.distanceToBrood; empiricalMeans.distanceToCenter; ...
    empiricalMeans.porTimeOnNest; empiricalMeans.interactionRate];
Ydatastddev = [nanstd(empiricalDistributions.activity); nanstd(empiricalDistributions.distanceToNestmates); ...
    nanstd(empiricalDistributions.distanceToBrood); nanstd(empiricalDistributions.distanceToCenter); ...
    nanstd(empiricalDistributions.porTimeOnNest); nanstd(empiricalDistributions.interactionRate)];
% cases C & D
% Ydata = [empiricalMeans.activity; empiricalMeans.distanceToNestmates; ...
%     empiricalMeans.distanceToBrood; empiricalMeans.distanceToCenter; ...
%     empiricalMeans.porTimeOnNest];
% Ydatastddev = [nanstd(empiricalDistributions.activity); nanstd(empiricalDistributions.distanceToNestmates); ...
%     nanstd(empiricalDistributions.distanceToBrood); nanstd(empiricalDistributions.distanceToCenter); ...
%     nanstd(empiricalDistributions.porTimeOnNest)];

% YdatastddevNormalized = [nanstd(empiricalDistributions.activity./empiricalMeans.activity);...
%     nanstd(empiricalDistributions.distanceToNestmates./empiricalMeans.distanceToNestmates); ...
%     nanstd(empiricalDistributions.distanceToBrood./empiricalMeans.distanceToBrood); ...
%     nanstd(empiricalDistributions.distanceToCenter./empiricalMeans.distanceToCenter); ...
%     nanstd(empiricalDistributions.porTimeOnNest./empiricalMeans.porTimeOnNest); ...
%     nanstd(empiricalDistributions.interactionRate./empiricalMeans.interactionRate)]
% Case E % ydata is difference between pre and postexposure means and
% ydatastddev is the std deviation of the difference
% Ydata = [empiricalMeanspre.activity; empiricalMeanspre.distanceToNestmates; ...
%     empiricalMeanspre.distanceToBrood; empiricalMeanspre.distanceToCenter; ...
%     empiricalMeanspre.porTimeOnNest; empiricalMeanspre.interactionRate]-...
%     [empiricalMeanspost.activity; empiricalMeanspost.distanceToNestmates; ...
%     empiricalMeanspost.distanceToBrood; empiricalMeanspost.distanceToCenter; ...
%     empiricalMeanspost.porTimeOnNest; empiricalMeanspost.interactionRate];
% Ydatastddev = [nanstd(empiricalDistributionspre.activity); nanstd(empiricalDistributionspre.distanceToNestmates); ...
%     nanstd(empiricalDistributionspre.distanceToBrood); nanstd(empiricalDistributionspre.distanceToCenter); ...
%     nanstd(empiricalDistributionspre.porTimeOnNest); nanstd(empiricalDistributionspre.interactionRate)];
% sigma = CI95percent.LB;
% sigma(3) = Ydatastddev(3);
% sigma = sigma';
sigma = Ydatastddev;
fval = 5;
fvalthreshold = 1; % estimated by inspection of typical min f(x) with 
%starting guess of coef = 0.06; min outputStochastic

% define motion parameters
load(['motionsParams', exposure_state, '.mat'],'motionParamsCurrent');

%% Minimize
options = optimset('Display','iter','PlotFcns',@optimplotfval,'tolfun',1e-3,'tolx',1e-3); %both tolfun & tolx must be satisfied
% tolfun estimated as the same order of magnitude as the natural stochastic standard deviation of the
% output from 100 repeated trials with the same coefficients:
% outputStochasticStddev from localSensitivity = 0.5 for
% environmentalStimuliWeight = 0.06
bestcoef = startcoefficientsguess;
bestfval = fval;
numTrials = 1;
store_coefficients = zeros(numTrials,length(startcoefficientsguess));
for i = 1:numTrials
    i
    coefficientsguess = startcoefficientsguess;
    fval = 5;
    % tolx set because to have 2 sig figs of accuracy with the coefficient values
%     while fval > fvalthreshold 
        % coefficients = fmincon(@(coefficients) objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, vis),coefficientsguess,'','','','',LB,UB,'',options) % default function tolerance is 1e-6
        [coefficients,fval,exitflag,output] = fminsearch(@(coefficients) objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, motionParamsCurrent,vis,paramcase),coefficientsguess,options)
        coefficientsguess = coefficients;
%    end
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

    function output = objectiveFunction(coefficients, colony, estimatedData, exposure_state, totalTimePoints, motionParamsCurrent,vis,paramcase)
    	nestSimulationData = simulationOutputSpatial(colony, estimatedData, exposure_state, totalTimePoints,motionParamsCurrent, vis,coefficients,paramcase);
        if paramcase == 3
            % Case E
            [means distributions]=calculateSummaryStatistics(nestSimulationData(:,tags==3,:),brood);
            % Cases B & D
%           [means distributions]=calculateSummaryStatistics(nestSimulationData,brood);
        else
            [means distributions]=calculateSummaryStatistics(nestSimulationData,brood);
        end

        activity = means.activity;
        distanceToNestmates = means.distanceToNestmates;
        distanceToBrood = means.distanceToBrood;
        distanceToCenter = means.distanceToCenter;
        porTimeOnNest = means.porTimeOnNest;
        interactionRate = means.interactionRate;
        % cases A & B
         Ycalc = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest; interactionRate];
        % cases C & D
%         Ycalc = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest];
%         sqdiffs = ((Ycalc-Ydata)).^2
%         sqpercentdiffs = ((Ycalc-Ydata)./Ydata).^2
        sqdiffs_over_sigmasq = ((Ycalc-Ydata)./sigma).^2;
%         sqpercentdiffs_over_sigmasq = ((Ycalc-Ydata)./Ydata./sigma).^2
        output = sum(sqdiffs_over_sigmasq);%weighted to by data standard deviation
    end
% save('pre-exposure_spatialOptimizer.mat')
end


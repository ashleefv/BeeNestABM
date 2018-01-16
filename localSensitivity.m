%% local Sensitivity function for bee pesticide project
% Local sensitivity of optimizerShellSpatial's objective function 
% (the individual terms of calculateSummaryStatistics) with respect to the
% coefficients
% Uses the "simulationOuput" function. Runs the model for numTrials number 
% of times for std dev calculations. For sensitivity = 0, the std dev calcs
% are done for 1 set of coefficients. For senstivity = 1, the local
% sensitivity is calculcated for a range of local perturbations in
% coefficients one at a time. For sensitivity = 2, the std dev calcs are
% done for user-defined changes in any of the coefficients (currently set
% to explore ranges of values of coefficients(1). 
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

clear all
sensitivity = 0; 
% if 0, just run the stddev calculations 
% if 1, run the sensitivity and the static stddev calculations for constant
% coefficients
% if 2, run through a range of coefficient values rather than true
% sensitivity

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
    % cases A & B: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics, pre-exposure
    % cases C & D: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics EXCEPT interactionRate, pre-exposure 
    % cases A & C: post-exposure estimate coef 4 and 5 distinctly with
    % paramcase 2
    % cases B & D: post-exposure estimate coef 5 distinctly and set 
    % coef 4 = to control/untreated with paramcase 3
if paramcase == 0
    % case 0; pre-exposure Cases A, B, C, & D 
    coefficients(1) = 0.09056; %COHORTenvironmentalStimuliWeight
elseif paramcase == 1 
    % case 1 is four coefficients for the pre exposure case
    coefficients(1) = 1; %COHORTenvironmentalStimuliWeight
    coefficients(2) = 1; %COHORTlambdaBrood 
    coefficients(3) = 0; %COHORTlambdaFullFood 
    coefficients(4) = 0; %COHORTlambdaEmptyFood
elseif paramcase == 2
    % case A & C: post-exposure estimate
    % coef 4 and 5 distinctly
    coefficients(1) = 0.05014; %COHORTenvironmentalStimuliWeight(4)
    coefficients(2) = 0.03489; %COHORTenvironmentalStimuliWeight(5)
elseif paramcase == 3
    % case B & D: post-exposure estimate
    % coef 5 distinctly and set coef 4 = to control/untreated
    coefficients(1) =   0.0519; %COHORTenvironmentalStimuliWeight(5) % 
end
 


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
totalTimePoints = 600; % 600 5min of simulation time % should take about 45 min if totalTimePoints = 7200
vis = 0;

nestSimulationData = zeros(totalTimePoints,numBees,5); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);


if strcmp(exposure_state,'pre')
    % from Plot Digitizer based on James' first blox plot emailed on 10/30
    % and saved as beeCItoDigitize.png
    CI95percent.activity = [0.74 0.83];
    CI95percent.distanceToNestmates = [0.081 0.085];
    CI95percent.distanceToBrood = [ NaN NaN]; %dummy values
    CI95percent.distanceToCenter = [0.055  0.062];
    CI95percent.porTimeOnNest = [0.56 0.68];
    CI95percent.interactionRate = [0.0050 0.0060];

    % from Plot Digitizer based on James' blox plot emailed on 11/30 for colony
    % 1 only and saved as beeCItoDigitizeColony1.png
    CI95percentColony1.activity = [0.61 0.84];
    CI95percentColony1.distanceToNestmates = [NaN NaN];
    CI95percentColony1.distanceToBrood = [ NaN NaN]; 
    CI95percentColony1.distanceToCenter = [0.053 0.067];
    CI95percentColony1.porTimeOnNest = [0.40 0.66];
    CI95percentColony1.interactionRate = [0.0051 0.0076];
elseif strcmp(exposure_state,'post')
    % from Plot Digitizer based on James' first blox plot emailed on 10/30
    % and saved as beeCItoDigitize.png
    CI95percent.activity = [0.45 0.57];
    CI95percent.distanceToNestmates = [0.085 0.090];
    CI95percent.distanceToBrood = [ NaN NaN]; %dummy values
    CI95percent.distanceToCenter = [0.063 0.070];
    CI95percent.porTimeOnNest = [0.36 0.50];
    CI95percent.interactionRate = [0.0043 0.0054];

    % from Plot Digitizer based on James' blox plot emailed on 11/30 for colony
    % 1 only and saved as beeCItoDigitizeColony1.png
    CI95percentColony1.activity = [0.36 0.59];
    CI95percentColony1.distanceToNestmates = [NaN NaN];
    CI95percentColony1.distanceToBrood = [ NaN NaN]; 
    CI95percentColony1.distanceToCenter = [0.056 0.073];
    CI95percentColony1.porTimeOnNest = [0.25 0.52];
    CI95percentColony1.interactionRate = [0.0037 0.0060];
end
CI95percent.mean = [mean(CI95percent.activity), mean(CI95percent.distanceToNestmates), mean(CI95percent.distanceToBrood), ...
    mean(CI95percent.distanceToCenter), mean(CI95percent.porTimeOnNest), mean( CI95percent.interactionRate)];
CI95percent.LB = [mean(CI95percent.activity)-CI95percent.activity(1), mean(CI95percent.distanceToNestmates)-CI95percent.distanceToNestmates(1), mean(CI95percent.distanceToBrood)-CI95percent.distanceToBrood(1), ...
    mean(CI95percent.distanceToCenter)-CI95percent.distanceToCenter(1), mean(CI95percent.porTimeOnNest)-CI95percent.porTimeOnNest(1), mean( CI95percent.interactionRate)-CI95percent.interactionRate(1)];
CI95percent.UB = [CI95percent.activity(2)-mean(CI95percent.activity), CI95percent.distanceToNestmates(2)-mean(CI95percent.distanceToNestmates), CI95percent.distanceToBrood(2)-mean(CI95percent.distanceToBrood), ...
    CI95percent.distanceToCenter(2)-mean(CI95percent.distanceToCenter), CI95percent.porTimeOnNest(2)-mean(CI95percent.porTimeOnNest), CI95percent.interactionRate(2)-mean( CI95percent.interactionRate)];

CI95percentColony1.mean = [mean(CI95percentColony1.activity), mean(CI95percentColony1.distanceToNestmates), mean(CI95percentColony1.distanceToBrood), ...
    mean(CI95percentColony1.distanceToCenter), mean(CI95percentColony1.porTimeOnNest), mean( CI95percentColony1.interactionRate)];
CI95percentColony1.LB = [mean(CI95percentColony1.activity)-CI95percentColony1.activity(1), mean(CI95percentColony1.distanceToNestmates)-CI95percentColony1.distanceToNestmates(1), mean(CI95percentColony1.distanceToBrood)-CI95percentColony1.distanceToBrood(1), ...
    mean(CI95percentColony1.distanceToCenter)-CI95percentColony1.distanceToCenter(1), mean(CI95percentColony1.porTimeOnNest)-CI95percentColony1.porTimeOnNest(1), mean( CI95percentColony1.interactionRate)-CI95percentColony1.interactionRate(1)];
CI95percentColony1.UB = [CI95percentColony1.activity(2)-mean(CI95percentColony1.activity), CI95percentColony1.distanceToNestmates(2)-mean(CI95percentColony1.distanceToNestmates), CI95percentColony1.distanceToBrood(2)-mean(CI95percentColony1.distanceToBrood), ...
    CI95percentColony1.distanceToCenter(2)-mean(CI95percentColony1.distanceToCenter), CI95percentColony1.porTimeOnNest(2)-mean(CI95percentColony1.porTimeOnNest), CI95percentColony1.interactionRate(2)-mean( CI95percentColony1.interactionRate)];
    
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
Ydata = [empiricalMeans.activity; empiricalMeans.distanceToNestmates; ...
    empiricalMeans.distanceToBrood; empiricalMeans.distanceToCenter; ...
    empiricalMeans.porTimeOnNest; empiricalMeans.interactionRate];
Ydatastddev = [nanstd(empiricalDistributions.activity); nanstd(empiricalDistributions.distanceToNestmates); ...
    nanstd(empiricalDistributions.distanceToBrood); nanstd(empiricalDistributions.distanceToCenter); ...
    nanstd(empiricalDistributions.porTimeOnNest); nanstd(empiricalDistributions.interactionRate)];
YdatastddevNormalized = [nanstd(empiricalDistributions.activity./empiricalMeans.activity);...
    nanstd(empiricalDistributions.distanceToNestmates./empiricalMeans.distanceToNestmates); ...
    nanstd(empiricalDistributions.distanceToBrood./empiricalMeans.distanceToBrood); ...
    nanstd(empiricalDistributions.distanceToCenter./empiricalMeans.distanceToCenter); ...
    nanstd(empiricalDistributions.porTimeOnNest./empiricalMeans.porTimeOnNest); ...
    nanstd(empiricalDistributions.interactionRate./empiricalMeans.interactionRate)];
sigma = Ydatastddev;
% sigma =CI95percent.LB';
% sigma(3) = Ydatastddev(3); % replace the value that is missing from the CI data with that from the experimental data

sigmahat = YdatastddevNormalized;
[A, Aindex] = sort(Ydatastddev);
[B, Bindex] = sort(CI95percent.LB);
[C, Cindex] = sort(CI95percentColony1.LB);
figure(1)
errorbar(Ydata,sigma)
hold on
errorbar(1:6,CI95percentColony1.mean,CI95percentColony1.LB,CI95percentColony1.UB,'k')
errorbar(1:6,CI95percent.mean,CI95percent.LB,CI95percent.UB,'g')
gca
set(gca,'XTickLabel',{' ','activity','distanceToNestmates','distanceToBrood','distanceToCenter','porTimeOnNest','interactionRate',' '})
axis([0 7 0 1.5])
hold off
legend('data stddev and mean','data 95% CI for Colony 1','data 95% CI for all 4 colonies')
% figure(2)
% errorbar(Ydata./Ydata,sigmahat)
% axis([0 7 0 1.5])
% gca
% set(gca,'XTickLabel',{' ','activity','distanceToNestmates','distanceToBrood','distanceToCenter','porTimeOnNest','interactionRate',' '})

%% define motion parameters
load(['motionsParams', exposure_state, '.mat'],'motionParamsCurrent');


 % Baseline stochastic variation
 numberTrials = 100; %number of time to repeat to assess stochastic variation 
 % in Ycalc values
 YcalcStochastic = zeros(6,numberTrials);
 for i = 1:numberTrials
     i
    nestSimulationData = simulationOutputSpatial(colony, estimatedData, exposure_state, totalTimePoints,motionParamsCurrent, vis,coefficients,paramcase);
    if paramcase == 3
        [means, ~]=calculateSummaryStatistics(nestSimulationData(:,tags==3,:),brood);
    else
        [means, ~]=calculateSummaryStatistics(nestSimulationData,brood);
    end
    activity = means.activity;
    distanceToNestmates = means.distanceToNestmates;
    distanceToBrood = means.distanceToBrood;
    distanceToCenter = means.distanceToCenter;
    porTimeOnNest = means.porTimeOnNest;
    interactionRate = means.interactionRate; 
    Ycalc = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest; interactionRate];
    YcalcStochastic(:,i) = Ycalc;
        sqdiffs = ((Ycalc-Ydata)).^2;
        sqpercentdiffs = ((Ycalc-Ydata)./Ydata).^2;
        sqpercentdiffs_over_sigmasq = ((Ycalc-Ydata)./Ydata./sigma).^2;
        sqdiffs_over_sigmasq = ((Ycalc-Ydata)./sigma).^2;
%         sqpercentdiffs_over_sigmahatsq = ((Ycalc-Ydata)./Ydata./sigmahat).^2
    outputStochastic(i) = sum(sqdiffs_over_sigmasq);
 end
 YcalcStochasticMean = mean(YcalcStochastic,2)
 YcalcStochasticStddev = std( YcalcStochastic')'
 outputStochasticMean = mean(outputStochastic)
 outputStochasticStddev = std( outputStochastic)
figure(1)
hold on
if numberTrials>1
    errorbar(YcalcStochasticMean,YcalcStochasticStddev,'ro')
else
    plot(YcalcStochasticMean,'o')
end
axis([0 7 0 1.5])
legend('data stddev and mean','data 95% CI for Colony 1','data 95% CI for all 4 colonies',['simulation results for ' num2str(numberTrials) ' runs'])
hold off

if sensitivity == 1
     vectorValues = 1:8;
 elseif sensitivity == 2
     % for range of coefficients
%      vectorValues = 0:0.1:0.9;
    vectorValues = 0.05:0.1:0.95;
 end
 if sensitivity >0
     % local sensitivity with increasing each coefficient one-at-a-time by
     % 5-fold decrease to 0.02
%      % 5-fold to 0.5 and then decreasing by 5-fold to 0.02
     for j = 1:length(vectorValues)
         if paramcase == 1
            coefficients = [0.1 0.1 0.1 0.1];
             if sensitivity == 1
             % for local sensitivity
                 if j <=4
                 coefficients(j) = 0.5;
                 elseif j >4
                     coefficients(j-4) = 0.02;
                 end
%                  coefficients(j) = 0.02;
             elseif sensitivity == 2
                % for range of coefficients
                coefficients(1) = vectorValues(j);
             end
         elseif paramcase == 2
             coefficients = [0.1 0.05];
             if sensitivity == 1
                 % for local sensitivity
                  sensitivitychanges = [.8 .9 1.1 1.2];
                  if j <= 4
                     coefficients(1) = coefficients(1)*sensitivitychanges(j);
                  else
                     coefficients(2) = coefficients(2)*sensitivitychanges(j-4);
                  end
             end             
         end
         
         for i = 1:numberTrials
            nestSimulationData = simulationOutputSpatial(colony, estimatedData, exposure_state, totalTimePoints,motionParamsCurrent, vis,coefficients,paramcase);
            [means distributions]=calculateSummaryStatistics(nestSimulationData,brood);
            activity = means.activity;
            distanceToNestmates = means.distanceToNestmates;
            distanceToBrood = means.distanceToBrood;
            distanceToCenter = means.distanceToCenter;
            porTimeOnNest = means.porTimeOnNest;
            interactionRate = means.interactionRate;
            Ycalc(:,i,j) = [activity; distanceToNestmates; distanceToBrood; distanceToCenter; porTimeOnNest; interactionRate];
            sqdiffs_over_sigmasq = ((Ycalc(:,i,j)-Ydata)./sigma).^2
            output(i,j) = sum(sqdiffs_over_sigmasq);
         end
        YcalcMean(:,j) = mean(Ycalc(:,:,j),2);
        for k = 1:5
            YcalcStddev(k,j) = std( Ycalc(k,:,j));
        end
        outputMean(j) = mean(output(:,j));
        outputStddev(j) = std( output(:,j));
     end
     YcalcMean
     YcalcStddev
     outputMean
     outputStddev
     if sensitivity ==2
        plot(vectorValues,outputMean)
        errorbar(vectorValues,outputMean,outputStddev)
     end
 end
% save('pre-exposure_spatialOptimizer.mat')
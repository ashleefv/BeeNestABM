%% Run simulation R times in parallel for run-to-run variation statistics
R = 8; %ran with 1000, took about 2 hours
%% Step 1: Defining the initial state & parameters
tic
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['data\allDataCol' num2str(colonyNumber) '.mat']);
brood = relabelBroodObject(colony.broodPre);
tags = colony.orTagTreat;
% number of time points
totalTimePoints = 7200;
vis = 0;
numBees = size(colony.preNest,2);

nestSimulationData = zeros(totalTimePoints,numBees,5); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);

exposure_state = 'pre'; 

activity = zeros(1,R);
distanceToNestmates = zeros(1,R);
distanceToBrood = zeros(1,R);

parpool(4)
parfor i = 1:R
     nestSimulationData = simulationOutput(colony, estimatedData, exposure_state, totalTimePoints, vis);
     [means distributions]=calculateSummaryStatistics(nestSimulationData,brood,tags);
     activity(i) = means.activity;
     distanceToNestmates(i) = means.distanceToNestmates;
     distanceToBrood(i) = means.distanceToBrood;
end
Mean.activity = mean(activity);
Std.activity = std(activity);
Mean.distanceToNestmates = mean(distanceToNestmates);
Std.distanceToNestmates = std(distanceToNestmates);
Mean.distanceToBrood = mean(distanceToBrood);
Std.distanceToBrood = std(distanceToBrood);
Mean
Std
delete(gcp)
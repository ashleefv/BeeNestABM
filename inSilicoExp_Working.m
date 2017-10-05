%Load initial data
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
estimatedData = load(['data/Essential_Info_Col_' num2str(colonyNumber) '.mat']);
brood = colony.broodPre;
treatmentList = colony.orTagTreat;
%%
treatments = {'full', 'sp', 'soc', 'control', 'fullNoSp', 'controlNoSp'}; %Simulation types -
%full is fully parameterized model (including treatment-specific effects
%
%  "sp" is sppontaneous motion parameters driven  differently between control
%and treatment
%
% "soc" - inverse of "sp", only socially-modulated switching parameters are
% affected
%
% "control" - all "treatment" group values are set equal to control values
%
% "fullNoSp" - including all treatment effects, but average parameters
% between on and off nest (wiping out spatial effects)
%
%  "controlNoSp" -

%% Set simulation parameters
ntrt = length(treatments); %How many treatments?
nrep = 3; %How many simulated replicates to run for each condition??
groups = 1:3;
samplesPerRun = 1000;
burnIn = 100; %How many timesteps to ignore at the beginning of the simulation?

totSamps = ntrt*nrep*numel(groups); %Total samples


%% Initialize output data structure
outData = array2table(nan(totSamps, 7));
outData.Properties.VariableNames = {'activity', 'distanceToNestmates', 'distanceToBrood', 'distanceToCenter', 'porTimeOnNest', 'simType', 'treatmentGroup'};

%% define motion parameters
%Add documentation about where these values come from


optimization = 1; %are we running optimization trials, or experimental trials? Optimization trials will use
%empirically derived parameters combined for all groups pre-exposure,
%otherwise group-specific post-exposure parameters are implemented

if optimization  == 1
    
    
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
    
elseif optimization == 0
    
    
    %On nest parameters
    AIU1 = .015; %Active to inactive unbumped, group 1
    AIU2 = .015; % group 2...
    AIU3 = .015; %Etc
    AIU4 = .027;
    
    IAU1 = .294; %Inactive to active, unbumped, group 1...
    IAU2 = .294;
    IAU3 = .294;
    IAU4 = .294;
    
    AIB1 = .01; %Active to inactive, bumped, group 1...
    AIB2 = .01;
    AIB3 = .01;
    AIB4 = .01;
    
    IAB1 = .32;
    IAB2 = .32;
    IAB3 = .32;
    IAB4 = .32;
    
    motionParamsOnNest = [AIU1 AIU2 AIU3 AIU4; IAU1 IAU2 IAU3 IAU4; AIB1 AIB2 AIB3 AIB4; IAB1 IAB2 IAB3 IAB4];
    
    %Off nest parameters
    AIU1 = .0305; %Active to inactive unbumped, group 1
    AIU2 = .0305; % group 2...
    AIU3 = .0305; %Etc
    AIU4 = .087;
    
    IAU1 = .2137; %Inactive to active, unbumped, group 1...
    IAU2 = .2137;
    IAU3 = .2137;
    IAU4 = .1047;
    
    AIB1 = .0206; %Active to inactive, bumped, group 1...
    AIB2 = .0206;
    AIB3 = .0206;
    AIB4 = .079;
    
    IAB1 = .25;
    IAB2 = .25;
    IAB3 = .25;
    IAB4 = .134;
    
    %For "full" model
    motionParamsOffNest = [AIU1 AIU2 AIU3 AIU4; IAU1 IAU2 IAU3 IAU4; AIB1 AIB2 AIB3 AIB4; IAB1 IAB2 IAB3 IAB4];
    
    %Combine on and off-nest parameters sets
    motionParamsFull = cat(3,motionParamsOnNest, motionParamsOffNest);
    
end
%%
%Change to add on means for each separate treatment group within colonies
zz = 1; %Indexing variable
tic
h = waitbar(0, 'Nest Simulation Progress');
for i = 1:ntrt
    %%
    %treatments = {'full', 'sp', 'soc', 'control', 'fullNoSp', 'controlNoSp'}; %Simulation types -
    treatment = treatments{i};
    
    motionParamsCurrent = motionParamsFull; %Initialize treatment-specific parameters
    
    if strcmp(treatment, 'full')
        
        disp('full model - no parameters adjustment');
        
    elseif strcmp(treatment, 'sp')
        
        disp('spontaneous motion model - setting treatments social effects equal to controls');
        motionParamsCurrent(3:4,4,:) = motionParamsCurrent(3:4,2,:) %Set bump parameters for treatment group equal to controls
        
    elseif strcmp(treatment, 'soc')
        
        disp('social switching model - setting treatments spontaneous effects equal to controls');
        motionParamsCurrent(1:2,4,:) = motionParamsCurrent(1:2,2,:) %Set spontaneouns parameters for treatment group equal to controls
        
    elseif strcmp(treatment, 'control')
        
        disp('control model - setting treatment parameters equal to controls for social and spontaneous effecs');
        motionParamsCurrent(1:4,4,:) = motionParamsCurrent(1:4,2,:) %Set both parameter sets for treatment group equal to controls
        
    elseif strcmp(treatment, 'fullNoSp')
        
        disp('no spatial effects full model - setting on and off nest parameter sets equal');
        motionParamsCurrent = mean(motionParamsCurrent,3); %average on and off nest parameters
        motionParamsCurrent = cat(3,motionParamsCurrent, motionParamsCurrent)
        
    elseif strcmp(treatment, 'controlNoSp')
        
        disp('no spatial effects control model - setting on and off nest parameter sets equal, then setting treatment effects equal to controls');
        motionParamsCurrent = mean(motionParamsCurrent,3); %average on and off nest parameters
        motionParamsCurrent = cat(3,motionParamsCurrent, motionParamsCurrent);
        motionParamsCurrent(1:4,4,:) = motionParamsCurrent(1:4,2,:)
    end
    
    %% Run simulations
    
    for j = 1:nrep
        %%
        if optimization == 1
            %If we're running parameter optimization, use brood data and initial
            %position from pre-exposure trials
            simDat = simulationOutputSpatial(colony, estimatedData, 'pre', samplesPerRun, motionParamsCurrent,1);
        elseif optimization == 0
            %If we're running an experiment, use post-exposure trial data
            simDat = simulationOutputSpatial(colony, estimatedData, 'post', samplesPerRun, motionParamsCurrent,1);
        end
        
        %% extract data separately for each group and append
        for k = 1:numel(groups)
            %%
            [means distributions] = calculateSummaryStatistics(simDat(burnIn:end,treatmentList == groups(k),:), brood);
            outData.activity(zz) = means.activity;
            outData.distanceToNestmates(zz) = means.distanceToNestmates;
            outData.distanceToBrood(zz) = means.distanceToBrood;
            outData.distanceToCenter(zz) = means.distanceToCenter;
            outData.porTimeOnNest(zz) = means.porTimeOnNest;
            outData.simType(zz) = i;
            outData.treatmentGroup(zz) = k;
            waitbar(zz/totSamps, h)
            zz = zz +1;
        end
    end
    clear motionParamsCurrent
    
end
toc
close(h)
%
writetable(outData, 'inSilicoExpOutputSep12.csv');
%
treatmentTable = cell2table(treatments');
writetable(treatmentTable, 'inSilicoTreatmentTypes.csv');
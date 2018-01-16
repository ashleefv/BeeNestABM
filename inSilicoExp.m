%Load initial data
%base = '/n/home04/jcrall/pesticideModel';
base = '/Users/james/Documents/pesticideProject';

runLabel = 'Jan9_master';
addpath(base);
%base = '/Users/james/Documents/pesticideProject';
cd(base);
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
estimatedData = load(['data/Essential_Info_Col_' num2str(colonyNumber) '.mat']);
brood = colony.broodPre;
treatmentList = colony.orTagTreat;
taglist = colony.tags; %Read out colony's list of marked bees

optimization = 0; %Optimiization(1) or experiment(0)?

%%
%treatments = {'full', 'sp', 'soc', 'control', 'fullNoSp', 'controlNoSp'}; %Simulation types -
%full is fully parameterized model (including treatment-specific effects
%
%  "sp" is sppontaneous motion parameters driven  differently between control
%and treatment
%
% "soc" - inverse of "sp", only socially-modulated switching parameters are
% affected
%
% 'act' - full activity switching model, but no group-specific
% environmental attraction
%
% "control" - all "treatment" group values are set equal to control values
%
% "fullNoSp" - including all treatment effects, but average parameters
% between on and off nest (wiping out spatial effects)
%
% 'ap' - include group-specific environmental attraction parameters, but
% equivalent mobility-switching parameters for all
%
%  "controlNoSp" -

treatments = {'control', 'sp', 'soc', 'act', 'ap', 'full'}; %Simulation types  for current run

%% Set simulation parameters
ntrt = length(treatments); %How many treatments?
nrep = 1000; %How many simulated replicates to run for each condition??
groups = 1:3;
ngroups = numel(groups);
samplesPerRun = 7200;
burnIn = 200; %How many timesteps to ignore at the beginning of the simulation?

totSamps = ntrt*nrep*numel(groups); %Total samples


%% define motion parameters
params = readtable('parameterEstimates.csv');
%rng(1); %Set random number generation seed for replicability
if optimization  == 1 %For setting motion parameters, are these being drawn from the empirical values (specific to
    %individual and colony), or by reshuffling the distribution of values by cohort? (non-empirical)?
    empirical = 1;
    pre = 1;
    motionParamsCurrent = generateMotionParamsObject(params, colonyNumber, taglist, pre, empirical);
    
elseif optimization == 0
    disp('running a simulation - will resample transition matrixes for each timestep');
    empirical = 0;
    pre = 0;
    motionParamsCurrent = generateMotionParamsObject(params, colonyNumber, taglist, pre, empirical); %Generate randomly sampled, treatment-specific parameters
    motionParamsControl = generateMotionParamsObject(params, colonyNumber, taglist, pre, 2); %Generate randomly sampled control data
    
end
% motionParamsCurrent
% subplot(2,1,1);
% boxplot(motionParamsCurrent.AtoI_Unbumped_OffNest, treatmentList == 3)
% subplot(2,1,2);
% boxplot(motionParamsControl.AtoI_Unbumped_OffNest, treatmentList == 3)

%%
parpool('local', 20);

%%


for i = 1:ntrt
    treatment = treatments{i};
    
    % Initialize separate output structures for each treatment groups
    g1Dat = zeros(nrep, 8);
    g2Dat = zeros(nrep, 8);
    g3Dat = zeros(nrep, 8);
    
    %% Run simulations
    
    
    parfor j = 1:nrep
        %%
        
        if optimization == 1
            %If we're running parameter optimization, use brood data and initial
            %position from pre-exposure trials
            %simDat = simulationOutputSpatial(colony, estimatedData, 'pre', samplesPerRun, motionParamsCurrent,0);
            
        elseif optimization == 0
            %If we're running an experiment, use post-exposure trial data, and
            %resample motion parameters for each run
            disp('running a simulation - will resample transition matrixes for each timestep');
            empirical = 0;
            pre = 0;
            motionParamsCurrent = generateMotionParamsObject(params, colonyNumber, taglist, pre, empirical); %Generate randomly sampled, treatment-specific parameters
            motionParamsControl = generateMotionParamsObject(params, colonyNumber, taglist, pre, 2); %Generate randomly sampled control data
            %%
            %treatments = {'full', 'sp', 'soc', 'control', 'fullNoSp', 'controlNoSp'}; %Simulation types -
            
            
            
            socVars = {'AtoI_Bumped_OnNest', 'ItoA_Bumped_OnNest', 'AtoI_Bumped_OffNest', 'ItoA_Bumped_OffNest'};
            sponVars = {'AtoI_Unbumped_OnNest', 'ItoA_Unbumped_OnNest', 'AtoI_Unbumped_OffNest', 'ItoA_Unbumped_OffNest'};
            
            onNestVars = {'AtoI_Bumped_OnNest', 'ItoA_Bumped_OnNest','AtoI_Unbumped_OnNest', 'ItoA_Unbumped_OnNest'};
            offNestVars = {'AtoI_Bumped_OffNest', 'ItoA_Bumped_OffNest', 'AtoI_Unbumped_OffNest', 'ItoA_Unbumped_OffNest'};
            
            if strcmp(treatment, 'full')
                
                disp('full model - no parameters adjustment');
                exposure_state = 'post'; %Use post-exposure attraction parameters
                
            elseif strcmp(treatment, 'ap')
                
                disp('attraction parameter model - no parameters adjustment');
                exposure_state = 'post'; %Use post-exposure attraction parameters
                motionParamsCurrent = motionParamsControl;
                
            elseif strcmp(treatment, 'sp')
                
                disp('spontaneous motion model - setting treatments social effects equal to controls');
                motionParamsCurrent(:,socVars) = motionParamsControl(:,socVars);
                exposure_state = 'pre'; %Use equal attraction for all
                
            elseif strcmp(treatment, 'soc')
                
                disp('social switching model - setting treatments spontaneous effects equal to controls');
                motionParamsCurrent(:,sponVars) = motionParamsControl(:,sponVars);
                exposure_state = 'pre';
                
            elseif strcmp(treatment, 'act')
                disp('activity switching model - setting treatments spontaneous effects equal to controls, equivalent environmental attraction across groups');
                exposure_state = 'pre';
                
                
            elseif strcmp(treatment, 'control')
                
                disp('control model - setting treatment parameters equal to controls for social and spontaneous effecs');
                motionParamsCurrent = motionParamsControl;
                exposure_state = 'pre';
                
            elseif strcmp(treatment, 'fullNoSp')
                
                disp('no spatial effects full model - setting on and off nest parameter sets equal');
                motionParamsMean = cat(3, table2array(motionParamsCurrent(:,onNestVars)), table2array(motionParamsCurrent(:,offNestVars)));
                motionParamsMean = mean(motionParamsMean, 3);
                motionParamsCurrent(:,onNestVars) = array2table(motionParamsMean);
                motionParamsCurrent(:,offNestVars) = array2table(motionParamsMean);
                exposure_state = 'post';
                
            elseif strcmp(treatment, 'controlNoSp')
                
                disp('no spatial effects control model - setting on and off nest parameter sets equal, then setting treatment effects equal to controls');
                motionParamsMean = cat(3, table2array(motionParamsControl(:,onNestVars)), table2array(motionParamsControl(:,offNestVars)));
                motionParamsMean = mean(motionParamsMean, 3);
                motionParamsCurrent(:,onNestVars) = array2table(motionParamsMean);
                motionParamsCurrent(:,offNestVars) = array2table(motionParamsMean);
                exposure_state = 'post';
            end
            
            simDat = simulationOutputSpatial(colony, estimatedData, exposure_state, samplesPerRun, motionParamsCurrent,0);
            
        end
        
        %% Initialize output data structure
        [means, distributions] = calculateSummaryStatistics(simDat(burnIn:end,:,:), brood);
        
        %outTemp = zeros(numel(groups), 8);
        
        %%
        k = 1; %Group 1
        dat1 = [ nanmean(distributions.activity(treatmentList == groups(k)))  nanmean(distributions.distanceToNestmates(treatmentList == groups(k)))...
            nanmean(distributions.distanceToBrood(treatmentList == groups(k)))  nanmean(distributions.distanceToCenter(treatmentList == groups(k)))...
            nanmean(distributions.porTimeOnNest(treatmentList == groups(k))) nanmean(distributions.interactionRate(treatmentList == groups(k))) i k];
        
        g1Dat(j,:) = dat1;
        
        k = 2; % Group 2
        dat2 = [ nanmean(distributions.activity(treatmentList == groups(k)))  nanmean(distributions.distanceToNestmates(treatmentList == groups(k)))...
            nanmean(distributions.distanceToBrood(treatmentList == groups(k)))  nanmean(distributions.distanceToCenter(treatmentList == groups(k)))...
            nanmean(distributions.porTimeOnNest(treatmentList == groups(k))) nanmean(distributions.interactionRate(treatmentList == groups(k))) i k];
        
        g2Dat(j,:) = dat2;
        
        k = 3; %Group 3
        dat3 = [ nanmean(distributions.activity(treatmentList == groups(k)))  nanmean(distributions.distanceToNestmates(treatmentList == groups(k)))...
            nanmean(distributions.distanceToBrood(treatmentList == groups(k)))  nanmean(distributions.distanceToCenter(treatmentList == groups(k)))...
            nanmean(distributions.porTimeOnNest(treatmentList == groups(k))) nanmean(distributions.interactionRate(treatmentList == groups(k))) i k];
        
        g3Dat(j,:) = dat3;
        
    end
    
    
    %pestParSave(strcat(pwd,'/inSilicoExpOutput/',  treatment,num2str(j),'.mat'), 'outTemp', 'variables');
    variables = {'activity', 'distanceToNestmates', 'distanceToBrood', 'distanceToCenter', 'porTimeOnNest', 'interactionRate','simType', 'treatmentGroup'};
    
    outData = array2table([g1Dat; g2Dat; g3Dat]);
    outData.Properties.VariableNames = variables;
    writetable(outData, strcat(pwd,'/inSilicoExpOutput/', runLabel, '_', treatment, '.csv'));
    
end
%% save list of treatment used

treatmentTable = cell2table(treatments');
writetable(treatmentTable, strcat(pwd,'/inSilicoExpOutput/', runLabel, '_inSilicoTreatmentTypes.csv'));


%% Load and save empircal data for comparison

for i = 1:4
    %%
    cd(base);
    colonyNumber = i; % 4 colonies numbered 1 - 4
    colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
    estimatedData = load(['data/Essential_Info_Col_' num2str(colonyNumber) '.mat']);
    brood = colony.broodPost;
    treatmentList = colony.orTagTreat;
    taglist = colony.tags; %Read out colony's list of marked bees
    nbees = numel(taglist);
    
    [means, distributions] = calculateSummaryStatistics(colony.postNest, brood);
    
    outData = zeros(nbees, 9);
    %%
    outData(:,1) = distributions.activity;
    outData(:,2) = distributions.distanceToNestmates;
    outData(:,3) = distributions.distanceToBrood;
    outData(:,4) = distributions.distanceToCenter;
    outData(:,5) = distributions.porTimeOnNest;
    outData(:,6) = distributions.interactionRate;
    outData(:,7) = repelem(0,nbees);
    outData(:,8) = treatmentList;
    outData(:,9) = repelem(i,nbees);
    outData = array2table(outData);
    outData.Properties.VariableNames = {'activity', 'distanceToNestmates', 'distanceToBrood', 'distanceToCenter', 'porTimeOnNest', 'interactionRate','simType', 'treatmentGroup', 'colony'};
    
    if i == 1
        outDataM = outData
    else
        outDataM = [outDataM; outData];
    end
end

writetable(outDataM, strcat(pwd,'/inSilicoExpOutput/', runLabel, '_empirical','.csv'));


function motionParams = generateMotionParamsObject(params, colonyNum, taglist, pre, empirical)
    
    %%
    %
    %  Inputs:
    %       params - empirical estimations (drawn from
    %         "parameterEstimates.csv", output by
    %         "alternateParameterEstimation.m"
    %
    %       colonyNum - which of 4 colonies to draw data from?
    %
    %       taglist - taglist supplied from "colony" object from empirical data - just for checking that ordering is correct
    %
    %       pre - draw on pre or post exposure parameter estimates?
    %
    %       empirical - set all values possible equal to empirical estimates (1),
    %           or randomly sampling with treatment category (0)?
    %
    %
    %  Outputs
    %
    %       motionParams - n x 8 table, where n = number of bees in this
    %       simulation, with columns for each of 8 mobility-switching
    %       parameters
    
    
    %%
    %
    varNames = {'AtoI_Unbumped_OnNest', 'ItoA_Unbumped_OnNest', 'AtoI_Bumped_OnNest', 'ItoA_Bumped_OnNest',...
        'AtoI_Unbumped_OffNest', 'ItoA_Unbumped_OffNest', 'AtoI_Bumped_OffNest', 'ItoA_Bumped_OffNest'};
    
    %Subset to correct
    colParams = params(params.colony == colonyNum & params.PreExposure == pre,:); %Subset to correct colony and time period (pre vs. post)
    
    %Check for correct ordering with respect to taglist
    check = sum(taglist ~= colParams.bee);
    if check > 0
        disp('ordering of taglist and parameter data dont agree - check input data!');
        return
    else
        disp('ordering of taglist and parameter data in agreement!');
    end
    
    %Subset data into queen and worker sets
    queenParams = colParams(colParams.queen == 1,:);
    workerParams = colParams(colParams.queen ~= 1,:);
    
    %Remove any rows from sampled data with missing values
    rmInd = sum(isnan(table2array(workerParams(:,varNames))),2) > 0;
    workerParams = workerParams(~rmInd,:);
    
    treatedInd = workerParams.group == 3;
    
    treated = workerParams(treatedInd,varNames);
    control = workerParams(~treatedInd,varNames);
    
    
    %     treated.AtoI_Unbumped_OnNest = workerParams.AtoI_Unbumped_OnNest(treatedInd); %This parameter gets treatment-specific values
    %     treated.ItoA_Unbumped_OnNest = workerParams.ItoA_Unbumped_OnNest; %This parameter isn't significantly affected by treatment (see 'plotParameterEffectsV2.R'), so gets the same pool of values as controls
    %     treated.AtoI_Bumped_OnNest = workerParams.AtoI_Bumped_OnNest;
    %     treated.ItoA_Bumped_OnNest = workerParams.ItoA_Bumped_OnNest;
    %     treated.AtoI_Unbumped_OffNest = workerParams.AtoI_Unbumped_OffNest(treatedInd);
    %     treated.ItoA_Unbumped_OffNest = workerParams.ItoA_Unbumped_OffNest(treatedInd);
    %     treated.AtoI_Bumped_OffNest = workerParams.AtoI_Bumped_OffNest(treatedInd);
    %     treated.ItoA_Bumped_OffNest = workerParams.ItoA_Bumped_OffNest(treatedInd);
    %
    %
    %     control.AtoI_Unbumped_OnNest = workerParams.AtoI_Unbumped_OnNest(~treatedInd); %This parameter gets treatment-specific values
    %     control.ItoA_Unbumped_OnNest = workerParams.ItoA_Unbumped_OnNest; %This parameter isn't significantly affected by treatment (see 'plotParameterEffectsV2.R'), so gets the same pool of values as controls
    %     control.AtoI_Bumped_OnNest = workerParams.AtoI_Bumped_OnNest;
    %     control.ItoA_Bumped_OnNest = workerParams.ItoA_Bumped_OnNest;
    %     control.AtoI_Unbumped_OffNest = workerParams.AtoI_Unbumped_OffNest(~treatedInd);
    %     control.ItoA_Unbumped_OffNest = workerParams.ItoA_Unbumped_OffNest(~treatedInd);
    %     control.AtoI_Bumped_OffNest = workerParams.AtoI_Bumped_OffNest(~treatedInd);
    %     control.ItoA_Bumped_OffNest = workerParams.ItoA_Bumped_OffNest(~treatedInd);
    
    motionParams = colParams(:,varNames);
    
    if ~empirical %if we're randomly assigning values, wipe out all empirical values first
        motionParams(:,:) = array2table(nan(size(motionParams)));
    end
    
    
    %Fill in missing values with randomly assigned, treatment-specific values
    
    varNames = motionParams.Properties.VariableNames
    %Double check that variables are in the same order for input and output
    %matrices
    check = sum(strcmp(control.Properties.VariableNames, varNames) == 0);
    if check > 0
        disp('mismatch in order of variable names for aggregated individual parameters and output parameter matrix - check data');
        return
    else
        disp('order of variables all good!');
    end
    
    
    %%
    
    if empirical ==1 %Option 1 - use empirical values for each specific bee, sample only missing values
        for j = 1:numel(motionParams(:,1))
            for i = 1:numel(varNames)
                %%
                if isnan(table2array(motionParams(j,i)))
                    %%
                    if colParams.group(j) == 3
                        sampInd = randsample(1:numel(treated(:,1)), 1);
                        motionParams(j,i) = treated(sampInd,i);
                    else
                        sampInd = randsample(1:numel(control(:,1)), 1);
                        motionParams(j,i) = control(sampInd,i);
                    end
                end
            end
        end
        
        
    elseif empirical == 0 %Option 2 - sample all values but from treatment-specific groups
        trtInd = colParams.group == 3;
        numTrt = sum(trtInd);
        numCtrl = sum(~trtInd);
        trtSample = randsample(1:numel(treated(:,1)), numTrt, 1);
        ctrlSample = randsample(1:numel(control(:,1)), numCtrl, 1);
        motionParams(trtInd,:) = treated(trtSample,:);
        motionParams(~trtInd,:) = control(ctrlSample,:);
        
    elseif empirical == 2 %Randomly sample everyone from controls
        nbees = numel(colParams.group);
        samp = randsample(1:numel(control(:,1)), nbees, 1);
        motionParams(:,:) = control(samp,:);
        
    end
    
    motionParams(1,:) = colParams(1,varNames); %In either scenario, set queen switching parameters equal to empirical
    

function motionParams = generateMotionParamsObject_general(params, orTagTreat)
    
    %%
    %
    %  Inputs:
    %       params - empirical estimations (drawn from
    %         "parameterEstimates.csv", output by
    %         "alternateParameterEstimation.m"
    %
    %       orTagTreat - vector identifying which bees are treated (0s for
    %       untreated, 3s for treated)
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
    

    
    
    workerParams = params(params.queen ~= 1 & params.group > 0,:);
    
    %Remove any rows from sampled data with missing values
    rmInd = sum(isnan(table2array(workerParams(:,varNames))),2) > 0;
    workerParams = workerParams(~rmInd,:);
    
    treatedInd = workerParams.group == 3 & workerParams.PreExposure == 0;
    
    treated = workerParams(treatedInd,varNames);
    control = workerParams(~treatedInd,varNames);
    
    %Calculate group-mean values from empirical data
    treatedMeans = nanmean(table2array(treated));
    controlMeans = nanmean(table2array(control));

    %Initialize empty motionParams object
    motionParams = nan(numel(orTagTreat), 8);
    
    
    %Allocate switching parameters
   for i = 1:numel(orTagTreat)
       
      if orTagTreat(i) == 3
          motionParams(i,:) = treatedMeans;
      else
          motionParams(i,:) = controlMeans;
      end
      
   end
    motionParams = array2table(motionParams);
    motionParams.Properties.VariableNames = varNames;

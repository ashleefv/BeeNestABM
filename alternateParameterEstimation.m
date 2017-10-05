cd /Users/james/Documents/pesticideProject


for pre = [0 1]
    for colonyNumber = 1:4
        %% Load data
        %colonyNumber = 1; % 4 colonies numbered 1 - 4
        colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
        %estimatedData = load(['data/Essential_Info_Col_' num2str(colonyNumber) '.mat']);
        brood = colony.broodPre;
        treatmentList = colony.orTagTreat;
        tags = colony.tags;
        %% Extract nest data
        %
        %Yields M x N x 7 nestData matrix, where sheets are x (1), y (2),
        %orientation(3), speed(4), activity(5), on nest (6), and bumped (7) states
        %for each frame
        if pre == 1
            nestData = colony.preNest;
        else
            
            nestData = colony.postNest;
        end
        
        nestData = calculateActivityMatrix(nestData);
        [porTimeOnNest onNest] = calculatePortionOfTimeOnNest(nestData, brood, .01);
        nestData(:,:,6) = onNest;
        distMat = calculatePairwiseDistanceMatrix(nestData);
        distThresh = 0.01;
        intMat = distMat < distThresh;
        bumped = nan(size(nestData(:,:,5)));
        idInd = logical(eye(size(intMat(:,:,1),1)));
        
        for i = 1:size(nestData,1)
            intMatC = intMat(:,:,i);
            intMatC(idInd) = 0;
            
            bumped(i,:) =  max(intMatC);
        end
        nestData(:,:,7) = bumped;
        %%
        subplot(2,1,1);
        imagesc(onNest);
        subplot(2,1,2);
        imagesc(~isnan(onNest));
        
        %% Estimate transition probabilities
        nbees = size(nestData,2)
        nframes = size(nestData,1);
        
        parameterData = array2table(nan(nbees, 15));
        parameterData.Properties.VariableNames = {'bee', 'colony', 'queen','group','PreExposure','framesOnNest', 'framesOffNest','AtoI_UnbumpedOnNest', 'ItoA_UnbumpedOnNest', 'AtoI_BumpedOnNest', 'ItoA_BumpedOnNest',...
            'AtoI_UnbumpedOffNest', 'ItoA_UnbumpedOffNest', 'AtoI_BumpedOffNest', 'ItoA_BumpedOffNest'};
        
        
        %% loop across tags to estimate motion parameters
        
        for i = 1:numel(tags)
            %% Specify metadata parameters
            activity = nestData(:,i,5);
            onNest = nestData(:,i,6);
            bumped = nestData(:,i,7);
            nanind = isnan(nestData(:,i,1));
            
            bee = tags(i);
            colony = colonyNumber;
            if i == 1
                queen = 1;
            else
                queen = 0;
            end
            group = treatmentList(i);
            
            %%
            %Calculate parameters for on nest
            nframes = size(activity,1);
            tmpInd = find(activity == 1 & onNest == 1 & bumped  == 0);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            AtoI_UnbumpedOnNest = 1- nanmean(activity(ind2));
            
            tmpInd = find(activity == 0 & onNest == 1 & bumped  == 0);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            ItoA_UnbumpedOnNest = nanmean(activity(ind2));
            
            tmpInd = find(activity == 1 & onNest == 1 & bumped  == 1);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            AtoI_BumpedOnNest = 1 - nanmean(activity(ind2));
            
            tmpInd = find(activity == 0 & onNest == 1 & bumped  == 1);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            ItoA_BumpedOnNest = nanmean(activity(ind2));
            %%
            %Calculate parameters for on nest
            nframes = size(activity,1);
            tmpInd = find(activity == 1 & onNest == 0 & bumped  == 0);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            AtoI_UnbumpedOffNest = 1- nanmean(activity(ind2));
            
            tmpInd = find(activity == 0 & onNest == 0 & bumped  == 0);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            ItoA_UnbumpedOffNest = nanmean(activity(ind2));
            
            tmpInd = find(activity == 1 & onNest == 0 & bumped  == 1);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            AtoI_BumpedOffNest = 1 - nanmean(activity(ind2));
            
            tmpInd = find(activity == 0 & onNest == 0 & bumped  == 1);
            ind2 = tmpInd+1;
            ind2 = ind2(ind2 <= nframes);
            ItoA_BumpedOffNest = nanmean(activity(ind2));
            %% write data
            %Metadata
            parameterData.bee(i) = bee;
            parameterData.queen(i) = queen;
            parameterData.group(i) = group;
            parameterData.colony(i) = colony;
            parameterData.framesOnNest(i) = sum(onNest == 1);
            parameterData.framesOffNest(i) = sum(onNest == 0);
            parameterData.PreExposure(i) = pre;
            %On nest parameters
            parameterData.AtoI_UnbumpedOnNest(i) = AtoI_UnbumpedOnNest;
            parameterData.ItoA_UnbumpedOnNest(i) = ItoA_UnbumpedOnNest;
            parameterData.AtoI_BumpedOnNest(i) = AtoI_BumpedOnNest;
            parameterData.ItoA_BumpedOnNest(i) = ItoA_BumpedOnNest;
            
            %Off nest parameters
            parameterData.AtoI_UnbumpedOffNest(i) = AtoI_UnbumpedOffNest;
            parameterData.ItoA_UnbumpedOffNest(i) = ItoA_UnbumpedOffNest;
            parameterData.AtoI_BumpedOffNest(i) = AtoI_BumpedOffNest;
            parameterData.ItoA_BumpedOffNest(i) = ItoA_BumpedOffNest;
            
        end
        
        %% Append to master data
        if colonyNumber == 1 & pre == 0
            parameterDataMaster = parameterData;
        else
            parameterDataMaster = [parameterDataMaster ; parameterData]
        end
    end
end

writetable(parameterDataMaster,'parameterEstimates.csv')
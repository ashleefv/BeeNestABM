function [means distributions] = calculateSummaryStatistics(data, brood)
    
    
    %%
    %Input
    %data - nestData-shaped object, in MxNx5 format, after passing through
    %"calculateActivityMatrix" OR if there are less then 5 "sheets", we run
    %
    %brood - Nx3 matrix with coordinates data on brood and nest structure map
    %
    %
    %Output:
    %means - struct with 3 values, one each for the mean activity, mean distance to
    %other bees, and mean distance to brood
    %
    %distributions - struct with 3 vectors, same values but with separate
    %estimates for each bee
    
    if size(data, 3) < 5
        data = calculateActivityMatrix(data);
    end
    %% Clean out zeros if there are any
    %Why are there zeros here?
    %Get index of zeros
    zeroInd = data(:,:,1) == 0;
    
    %Extract and replace x data
    x = data(:,:,1);
    x(zeroInd) = NaN;
    
    %Extract and replace ydata
    y = data(:,:,2);
    y(zeroInd) = NaN;
    
    %Write back into "data"
    data(:,:,1) = x;
    data(:,:,2) = y;
    
    %% Crate empty structure arrays
    means = struct();
    distributions = struct();
    
    %% Calculate Activity level
    act = data(:,:,5);
    act(zeroInd) = NaN;
    activity = nanmean(act);
    distributions.activity = activity;
    means.activity = nanmean(activity);
    
    %% Distance to nestmates
    distMat = calculatePairwiseDistanceMatrix(data(:,:,1:2));
    meanDistanceToOtherBees = nanmean(nanmean(distMat,3));
    means.distanceToNestmates = nanmean(meanDistanceToOtherBees);
    distributions.distanceToNestmates = meanDistanceToOtherBees;
    
    %% distance To nest structures
    
    broodDist = nan(size(data(:,:,1)));
    minBroodDist = broodDist;
    for i = 1:size(data,1)
        %%
        dists = pdist2(permute(data(i,:,1:2), [2 3 1]), brood(:,1:2));
        broodDist(i,:) = nanmean(dists,2);
        minBroodDist(i,:) = nanmin(dists');
    end
    
    distanceToBrood = nanmean(broodDist);
    means.distanceToBrood = nanmean(distanceToBrood);
    distributions.distanceToBrood = distanceToBrood;
    
    %% Distance to Nest Center
    xc = nanmean(nanmean(data(:,:,1))); %Calculate x center
    yc = nanmean(nanmean(data(:,:,2))); %Calculate y center
    
    totDist = sqrt((data(:,:,1) - xc).^2 + (data(:,:,2) - yc).^2);
    distanceToCenter = nanmean(totDist);
     
    distributions.distanceToCenter = distanceToCenter;
    means.distanceToCenter = nanmean(distanceToCenter);
    %% Portion of time on brood
%     thresh = 0.01; %1 cm (i.e. 1 body length)
%     onNest = double(minBroodDist < thresh);
%     onNest(isnan(broodDist)) = NaN;
%     porTimeOnNest = nanmean(onNest);

% Second try
    [porTimeOnNest onNest] = calculatePortionOfTimeOnNest(data, brood, 0.01);
    
    means.porTimeOnNest = nanmean(porTimeOnNest);
    distributions.porTimeOnNest = porTimeOnNest;